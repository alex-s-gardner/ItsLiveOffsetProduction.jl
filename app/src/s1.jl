# Sentinel-1 IW SLC pairs: the mosaic geometry geogrid is handed, the reference's bursts merged into the
# mosaic, the secondary resampled onto it burst by burst, and the product's metadata.
#
# Each secondary burst is coregistered onto its reference burst, as `s1_isce3.process_slc` does with
# COMPASS: here through an offset field solved on a lattice from the two orbits and a DEM, and an eight-tap
# sinc over the deramped secondary burst. The bursts are merged as ISCE2's topsApp merges them
# (`TopsProc/runMergeBursts.py`, `VRTManager.py`): each burst at its own sensing time, over the valid
# window its reference and coregistered secondary share, a later burst or subswath overwriting an earlier
# one where they overlap.

const S1Product = SLCDatasets.Sentinel1Product
const S1_SWATHS = 1:3

"""
    s1_polarization(name) -> String

The like-polarized channel the pipeline correlates, lowercase: `hh` for an `SH`/`DH` product, `vv` for an
`SV`/`DV` one.
"""
function s1_polarization(name::String)
    m = match(r"_1S(S|D)(H|V)_", name)
    m === nothing && throw(ArgumentError("cannot read a polarization out of " * name))
    return something(m.captures[2]) == "H" ? "hh" : "vv"
end

"""
    s1_orbit_file(name, orbit_dir, scratch) -> String

The precise orbit file for product `name`: the one in `orbit_dir` whose validity covers the product's
sensing start, when `orbit_dir` is given, else the one [`fetch_s1_orbit`](@ref) downloads into `scratch`.
"""
function s1_orbit_file(name::String, orbit_dir::String, scratch::String)
    isempty(orbit_dir) && return fetch_s1_orbit(name, scratch)
    _, start, _ = s1_name_times(name)
    for f in sort(readdir(orbit_dir))
        m = match(r"_V(\d{8}T\d{6})_(\d{8}T\d{6})\.EOF$", f)
        m === nothing && continue
        a, b = _ccsds(String(something(m.captures[1]))), _ccsds(String(something(m.captures[2])))
        a <= start <= b && return joinpath(orbit_dir, f)
    end
    throw(ArgumentError("no orbit file in " * orbit_dir * " covers " * name))
end

"""
    is_burst_list(input) -> Bool

Whether a burst pair's input is its bursts' extractor URLs rather than a product.
"""
is_burst_list(input::String) = SLCDatasets.is_asf_burst_url(strip(first(split(input, ','))))

"""
    burst_list(input, dir) -> Vector{SLCDatasets.AsfBurst}

The bursts a comma-separated list of extractor URLs names, each raster downloaded into `dir` when it is
first read.
"""
burst_list(input::String, dir::Union{Nothing,String}) =
    [SLCDatasets.AsfBurst(strip(u); dir) for u in split(input, ','; keepempty = false)]

# A scene's product, opened for `swaths`. Extracted bursts are downloaded, and a zip's rasters or bursts
# inflated, into `scratch`.
function open_s1(scene::SceneSource, swaths::Vector{Int}, orbit_dir::String, scratch::String)
    orbit = s1_orbit_file(scene.name, orbit_dir, scratch)
    is_burst_list(scene.band) && return SLCDatasets.Sentinel1Product(burst_list(scene.band, scratch); orbit, swaths)
    is_zip_burst_list(scene.band) &&
        return SLCDatasets.Sentinel1Product(zip_burst_list(scene.band, scratch); orbit, swaths)
    return SLCDatasets.Sentinel1Product(scene.band; orbit, polarization = s1_polarization(scene.name), swaths,
                                        dir = scratch)
end

"""
    product_swaths(path, polarization) -> Vector{Int}

The subswaths a product holds bursts of in `polarization`: all three for a full SLC, the bursts' own for a
burst product (`s1_isce3.py:58` processes the subswaths its bursts lie in). `burst2safe` writes an
annotation for every subswath, listing no bursts for one it holds none of.
"""
function product_swaths(path::String, polarization::String)
    return SLCDatasets.with_annotations(path) do r
        sws = Int[]
        for n in r.names
            m = match(r"-iw(\d)-slc-([a-z]{2})-", n)
            m === nothing && continue
            something(m.captures[2]) == polarization || continue
            sw = parse(Int, something(m.captures[1]))
            sw in sws && continue
            occursin("<burstList count=\"0\"", SLCDatasets.read_annotation_entry(r, sw, polarization)) && continue
            push!(sws, sw)
        end
        isempty(sws) && throw(ArgumentError(path * " has no annotation in polarization " * uppercase(polarization)))
        sort!(sws)
    end
end

# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------

"""
    MosaicBox(p, swaths)

The merged grid the subswaths are laid on, ISCE2's `mergeBox` (`TopsProc/runMergeBursts.py:55-82`): line 0
at the earliest subswath's first burst, sample 0 at the nearest subswath's starting range, and extents
reaching the latest burst's last line and the farthest subswath's last sample, at the earliest subswath's
line interval and the nearest's sample spacing.
"""
struct MosaicBox
    sensing_start::SLCDatasets.UtcTime
    starting_range::Float64
    dt::Float64
    dr::Float64
    nlines::Int
    nsamples::Int
end

function MosaicBox(p::S1Product, swaths)
    ann = [SLCDatasets.annotation(p, sw) for sw in swaths]
    t0 = first(ann[1].burst_start)
    top = argmin(i -> SLCDatasets.seconds_between(t0, first(ann[i].burst_start)), eachindex(ann))
    left = argmin(i -> ann[i].starting_range, eachindex(ann))
    start = first(ann[top].burst_start)
    dt = ann[top].azimuth_time_interval
    dr = ann[left].range_pixel_spacing
    near = ann[left].starting_range
    stop = maximum(a -> SLCDatasets.seconds_between(start, last(a.burst_start)) + (a.lines_per_burst - 1) * a.azimuth_time_interval, ann)
    far = maximum(a -> a.starting_range + (a.samples_per_burst - 1) * a.range_pixel_spacing, ann)
    return MosaicBox(start, near, dt, dr, round(Int, stop / dt + 1), round(Int, (far - near) / dr + 1))
end

"""
    s1_mosaic(p, swaths) -> (RadarCoordinate, UtcTime)

The merged radar grid geogrid is handed for product `p`: the [`MosaicBox`](@ref) with the PRF, wavelength
and orbit of the near subswath. Returns the absolute sensing start beside it, which the pair's interval is
taken from.
"""
function s1_mosaic(p::S1Product, swaths)
    sws = collect(swaths)
    box = MosaicBox(p, sws)
    ilo = argmin(i -> SLCDatasets.annotation(p, sws[i]).starting_range, eachindex(sws))
    lo = SLCDatasets.annotation(p, sws[ilo])
    base = SLCDatasets.merge_bursts(collect(SLCDatasets.bursts(p, sws[ilo])))
    c = radar_coordinate(base)
    base_start = first(base.backend.annotation.burst_start)
    kwargs = (; orbit = c.orbit, starting_range = box.starting_range, dr = box.dr,
              # `c.sensing_start` is seconds on the orbit's clock at the near subswath's first burst.
              sensing_start = c.sensing_start + SLCDatasets.seconds_between(base_start, box.sensing_start),
              prf = 1 / box.dt, nsamples = box.nsamples, nlines = box.nlines, look_side = c.look_side,
              wavelength = lo.wavelength, orbit_epoch_offset = c.orbit_epoch_offset)
    coord = ImagePairGeometry.RadarCoordinate(; kwargs..., incidence_angle = ImagePairGeometry.incidence_angle(; kwargs...))
    return coord, box.sensing_start
end

"""
    s1_setup(c, orbit_dir, scratch) -> RadarSetup

`setup` for a Sentinel-1 pair: the reference's mosaic geometry and the interval between the two
acquisitions' sensing starts — the secondary contributes nothing else, as `process_slc` copies the
reference's metadata and replaces only the sensing times. A full SLC pair mosaics all three subswaths, a
burst pair the subswaths the reference's bursts lie in.
"""
function s1_setup(c::Job, orbit_dir::String, scratch::String)
    early, late = acquisition_order(c)
    sws = c.platform == "S1-SLC" ? collect(S1_SWATHS) :
          is_burst_list(c.reference.band) ? sort!(unique([b.swath for b in burst_list(c.reference.band, nothing)])) :
          is_zip_burst_list(c.reference.band) ?
          sort!(unique([b.swath for b in zip_burst_list(c.reference.band, nothing)])) :
          product_swaths(c.reference.band, s1_polarization(c.reference.name))
    # The two products open together. Their rasters then download or inflate while the geometry is solved
    # and the DEM read — the reference's first, so its mosaic can be filled while the secondary's arrive.
    tl = Threads.@spawn open_s1(late, sws, orbit_dir, scratch)
    rp = open_s1(early, sws, orbit_dir, scratch)
    first_fetch = Threads.@spawn SLCDatasets.prefetch_rasters(rp)
    sp = _await(tl)::S1Product
    pending = Task[first_fetch, Threads.@spawn (wait(first_fetch); SLCDatasets.prefetch_rasters(sp))]
    ref, ref_start = s1_mosaic(rp, sws)
    _, sec_start = s1_mosaic(sp, sws)
    pair = CoregisteredPair(ref; dt = SLCDatasets.seconds_between(ref_start, sec_start))
    b = footprint_bounds(IdentityTransform(), pair.coordinate)
    region = parameter_region((b.X[1] + b.X[2]) / 2, (b.Y[1] + b.Y[2]) / 2)
    grid = parameter_grid(region)
    nd = with_dataset(band_nodata, region.paths[4])
    g = FGP.native_crs(region.epsg)
    g isa FGP.UTMCRS && return _radar_setup(c, early, late, rp, sp, pair, region, grid, nd, g; pending)
    g isa FGP.PolarStereographicCRS && return _radar_setup(c, early, late, rp, sp, pair, region, grid, nd, g; pending)
    throw(ArgumentError("a product grid is UTM or polar stereographic, not EPSG:" * string(region.epsg)))
end

# ---------------------------------------------------------------------------
# The merge
# ---------------------------------------------------------------------------

"""
    valid_window(a, i) -> NTuple{4,Int}

Burst `i`'s valid window in annotation `a`, zero-based first and last line and sample.

The window is ISCE2's (`Sensor/TOPS/Sentinel1.py:662-676`), which keeps every valid line but ends one
sample before the annotation's last valid sample: `numValidSamples = lastValidSample - firstValidSample`
with `lastValidSample` the annotation's zero-based inclusive index, read as a half-open range.
"""
function valid_window(a, i::Int)
    l, s = SLCDatasets.valid_lines(a, i), SLCDatasets.valid_samples(a, i)
    return (first(l) - 1, last(l) - 1, first(s) - 1, last(s) - 2)
end

"""
    coregistered_window(rwin, swin, dl, ds) -> NTuple{4,Int}

The window on the reference burst that the reference and its coregistered secondary share, zero-based
first and last line and sample: the part of the reference's [`valid_window`](@ref) `rwin` whose resampled
pixels take all eight sinc taps of [`resample!`](@ref) from inside the secondary's valid window `swin`.
This is the guarantee topsApp's `adjustValidSampleLine` (`TopsProc/runFineResamp.py:214-273`) trims for,
tested directly on the offsets `dl`, `ds`.

Each edge moves inward until its whole row or column is inside; the rows and columns behind it then are
too, as a position advances with its index. Rows are tested along all of `rwin`'s samples, columns along
the rows kept.
"""
function coregistered_window(rwin, swin, dl, ds)
    rfl, rll, rfs, rls = rwin
    sfl, sll, sfs, sls = swin
    # One-based source position `y` reads one-based rows `floor(y) - 3:floor(y) + 4`, inside zero-based
    # `sfl:sll` when `sfl + 4 <= y < sll - 2`. Columns likewise.
    inside(lo, hi, first, last) = lo >= first + 4 && hi < last - 2
    fl, ll = _shrink(rfl, rll) do l
        ymin, ymax = support_bounds(dl, ds, l:l, rfs:rls)
        inside(ymin, ymax, sfl, sll)
    end
    fl <= ll || return (fl, ll, rfs, rls)
    fs, ls = _shrink(rfs, rls) do s
        _, _, xmin, xmax = support_bounds(dl, ds, fl:ll, s:s)
        inside(xmin, xmax, sfs, sls)
    end
    return (fl, ll, fs, ls)
end

# `a:b` narrowed from each end until `keep` holds at both ends.
function _shrink(keep, a::Int, b::Int)
    while a <= b && !keep(a)
        a += 1
    end
    while a <= b && !keep(b)
        b -= 1
    end
    return a, b
end

"""
    MergePiece

One burst's part of the mosaic, `VRTConstructor.addBurst` with `validOnly` (`TopsProc/VRTManager.py:165-205`):
the burst's valid window, `lines` x `samples` zero-based on its own grid, lands at `mosaic_rows` x
`mosaic_cols`, offset from those indices by the burst's sensing time and starting range on the
[`MosaicBox`](@ref). `swath` and `burst` name it.
"""
struct MergePiece
    swath::Int
    burst::Int
    lines::UnitRange{Int}
    samples::UnitRange{Int}
    mosaic_rows::UnitRange{Int}
    mosaic_cols::UnitRange{Int}
end

"""
    merge_pieces(rp, bursts) -> (Vector{MergePiece}, MosaicBox)

Every burst of the reference's subswaths in the order `mergeBurstsVirtual` adds them, subswath by subswath
and burst by burst, so that a later piece overwrites an earlier one where they overlap. The windows are the
reference's and coregistered secondary's shared ones ([`coregistered_window`](@ref)) over the secondary's
`bursts` ([`solve_secondary_bursts`](@ref)); both mosaics use the same pieces, so the reference and the
secondary carry the same burst in every pixel.
"""
function merge_pieces(rp::S1Product, bursts)
    box = MosaicBox(rp, rp.swaths)
    pieces = MergePiece[]
    for sw in rp.swaths
        ra = SLCDatasets.annotation(rp, sw)
        xoff = round(Int, (ra.starting_range - box.starting_range) / box.dr)
        for i in 1:SLCDatasets.nbursts(ra)
            yoff = round(Int, SLCDatasets.seconds_between(box.sensing_start, ra.burst_start[i]) / box.dt)
            p = bursts[(sw, i)]
            fl, ll, fs, ls = coregistered_window(valid_window(ra, i), p.svalid, p.dl, p.ds)
            (fl <= ll && fs <= ls) || continue
            push!(pieces, MergePiece(sw, i, fl:ll, fs:ls, (yoff + fl + 1):(yoff + ll + 1), (xoff + fs + 1):(xoff + ls + 1)))
        end
    end
    return pieces, box
end

"""
    BurstMosaic <: AbstractMatrix{Float32}

Bursts merged as topsApp merges them, read a window at a time: each [`MergePiece`](@ref) in order, a later
one overwriting an earlier one, as GDAL renders the sources of a VRT band. `read!(dest, burst, lines,
samples)` writes the amplitude of burst `burst`'s zero-based `lines` x `samples` into `dest`.
"""
struct BurstMosaic{B} <: AbstractMatrix{Float32}
    bursts::B
    pieces::Vector{MergePiece}
    dims::Tuple{Int,Int}
end

Base.size(m::BurstMosaic) = m.dims

function Base.getindex(m::BurstMosaic, rows::AbstractUnitRange, cols::AbstractUnitRange)
    checkbounds(m, rows, cols)
    out = zeros(Float32, length(rows), length(cols))
    # The part `inner` of range `outer` re-expressed on the range `target` that `outer` maps to.
    shift(inner, outer, target) = (inner .- first(outer)) .+ first(target)
    for q in m.pieces
        mr, mc = intersect(rows, q.mosaic_rows), intersect(cols, q.mosaic_cols)
        (isempty(mr) || isempty(mc)) && continue
        dest = view(out, shift(mr, rows, 1:length(rows)), shift(mc, cols, 1:length(cols)))
        read_burst!(dest, m.bursts, q, shift(mr, q.mosaic_rows, q.lines), shift(mc, q.mosaic_cols, q.samples))
    end
    return out
end

@_window_indexing BurstMosaic

"""
    ReferenceBursts

The reference's bursts, each subswath's raster stacking its bursts `lines_per_burst` rows apart.
"""
struct ReferenceBursts{R}
    rasters::Dict{Int,R}
    lines_per_burst::Dict{Int,Int}
end

function read_burst!(dest, r::ReferenceBursts, q::MergePiece, lines::UnitRange{Int}, samples::UnitRange{Int})
    rows = ((q.burst - 1) * r.lines_per_burst[q.swath] + 1) .+ lines
    dest .= abs.(Matrix{Complex{Int16}}(r.rasters[q.swath][rows, first(samples) + 1:last(samples) + 1]))
    return dest
end

function reference_mosaic(rp::S1Product, pieces, box::MosaicBox, raster_of)
    rasters = Dict(sw => raster_of(sw) for sw in rp.swaths)
    lpb = Dict(sw => SLCDatasets.annotation(rp, sw).lines_per_burst for sw in rp.swaths)
    return BurstMosaic(ReferenceBursts(rasters, lpb), pieces, (box.nlines, box.nsamples))
end

# ---------------------------------------------------------------------------
# The TOPS carrier
# ---------------------------------------------------------------------------

const CLIGHT = 299792458.0

# A `RangePolynomial` at slant `range`: its coefficients are in powers of two-way slant-range time.
range_poly(p::SLCDatasets.RangePolynomial, range::Float64) = evalpoly(2 * (range - p.r0) / CLIGHT, p.coeffs)

"""
    TopsCarrier

One burst's azimuth carrier phase, `s1reader.az_carrier_components`: azimuth time from the middle line
index, a reference time that is a difference of Doppler-over-FM-rate terms, the quadratic
`pi kt (eta - eta_ref)^2`, and the linear Doppler term `ResampSlc` applies beside it.
"""
struct TopsCarrier
    fm::SLCDatasets.RangePolynomial
    dc::SLCDatasets.RangePolynomial
    ks::Float64
    dt::Float64
    r0::Float64
    dr::Float64
    mid_line::Int
    eta_ref0::Float64
end

function TopsCarrier(slc, coord, lines_per_burst::Int)
    dp = SLCDatasets.deramp_parameters(slc)
    _, v = ImagePairGeometry.interpolate(coord.orbit, ImagePairGeometry.orbit_midtime(coord))
    # The annotation gives the steering rate in degrees per second; `s1reader` holds it in radians.
    ks_rad = dp.azimuth_steering_rate * (abs(dp.azimuth_steering_rate) > 0.1 ? pi / 180 : 1.0)
    ks = ks_rad * 2 * sqrt(v[1] * v[1] + v[2] * v[2] + v[3] * v[3]) / dp.wavelength
    r0 = dp.starting_range
    eta_ref0 = range_poly(dp.doppler_centroid, r0) / range_poly(dp.azimuth_fm_rate, r0)
    return TopsCarrier(dp.azimuth_fm_rate, dp.doppler_centroid, ks, dp.azimuth_time_interval, r0,
                       dp.range_pixel_spacing, lines_per_burst ÷ 2, eta_ref0)
end

function (c::TopsCarrier)(line::Int, sample::Int)
    r = c.r0 + sample * c.dr
    ka = range_poly(c.fm, r)
    fdc = range_poly(c.dc, r)
    kt = c.ks / (1 - c.ks / ka)
    eta = (line - c.mid_line) * c.dt
    de = eta - (c.eta_ref0 - fdc / ka)
    return pi * kt * de * de + 2pi * fdc * de
end

"""
    CarrierFit

The carrier as the polynomial ISCE3 evaluates (`get_az_carrier_poly`): sampled every 50 lines and 500
samples and least-squares fitted in normalized range and azimuth time over the terms `x^j y^i` with
`i <= 5`, `j <= 3`, `i + j <= 5`, against the coordinates of the pixel one line and one sample on.
"""
struct CarrierFit
    table::Matrix{Float64}
    xmin::Float64
    xnorm::Float64
    ymin::Float64
    ynorm::Float64
    r0::Float64
    dr::Float64
    dt::Float64
end

const CARRIER_TERMS = [(i, j) for i in 0:5 for j in 0:3 if i + j <= 5]

function CarrierFit(c::TopsCarrier, lines::Int, samples::Int; ystep::Int = 50, xstep::Int = 500)
    xs = 0:xstep:(samples - 1)
    ys = 0:ystep:(lines - 1)
    n = length(xs) * length(ys)
    rg = [c.r0 + (x + 1) * c.dr for x in xs]
    az = [(y + 1) * c.dt for y in ys]
    xmin, xnorm = minimum(rg), max(maximum(rg) - minimum(rg), 1.0)
    ymin, ynorm = minimum(az), max(maximum(az) - minimum(az), 1.0)
    A = Matrix{Float64}(undef, n, length(CARRIER_TERMS))
    z = Vector{Float64}(undef, n)
    k = 0
    for (jy, y) in enumerate(ys), (jx, x) in enumerate(xs)
        k += 1
        xn = (rg[jx] - xmin) / xnorm
        yn = (az[jy] - ymin) / ynorm
        for (t, (i, j)) in enumerate(CARRIER_TERMS)
            A[k, t] = xn^j * yn^i
        end
        z[k] = c(y, x)
    end
    coef = A \ z
    table = zeros(Float64, 6, 4)
    for (t, (i, j)) in enumerate(CARRIER_TERMS)
        table[i + 1, j + 1] = coef[t]
    end
    return CarrierFit(table, xmin, xnorm, ymin, ynorm, c.r0, c.dr, c.dt)
end

# The azimuth polynomial's coefficients at one range sample: the fit is separable, so the range half is
# evaluated once per column.
function carrier_column(f::CarrierFit, sample::Int)
    xn = ((f.r0 + sample * f.dr) - f.xmin) / f.xnorm
    T = f.table
    return ntuple(i -> T[i, 1] + xn * (T[i, 2] + xn * (T[i, 3] + xn * T[i, 4])), Val(6))
end

function carrier_phase(b::NTuple{6,Float64}, f::CarrierFit, line::Int)
    yn = (line * f.dt - f.ymin) / f.ynorm
    return ((((b[6] * yn + b[5]) * yn + b[4]) * yn + b[3]) * yn + b[2]) * yn + b[1]
end

# ---------------------------------------------------------------------------
# The secondary on the reference grid
# ---------------------------------------------------------------------------

"""
    SecondaryBurst

One secondary burst resampled onto its reference burst: its rows in the stacked raster, its carrier, the
offset field onto the reference burst, its own valid rectangle (zero-based first/last line and sample) and
its own `(lines, samples)`, which bounds what its source offers.
"""
struct SecondaryBurst
    source_rows::UnitRange{Int}
    carrier::CarrierFit
    dl::LatticeOffsets
    ds::LatticeOffsets
    svalid::NTuple{4,Int}
    extent::Tuple{Int,Int}
end

"""
    solve_secondary_bursts(rp, sp, dem) -> Dict{Tuple{Int,Int},SecondaryBurst}

Every secondary burst's [`SecondaryBurst`](@ref), keyed by subswath and burst, from the annotations, the
orbits and `dem` alone: no raster is read.
"""
function solve_secondary_bursts(rp::S1Product, sp::S1Product, dem::D) where {D}
    for sw in rp.swaths
        SLCDatasets.nbursts(SLCDatasets.annotation(sp, sw)) >= SLCDatasets.nbursts(SLCDatasets.annotation(rp, sw)) ||
            throw(ArgumentError("the secondary's IW" * string(sw) * " has fewer bursts than the reference's"))
    end
    jobs = [(sw, i) for sw in rp.swaths for i in 1:SLCDatasets.nbursts(SLCDatasets.annotation(rp, sw))]
    tasks = map(jobs) do (sw, i)
        Threads.@spawn begin
            a, sa = SLCDatasets.annotation(rp, sw), SLCDatasets.annotation(sp, sw)
            # The source burst's shape is the secondary's own, which need not be the reference's: it sets
            # where burst `i` starts in the stacked raster, the carrier's middle line, and every bound
            # `ResampSlc` takes off its source grid.
            slpb, sspb = sa.lines_per_burst, sa.samples_per_burst
            ss = SLCDatasets.bursts(sp, sw)[i]
            cs = radar_coordinate(ss)
            carrier = CarrierFit(TopsCarrier(ss, cs, slpb), slpb, sspb)
            dl, ds = offset_lattice(radar_coordinate(SLCDatasets.bursts(rp, sw)[i]), cs, dem,
                                    a.lines_per_burst, a.samples_per_burst)
            SecondaryBurst(((i - 1) * slpb + 1):(i * slpb), carrier, dl, ds, valid_window(sa, i), (slpb, sspb))
        end
    end
    return Dict{Tuple{Int,Int},SecondaryBurst}(j => fetch(t)::SecondaryBurst for (j, t) in zip(jobs, tasks))
end

"""
    SecondaryBursts

The secondary's bursts, each resampled onto its reference burst on demand, from its subswath's raster.
"""
struct SecondaryBursts{R}
    sources::Dict{Int,R}
    bursts::Dict{Tuple{Int,Int},SecondaryBurst}
end

"""
    deramped_burst!(buf, raster, rows, cols, carrier, origin) -> view

The samples `raster[rows, cols]` with the burst's azimuth carrier removed, written into `buf`. `origin` is
the burst-local zero-based `(line, sample)` of the band's first element, which is what the carrier is a
function of. The phase stays `Float64` up to the trig.
"""
function deramped_burst!(buf::Matrix{ComplexF32}, raster, rows::UnitRange{Int}, cols::UnitRange{Int},
                         carrier::CarrierFit, origin::NTuple{2,Int})
    nr, nc = length(rows), length(cols)
    out = view(buf, 1:nr, 1:nc)
    band = Matrix{Complex{Int16}}(raster[rows, cols])
    l0, s0 = origin
    Threads.@threads for j in 1:nc
        b = carrier_column(carrier, s0 + j - 1)
        for i in 1:nr
            out[i, j] = ComplexF32(band[i, j]) * ComplexF32(cis(-carrier_phase(b, carrier, l0 + i - 1)))
        end
    end
    return out
end

# Burst `q`'s resampled amplitude over its zero-based output `lines` x `samples`, into `dest`. The band read
# is the support those outputs reach, clamped to the burst. The pieces' windows keep every tap inside the
# secondary's valid window ([`coregistered_window`](@ref)), so the samples outside it are never read.
function read_burst!(dest::AbstractMatrix{Float32}, r::SecondaryBursts, q::MergePiece, lines::UnitRange{Int},
                     samples::UnitRange{Int})
    p = r.bursts[(q.swath, q.burst)]
    slpb, sspb = p.extent
    ymin, ymax, xmin, xmax = support_bounds(p.dl, p.ds, lines, samples)
    brows = clamp(floor(Int, ymin) - 3, 1, slpb):clamp(ceil(Int, ymax) + 4, 1, slpb)
    bcols = clamp(floor(Int, xmin) - 3, 1, sspb):clamp(ceil(Int, xmax) + 4, 1, sspb)
    (isempty(brows) || isempty(bcols)) && return dest
    origin = (first(brows) - 1, first(bcols) - 1)
    srows = (first(p.source_rows) - 1 + first(brows)):(first(p.source_rows) - 1 + last(brows))
    scratch = Matrix{ComplexF32}(undef, length(brows), length(bcols))
    src = deramped_burst!(scratch, r.sources[q.swath], srows, bcols, p.carrier, origin)
    return resample!(dest, src, p.dl, p.ds, lines, samples; origin, extent = p.extent, doppler = p.extent)
end

function secondary_mosaic(rp::S1Product, bursts, pieces, box::MosaicBox, raster_of)
    sources = Dict(sw => raster_of(sw) for sw in rp.swaths)
    return BurstMosaic(SecondaryBursts(sources, bursts), pieces, (box.nlines, box.nsamples))
end

# A product's subswath measurement raster: a file of a `.SAFE`, a zip's member inflated, or assembled bursts.
function _swath_raster(p::S1Product, sw::Int)
    return SLCDatasets.burst_raster(SLCDatasets.bursts(p, sw)[1].backend).raster
end

"""
    s1_imagery(s, dem; dir, cache_io) -> Imagery

The correlator's two images for a Sentinel-1 SLC pair, each filled into a tile cache under `dir`: the
reference's mosaic, and the secondary resampled onto it. A product read from a `.SAFE` or a zip and one
assembled from extracted bursts read their rasters through different types, so each has its own chain.
"""
function s1_imagery(s::RadarSetup, dem::D; dir::String, cache_io::Symbol) where {D}
    rp, sp = s.early_product, s.late_product
    # Decided without reading a raster, so the reference's mosaic is filled while the secondary's rasters
    # are still being fetched; only the two products' annotations, orbits and the DEM come before it.
    assembled = !isempty(rp.bursts)
    assembled == !isempty(sp.bursts) || throw(ArgumentError(
        "the two Sentinel-1 products must be given alike: both as products, or both as bursts"))
    assembled && return _s1_imagery(s, dem, sw -> _swath_raster(rp, sw)::SLCDatasets.BurstStack,
                                    sw -> _swath_raster(sp, sw)::SLCDatasets.BurstStack; dir, cache_io)
    return _s1_imagery(s, dem, sw -> _swath_raster(rp, sw)::SLCDatasets.MeasurementRaster,
                       sw -> _swath_raster(sp, sw)::SLCDatasets.MeasurementRaster; dir, cache_io)
end

function _s1_imagery(s::RadarSetup, dem::D, rraster::F1, sraster::F2; dir::String, cache_io::Symbol) where {D,F1,F2}
    g = s.geometry
    lx, ly = permutedims(g.location_x), permutedims(g.location_y)
    sentinel = Int32(g.nodata.output)
    declined = _sentinel_declined(lx, ly, sentinel)
    rp, sp = s.early_product, s.late_product
    # The secondary's offsets set both mosaics' windows, so they are solved before either is read.
    bursts = solve_secondary_bursts(rp, sp, dem)
    pieces, box = merge_pieces(rp, bursts)
    ref = reference_mosaic(rp, pieces, box, rraster)
    co = s.pair.coordinate
    (co.nlines, co.nsamples) == size(ref) || throw(DimensionMismatch(
        "the geogrid was built on a different mosaic than the imagery"))
    early = fill_tiles!(TileCache(ref; tile = 1024, dir, io = cache_io))
    _accumulate_zeros!(declined, early, lx, ly, sentinel, NODATA_SAMPLE_OFFSET)
    sec = secondary_mosaic(rp, bursts, pieces, box, sraster)
    size(sec) == size(ref) || throw(DimensionMismatch("the reference and secondary mosaics differ in size"))
    late = fill_tiles!(TileCache(sec; tile = 1024, dir, io = cache_io))
    _accumulate_zeros!(declined, late, lx, ly, sentinel, NODATA_SAMPLE_OFFSET)
    return Imagery(early, late, declined)
end

# ---------------------------------------------------------------------------
# Packaging
# ---------------------------------------------------------------------------

# A `UtcTime` as a `DateTime`, its fraction rounded to the millisecond.
s1_datetime(t::SLCDatasets.UtcTime) = t.datetime + Dates.Millisecond(round(Int, t.seconds * 1000))

# A product's id: its file name without `.SAFE` or `.zip`.
function _s1_name(p::S1Product)
    n = basename(rstrip(p.path, '/'))
    endswith(n, ".SAFE") && return String(n[1:(end - 5)])
    endswith(n, ".zip") && return String(n[1:(end - 4)])
    return String(n)
end

"""
    s1_img_pair_info(s, roi_valid, lon, lat) -> ItsLiveOffsetProduction.ImagePairInfo

`IMG_INFO_DICT` for a Sentinel-1 pair: the satellites, orbits, datatakes and product ids out of the two
product ids (`S1A_IW_SLC__1SSH_<start>_<stop>_<orbit>_<datatake>_<id>`), the pass directions out of the
annotations, and the acquisition times as the two mosaics' sensing starts.
"""
function s1_img_pair_info(s::RadarSetup, roi_valid::Float64, lon::Float64, lat::Float64)
    rp, sp = s.early_product, s.late_product
    t1 = MosaicBox(rp, rp.swaths).sensing_start
    t2 = MosaicBox(sp, rp.swaths).sensing_start
    a1, a2 = SLCDatasets.annotation(rp, first(rp.swaths)), SLCDatasets.annotation(sp, first(rp.swaths))
    n1, n2 = _s1_name(rp), _s1_name(sp)
    f1, f2 = split(n1, '_'), split(n2, '_')
    (length(f1) == 10 && length(f2) == 10) || throw(ArgumentError(
        n1 * " or " * n2 * " is not a Sentinel-1 product id of ten underscore-delimited fields"))
    # Sized for its twelve entries, so the attributes iterate in the order a twelve-pair `Dict` gives.
    extra = sizehint!(Dict{String,AttrValue}(), 12)
    extra["id_img1"] = n1
    extra["id_img2"] = n2
    extra["sensor_img1"] = "C"
    extra["sensor_img2"] = "C"
    extra["flight_direction_img1"] = lowercase(a1.pass_direction)
    extra["flight_direction_img2"] = lowercase(a2.pass_direction)
    extra["absolute_orbit_number_img1"] = String(f1[8])
    extra["absolute_orbit_number_img2"] = String(f2[8])
    extra["mission_data_take_ID_img1"] = String(f1[9])
    extra["mission_data_take_ID_img2"] = String(f2[9])
    extra["product_unique_ID_img1"] = String(f1[10])
    extra["product_unique_ID_img2"] = String(f2[10])
    return ItsLiveOffsetProduction.ImagePairInfo(s1_datetime(t1), s1_datetime(t2), "S", "S",
                                                 String(f1[1][2:end]), String(f2[1][2:end]), roi_valid,
                                                 round(lat; digits = 2), round(lon; digits = 2), extra)
end

"""
    s1_swath_bias(s, c) -> Union{ItsLiveOffsetProduction.SwathOffsetBias,Nothing}

The subswath offset-bias correction's inputs for a full-SLC pair of two different spacecraft, `nothing` for
one spacecraft, where there is no antenna difference to remove. A burst pair has no correction.
"""
function s1_swath_bias(s::RadarSetup, c::Correlation)
    bd = SLCDatasets.subswath_borders(s.early_product, s.late_product)
    bd.same_platform && return nothing
    lx, ly = ItsLiveOffsetProduction.image_location(s.geometry)
    return ItsLiveOffsetProduction.SwathOffsetBias(lx, ly, c.grid_spacing_x, bd.ncols, Float64(bd.border12),
                                                   Float64(bd.border23), bd.same_platform, bd.reference_platform)
end

_radar_info(s::RadarSetup{P,G,S1Product}, c::Correlation, roi::Float64, lon::Float64, lat::Float64) where {P,G} =
    (s1_img_pair_info(s, roi, lon, lat), s.job.platform == "S1-SLC" ? s1_swath_bias(s, c) : nothing)
