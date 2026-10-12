# Sentinel-1 IW SLC pairs: the mosaic geometry geogrid is handed, the reference's bursts merged into the
# `merge_swaths` mosaic, the secondary resampled onto it burst by burst, and the product's metadata.
#
# The pipeline being reproduced is `s1_isce3.process_slc`: COMPASS coregisters each secondary burst onto
# its reference burst, `merge_bursts_in_swath` and `merge_swaths` stack the bursts and subswaths, and
# `read_slc_gdal` takes the amplitude. Here the coregistration is an offset field solved on a lattice from
# the two orbits and a DEM, and the resample an eight-tap sinc over the deramped secondary burst.

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

# A scene's product, opened for `swaths`. Extracted bursts are downloaded, and a zip's rasters inflated,
# into `scratch`.
function open_s1(scene::SceneSource, swaths::Vector{Int}, orbit_dir::String, scratch::String)
    orbit = s1_orbit_file(scene.name, orbit_dir, scratch)
    is_burst_list(scene.band) && return SLCDatasets.Sentinel1Product(burst_list(scene.band, scratch); orbit, swaths)
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
    s1_mosaic(p, swaths) -> (RadarCoordinate, UtcTime)

The merged radar grid geogrid is handed for product `p` (`loadMetadataSlc`): the range origin, spacing,
PRF, wavelength and orbit of the near subswath; the sample count to the far subswath's far edge plus the
floored range offset to it; and the line count `merge_swaths` lays out, which adds a whole merged subswath
length to the last burst's start. Returns the absolute sensing start beside it, which the pair's interval
is taken from.
"""
function s1_mosaic(p::S1Product, swaths)
    sws = collect(swaths)
    ann = [SLCDatasets.annotation(p, sw) for sw in sws]
    shape = [(a.lines_per_burst, a.samples_per_burst) for a in ann]
    ilo = argmin(i -> ann[i].starting_range, eachindex(ann))
    ihi = argmax(i -> ann[i].starting_range, eachindex(ann))
    lo, hi = ann[ilo], ann[ihi]
    dr = lo.range_pixel_spacing
    adt = lo.azimuth_time_interval
    swath_lines = [1 + round(Int, (SLCDatasets.seconds_between(first(a.burst_start), last(a.burst_start)) +
                                   (shape[k][1] - 1) * a.azimuth_time_interval) / a.azimuth_time_interval)
                   for (k, a) in enumerate(ann)]
    start = first(ann[1].burst_start)
    for a in ann
        SLCDatasets.seconds_between(start, first(a.burst_start)) < 0 && (start = first(a.burst_start))
    end
    span = -Inf
    for (a, n) in zip(ann, swath_lines)
        span = max(span, SLCDatasets.seconds_between(start, last(a.burst_start)) + (n - 1) * a.azimuth_time_interval)
    end
    nlines = 1 + round(Int, span / adt)
    nsamples = shape[ihi][2] + floor(Int, (hi.starting_range - lo.starting_range) / dr)
    base = SLCDatasets.merge_bursts(collect(SLCDatasets.bursts(p, sws[ilo])))
    c = radar_coordinate(base)
    base_start = first(base.backend.annotation.burst_start)
    kwargs = (; orbit = c.orbit, starting_range = lo.starting_range, dr,
              # `c.sensing_start` is seconds on the orbit's clock at the near subswath's first burst.
              sensing_start = c.sensing_start + SLCDatasets.seconds_between(base_start, start),
              prf = 1 / adt, nsamples, nlines, look_side = c.look_side, wavelength = lo.wavelength,
              orbit_epoch_offset = c.orbit_epoch_offset)
    coord = ImagePairGeometry.RadarCoordinate(; kwargs..., incidence_angle = ImagePairGeometry.incidence_angle(; kwargs...))
    return coord, start
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
# The reference mosaic
# ---------------------------------------------------------------------------

# `bursts[-1].last_valid_line` of a subswath, 0-based.
_lvl_last(a) = a.last_valid_line[end] - 1

"""
    SwathAmplitude <: AbstractMatrix{Float32}

One subswath's bursts merged as `merge_bursts_in_swath` merges them, read a window at a time: the azimuth
seam halfway through each overlap, the first burst's valid samples applied to every burst. `pieces` is
the merge as `(merged rows, raster rows, columns)` placements, applied in order, so that a line two
bursts both write keeps the later one.
"""
struct SwathAmplitude{R} <: AbstractMatrix{Float32}
    raster::R
    pieces::Vector{NTuple{3,UnitRange{Int}}}
    dims::Tuple{Int,Int}
end

function SwathAmplitude(p::S1Product, swath::Int, raster::R) where {R}
    a = SLCDatasets.annotation(p, swath)
    n = SLCDatasets.nbursts(a)
    lpb, spb = a.lines_per_burst, a.samples_per_burst
    dt = a.azimuth_time_interval
    burst_rows(i, rows) = ((i - 1) * lpb) .+ rows
    n == 1 && return SwathAmplitude(raster, [(1:lpb, burst_rows(1, 1:lpb), 1:spb)], (lpb, spb))
    fvl = a.first_valid_line .- 1
    lvl = a.last_valid_line .- 1
    fvs = a.first_valid_sample .- 1
    lvs = a.last_valid_sample .- 1
    lims = [(s = round(Int, (SLCDatasets.seconds_between(first(a.burst_start), a.burst_start[i]) + fvl[i] * dt) / dt);
             (s, s + (lvl[i] - fvl[i]) + 1)) for i in 1:n]
    nlines = 1 + round(Int, (SLCDatasets.seconds_between(first(a.burst_start), last(a.burst_start)) +
                             (lpb - 1) * dt) / dt)
    pieces = map(1:n) do i
        prev = i > 1 ? fld(lims[i - 1][2] - lims[i][1], 2) : 0
        nxt = i < n ? fld(lims[i][2] - lims[i + 1][1], 2) : 0
        bstart, bend = fvl[i] + prev, 1 + lvl[i] - nxt
        mstart, mend = lims[i][1] + prev, lims[i][2] - nxt
        ((mstart + 1):mend, burst_rows(i, (bstart + 1):bend), (fvs[i] + 1):lvs[i])
    end
    return SwathAmplitude(raster, pieces, (nlines, spb))
end

Base.size(a::SwathAmplitude) = a.dims

function Base.getindex(a::SwathAmplitude, rows::AbstractUnitRange, cols::AbstractUnitRange)
    checkbounds(a, rows, cols)
    out = zeros(Float32, length(rows), length(cols))
    for (mr, rr, cc) in a.pieces
        wr, wc = intersect(rows, mr), intersect(cols, cc)
        (isempty(wr) || isempty(wc)) && continue
        src = rr[(first(wr) - first(mr) + 1):(last(wr) - first(mr) + 1)]
        out[wr .- (first(rows) - 1), wc .- (first(cols) - 1)] = abs.(Matrix{Complex{Int16}}(a.raster[src, wc]))
    end
    return out
end

@_window_indexing SwathAmplitude

"""
    SwathPlacement

Which window of one merged subswath goes where in the `merge_swaths` mosaic: the azimuth window is
cropped and moved to the subswath's azimuth offset, the range window cropped and kept at its own column
indices (`s1_isce3.py:487-500`).
"""
struct SwathPlacement
    swath::Int
    mosaic_rows::UnitRange{Int}
    mosaic_cols::UnitRange{Int}
    swath_rows::UnitRange{Int}
    swath_cols::UnitRange{Int}
end

"""
    mosaic_layout(p, swaths, shapes) -> (Vector{SwathPlacement}, (rows, cols))

Where each merged subswath of shape `shapes[k]` lands in the mosaic, and the mosaic's size. Every subswath
but the far one gives up 64 samples at its far edge, `invalid_pixel_buffer`.
"""
function mosaic_layout(p::S1Product, swaths, shapes::Vector{Tuple{Int,Int}})
    sws = collect(swaths)
    ann = [SLCDatasets.annotation(p, sw) for sw in sws]
    dr = first(ann).range_pixel_spacing
    dt = first(ann).azimuth_time_interval
    starts = [first(a.burst_start) for a in ann]
    sensing_start = starts[1]
    for s in starts
        SLCDatasets.seconds_between(sensing_start, s) < 0 && (sensing_start = s)
    end
    span = -Inf
    for (a, sh) in zip(ann, shapes)
        span = max(span, SLCDatasets.seconds_between(sensing_start, last(a.burst_start)) + (sh[1] - 1) * a.azimuth_time_interval)
    end
    total_az = 1 + round(Int, span / dt)
    rng_offsets = [sw == first(sws) ? 0 :
                   floor(Int, (SLCDatasets.annotation(p, sw).starting_range - first(ann).starting_range) / dr)
                   for sw in sws]
    total_rng = last(shapes)[2] + floor(Int, (last(ann).starting_range - first(ann).starting_range) / dr)
    places = map(eachindex(sws)) do k
        a, nrows = ann[k], shapes[k][1]
        fvl, lvl = a.first_valid_line[1] - 1, a.last_valid_line[1] - 1
        fvs, lvs = a.first_valid_sample[1] - 1, a.last_valid_sample[1] - 1
        az_offset = floor(Int, SLCDatasets.seconds_between(sensing_start, starts[k]) / dt)
        rng_offset = rng_offsets[k] + fvs
        rng_end = rng_offset + (lvs - fvs)
        buffer = sws[k] == maximum(sws) ? 0 : 64
        slc_az_end = SLCDatasets.nbursts(a) > 1 ? nrows - _lvl_last(a) : lvl
        merged_az_end = SLCDatasets.nbursts(a) > 1 ? az_offset + nrows - _lvl_last(a) - fvl : az_offset + (lvl - fvl)
        SwathPlacement(sws[k], (az_offset + 1):merged_az_end, (rng_offset + 1):(rng_end - buffer),
                       (fvl + 1):slc_az_end, (fvs + 1):(lvs - buffer))
    end
    return (places, (total_az, total_rng))
end

"""
    S1Mosaic <: AbstractMatrix{Float32}

Merged subswaths stacked as `merge_swaths` stacks them, read a window at a time, first-come in subswath
order: a pixel already non-zero is not overwritten.
"""
struct S1Mosaic{S} <: AbstractMatrix{Float32}
    swaths::Vector{S}
    places::Vector{SwathPlacement}
    dims::Tuple{Int,Int}
end

Base.size(m::S1Mosaic) = m.dims

function Base.getindex(m::S1Mosaic, rows::AbstractUnitRange, cols::AbstractUnitRange)
    checkbounds(m, rows, cols)
    out = zeros(Float32, length(rows), length(cols))
    for (k, q) in enumerate(m.places)
        mr = intersect(rows, q.mosaic_rows)
        isempty(mr) && continue
        mc = intersect(cols, q.mosaic_cols)
        isempty(mc) && continue
        r0 = first(q.swath_rows) + (first(mr) - first(q.mosaic_rows))
        c0 = first(q.swath_cols) + (first(mc) - first(q.mosaic_cols))
        got = m.swaths[k][r0:(r0 + length(mr) - 1), c0:(c0 + length(mc) - 1)]
        i0, j0 = first(mr) - first(rows), first(mc) - first(cols)
        for j in axes(got, 2), i in axes(got, 1)
            v = got[i, j]
            (out[i0 + i, j0 + j] == 0 && v != 0) && (out[i0 + i, j0 + j] = v)
        end
    end
    return out
end

@_window_indexing S1Mosaic

function reference_mosaic(p::S1Product, raster_of)
    sws = p.swaths
    lazy = [SwathAmplitude(p, sw, raster_of(sw)) for sw in sws]
    places, dims = mosaic_layout(p, sws, [size(x) for x in lazy])
    return S1Mosaic(lazy, places, dims)
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
    BurstPlacement

One secondary burst's part of the merged subswath: the burst, its rows in the stacked raster, its carrier,
the offset field onto the reference burst, its own valid rectangle (zero-based first/last line and
sample), the reference burst's zero-based output `lines` and `samples` it covers, the merged subswath's
one-based rows and columns they land on — positionally, so `mosaic_rows[t]` carries `lines[t]` — and the
secondary burst's own `(lines, samples)`, which bounds what its source offers.
"""
struct BurstPlacement
    burst::Int
    source_rows::UnitRange{Int}
    carrier::CarrierFit
    dl::LatticeOffsets
    ds::LatticeOffsets
    svalid::NTuple{4,Int}
    lines::UnitRange{Int}
    samples::UnitRange{Int}
    mosaic_rows::UnitRange{Int}
    mosaic_cols::UnitRange{Int}
    extent::Tuple{Int,Int}
end

"""
    ResampledSwath <: AbstractMatrix{Float32}

The secondary's amplitude on the reference's merged subswath grid, resampled a window at a time. The
placement is the reference's own seam arithmetic applied to the secondary, as `merge_bursts_in_swath`
applies it. Each burst's carrier and offset field are built once, at construction.
"""
struct ResampledSwath{R} <: AbstractMatrix{Float32}
    source::R
    places::Vector{BurstPlacement}
    dims::Tuple{Int,Int}
end

function ResampledSwath(rp::S1Product, sp::S1Product, swath::Int, dem::D, sraster::R) where {D,R}
    a = SLCDatasets.annotation(rp, swath)
    sa = SLCDatasets.annotation(sp, swath)
    n = SLCDatasets.nbursts(a)
    lpb, spb = a.lines_per_burst, a.samples_per_burst
    # The source burst's shape is the secondary's own, which need not be the reference's: it sets where
    # burst `i` starts in the stacked raster, the carrier's middle line, and every bound `ResampSlc` takes
    # off its source grid (`s1_resample.py` builds the radar grid, Doppler LUT and carrier from the
    # secondary burst).
    slpb, sspb = sa.lines_per_burst, sa.samples_per_burst
    dt = a.azimuth_time_interval
    fvl = a.first_valid_line .- 1
    lvl = a.last_valid_line .- 1
    fvs = a.first_valid_sample .- 1
    lvs = a.last_valid_sample .- 1
    lims = [(s = round(Int, (SLCDatasets.seconds_between(first(a.burst_start), a.burst_start[i]) + fvl[i] * dt) / dt);
             (s, s + (lvl[i] - fvl[i]) + 1)) for i in 1:n]
    nlines = n == 1 ? lpb :
             1 + round(Int, (SLCDatasets.seconds_between(first(a.burst_start), last(a.burst_start)) + (lpb - 1) * dt) / dt)
    sb = collect(SLCDatasets.bursts(sp, swath))
    rb = collect(SLCDatasets.bursts(rp, swath))
    length(sb) >= n || throw(ArgumentError("the secondary's IW" * string(swath) * " has fewer bursts than the reference's"))
    tasks = map(1:n) do i
        Threads.@spawn begin
            ss = sb[i]
            cs = radar_coordinate(ss)
            carrier = CarrierFit(TopsCarrier(ss, cs, slpb), slpb, sspb)
            dl, ds = offset_lattice(radar_coordinate(rb[i]), cs, dem, lpb, spb)
            prev = i > 1 ? fld(lims[i - 1][2] - lims[i][1], 2) : 0
            nxt = i < n ? fld(lims[i][2] - lims[i + 1][1], 2) : 0
            bstart, bend = n == 1 ? (0, lpb) : (fvl[i] + prev, 1 + lvl[i] - nxt)
            mstart, mend = n == 1 ? (0, lpb) : (lims[i][1] + prev, lims[i][2] - nxt)
            cols = n == 1 ? (1:spb) : ((fvs[i] + 1):min(lvs[i], spb))
            BurstPlacement(i, ((i - 1) * slpb + 1):(i * slpb), carrier, dl, ds,
                           (sa.first_valid_line[i] - 1, sa.last_valid_line[i] - 1,
                            sa.first_valid_sample[i] - 1, sa.last_valid_sample[i] - 1),
                           bstart:(bend - 1), (first(cols) - 1):(last(cols) - 1),
                           (mstart + 1):mend, cols, (slpb, sspb))
        end
    end
    places = BurstPlacement[fetch(t)::BurstPlacement for t in tasks]
    return ResampledSwath(sraster, places, (nlines, spb))
end

Base.size(r::ResampledSwath) = r.dims

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

# The burst's samples zeroed outside its valid rectangle, as `slc_to_vrt_file` exposes the burst:
# `origin` is the burst-local zero-based `(line, sample)` of `A[1, 1]`, the bounds zero-based inclusive.
function zero_outside_valid!(A::AbstractMatrix{ComplexF32}, valid::NTuple{4,Int}, origin::NTuple{2,Int})
    fvl, lvl, fvs, lvs = valid
    nr, nc = size(A)
    oy, ox = origin
    z = zero(ComplexF32)
    r0, r1 = clamp(fvl - oy + 1, 1, nr + 1), clamp(lvl - oy + 1, 0, nr)
    c0, c1 = clamp(fvs - ox + 1, 1, nc + 1), clamp(lvs - ox + 1, 0, nc)
    r0 > 1 && fill!(view(A, 1:(r0 - 1), :), z)
    r1 < nr && fill!(view(A, (r1 + 1):nr, :), z)
    if r0 <= r1
        c0 > 1 && fill!(view(A, r0:r1, 1:(c0 - 1)), z)
        c1 < nc && fill!(view(A, r0:r1, (c1 + 1):nc), z)
    end
    return A
end

# One burst's contribution to a window over the burst-local output indices `lines` x `samples`, into
# `dest`. The band read is the support those outputs reach, clamped to the burst: a pixel whose support
# leaves the burst is rejected by the resample and never read.
function resample_piece!(dest::AbstractMatrix{Float32}, r::ResampledSwath, p::BurstPlacement, lines::UnitRange{Int},
                         samples::UnitRange{Int}, scratch::Matrix{ComplexF32})
    slpb, sspb = p.extent
    ymin, ymax, xmin, xmax = support_bounds(p.dl, p.ds, lines, samples)
    brows = clamp(floor(Int, ymin) - 3, 1, slpb):clamp(ceil(Int, ymax) + 4, 1, slpb)
    bcols = clamp(floor(Int, xmin) - 3, 1, sspb):clamp(ceil(Int, xmax) + 4, 1, sspb)
    (isempty(brows) || isempty(bcols)) && return dest
    origin = (first(brows) - 1, first(bcols) - 1)
    srows = (first(p.source_rows) - 1 + first(brows)):(first(p.source_rows) - 1 + last(brows))
    src = deramped_burst!(scratch, r.source, srows, bcols, p.carrier, origin)
    zero_outside_valid!(src, p.svalid, origin)
    return resample!(dest, src, p.dl, p.ds, lines, samples; origin, extent = p.extent, doppler = p.extent)
end

function Base.getindex(r::ResampledSwath, rows::AbstractUnitRange, cols::AbstractUnitRange)
    checkbounds(r, rows, cols)
    out = zeros(Float32, length(rows), length(cols))
    slpb, sspb = first(r.places).extent
    scratch = Matrix{ComplexF32}(undef, slpb, sspb)
    for p in r.places
        mr = intersect(rows, p.mosaic_rows)
        isempty(mr) && continue
        mc = intersect(cols, p.mosaic_cols)
        isempty(mc) && continue
        t0 = first(mr) - first(p.mosaic_rows)
        u0 = first(mc) - first(p.mosaic_cols)
        lines = (first(p.lines) + t0):(first(p.lines) + t0 + length(mr) - 1)
        samples = (first(p.samples) + u0):(first(p.samples) + u0 + length(mc) - 1)
        dest = view(out, (first(mr) - first(rows) + 1):(last(mr) - first(rows) + 1),
                    (first(mc) - first(cols) + 1):(last(mc) - first(cols) + 1))
        resample_piece!(dest, r, p, lines, samples, scratch)
    end
    return out
end

@_window_indexing ResampledSwath

function secondary_mosaic(rp::S1Product, sp::S1Product, dem::D, raster_of) where {D}
    sws = rp.swaths
    lazy = [ResampledSwath(rp, sp, sw, dem, raster_of(sw)) for sw in sws]
    places, dims = mosaic_layout(rp, sws, [size(r) for r in lazy])
    return S1Mosaic(lazy, places, dims)
end

# A product's subswath measurement raster: a file of a `.SAFE`, a zip's member inflated, or extracted bursts.
function _swath_raster(p::S1Product, sw::Int)
    return SLCDatasets.burst_raster(first(collect(SLCDatasets.bursts(p, sw))).backend).raster
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
    # are still being fetched.
    assembled = !isempty(rp.bursts)
    assembled == !isempty(sp.bursts) || throw(ArgumentError(
        "the two Sentinel-1 products must be given alike: both as products, or both as extractor URLs"))
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
    ref = reference_mosaic(s.early_product, rraster)
    co = s.pair.coordinate
    (co.nlines, co.nsamples) == size(ref) || throw(DimensionMismatch(
        "the geogrid was built on a different mosaic than the imagery"))
    early = fill_tiles!(TileCache(ref; tile = 1024, dir, io = cache_io))
    _accumulate_zeros!(declined, early, lx, ly, sentinel, NODATA_SAMPLE_OFFSET)
    sec = secondary_mosaic(s.early_product, s.late_product, dem, sraster)
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
    _, t1 = s1_mosaic(rp, rp.swaths)
    _, t2 = s1_mosaic(sp, rp.swaths)
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
