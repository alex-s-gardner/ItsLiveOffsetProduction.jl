# NISAR L1 (RSLC) pairs: the radar-geometry setup, the terrain heights the coregistration solves against,
# and the secondary resampled onto the reference's radar grid.

const RslcProduct = SLCDatasets.SLC{SLCDatasets.NisarBackend}

"""
    open_rslc(uri) -> RslcProduct

The RSLC at `uri` — a local path, an `s3://` URI, or an `https://` URL read with an Earthdata Login token
from `EARTHDATA_TOKEN` or `~/.netrc`.
"""
function open_rslc(uri::String)
    src = SLCDatasets.source_for(uri)
    src isa SLCDatasets.LocalFile && return open_rslc(src)
    src isa SLCDatasets.RemoteHTTP && return open_rslc(src)
    src isa SLCDatasets.RemoteS3 && return open_rslc(src)
    throw(ArgumentError("no source for " * uri))
end

open_rslc(src::SLCDatasets.AbstractSLCSource) =
    require_frequency_a(SLCDatasets._open_nisar(src; frequency = nothing)::RslcProduct)

_look(side) = side == SLCDatasets.LookLeft ? ImagePairGeometry.LookLeft : ImagePairGeometry.LookRight

"""
    radar_coordinate(s::SLCDatasets.SLC) -> RadarCoordinate

The acquisition as a radar coordinate, with its scene-centre incidence angle — what
`RadarCoordinate(::SLCDatasets.SLC)` builds, with the orbit interpolant fixed to `Orbit`.
"""
function radar_coordinate(s::SLCDatasets.SLC)
    g = s.geometry
    sv = SLCDatasets.orbit(s)
    g.epoch == sv.epoch || throw(ArgumentError(
        "an acquisition's azimuth times and its state vectors are measured from different epochs"))
    first(sv.time) <= g.sensing_start && g.sensing_stop <= last(sv.time) || throw(ArgumentError(
        "an acquisition's state vectors do not span it"))
    orb = ImagePairGeometry.Orbit(sv)
    kwargs = (; orbit = orb, starting_range = g.starting_range, dr = g.range_pixel_spacing,
              sensing_start = g.sensing_start, prf = g.prf, nsamples = g.nsamples, nlines = g.nlines,
              look_side = _look(g.look_side), wavelength = g.wavelength, orbit_epoch_offset = 0.0)
    return ImagePairGeometry.RadarCoordinate(; kwargs..., incidence_angle = ImagePairGeometry.incidence_angle(; kwargs...))
end

"""
    RadarSetup

[`Setup`](@ref) for a radar-geometry pair (NISAR L1, Sentinel-1): the scenes in acquisition order and their
products, the pair — the reference's radar coordinate, for NISAR the secondary's for the resample, and the
interval — the parameter
region, the output grid and the window on it, the parameter rasters over it, and the geogrid result.
"""
struct RadarSetup{P,G,R}
    job::Job
    early::SceneSource
    late::SceneSource
    early_product::R
    late_product::R
    # The NISAR polarization both products are read at; empty for Sentinel-1, whose product holds its own.
    polarization::String
    pair::P
    region::Region
    grid::Grid
    window::GridWindow
    inputs::Inputs
    geometry::G
    # Rasters being fetched in the background while the geometry is solved; waited for before the imagery.
    pending::Vector{Task}
end

"""
    radar_setup(c::Job) -> RadarSetup

`setup` for a NISAR L1 pair (`loadMetadataRslc`, `testGeogrid.py:240-260`): the geometry comes from the two
products rather than from a map raster. The grid transform maps the parameter grid to **geodetic degrees**,
because a radar footprint is solved for with `rdr2geo` rather than transformed, and the nodata sentinel is
read from the `vx` raster rather than the DEM (`geogridRadar.cpp:509-512`).
"""
function radar_setup(c::Job)
    early, late = acquisition_order(c)
    re, rl = open_rslc(early.band), open_rslc(late.band)
    ref, sec = c.reference.name == early.name ? (re, rl) : (rl, re)
    pol = nisar_polarization(ref.backend, sec.backend)
    dt = SLCDatasets.repeat_interval(re, rl)
    dt > 0 || throw(ArgumentError("the secondary acquisition does not start after the reference"))
    pair = CoregisteredPair(radar_coordinate(re); dt, secondary = radar_coordinate(rl))
    b = footprint_bounds(IdentityTransform(), pair.coordinate)
    region = parameter_region((b.X[1] + b.X[2]) / 2, (b.Y[1] + b.Y[2]) / 2)
    grid = parameter_grid(region)
    nd = with_dataset(band_nodata, region.paths[4])
    g = FGP.native_crs(region.epsg)
    g isa FGP.UTMCRS && return _radar_setup(c, early, late, re, rl, pair, region, grid, nd, g; polarization = pol)
    g isa FGP.PolarStereographicCRS && return _radar_setup(c, early, late, re, rl, pair, region, grid, nd, g; polarization = pol)
    throw(ArgumentError("a product grid is UTM or polar stereographic, not EPSG:" * string(region.epsg)))
end

function _radar_setup(c::Job, early::SceneSource, late::SceneSource, re::R, rl::R, pair,
                      region::Region, grid::Grid, nd::Union{Nothing,Float64}, g::FGP.NativeCRS;
                      pending::Vector{Task} = Task[], polarization::String = "") where {R}
    f = FGP.Transformation(g, FGP.GeographicCRS())
    tf = TransformPair(f, inv(f))
    window = grid_window(grid, footprint_bounds(tf, pair.coordinate))
    inputs = geometry_inputs(region, window)
    geom = pairgeometry_blocked(grid, pair, InMemoryInputs(inputs, window); transform = () -> tf, window,
                                nodata = nodata_from(nd === nothing ? 0.0 : nd))
    return RadarSetup(c, early, late, re, rl, polarization, pair, region, grid, window, inputs, geom, pending)
end

dy_sign(::RadarSetup) = -DY_SIGN

# The footprint's midpoint in geodetic degrees, which the region lookup and the product's centroid take.
function radar_centroid(s::RadarSetup)
    b = footprint_bounds(IdentityTransform(), s.pair.coordinate)
    return ((b.X[1] + b.X[2]) / 2, (b.Y[1] + b.Y[2]) / 2)
end

# ---------------------------------------------------------------------------
# Terrain heights
# ---------------------------------------------------------------------------

"""
    DemSampler

Terrain height as `(lon_degrees, lat_degrees) -> height`, bilinear on a DEM window held in memory:
`to_grid` maps a point into the DEM's CRS, and the sample at zero-based index `i` sits at
`origin + (i + shift) * spacing` — `shift` is `0.5` for cell centres, `0` for cell corners.
"""
struct DemSampler{F}
    z::Matrix{Float64}
    to_grid::F
    origin::NTuple{2,Float64}
    spacing::NTuple{2,Float64}
    shift::Float64
end

function (d::DemSampler)(lon_d::Float64, lat_d::Float64)
    x, y = d.to_grid(lon_d, lat_d)
    px = (x - d.origin[1]) / d.spacing[1] - d.shift
    py = (y - d.origin[2]) / d.spacing[2] - d.shift
    nx, ny = size(d.z)
    i = clamp(floor(Int, px), 0, nx - 2)
    j = clamp(floor(Int, py), 0, ny - 2)
    fx, fy = px - i, py - j
    z = d.z
    return (1 - fx) * (1 - fy) * z[i + 1, j + 1] + fx * (1 - fy) * z[i + 2, j + 1] +
           (1 - fx) * fy * z[i + 1, j + 2] + fx * fy * z[i + 2, j + 2]
end

# Longitude and latitude unchanged, for a DEM already in geodetic degrees.
struct LonLat end
(::LonLat)(lon::Float64, lat::Float64) = (lon, lat)

"""
    geographic_dem(path) -> DemSampler

A whole EPSG:4326 DEM of heights above the ellipsoid — the `dem.tif` the reference container downloads for
a pair — sampled with its values at cell corners, as the harness's `dem_sampler` reads it.
"""
function geographic_dem(path::String)
    return with_dataset(gdal_path(path)) do d
        projected_wkt_is_geographic(projection_wkt(d)) || throw(ArgumentError(
            path * " is not on an EPSG:4326 grid; --dem takes a geographic DEM"))
        gt = geotransform(d)
        nx, ny = rastersize(d)
        DemSampler(read_window(d, Float64, 0, 0, nx, ny), LonLat(), (gt[1], gt[4]), (gt[2], gt[6]), 0.0)
    end
end

projected_wkt_is_geographic(wkt::String) = startswith(wkt, "GEOGCS") || startswith(wkt, "GEOGCRS")

"""
    EllipsoidalDem(dem, geoid)

Heights above the WGS84 ellipsoid from an orthometric `dem` and the `geoid` undulation, both
[`DemSampler`](@ref)s, as `(lon_degrees, lat_degrees) -> height`.
"""
struct EllipsoidalDem{F}
    dem::DemSampler{F}
    geoid::DemSampler{LonLat}
end

(d::EllipsoidalDem)(lon_d::Float64, lat_d::Float64) = d.dem(lon_d, lat_d) + d.geoid(lon_d, lat_d)

# EGM2008 on a 2.5' grid, PROJ's copy of the NGA model.
const EGM2008_GRID = "/vsicurl/https://cdn.proj.org/us_nga_egm08_25.tif"

"""
    parameter_dem(s::RadarSetup; margin_px = 8) -> EllipsoidalDem

Terrain heights above the ellipsoid over the pair's footprint, which is what `rdr2geo` solves against: the
parameter region's DEM, whose heights are above the EGM2008 geoid, plus EGM2008's undulation.

The DEM window is the footprint's geodetic bounds — their edges sampled densely enough that the projected
box contains the footprint — projected into the region and grown by `margin_px` cells, read once; the
geoid window is the geodetic bounds grown by `margin_px` cells. Both are bilinear with values at cell
centres.
"""
function parameter_dem(s::RadarSetup; margin_px::Int = 8)
    g = FGP.native_crs(s.region.epsg)
    g isa FGP.UTMCRS && return parameter_dem(s, FGP.Transformation(FGP.GeographicCRS(), g); margin_px)
    g isa FGP.PolarStereographicCRS && return parameter_dem(s, FGP.Transformation(FGP.GeographicCRS(), g); margin_px)
    throw(ArgumentError("a product grid is UTM or polar stereographic, not EPSG:" * string(s.region.epsg)))
end

function parameter_dem(s::RadarSetup, to_grid::F; margin_px::Int) where {F}
    b = footprint_bounds(IdentityTransform(), s.pair.coordinate)
    n = 33
    lons = vcat(collect(range(b.X[1], b.X[2]; length = n)), collect(range(b.X[1], b.X[2]; length = n)),
                fill(b.X[1], n), fill(b.X[2], n))
    lats = vcat(fill(b.Y[1], n), fill(b.Y[2], n), collect(range(b.Y[1], b.Y[2]; length = n)),
                collect(range(b.Y[1], b.Y[2]; length = n)))
    xs, ys = to_grid(lons, lats)
    geoid = with_dataset(EGM2008_GRID) do d
        gt = geotransform(d)
        w, h = rastersize(d)
        c0 = clamp(floor(Int, (b.X[1] - gt[1]) / gt[2]) - margin_px, 0, w - 1)
        c1 = clamp(ceil(Int, (b.X[2] - gt[1]) / gt[2]) + margin_px, 0, w - 1)
        r0 = clamp(floor(Int, (b.Y[2] - gt[4]) / gt[6]) - margin_px, 0, h - 1)
        r1 = clamp(ceil(Int, (b.Y[1] - gt[4]) / gt[6]) + margin_px, 0, h - 1)
        z = read_window(d, Float64, c0, r0, c1 - c0 + 1, r1 - r0 + 1)
        DemSampler(z, LonLat(), (gt[1] + c0 * gt[2], gt[4] + r0 * gt[6]), (gt[2], gt[6]), 0.5)
    end
    dem = with_dataset(dem_path(s.region)) do d
        gt = geotransform(d)
        w, h = rastersize(d)
        c0 = clamp(floor(Int, (minimum(xs) - gt[1]) / gt[2]) - margin_px, 0, w - 1)
        c1 = clamp(ceil(Int, (maximum(xs) - gt[1]) / gt[2]) + margin_px, 0, w - 1)
        r0 = clamp(floor(Int, (maximum(ys) - gt[4]) / gt[6]) - margin_px, 0, h - 1)
        r1 = clamp(ceil(Int, (minimum(ys) - gt[4]) / gt[6]) + margin_px, 0, h - 1)
        z = read_window(d, Float64, c0, r0, c1 - c0 + 1, r1 - r0 + 1)
        DemSampler(z, to_grid, (gt[1] + c0 * gt[2], gt[4] + r0 * gt[6]), (gt[2], gt[6]), 0.5)
    end
    return EllipsoidalDem(dem, geoid)
end

# ---------------------------------------------------------------------------
# The offset field
# ---------------------------------------------------------------------------

"""
    coregistration_offset(cr, cs, line, sample, dem) -> (dline, dsample)

Where the secondary images the ground point the reference images at zero-based `(line, sample)`, as an
offset in the reference's own pixels: `rdr2geo` on the reference, solving for the terrain against `dem`
as isce3's `Rdr2Geo` does, then `geo2rdr` on the secondary. The secondary line is read from the orbit's
clock, and the offsets from the solved time and range rather than through the geogrid's rounded indices.
"""
function coregistration_offset(cr, cs, line::Int, sample::Int, dem::D) where {D}
    el = ImagePairGeometry.Ellipsoid()
    az = cr.sensing_start + line / cr.prf + cr.orbit_epoch_offset
    rg = cr.starting_range + sample * cr.dr
    # `rdr2geo` hands a height source radians; the DEM takes degrees.
    height = (lon, lat) -> dem(lon / ImagePairGeometry.DEG2RAD, lat / ImagePairGeometry.DEG2RAD)
    llh = ImagePairGeometry.rdr2geo(cr.orbit, el, az, rg; height, wavelength = cr.wavelength, side = cr.look_side)
    xyz = ImagePairGeometry.lonlat_to_xyz(el, llh)
    t0 = ImagePairGeometry.orbit_midtime(cs)
    pm, vm = ImagePairGeometry.interpolate(cs.orbit, t0)
    p = ImagePairGeometry.geo2rdr(cs.orbit, xyz, t0, t0, pm, vm)
    return ((p.orbittime - cs.orbit_epoch_offset - cs.sensing_start) * cs.prf - line,
            (p.range - cs.starting_range) / cs.dr - sample)
end

# The lattice the offset field is solved on, in lines and samples.
const LSTEP, SSTEP = 64, 512

"""
    LatticeOffsets

One offset component on the [`offset_lattice`](@ref) grid, bilinear between nodes and clamped to the edge
cell outside it. Inside one cell it is linear in the line index at a fixed sample, which is what lets the
resampler advance it a cell at a time.
"""
struct LatticeOffsets
    A::Matrix{Float64}
    nl::Int
    ns::Int
end

function (f::LatticeOffsets)(l::Int, s::Int)
    fi = l / LSTEP + 1
    fj = s / SSTEP + 1
    i = clamp(floor(Int, fi), 1, f.nl - 1)
    j = clamp(floor(Int, fj), 1, f.ns - 1)
    u, v = fi - i, fj - j
    A = f.A
    return (1 - u) * (1 - v) * A[i, j] + u * (1 - v) * A[i + 1, j] + (1 - u) * v * A[i, j + 1] + u * v * A[i + 1, j + 1]
end

"""
    offset_lattice(cr, cs, dem, nlines, nsamples) -> (dl, ds)

The offset field solved at a node every [`LSTEP`](@ref) lines and [`SSTEP`](@ref) samples over an
`nlines x nsamples` grid and one step past it, as two [`LatticeOffsets`](@ref). Exact to about 1e-5 px on
this field, against a per-pixel solve.
"""
function offset_lattice(cr, cs, dem::D, nlines::Int, nsamples::Int) where {D}
    ls = 0:LSTEP:(nlines + LSTEP)
    ss = 0:SSTEP:(nsamples + SSTEP)
    DL = Matrix{Float64}(undef, length(ls), length(ss))
    DS = similar(DL)
    Threads.@threads for jj in eachindex(ss)
        for ii in eachindex(ls)
            DL[ii, jj], DS[ii, jj] = coregistration_offset(cr, cs, ls[ii], ss[jj], dem)
        end
    end
    return (LatticeOffsets(DL, length(ls), length(ss)), LatticeOffsets(DS, length(ls), length(ss)))
end

_next_mark(v::Int, step::Int) = (fld(v, step) + 1) * step

# The extremes of the source positions `lines x samples` reach. Inside one lattice cell the position is
# bilinear, so its extremes over a window are at the corners of the cells the window meets: the window's
# edges and the lattice lines inside it.
function support_bounds(dl::LatticeOffsets, ds::LatticeOffsets, lines::UnitRange{Int}, samples::UnitRange{Int})
    ymin = xmin = Inf
    ymax = xmax = -Inf
    s = first(samples)
    while true
        l = first(lines)
        while true
            y = l + dl(l, s) + 1.0
            x = s + ds(l, s) + 1.0
            ymin = min(ymin, y); ymax = max(ymax, y)
            xmin = min(xmin, x); xmax = max(xmax, x)
            l == last(lines) && break
            l = min(_next_mark(l, LSTEP), last(lines))
        end
        s == last(samples) && break
        s = min(_next_mark(s, SSTEP), last(samples))
    end
    return (ymin, ymax, xmin, xmax)
end

# ---------------------------------------------------------------------------
# The interpolation kernel
# ---------------------------------------------------------------------------

# The eight taps of a Hann-windowed sinc at fractional offset `f`, for taps at −3 through +4:
# `Sinc2dInterpolator`'s kernel, ISCE's `sinc_coef` at `pedestal = 0` and `ns = 8`.
function sinc8(f::Float64)
    return ntuple(8) do k
        x = (k - 4) - f
        s = x == 0 ? 1.0 : sinpi(x) / (pi * x)
        s * (0.5 + 0.5 * cospi(x / 4))
    end
end

# `sinc8` tabulated at `SINC_SUBDIVISIONS + 1` fractional offsets, each column normalized to sum to one.
# 2048 subdivisions bound the position quantization at 0.00024 px.
const SINC_SUBDIVISIONS = 2048
const SINC_TAPS = let n = SINC_SUBDIVISIONS
    T = Matrix{Float64}(undef, 8, n + 1)
    for k in 0:n
        t = sinc8(k / n)
        s = sum(t)
        for m in 1:8
            T[m, k + 1] = t[m] / s
        end
    end
    T
end

# `floor` and the tap bin to an `Int` without the conversion's range check, which is 20% of the kernel:
# the positions are finite pixel indices by the time they reach it, because `support_bounds` evaluates the
# field at the corners of every cell a window meets and `floor`s those first.
_ifloor(x::Float64) = unsafe_trunc(Int, floor(x))
_tap_bin(f::Float64) = unsafe_trunc(Int, round(f * SINC_SUBDIVISIONS)) + 1
sinc8_taps(f::Float64) = (k = _tap_bin(f); ntuple(m -> SINC_TAPS[m, k], 8))

"""
    resample!(out, src, dl, ds, lines, samples; origin, extent) -> out

The amplitude of `src` read at each output pixel's own position in the secondary, through an eight-tap
sinc. `lines` and `samples` are the zero-based output indices; `src[1, 1]` is the secondary's zero-based
`origin`, and `extent` its full size. A pixel whose support leaves the secondary is left at zero, as
`ResampSlc` has no samples there; one whose support leaves `src` is a band sized too small, and throws.

`doppler` is `ResampSlc`'s Doppler-LUT containment test, evaluated at the output pixel against the source
burst's own `(lines, samples)` (`s1_reader.py:131-168`): an output line past `doppler[1]` or sample at or
past `doppler[2] - 1` is skipped. A NISAR RSLC has no such test.
"""
function resample!(out::AbstractMatrix{Float32}, src::AbstractMatrix{ComplexF32}, dl::LatticeOffsets,
                   ds::LatticeOffsets, lines::UnitRange{Int}, samples::UnitRange{Int}; origin::NTuple{2,Int},
                   extent::NTuple{2,Int}, doppler::NTuple{2,Int} = (typemax(Int), typemax(Int)))
    oy, ox = origin
    ey, ex = extent
    ny, nx = size(src)
    Threads.@threads for jj in eachindex(samples)
        s = samples[jj]
        cell = typemin(Int)
        l_cell = 0
        pl = ql = ps = qs = 0.0
        for (ii, l) in enumerate(lines)
            (l > doppler[1] || s > doppler[2] - 1) && continue
            # The lattice field is advanced per cell: two evaluations every `LSTEP` lines.
            k = fld(l, LSTEP)
            if k != cell
                l_cell = k * LSTEP
                pl = dl(l_cell, s)
                ps = ds(l_cell, s)
                ql = (dl(l_cell + LSTEP, s) - pl) / LSTEP
                qs = (ds(l_cell + LSTEP, s) - ps) / LSTEP
                cell = k
            end
            d = l - l_cell
            y = l + (pl + ql * d) + 1.0
            x = s + (ps + qs * d) + 1.0
            iy, ix = _ifloor(y), _ifloor(x)
            (iy - 4 < 1 || iy + 4 > ey || ix - 4 < 1 || ix + 4 > ex) && continue
            by, bx = iy - oy, ix - ox
            (by - 3 < 1 || by + 4 > ny || bx - 3 < 1 || bx + 4 > nx) && throw(ArgumentError(
                "the resampler's source band was sized too small for its window"))
            ty = sinc8_taps(y - iy)
            tx = sinc8_taps(x - ix)
            # The guard above admits exactly rows `by - 3:by + 4` and columns `bx - 3:bx + 4`, the accesses
            # below; the bounds check is 42% of this kernel, which the compiler cannot elide for an index
            # derived from `floor`.
            acc = ComplexF64(0)
            @inbounds for m in 1:8
                row = ComplexF64(0)
                @simd for k8 in 1:8
                    row += tx[k8] * ComplexF64(src[by - 4 + m, bx - 4 + k8])
                end
                acc += ty[m] * row
            end
            out[ii, jj] = Float32(abs(acc))
        end
    end
    return out
end

# ---------------------------------------------------------------------------
# The secondary on the reference grid
# ---------------------------------------------------------------------------

"""
    DopplerCentroid

An RSLC's Doppler centroid LUT (`metadata/processingInformation/parameters/frequencyA/dopplerCentroid`), as
`(time, range) -> Hz`, bilinear between its nodes and clamped at its edges. `time` is seconds on the swath's
`zeroDopplerTime` clock and `range` slant range in metres; `lut[i, j]` is at `range[i]`, `time[j]`.
"""
struct DopplerCentroid
    time::Vector{Float64}
    range::Vector{Float64}
    lut::Matrix{Float64}
end

function DopplerCentroid(s::RslcProduct)
    b = s.backend
    p = "science/" * b.band * "/" * b.product_type * "/metadata/processingInformation/parameters/frequency" *
        b.frequency * "/"
    f = ChunkManifests.h5open_typed(b.path; transport = b.transport)
    try
        t = ChunkManifests.h5_read(f, p * "zeroDopplerTime", Float64)
        r = ChunkManifests.h5_read(f, p * "slantRange", Float64)
        # Row-major `(time, range)`: the range index runs fastest.
        lut = ChunkManifests.h5_read(f, p * "dopplerCentroid", Float64)
        length(lut) == length(t) * length(r) || throw(ArgumentError(
            b.path * ": the Doppler LUT's size is not its time axis by its range axis"))
        return DopplerCentroid(t, r, reshape(lut, length(r), length(t)))
    finally
        close(f)
    end
end

function _node(ax::Vector{Float64}, v::Float64)
    k = clamp(searchsortedlast(ax, v), 1, length(ax) - 1)
    return k, clamp((v - ax[k]) / (ax[k + 1] - ax[k]), 0.0, 1.0)
end

function (d::DopplerCentroid)(time::Float64, range::Float64)
    i, u = _node(d.range, range)
    j, w = _node(d.time, time)
    A = d.lut
    return (1 - u) * (1 - w) * A[i, j] + u * (1 - w) * A[i + 1, j] + (1 - u) * w * A[i, j + 1] + u * w * A[i + 1, j + 1]
end

"""
    ResampledRSLC(reference, secondary, dem) <: AbstractMatrix{Float32}

The secondary RSLC's amplitude on the reference's radar grid, resampled a window at a time through the
[`offset_lattice`](@ref) field and an eight-tap sinc.

**The azimuth spectrum is brought to baseband first.** An RSLC is focused to zero Doppler but keeps its
Doppler centroid, and the eight-tap kernel is accurate only near zero frequency; so each source column is
multiplied by `exp(-2πi f_dc t)` at its own slant range before interpolating, which is what `ResampSlc` does
with the Doppler LUT it is handed. The remodulation after is a unit phasor and cannot change an amplitude.
"""
struct ResampledRSLC{S} <: AbstractMatrix{Float32}
    source::S
    dl::LatticeOffsets
    ds::LatticeOffsets
    doppler::DopplerCentroid
    # The secondary's first-line time, line spacing, near range and range spacing, on the clocks
    # `doppler` is indexed by.
    t0::Float64
    dt::Float64
    r0::Float64
    dr::Float64
    dims::Tuple{Int,Int}
    extent::Tuple{Int,Int}
end

function ResampledRSLC(s::RadarSetup, dem::D; dir::String) where {D}
    # Each window read is the support of a shifted output tile, so neighbouring windows share chunks; a
    # remote granule's chunks are cached under `dir` as they arrive, so each is fetched once.
    px = nisar_pixels(s.late_product, s.polarization)
    src = permutedims(ChunkManifests.cached(parent(px), dir), (2, 1))
    cr, cs = s.pair.coordinate, s.pair.secondary_offset_coordinate
    dims = (s.early_product.geometry.nlines, s.early_product.geometry.nsamples)
    dl, ds = offset_lattice(cr, cs, dem, dims[1], dims[2])
    return ResampledRSLC(src, dl, ds, DopplerCentroid(s.late_product), cs.sensing_start, 1 / cs.prf,
                         cs.starting_range, cs.dr, dims, size(src))
end

Base.size(r::ResampledRSLC) = r.dims

function Base.getindex(r::ResampledRSLC, rows::AbstractUnitRange, cols::AbstractUnitRange)
    checkbounds(r, rows, cols)
    # One continuous swath: the grid's zero-based line and sample are the output index less one.
    lines = (first(rows) - 1):(last(rows) - 1)
    samples = (first(cols) - 1):(last(cols) - 1)
    ymin, ymax, xmin, xmax = support_bounds(r.dl, r.ds, lines, samples)
    ny, nx = r.extent
    brows = clamp(floor(Int, ymin) - 3, 1, ny):clamp(ceil(Int, ymax) + 4, 1, ny)
    bcols = clamp(floor(Int, xmin) - 3, 1, nx):clamp(ceil(Int, xmax) + 4, 1, nx)
    out = zeros(Float32, length(rows), length(cols))
    # A window whose support lies wholly outside the secondary reads nothing; every pixel of it is rejected.
    (isempty(brows) || isempty(bcols)) && return out
    src = Matrix{ComplexF32}(r.source[brows, bcols])
    _demodulate!(src, r, brows, bcols)
    return resample!(out, src, r.dl, r.ds, lines, samples; origin = (first(brows) - 1, first(bcols) - 1),
                     extent = r.extent)
end

@_window_indexing ResampledRSLC

# `src` holds source lines `brows` and samples `bcols`; each column is brought to baseband at its own slant
# range. The phase is taken from the absolute line index, so two windows agree where they overlap, and the
# centroid at the band's middle line, since it moves by under 2% of the line rate over a scene.
function _demodulate!(src::Matrix{ComplexF32}, r::ResampledRSLC, brows::UnitRange{Int}, bcols::UnitRange{Int})
    tmid = r.t0 + (first(brows) + last(brows) - 2) / 2 * r.dt
    Threads.@threads for j in eachindex(bcols)
        cycles = r.doppler(tmid, r.r0 + (bcols[j] - 1) * r.dr) * r.dt
        for (i, n) in pairs(brows)
            src[i, j] *= cis(-2pi * rem((n - 1) * cycles, 1.0))
        end
    end
    return src
end

# ---------------------------------------------------------------------------
# Imagery and packaging
# ---------------------------------------------------------------------------

"""
    radar_reference(s; dir, cache_io) -> TileFile

The reference's own amplitude for a NISAR L1 pair, which defines the grid, filled into a tile cache under
`dir`. It needs no terrain, so it is read while the DEM is.
"""
function radar_reference(s::RadarSetup; dir::String, cache_io::Symbol)
    ref = SLCDatasets.amplitude(nisar_pixels(s.early_product, s.polarization))
    return fill_tiles!(TileCache(LazyWindow(ref, 0, 0, size(ref)); tile = 1024, dir, io = cache_io);
                       inflight = remote_tiles_in_flight())
end

"""
    radar_imagery(s, dem, early; dir, cache_io) -> Imagery

The correlator's two images for a NISAR L1 pair: the reference's amplitude `early`
([`radar_reference`](@ref)), and the secondary resampled onto it with terrain heights from `dem`, filled
into a tile cache under `dir`. Both share the geometry's grid, so the location bands apply to them with no
shift.
"""
function radar_imagery(s::RadarSetup, dem::D, early::TileFile; dir::String, cache_io::Symbol) where {D}
    g = s.geometry
    lx, ly = permutedims(g.location_x), permutedims(g.location_y)
    sentinel = Int32(g.nodata.output)
    declined = _sentinel_declined(lx, ly, sentinel)
    _accumulate_zeros!(declined, early, lx, ly, sentinel, NODATA_SAMPLE_OFFSET)
    sec = ResampledRSLC(s, dem; dir)
    size(sec) == size(early) || throw(DimensionMismatch(
        "the reference grid and the resampled secondary are not the same grid"))
    # A secondary tile is resampled as well as read, so half as many are kept in flight.
    late = fill_tiles!(TileCache(sec; tile = 1024, dir, io = cache_io); inflight = remote_tiles_in_flight() ÷ 2)
    _accumulate_zeros!(declined, late, lx, ly, sentinel, NODATA_SAMPLE_OFFSET)
    return Imagery(early, late, declined)
end

# The product's image-pair metadata and, for a Sentinel-1 pair of two spacecraft, its subswath bias inputs.
_radar_info(s::RadarSetup{P,G,RslcProduct}, c::Correlation, roi::Float64, lon::Float64, lat::Float64) where {P,G} =
    (nisar_img_pair_info(s.early.name, s.late.name, s.early_product.identification,
                         s.late_product.identification, roi, lon, lat), nothing)

"""
    package(s::RadarSetup, c, output)

Write a NISAR L1 or Sentinel-1 product: the `:radar` schema, in `"radar, map"` coordinates, with the pair's interval.
`pixel_size_x` is the ground-range spacing at the scene centre, `dr / sin(incidence_angle)`, which matches
the reference's `rangePixelSize`; nothing reads `pixel_size_y` for this schema, so it is the same value.
"""
function package(s::RadarSetup, c::Correlation, output::String)
    g = FGP.native_crs(s.region.epsg)
    g isa FGP.UTMCRS && return package(s, c, output, FGP.Transformation(g, FGP.GeographicCRS()))
    g isa FGP.PolarStereographicCRS && return package(s, c, output, FGP.Transformation(g, FGP.GeographicCRS()))
    throw(ArgumentError("a product grid is UTM or polar stereographic, not EPSG:" * string(s.region.epsg)))
end

function package(s::RadarSetup, c::Correlation, output::String, lonlat::F) where {F}
    gt = ImagePairGeometry.window_geotransform(s.grid, s.window)
    nx, ny = size(s.window)
    x = collect(gt[1] + gt[2] / 2 .+ (0:(nx - 1)) .* gt[2])
    y = collect(gt[4] + gt[6] / 2 .+ (0:(ny - 1)) .* gt[6])
    coord = s.pair.coordinate
    px = coord.dr / sin(coord.incidence_angle)
    georef = ItsLiveOffsetProduction.ItsLiveGeoref(x, y, ItsLiveOffsetProduction.cf_grid_mapping(s.region.epsg),
                                                   lonlat, px, px)
    roi = ItsLiveOffsetProduction.roi_valid_percentage(c.chip_size, c.search)
    lon, lat = radar_centroid(s)
    info, swath_bias = _radar_info(s, c, roi, lon, lat)
    coeffs = ItsLiveOffsetProduction.coefficients(s.geometry)
    refv, ssm = ItsLiveOffsetProduction.reference_velocity(s.inputs)
    input = ItsLiveOffsetProduction.ItsLiveInput(
        :radar, DETECTION_METHOD, "radar, map", AUTORIFT_VERSION,
        replace(PARAMETER_SHAPEFILE, "/vsicurl/" => ""), product_source(info), s.pair.dt,
        georef, info, coeffs, refv, ssm, swath_bias, c.dx, c.dy, c.chip_size, c.scale, c.interpolated)
    isempty(dirname(output)) || mkpath(dirname(output))
    ItsLiveOffsetProduction.write_product(output, input)
    return output
end
