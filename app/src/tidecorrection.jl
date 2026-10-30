# The vertical motion of floating ice between a radar pair's two acquisitions, removed from its offsets:
# the ocean tide from CATS2008_v2023 with the model's own grounding-zone flexure, and the
# inverse-barometer response to ERA5 mean sea level pressure, after Wuite et al. (2026, Remote Sens.
# Environ. 332, 115092, eq. 10).
#
# In radar geometry (Sentinel-1, NISAR L1) the correlator's `dx` is where a feature sits in the late image
# minus where it sits in the early one, in range samples of the early acquisition's grid (`correlate`
# passes the images to AutoRIFT late first). A surface raised by Δz between the two acquisitions is
# Δz·cos θ nearer the radar, θ the incidence angle against the ellipsoid normal, so `dx` holds
# −Δz·cos θ / dr on top of the ice's horizontal motion; [`range_offset_correction`](@ref) is the
# Δz·cos θ / dr added back. In map geometry (NISAR L2) the same rise appears as horizontal motion toward
# the sensor; see the `Setup` method of [`tide_correct`](@ref).
#
# The correction is applied only on the parameter region's floating-ice mask, so grounded and stable
# surfaces, and with them the stable-surface statistics, are unchanged. Across the grounding zone it is
# scaled by the model's flexure and by [`grounding_taper`](@ref), which brings it to zero at the ITS_LIVE
# grounding line wherever the model's coarser grounding line lies inland of it.

"""
    TIDE_MODEL

The CATS2008_v2023 height constants and flexure as one cloud-optimized GeoTIFF, built from the
model's netCDF by `tools/tidal/make_cats_cog.jl`: bands `hRe_<c>`, `hIm_<c>` for each constituent of
[`TIDE_CONSTITUENTS`](@ref) and `flexure`, in metres on the model's 2 km polar stereographic grid.
"""
const TIDE_MODEL = "/vsicurl/https://its-live-data.s3.amazonaws.com/tidal/CATS2008_v2023_z.tif"

const TIDE_BANDS = (ntuple(i -> (isodd(i) ? "hRe_" : "hIm_") * TIDE_CONSTITUENTS[(i + 1) ÷ 2], 20)..., "flexure")
const FLEXURE_BAND = length(TIDE_BANDS)

# The parameter region the tide model covers: Antarctic polar stereographic.
const ANTARCTIC_EPSG = 3031

# The model's projection, for the metres-based copy of its grid.
const CATS_PROJECTION = FGP.LonLatToPolarStereographic(; lat_ts = -71.0, lon_0 = -70.0)

# ERA5 hourly mean sea level pressure (Pa) on its 0.25° grid, one Zarr chunk per hour, from Google's
# public Analysis-Ready, Cloud-Optimized copy; the hour is appended as the time index.
const ERA5_MSL = "ZARR:\"/vsigs/gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3\":/mean_sea_level_pressure:"
const ERA5_EPOCH = DateTime(1900, 1, 1)
const ERA5_NLON = 1440
const ERA5_STEP = 0.25

# The inverse-barometer response, metres of sea surface per pascal: 1 / (ρ g) with the sea water density
# and gravity of Wuite et al. (2026).
const IB_M_PER_PA = 1 / (1030.0 * 9.81)

"""
    TideConstants

The model's constants over a window of its grid: `bands[x, y, b]` in [`TIDE_BANDS`](@ref) order, with the
window's first cell centre at `origin` and cells `spacing` apart, both in metres.
"""
struct TideConstants
    bands::Array{Float32,3}
    origin::NTuple{2,Float64}
    spacing::NTuple{2,Float64}
end

"""
    tide_constants(path, xs, ys) -> TideConstants

The window of the model at `path` covering the model-grid points `xs`, `ys` (metres) with one cell to
spare, all bands. Throws when the file's bands are not [`TIDE_BANDS`](@ref) or a point lies outside it.
"""
function tide_constants(path::String, xs::Vector{Float64}, ys::Vector{Float64})
    return with_dataset(path) do d
        band_count(d) == length(TIDE_BANDS) || throw(ArgumentError(
            path * " has " * string(band_count(d)) * " bands; a tide model has " * string(length(TIDE_BANDS))))
        for b in eachindex(TIDE_BANDS)
            band_description(d, b) == TIDE_BANDS[b] || throw(ArgumentError(
                path * " band " * string(b) * " is " * band_description(d, b) * ", not " * TIDE_BANDS[b]))
        end
        gt = geotransform(d)
        w, h = rastersize(d)
        c0 = floor(Int, (minimum(xs) - gt[1]) / gt[2] - 0.5) - 1
        c1 = ceil(Int, (maximum(xs) - gt[1]) / gt[2] - 0.5) + 1
        r0 = floor(Int, (maximum(ys) - gt[4]) / gt[6] - 0.5) - 1
        r1 = ceil(Int, (minimum(ys) - gt[4]) / gt[6] - 0.5) + 1
        (c0 >= 0 && r0 >= 0 && c1 < w && r1 < h) || throw(ArgumentError(
            "the floating ice lies outside the tide model " * path))
        bands = Array{Float32,3}(undef, c1 - c0 + 1, r1 - r0 + 1, length(TIDE_BANDS))
        for b in eachindex(TIDE_BANDS)
            bands[:, :, b] = read_window(d, Float32, c0, r0, c1 - c0 + 1, r1 - r0 + 1; band = b)
        end
        TideConstants(bands, (gt[1] + (c0 + 0.5) * gt[2], gt[4] + (r0 + 0.5) * gt[6]), (gt[2], gt[6]))
    end
end

"""
    tide_point(k, x, y) -> (z, flexure)

The ten complex constants at model-grid point `(x, y)`, each scaled by the flexure, and the flexure:
bilinear between cell centres, as pyTMD interpolates a TMD3 model whose constants it has multiplied by
the flexure first. The imaginary part is the file's `hIm` negated, pyTMD's convention for TMD3.
"""
function tide_point(k::TideConstants, x::Float64, y::Float64)
    px = (x - k.origin[1]) / k.spacing[1]
    py = (y - k.origin[2]) / k.spacing[2]
    i, j = floor(Int, px), floor(Int, py)
    fx, fy = px - i, py - j
    nx, ny = size(k.bands, 1), size(k.bands, 2)
    (0 <= i < nx - 1 && 0 <= j < ny - 1) || throw(ArgumentError("a point lies outside the tide model window"))
    w = ((1 - fx) * (1 - fy), fx * (1 - fy), (1 - fx) * fy, fx * fy)
    cells = ((i + 1, j + 1), (i + 2, j + 1), (i + 1, j + 2), (i + 2, j + 2))
    B = k.bands
    fl = ntuple(q -> Float64(B[cells[q][1], cells[q][2], FLEXURE_BAND]), Val(4))
    z = ntuple(Val(10)) do n
        re, im = 0.0, 0.0
        for q in 1:4
            c = cells[q]
            re += w[q] * fl[q] * B[c[1], c[2], 2n - 1]
            im -= w[q] * fl[q] * B[c[1], c[2], 2n]
        end
        complex(re, im)
    end
    flexure = sum(ntuple(q -> w[q] * fl[q], Val(4)))
    return z, flexure
end

# GDAL's Zarr driver deadlocks in `AdviseRead` (GDAL 3.13.3) when two threads read the store at once.
const ERA5_LOCK = ReentrantLock()

"""
    era5_pressure_rows(hour, row0, row1) -> Union{Matrix{Float32},Nothing}

ERA5 mean sea level pressure (Pa) at `hour` (hours since 1900-01-01T00 UTC), rows `row0`–`row1` (zero-based
from 90°N) of every longitude, indexed `[lon, lat]`; `nothing` when ERA5 does not yet reach that hour,
which the archive shows either by having no such index or by an hour of no values.
"""
function era5_pressure_rows(hour::Int, row0::Int, row1::Int)
    p = lock(ERA5_LOCK) do
        try
            with_dataset(d -> read_window(d, Float32, 0, row0, ERA5_NLON, row1 - row0 + 1), ERA5_MSL * string(hour))
        catch e
            # An hour past the end of the record is an index the store does not have.
            e isa ErrorException && (m = e.msg; m isa String && occursin("out of bounds", m)) && return nothing
            rethrow()
        end
    end
    (p === nothing || any(isnan, p)) && return nothing
    return p
end

# The two hours bracketing `t`, and the weight of the later.
function _era5_hours(t::DateTime)
    ms = Dates.value(t - ERA5_EPOCH)
    h0 = fld(ms, 3_600_000)
    return h0, (ms - h0 * 3_600_000) / 3_600_000
end

# Pressure at `(lon, lat)` from rows starting `row0` rows south of 90°N: bilinear between the grid's
# points, which lie on whole multiples of 0.25° of longitude from 0° and of latitude from 90°N.
function _era5_at(p::Matrix{Float32}, row0::Int, lon::Float64, lat::Float64)
    px = mod(lon, 360.0) / ERA5_STEP
    py = (90.0 - lat) / ERA5_STEP - row0
    i, j = floor(Int, px), floor(Int, py)
    fx, fy = px - i, py - j
    i0, i1 = mod(i, ERA5_NLON) + 1, mod(i + 1, ERA5_NLON) + 1
    (0 <= j < size(p, 2) - 1) || throw(ArgumentError("a point lies outside the ERA5 rows read"))
    return (1 - fx) * (1 - fy) * p[i0, j + 1] + fx * (1 - fy) * p[i1, j + 1] +
           (1 - fx) * fy * p[i0, j + 2] + fx * fy * p[i1, j + 2]
end

"""
    era5_pressure(t, lons, lats) -> Union{Vector{Float64},Nothing}

Mean sea level pressure (Pa) at each `(lons[k], lats[k])` at time `t`, linear in time between ERA5's hours;
`nothing` when ERA5 does not yet reach `t`.
"""
function era5_pressure(t::DateTime, lons::Vector{Float64}, lats::Vector{Float64})
    h0, f = _era5_hours(t)
    row0 = clamp(floor(Int, (90.0 - maximum(lats)) / ERA5_STEP) - 1, 0, 720)
    row1 = clamp(ceil(Int, (90.0 - minimum(lats)) / ERA5_STEP) + 1, 0, 720)
    p0 = era5_pressure_rows(h0, row0, row1)
    p0 === nothing && return nothing
    p1 = era5_pressure_rows(h0 + 1, row0, row1)
    p1 === nothing && return nothing
    return [(1 - f) * _era5_at(p0, row0, lons[k], lats[k]) + f * _era5_at(p1, row0, lons[k], lats[k])
            for k in eachindex(lons)]
end

"""
    range_offset_correction(dz, cosθ, dr) -> Float64

The range offset, in samples of spacing `dr` metres, that adds back a vertical displacement `dz` (late
minus early, metres, positive up) seen at incidence angle `θ`: raising the surface by `dz` shortens the
slant range by `dz·cos θ`, which the correlator reads as a `dx` of `−dz·cos θ / dr`.
"""
range_offset_correction(dz::Float64, cosθ::Float64, dr::Float64) = dz * cosθ / dr

const SV3 = ImagePairGeometry.SVector{3,Float64}

"""
    geocoded_shift(dz, los, vel, east, north) -> SVector{3,Float64}

Where an image geocoded onto a fixed surface places a point that rose by `dz` metres, relative to where
it is, to first order in `dz`: the point `δ` on the surface's tangent plane that images at the raised
point's range and zero-Doppler time. Raising the point by `dz·up` changes its range by `dz·(los·up)` and
its Doppler `(P − S)·v` by `dz·(up·v)`, so `δ` solves `δ·los = dz·(up·los)` and `δ·v = dz·(up·v)` within
the plane spanned by `east` and `north`.

`los` is the unit vector from the point toward the sensor, `vel` the sensor's velocity at the point's
zero-Doppler time, and `east`, `north` the unit vectors of the local horizontal, all in one frame (ECEF);
`up` completes them. The Doppler condition moves `δ` along track by `dz·(up·v̂)/sin²θ`, a few millimetres
per metre of `dz` against the `dz·cot θ` across it.
"""
function geocoded_shift(dz::Float64, los::SV3, vel::SV3, east::SV3, north::SV3)
    up = ImagePairGeometry.cross3(east, north)
    dot = ImagePairGeometry.dot3
    a11, a12 = dot(east, vel), dot(north, vel)
    a21, a22 = dot(east, los), dot(north, los)
    b1, b2 = dz * dot(up, vel), dz * dot(up, los)
    det = a11 * a22 - a12 * a21
    abs(det) > 1.0e-6 * sqrt(sum(abs2, vel)) || throw(ArgumentError(
        "the line of sight lies along the velocity or the vertical; the shift is undefined"))
    a = (b1 * a22 - a12 * b2) / det
    b = (a11 * b2 - a21 * b1) / det
    return a .* east .+ b .* north
end

# The geodetic east and north unit vectors at longitude `lon` and latitude `lat`, radians, in ECEF.
_east(lon::Float64) = SV3(-sin(lon), cos(lon), 0.0)
_north(lon::Float64, lat::Float64) = SV3(-sin(lat) * cos(lon), -sin(lat) * sin(lon), cos(lat))

# The distance over which the correction rises from zero at the ITS_LIVE grounding line: one cell of the
# tide model's grid, the scale below which the model cannot place its own grounding line.
const GROUNDING_TAPER_LENGTH = 2000.0

"""
    grounding_taper(d) -> Float64

The weight of the correction at distance `d` metres from the nearest grounded pixel: `3t² − 2t³` with
`t = d / GROUNDING_TAPER_LENGTH`, clamped to [0, 1] — zero at the grounding line, one beyond the taper
length, with zero slope at both ends.
"""
function grounding_taper(d::Float64)
    t = clamp(d / GROUNDING_TAPER_LENGTH, 0.0, 1.0)
    return t * t * (3 - 2t)
end

# A squared distance that stands for "no feature": finite, so the lower envelope's intersections stay
# finite, and far larger than any squared distance on a grid of 10⁵ pixels a side.
const _NO_FEATURE = 1.0e20

# The squared distance from each element of `f`'s index to the nearest feature along one line, where
# `f[q]` is the squared distance already found to `q` (0 at a feature): the lower envelope of parabolas
# of Felzenszwalb and Huttenlocher (2012, Theory Comput. 8, 415), exact.
function _sqdist_1d!(out::AbstractVector{Float64}, f::AbstractVector{Float64}, v::Vector{Int},
                     z::Vector{Float64})
    n = length(f)
    k = 1
    v[1] = 1
    z[1] = -Inf
    z[2] = Inf
    for q in 2:n
        s = ((f[q] + q^2) - (f[v[k]] + v[k]^2)) / (2q - 2v[k])
        while s <= z[k]
            k -= 1
            s = ((f[q] + q^2) - (f[v[k]] + v[k]^2)) / (2q - 2v[k])
        end
        k += 1
        v[k] = q
        z[k] = s
        z[k + 1] = Inf
    end
    k = 1
    for q in 1:n
        while z[k + 1] < q
            k += 1
        end
        out[q] = (q - v[k])^2 + f[v[k]]
    end
    return out
end

"""
    distance_to(feature::AbstractMatrix{Bool}, spacing) -> Matrix{Float64}

The Euclidean distance from every pixel to the nearest `true` pixel, for square pixels `spacing` apart
(so 0 on a feature and `spacing` beside one); about `1e10·spacing` where `feature` has none.
"""
function distance_to(feature::AbstractMatrix{Bool}, spacing::Float64)
    nx, ny = size(feature)
    d = [feature[i, j] ? 0.0 : _NO_FEATURE for i in 1:nx, j in 1:ny]
    n = max(nx, ny)
    v, z, f, out = Vector{Int}(undef, n), Vector{Float64}(undef, n + 1), Vector{Float64}(undef, n),
                   Vector{Float64}(undef, n)
    for j in 1:ny
        copyto!(f, 1, view(d, :, j), 1, nx)
        _sqdist_1d!(view(out, 1:nx), view(f, 1:nx), v, z)
        d[:, j] .= view(out, 1:nx)
    end
    for i in 1:nx
        copyto!(f, 1, view(d, i, :), 1, ny)
        _sqdist_1d!(view(out, 1:ny), view(f, 1:ny), v, z)
        d[i, :] .= view(out, 1:ny)
    end
    return sqrt.(d) .* spacing
end

"""
    grounding_weights(floating, ocean, offset, window_size, spacing) -> Matrix{Float64}

The taper weight of each pixel of a product window, indexed `[x, y]`: [`grounding_taper`](@ref) of the
distance to the nearest grounded pixel on floating ice, zero elsewhere. `floating` and `ocean` are the
two masks over the window grown by a margin, `offset` how far that margin reaches before the window's
first column and row. Every pixel that is neither floating ice nor ocean is grounded — grounded ice and
rock — so an ice front facing the open ocean is not tapered.
"""
function grounding_weights(floating::AbstractMatrix{Bool}, ocean::AbstractMatrix{Bool}, offset::NTuple{2,Int},
                           window_size::NTuple{2,Int}, spacing::Float64)
    size(floating) == size(ocean) || throw(DimensionMismatch("the floating-ice and ocean masks differ in size"))
    d = distance_to(.!(floating .| ocean), spacing)
    nx, ny = window_size
    return [floating[i + offset[1], j + offset[2]] ? grounding_taper(d[i + offset[1], j + offset[2]]) : 0.0
            for i in 1:nx, j in 1:ny]
end

"""
    cos_incidence(coord) -> (samples, cosθ)

The cosine of the incidence angle on the ellipsoid at 65 range samples spanning `coord`'s swath, at
mid-azimuth: it varies by about a quarter across a Sentinel-1 IW swath, and negligibly along it.
"""
function cos_incidence(coord)
    samples = collect(range(0.0, coord.nsamples - 1.0; length = 65))
    az = coord.sensing_start + (coord.nlines / 2) / coord.prf + coord.orbit_epoch_offset
    c = [cos(ImagePairGeometry.incidence_angle(coord.orbit, ImagePairGeometry.Ellipsoid(), az,
                                               coord.starting_range + s * coord.dr;
                                               wavelength = coord.wavelength, side = coord.look_side,
                                               zrange = (0.0,)))
         for s in samples]
    return samples, c
end

function _interp_linear(xs::Vector{Float64}, ys::Vector{Float64}, x::Float64)
    k = clamp(searchsortedlast(xs, x), 1, length(xs) - 1)
    f = (x - xs[k]) / (xs[k + 1] - xs[k])
    return (1 - f) * ys[k] + f * ys[k + 1]
end

"""
    tide_correct(s::RadarSetup, x, y, lonlat, dx, info, model) -> Matrix{Float32}

`dx` with the vertical motion of floating ice between the two acquisitions removed, and the corrections
recorded in `info.extra`; `dx` itself when the region is not Antarctic or the pair covers no floating
ice. `x`, `y` are the product grid's cell centres and `lonlat` its transform to longitude and latitude;
`model` is the tide model ([`TIDE_MODEL`](@ref)).

The tide and pressure are each scene's at its mid-time; across a scene they change by millimetres.
"""
function tide_correct(s::RadarSetup, x::Vector{Float64}, y::Vector{Float64}, lonlat::F, dx::Matrix{Float32},
                      info::ItsLiveOffsetProduction.ImagePairInfo, model::String) where {F}
    idx, w = floating_points(s.region, s.window, s.grid, dx)
    isempty(idx) && return dx
    lx, _ = ItsLiveOffsetProduction.image_location(s.geometry)
    idx = filter(I -> isfinite(lx[I]), idx)
    isempty(idx) && return dx
    coord = s.pair.coordinate
    half = Dates.Millisecond(round(Int, 1000 * coord.nlines / coord.prf / 2))
    t1, t2 = info.acquisition_date_img1 + half, info.acquisition_date_img2 + half
    dz, _, _, ib = vertical_motion(idx, w, x, y, lonlat, t1, t2, model)
    samples, cosθ = cos_incidence(coord)
    out = copy(dx)
    for (q, I) in enumerate(idx)
        out[I] = dx[I] + Float32(range_offset_correction(dz[q], _interp_linear(samples, cosθ, lx[I]), coord.dr))
    end
    _record!(info, dz, ib, t2, "range offsets")
    return out
end

"""
    floating_points(region, window, grid, dx) -> (idx, w)

The product-grid indices `[row, col]` that are floating ice with a measured `dx`, and the taper weight of
every pixel ([`grounding_weights`](@ref)), for the Antarctic polar stereographic region; empty for any
other. The southern UTM regions end at 56°S, outside the tide model.
"""
function floating_points(region::Region, window::GridWindow, grid::Grid, dx::Matrix{Float32})
    idx = CartesianIndex{2}[]
    region.epsg == ANTARCTIC_EPSG || return idx, zeros(size(dx))
    gt = grid.geotransform
    gt[2] == -gt[6] || throw(ArgumentError("the product grid's pixels are not square"))
    margin = ceil(Int, GROUNDING_TAPER_LENGTH / gt[2]) + 1
    floating, off = mask_window(region.floating, window, margin)
    ocean, off_ocean = mask_window(region.ocean, window, margin)
    off == off_ocean || throw(DimensionMismatch("the floating-ice and ocean masks are not on one grid"))
    w = permutedims(grounding_weights(floating, ocean, off, size(window), gt[2]))
    size(w) == size(dx) || throw(DimensionMismatch("the floating-ice window is not the product grid"))
    for I in CartesianIndices(dx)
        w[I] > 0 && isfinite(dx[I]) && push!(idx, I)
    end
    return idx, w
end

"""
    vertical_motion(idx, w, x, y, lonlat, t1, t2, model) -> (dz, lons, lats, ib)

The tapered vertical motion `dz` (metres, late minus early) of the floating-ice points `idx` between `t1`
and `t2`: ocean tide plus inverse barometer, scaled by the model's flexure and the taper `w`. Also each
point's longitude and latitude (degrees) and whether the inverse barometer was included.
"""
function vertical_motion(idx::Vector{CartesianIndex{2}}, w::Matrix{Float64}, x::Vector{Float64},
                         y::Vector{Float64}, lonlat::F, t1::DateTime, t2::DateTime, model::String) where {F}
    n = length(idx)
    lons, lats = Vector{Float64}(undef, n), Vector{Float64}(undef, n)
    mx, my = Vector{Float64}(undef, n), Vector{Float64}(undef, n)
    for (k, I) in enumerate(idx)
        lon, lat = lonlat(x[I[2]], y[I[1]])
        lons[k], lats[k] = lon, lat
        mx[k], my[k] = CATS_PROJECTION(lon, lat)
    end
    p1 = era5_pressure(t1, lons, lats)
    p2 = p1 === nothing ? nothing : era5_pressure(t2, lons, lats)
    k = tide_constants(model, mx, my)
    pred1, pred2 = TidePrediction(t1), TidePrediction(t2)
    dz = Vector{Float64}(undef, n)
    for (q, I) in enumerate(idx)
        z, flexure = tide_point(k, mx[q], my[q])
        d = tide_height(z, pred2) - tide_height(z, pred1)
        if p1 !== nothing && p2 !== nothing
            d -= flexure * (p2[q] - p1[q]) * IB_M_PER_PA
        end
        dz[q] = d * w[I]
    end
    return dz, lons, lats, p1 !== nothing && p2 !== nothing
end

function _record!(info::ItsLiveOffsetProduction.ImagePairInfo, dz::Vector{Float64}, ib::Bool, t::DateTime,
                  removed_from::String)
    info.extra["tide_correction"] = "CATS2008_v2023 ocean tide (doi:10.15784/601772) with its flexure, " *
                                    "tapered to zero over 2 km at the grounding line, on floating ice, " *
                                    "removed from the " * removed_from
    info.extra["inverse_barometer_correction"] = ib ?
        "ERA5 mean sea level pressure (ARCO-ERA5), -9.9 mm/hPa, with the tide" :
        "none: ERA5 does not yet reach " * _hour_string(t)
    m = round(_median(dz); digits = 3)
    info.extra["tide_correction_dz_median"] = m
    say("tide       " * string(length(dz)) * " floating points, median dz " * string(m) *
        " m" * (ib ? "" : ", no inverse barometer: ERA5 does not reach the pair"))
    return nothing
end

"""
    tide_correct(s::Setup, x, y, lonlat, dx, dy, info, model) -> (dx, dy)

A NISAR L2 (GSLC) pair's map offsets, `dx` east and `dy` north in pixels of the images (the GSLC grid's
spacing, not the product grid's), with the vertical motion of floating ice between the two acquisitions
removed, and the corrections recorded in `info.extra`; the offsets themselves for an optical pair, a
region that is not Antarctic, or a pair with no floating ice.

A GSLC is geocoded onto the DEM, which does not move with the tide, so a surface that rose by Δz between
the two acquisitions appears in the late image Δz·cot θ along the ground toward the sensor
([`geocoded_shift`](@ref)); the correction subtracts that apparent motion, projected onto the product grid.
The line of sight and the sensor's velocity at each point come from the early GSLC's orbit.
"""
function tide_correct(s::Setup, x::Vector{Float64}, y::Vector{Float64}, lonlat::F, dx::Matrix{Float32},
                      dy::Matrix{Float32}, info::ItsLiveOffsetProduction.ImagePairInfo, model::String) where {F}
    eg = s.early_gslc
    eg === nothing && return dx, dy
    idx, w = floating_points(s.region, s.window, s.grid, dx)
    isempty(idx) && return dx, dy
    idx = filter(I -> isfinite(dy[I]), idx)
    isempty(idx) && return dx, dy
    # The correlator measured in image pixels; the images are the GSLCs on their own map grid.
    spacing = s.pair.coordinate.spacing
    t1, t2 = _mid_time(eg.product.identification), _mid_time(s.late_gslc.product.identification)
    dz, lons, lats, ib = vertical_motion(idx, w, x, y, lonlat, t1, t2, model)

    sv = SLCDatasets.read_orbit(eg.product.backend)
    orb = ImagePairGeometry.Orbit(sv)
    el = ImagePairGeometry.Ellipsoid()
    to_grid = FGP.Transformation(FGP.GeographicCRS(), FGP.native_crs(s.region.epsg))
    t0 = (first(sv.time) + last(sv.time)) / 2
    pos0, vel0 = ImagePairGeometry.interpolate(orb, t0)
    out_x, out_y = copy(dx), copy(dy)
    for (q, I) in enumerate(idx)
        lon, lat = deg2rad(lons[q]), deg2rad(lats[q])
        p = ImagePairGeometry.lonlat_to_xyz(el, SV3(lon, lat, 0.0))
        r = ImagePairGeometry.geo2rdr(orb, p, t0, t0, pos0, vel0)
        los = -r.look / sqrt(sum(abs2, r.look))
        δ = geocoded_shift(dz[q], los, r.velocity, _east(lon), _north(lon, lat))
        # The shift on the product grid: where the moved point lands, minus where the point is.
        llh = ImagePairGeometry.xyz_to_lonlat(el, p .+ δ)
        x1, y1 = to_grid(rad2deg(llh[1]), rad2deg(llh[2]))
        x0, y0 = to_grid(lons[q], lats[q])
        mx, my = map_offset_pixels(x1 - x0, y1 - y0, spacing)
        out_x[I] = dx[I] - Float32(mx)
        out_y[I] = dy[I] - Float32(my)
    end
    _record!(info, dz, ib, t2, "map offsets")
    return out_x, out_y
end

# The midpoint of a NISAR scene's zero-Doppler start and end.
_mid_time(id::SLCDatasets.Identification) =
    (t = nisar_datetime(id.start_time); t + (nisar_datetime(id.stop_time) - t) ÷ 2)

# The median of a non-empty vector holding no NaN.
_median(v::Vector{Float64}) = (s = sort(v); m = length(s); isodd(m) ? s[(m + 1) ÷ 2] : (s[m ÷ 2] + s[m ÷ 2 + 1]) / 2)

"""
    map_offset_pixels(ex, ny, spacing) -> (px, py)

A map displacement `ex` metres along grid x and `ny` metres along grid y as the correlator's offsets: `px`
in columns, and `py` positive north in rows of height `|spacing[2]|`, the sign the map path's `dy`
carries after [`DY_SIGN`](@ref). `spacing` is the images' own pixel spacing, `(dx, dy)` with `dy`
negative for a north-up grid.
"""
map_offset_pixels(ex::Float64, ny::Float64, spacing::NTuple{2,Float64}) = (ex / spacing[1], ny / abs(spacing[2]))

# `yyyy-mm-ddTHH`, field by field: `Dates.format` is not available to a trimmed program.
_hour_string(t::DateTime) = string(Dates.year(t); pad = 4) * "-" * string(Dates.month(t); pad = 2) * "-" *
                            string(Dates.day(t); pad = 2) * "T" * string(Dates.hour(t); pad = 2)

