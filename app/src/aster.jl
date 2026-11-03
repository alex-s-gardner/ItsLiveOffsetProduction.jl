# ASTER L1T: the scene ID and the band the chain correlates, found through NASA's CMR catalog and read from
# LP DAAC; the pass and metadata the product records; and the corrections an ASTER pair needs that a Landsat
# or Sentinel-2 pair does not.
#
# The correlated band is VNIR band 3N (0.78–0.86 µm, 15 m, nadir). Over snow it saturates far less often than
# bands 1 and 2, but it does saturate at its normal gain, and a saturated area's edge is fixed in the image:
# where both scenes saturate, the correlator matches those edges and reads moving ice as stationary. Saturated
# samples and their neighbours are therefore invalid to the correlator ([`UnsaturatedMask`](@ref)).
#
# The L1T product's geolocation carries two errors the constant stable-surface shift every optical product
# receives does not remove: a smooth distortion across the scene, and the satellite's pitch and roll jitter,
# which shifts whole image lines and appears as sinusoids along and across the track. Both are fitted to the
# offsets over the stable surface and removed before the product's own stable shift ([`aster_correct`](@ref)).

const CMR_LPCLOUD =
    "https://cmr.earthdata.nasa.gov/search/granules.umm_json?provider=LPCLOUD&short_name=AST_L1T&version=004&page_size=2000"
# `AST_L1T_VVVMMDDYYYYHHMMSS_YYYYMMDDHHMMSS`: the processing version (`003`, `004`), the acquisition, then
# the production time.
const ASTER_ID = r"^AST_L1T_\d{3}(\d{2})(\d{2})(\d{4})(\d{2})(\d{2})(\d{2})_\d{14}$"
const ASTER_BAND_FILE = r"^(AST_L1T_\d{17}_\d{14})_VNIR_B03N\.tif$"

# The value a saturated band-3N sample holds, and how far around a saturated sample the correlator is
# kept off it: the resampling that put L1T on its grid spreads a saturated area's edge over about a pixel.
const ASTER_SATURATED = 0xff
const ASTER_SATURATION_GROW = 2
# The correlation a measurement must reach. ASTER pairs that lose their surface texture keep matches near
# zero motion with low correlation; below this the share of measurements more than 30 % from the motion of
# fast ice rises steeply.
const ASTER_MIN_CORRELATION = 0.3f0

is_aster_id(s::AbstractString) = occursin(ASTER_ID, s)

"""
    aster_acquisition(name) -> DateTime

The acquisition time an ASTER L1T granule's name carries, to the second.
"""
function aster_acquisition(name::AbstractString)
    m = match(ASTER_ID, name)
    m === nothing && throw(ArgumentError("\"" * name * "\" is not an ASTER L1T granule name"))
    c = [parse(Int, something(x)) for x in m.captures]
    return DateTime(c[3], c[1], c[2], c[4], c[5], c[6])
end

aster_date(name::AbstractString) = yyyymmdd(aster_acquisition(name))

"""
    aster_url(id) -> String

The band-3N GeoTIFF of ASTER L1T granule `id`, as LP DAAC serves it and CMR lists it.
"""
aster_url(id::String) = granule_link(only(cmr_lookup([id]; base = CMR_LPCLOUD, catalog = "LP DAAC's catalog (CMR)")),
                                     "_VNIR_B03N.tif")

"""
    aster_scene_source(input) -> SceneSource

The ASTER scene `input` names: a granule ID, resolved through CMR to its band-3N file at LP DAAC, or that
file itself as a path or URL.
"""
function aster_scene_source(input::String)
    is_aster_id(input) && return SceneSource(input, input, aster_url(input), "", "")
    m = match(ASTER_BAND_FILE, basename(input))
    m === nothing && throw(ArgumentError("\"" * input * "\" is neither an ASTER L1T granule ID nor its " *
                                         "<granule>_VNIR_B03N.tif, the band the chain correlates"))
    return SceneSource(input, String(something(m.captures[1])), input, "", "")
end

"""
    aster_identification(name) -> OpticalDatasets.Identification

An ASTER scene's identity: mission `"T"` (Terra), no satellite number (there was one Terra), sensor
`"ASTER"`, product level `"L1T"`, and the acquisition time the granule name carries. It carries none of
Landsat's path, row and collection.
"""
aster_identification(name::String) =
    OpticalDatasets.Identification(name, "T", "", "ASTER", "L1T", aster_acquisition(name), nothing, nothing, nothing,
                                   nothing, nothing)

# ---------------------------------------------------------------------------
# Reading the band
# ---------------------------------------------------------------------------

"""
    earthdata_transport(uri) -> ChunkManifests.LeafTransport

The transport an `https://` object on an Earthdata host is read through: an Earthdata Login bearer token,
from `EARTHDATA_TOKEN` or `~/.netrc` (`NETRC`), and the redirect to the pre-signed URL followed once.
"""
function earthdata_transport(uri::String)
    t = SLCDatasets.source_transport(SLCDatasets.RemoteHTTP(uri))
    t === nothing && throw(ArgumentError(uri * " is on an Earthdata host and no Earthdata Login credentials " *
                                         "were found; set EARTHDATA_TOKEN, or EARTHDATA_USERNAME and " *
                                         "EARTHDATA_PASSWORD, or add urs.earthdata.nasa.gov to ~/.netrc"))
    return t::ChunkManifests.HTTPTransport
end

"""
    stage_aster_band(s::SceneSource, dir) -> SceneSource

[`stage_band`](@ref) for an ASTER band, read with Earthdata Login ([`earthdata_transport`](@ref)). GDAL reads
the band too — its footprint and, for a pair in two UTM zones, the warp — and GDAL does not follow LP DAAC's
signed redirect with a bearer token, so the band is read from disk. One band-3N file is about 15 MB.
"""
stage_aster_band(s::SceneSource, dir::String) =
    is_http(s.band) ? stage_band(s, dir, earthdata_transport(s.band), "_VNIR_B03N.tif") : s

"""
    UnsaturatedMask(parent)

Where `parent` holds data the correlator may use: not zero (no data), and not within
[`ASTER_SATURATION_GROW`](@ref) samples of a saturated one ([`ASTER_SATURATED`](@ref)). Read a window at a
time, each window's margin read once.
"""
struct UnsaturatedMask{A<:AbstractMatrix} <: AbstractMatrix{Bool}
    parent::A
end

Base.size(m::UnsaturatedMask) = size(m.parent)
function Base.getindex(m::UnsaturatedMask, rows::AbstractUnitRange, cols::AbstractUnitRange)
    g = ASTER_SATURATION_GROW
    nr, nc = size(m.parent)
    R, C = _grown(rows, g, nr), _grown(cols, g, nc)
    v = m.parent[R, C]
    # Valid where there is data, then each saturated sample clears the box `g` around it.
    out = map(!=(0), view(v, (first(rows) - first(R) + 1):(last(rows) - first(R) + 1),
                         (first(cols) - first(C) + 1):(last(cols) - first(C) + 1)))::Matrix{Bool}
    for k in CartesianIndices(v)
        v[k] == ASTER_SATURATED || continue
        i, j = k[1] + first(R) - 1, k[2] + first(C) - 1
        for jj in max(j - g, first(cols)):min(j + g, last(cols)), ii in max(i - g, first(rows)):min(i + g, last(rows))
            out[ii - first(rows) + 1, jj - first(cols) + 1] = false
        end
    end
    return out
end
@_window_indexing UnsaturatedMask

# ---------------------------------------------------------------------------
# The along-track and stable-surface correction
# ---------------------------------------------------------------------------

# The pitch jitter's wavelength along the track, and the roll jitter's across it, are searched for in
# these ranges (metres): the scene's own peak is used, since it varies a little with the orbit.
const ASTER_ALONG_JITTER = (4000.0, 5600.0)
const ASTER_CROSS_JITTER = (5500.0, 9500.0)
# Profiles of the stable-surface offsets are binned at the product's grid spacing.
const ASTER_PROFILE_BIN = 120.0
# Fewer stable cells than this, and the scene keeps the product's constant stable shift alone.
const ASTER_MIN_STABLE = 2000

"""
    aster_track_frame(band, epsg) -> (time, detector)

Unit vectors, in coordinates of `epsg`, of the directions along which an ASTER scene's image-line
acquisition time and its detector position increase: the normals of the footprint's leading edge (an image
line) and of its side edge (one detector over time). Earth's rotation skews the two by about 3°. The
footprint's corners are the band's valid samples furthest along ±(x + y) and ±(x − y), read from an
eighth of its samples in each direction, which places the edges' directions within about 0.3°.
"""
function aster_track_frame(band::String, epsg::Int)
    corners, src = with_dataset(gdal_path(band)) do d
        gt = geotransform(d)
        nx, ny = rastersize(d)
        step = 8
        best = fill((-Inf, 0.0, 0.0), 4)
        for r0 in 0:(step * 64):(ny - 1)
            h = min(step * 64, ny - r0)
            w = read_window(d, UInt8, 0, r0, nx, h)
            for r in 1:step:h, c in 1:step:nx
                w[c, r] == 0 && continue
                x = gt[1] + (c - 0.5) * gt[2]; y = gt[4] + (r0 + r - 0.5) * gt[6]
                for (k, v) in enumerate((x + y, x - y, -x - y, -x + y))
                    v > best[k][1] && (best[k] = (v, x, y))
                end
            end
        end
        [(b[2], b[3]) for b in best], projected_epsg(projection_wkt(d))
    end
    p = src == epsg ? corners : with_transform(ct -> [transform_point(ct, c[1], c[2]) for c in corners], src, epsg)
    edges = [(p[mod1(i + 1, 4)][1] - p[i][1], p[mod1(i + 1, 4)][2] - p[i][2]) for i in 1:4]
    side = edges[argmax([abs(e[2]) / hypot(e...) for e in edges])]
    lead = edges[argmax([abs(e[1]) / hypot(e...) for e in edges])]
    side = side[2] > 0 ? side : .-side
    lead = lead[1] > 0 ? lead : .-lead
    unit(v) = v ./ hypot(v...)
    return unit((-lead[2], lead[1])), unit((side[2], -side[1]))
end

"""
    aster_correct(dx, dy, ssm, x, y, frame) -> (dx, dy)

`dx`, `dy` (pixels, on the product grid at cell centres `x`, `y`) with the ASTER L1T distortion removed, as
fitted over the stable surface `ssm`: per offset component, a quadratic in map position, cubics in the
along-track (time) and cross-track (detector) positions of `frame`, and sinusoids — with an amplitude that
varies linearly along their own axis — at the scene's along-track jitter wavelength and its first harmonic
and at its cross-track jitter wavelength. The fit is a Tukey-biweight robust least squares, so a stable cell
the correlator got wrong does not pull it. The constant term is left in the offsets: the product's own
stable-surface shift removes it, and records it as it does for every optical pair.

A scene with fewer than [`ASTER_MIN_STABLE`](@ref) stable cells is returned as it is.
"""
function aster_correct(dx::Matrix{Float32}, dy::Matrix{Float32}, ssm::AbstractMatrix{Bool},
                       x::Vector{Float64}, y::Vector{Float64}, frame)
    nr, nc = size(dx)
    (length(y), length(x)) == (nr, nc) || throw(DimensionMismatch("the offsets and the grid differ in size"))
    fit = ssm .& isfinite.(dx) .& isfinite.(dy)
    count(fit) < ASTER_MIN_STABLE && return dx, dy
    tv, dv = frame
    X = repeat(x', nr); Y = repeat(y, 1, nc)
    S = tv[1] .* X .+ tv[2] .* Y           # along-track (acquisition time) position, m
    T = dv[1] .* X .+ dv[2] .* Y           # cross-track (detector) position, m
    ddx, ddy = Float64.(dx), Float64.(dy)
    # The jitter wavelengths, from the cross-track offset after the quadratic. The detector direction is
    # within a few tens of degrees of east–west or north–south, so one map component carries most of it.
    cross = abs(dv[1]) >= abs(dv[2]) ? ddx : ddy
    Gq = _quadratic_columns(X, Y, fit)
    rq = cross .- reshape(_robust_fit(Gq, vec(cross), vec(fit)), nr, nc)
    ls = _peak_wavelength(S, rq, fit, ASTER_ALONG_JITTER)
    lt = _peak_wavelength(T, rq, fit, ASTER_CROSS_JITTER)
    G = hcat(Gq, _track_columns(S, T, fit, ls, lt))
    fx = reshape(_robust_fit(G, vec(ddx), vec(fit)), nr, nc)
    fy = reshape(_robust_fit(G, vec(ddy), vec(fit)), nr, nc)
    # Keep the constant: subtract the fit minus its value at the stable cells' median.
    cx = _median_nonan(fx[fit]); cy = _median_nonan(fy[fit])
    return Float32.(ddx .- (fx .- cx)), Float32.(ddy .- (fy .- cy))
end

_median_nonan(v) = ItsLiveOffsetProduction._median_nonan(v)
_mean_over(a, fit) = sum(a[fit]) / count(fit)

# A quadratic in map position, scaled to 10 km about the fitted cells' centre.
function _quadratic_columns(X, Y, fit::BitMatrix)
    x = vec((X .- _mean_over(X, fit)) ./ 1e4); y = vec((Y .- _mean_over(Y, fit)) ./ 1e4)
    return hcat(ones(length(x)), x, y, x .* x, x .* y, y .* y)
end

# Cubics in the standardized along- and cross-track positions, and the jitter sinusoids with an amplitude
# that varies linearly along their own axis: the along-track wavelength `ls` and its first harmonic, and the
# cross-track wavelength `lt`. A wavelength that could not be found (NaN) contributes no columns.
function _track_columns(S, T, fit::BitMatrix, ls::Float64, lt::Float64)
    sd(a, mu) = sqrt(sum((a[fit] .- mu) .^ 2) / count(fit)) + 1.0
    ms, mt = _mean_over(S, fit), _mean_over(T, fit)
    s = vec((S .- ms) ./ sd(S, ms)); t = vec((T .- mt) ./ sd(T, mt))
    # Filled column by column into one matrix: `reduce(hcat, ...)` over a vector of columns cannot be
    # compiled by `juliac --trim`.
    n = 2 + (isfinite(ls) ? 8 : 0) + (isfinite(lt) ? 4 : 0)
    G = Matrix{Float64}(undef, length(s), n)
    G[:, 1] .= s .^ 3; G[:, 2] .= t .^ 3
    k = 2
    if isfinite(ls)
        for h in 1:2, f in (cos, sin)
            G[:, k + 1] .= f.(2π .* h .* vec(S) ./ ls); G[:, k + 2] .= G[:, k + 1] .* s
            k += 2
        end
    end
    if isfinite(lt)
        for f in (cos, sin)
            G[:, k + 1] .= f.(2π .* vec(T) ./ lt); G[:, k + 2] .= G[:, k + 1] .* t
            k += 2
        end
    end
    return G
end

# Tukey-biweight robust least squares of `v` on `G` over the rows `fit`, evaluated at every row.
function _robust_fit(G::Matrix{Float64}, v::Vector{Float64}, fit::AbstractVector{Bool}; iters::Int = 6)
    A = G[fit, :]; b = v[fit]
    w = ones(length(b))
    coef = zeros(size(G, 2))
    for _ in 1:iters
        coef = (A .* w) \ (b .* w)
        r = b .- A * coef
        med = _median_nonan(r)
        s = 1.4826 * _median_nonan(abs.(r .- med)) + 1e-9
        u = r ./ (4.685 * s)
        w = [abs(q) < 1 ? (1 - q * q)^2 : 0.0 for q in u]
    end
    return G * coef
end

# The wavelength (m) of the strongest component, within `band`, of `v`'s median profile along `coord` over
# the cells `fit`; NaN when the profile is too short or gappy to say.
function _peak_wavelength(coord::Matrix{Float64}, v::Matrix{Float64}, fit::BitMatrix, band::NTuple{2,Float64})
    c0 = minimum(coord[fit])
    nb = floor(Int, (maximum(coord[fit]) - c0) / ASTER_PROFILE_BIN) + 1
    bins = [Float64[] for _ in 1:nb]
    for k in eachindex(coord)
        fit[k] && push!(bins[floor(Int, (coord[k] - c0) / ASTER_PROFILE_BIN) + 1], v[k])
    end
    prof = [length(b) >= 5 ? _median_nonan(b) : NaN for b in bins]
    good = findall(isfinite, prof)
    (length(good) < 64 || length(good) < nb ÷ 2) && return NaN
    f = similar(prof)
    for i in 1:nb   # linear between the nearest measured bins
        k = searchsortedfirst(good, i)
        f[i] = k == 1 ? prof[good[1]] : k > length(good) ? prof[good[end]] :
               prof[good[k - 1]] + (prof[good[k]] - prof[good[k - 1]]) * (i - good[k - 1]) / (good[k] - good[k - 1])
    end
    xs = collect(1.0:nb); xm = sum(xs) / nb; fm = sum(f) / nb
    slope = sum((xs .- xm) .* (f .- fm)) / sum((xs .- xm) .^ 2)
    f .-= fm .+ slope .* (xs .- xm)
    hann = [0.5 - 0.5 * cos(2π * (i - 1) / (nb - 1)) for i in 1:nb]
    # A direct Fourier sum at the candidate wavelengths: a profile is a few hundred bins.
    best, arg = -Inf, NaN
    for λ in range(band[1], band[2]; step = 20.0)
        k = ASTER_PROFILE_BIN / λ
        c = 0.0; sn = 0.0
        for i in 1:nb
            c += f[i] * hann[i] * cos(2π * k * (i - 1)); sn += f[i] * hann[i] * sin(2π * k * (i - 1))
        end
        p = c * c + sn * sn
        p > best && ((best, arg) = (p, λ))
    end
    return arg
end

"""
    aster_flight_direction(scene) -> String

`"ascending"` or `"descending"`: which way along the side edge of the scene's footprint the satellite moved.
Terra crosses the equator southbound in daylight, so a daytime scene is descending and a nighttime one
ascending; the side edge is the direction of travel up to its sign, and the sign is the one that agrees with
the local solar time of the acquisition — before 18:00 and after 06:00 descending.
"""
function aster_flight_direction(scene::SceneSource)
    t = aster_acquisition(scene.name)
    lon = with_dataset(gdal_path(scene.band)) do d
        gt = geotransform(d); nx, ny = rastersize(d)
        x = gt[1] + nx * gt[2] / 2; y = gt[4] + ny * gt[6] / 2
        with_transform(ct -> transform_point(ct, x, y), projected_epsg(projection_wkt(d)), 4326)[1]
    end
    solar = mod(Dates.hour(t) + Dates.minute(t) / 60 + lon / 15, 24)
    return 6 <= solar < 18 ? "descending" : "ascending"
end

"""
    aster_fit_mask(ssm, refv, dx) -> BitMatrix

The cells [`aster_correct`](@ref) is fitted over, chosen as the product's stable-surface correction chooses:
the stable-surface mask `ssm` where it holds a measurement, and otherwise the cells whose reference speed is
in its slowest quarter (`ItsLiveOffsetProduction.slowest_quarter`).
"""
aster_fit_mask(ssm::AbstractMatrix{Bool}, refv::ItsLiveOffsetProduction.ReferenceVelocity, dx::AbstractMatrix{Float32}) =
    any(ssm .& isfinite.(dx)) ? BitMatrix(ssm) : ItsLiveOffsetProduction.slowest_quarter(refv)
