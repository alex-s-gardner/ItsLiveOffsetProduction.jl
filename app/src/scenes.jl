# Scene names, bands, filters and footprints for the optical sensors, and the dates and centroids the
# geometry needs.

"""
    scene_band(platform) -> Symbol

Which band the reference correlates: `:pan` for Landsat 7/8/9 and Sentinel-2, `:green` for Landsat 4/5
(`process.py:77-89`).
"""
function scene_band(platform::AbstractString)
    platform in ("L4", "L5") && return :green
    platform in ("L7", "L8", "L9", "S2") && return :pan
    throw(ArgumentError("no band rule for platform \"" * String(platform) * "\""))
end

# The band file a Landsat scene's correlated band is: green is B2 on TM, the panchromatic band B8.
landsat_band_file(platform::AbstractString) = scene_band(platform) === :green ? "B2" : "B8"

"""
    acquisition_order(c::Job) -> (early, late)

`c`'s two scenes ordered by acquisition date, which is the order the pipeline correlates them in
whatever order the job listed them. Two of the golden jobs name the later acquisition as the
reference, and their products report the earlier one as `id_img1` with a positive `date_dt`; following
the job's order instead negates the offsets and inflates the search range.
"""
function acquisition_order(c::Job)
    a, b = c.reference, c.secondary
    da, db = scene_date(c.platform, a.name), scene_date(c.platform, b.name)
    return da <= db ? (a, b) : (b, a)
end

# Where the acquisition date sits in a scene name: Landsat field 4 (`LC08_L1TP_009011_20200703_...`),
# Sentinel-2 the first eight characters of field 3 (`S2B_MSIL1C_20200612T150759_...`), NISAR those of
# field 12, as the reference reads each.
function scene_date(platform::AbstractString, name::AbstractString)
    parts = split(name, '_')
    startswith(platform, "L") && return String(parts[4])
    platform == "S2" && return String(parts[3][1:8])
    # NISAR names it as field 12, `20251028T235201`.
    startswith(platform, "NISAR") && return String(parts[12][1:8])
    # Sentinel-1: `S1A_IW_SLC__1SSH_20150828T162412_…`, field 6 after the empty one `SLC__` makes.
    startswith(platform, "S1") && return String(parts[6][1:8])
    throw(ArgumentError("no date rule for platform \"" * String(platform) * "\""))
end

"""
    geogrid_seconds(c::Job) -> Float64

The interval the geogrid is given: whole calendar days between the two scene dates, in seconds.
`testGeogridOptical.py:161-165` differences two `datetime.date`s built from the names, so the time of day
is discarded; the product's `date_dt` is a different quantity, and using it here moves offsets that sit
near a rounding boundary.
"""
function geogrid_seconds(c::Job)
    early, late = acquisition_order(c)
    d(s::SceneSource) = parse_yyyymmdd(scene_date(c.platform, s.name))
    return Float64(Dates.value(d(late) - d(early))) * 86400.0
end

"""
    CorrelatorFilter

The filter `runAutorift` applies to the cropped pair (`testautoRIFT.py:293-308,533-540`): a reduction over
the pair's scenes, most stringent first. An `L[EO]07_` scene acquired on or after 2003-05-31, after the
Scan Line Corrector failed, makes it `WALLIS_FILL`; an `LT0[45]_` scene makes it `NO_FILTER`, the
reference's commented-out `fft` branch, because those scenes reach the correlator carrying their native
filter; anything else is `HIGHPASS`. The width is 5 for every optical sensor.
"""
@enum CorrelatorFilter HIGHPASS WALLIS_FILL NO_FILTER

const FILTER_WIDTH = 5

# The high-pass width: 21 for Sentinel-1, 5 for every other sensor (`testautoRIFT.py:526-528`).
filter_width(c::Job) = startswith(c.platform, "S1") ? 21 : FILTER_WIDTH
const WALLIS_STD_CUTOFF = 0.25

function correlator_filter(c::Job)
    methods = (_scene_method(c.reference.name), _scene_method(c.secondary.name))
    :wallis_fill in methods && return WALLIS_FILL
    :fft in methods && return NO_FILTER
    return HIGHPASS
end

function _scene_method(name::AbstractString)
    if occursin(r"^L[EO]07_", name)
        acquired = parse_yyyymmdd(split(name, '_')[4])
        return acquired >= Date(2003, 5, 31) ? :wallis_fill : :hps
    end
    occursin(r"^LT0[45]_", name) && return :fft
    return :hps
end

"""
    native_filter(c::Job, name) -> Symbol

The filter `process.py` applies to `name`'s *native* scene, or `:none` (`process.py:471-479`,
`apply_landsat_filtering`). Gated on the pair — only when `min(reference_platform, secondary_platform)` is
L4, L5 or L7 — and then chosen per scene: `:fft` (Wallis at width 5, then the band-reject
`AutoRIFT.Destripe` reproduces) for L4 and L5, `:wallis_fill` (`AutoRIFT.WallisGapfill(5, 0.25)`) for L7
and L8.
"""
function native_filter(c::Job, name::AbstractString)
    p1, p2 = _landsat_platform(c.reference.name), _landsat_platform(c.secondary.name)
    (isempty(p1) || isempty(p2)) && return :none
    min(p1, p2) in ("L4", "L5", "L7") || return :none
    own = _landsat_platform(name)
    own in ("L4", "L5") && return :fft
    own in ("L7", "L8") && return :wallis_fill
    throw(ArgumentError("no native filter rule for \"" * String(name) * "\"; apply_landsat_filtering " *
                        "dispatches on L4, L5, L7 and L8 only"))
end

# `LT05_...` gives "L5", the platform code `get_platform` returns — the string comparison `min` relies on
# orders L4 < L5 < L7 < L8 < L9 — or "" for a scene that is not Landsat.
function _landsat_platform(name::AbstractString)
    m = match(r"^L[TEOC]0?(\d)_", name)
    return m === nothing ? "" : "L" * String(something(m.captures[1]))
end

# ---------------------------------------------------------------------------
# Footprints
# ---------------------------------------------------------------------------

const Footprint = ImageFootprint{Float64,GFT.EPSG{1}}

"""
    scene_footprint(uri) -> Footprint

The image's footprint — origin, signed spacing, size — and the EPSG code of its projected CRS, read from
its header.
"""
function scene_footprint(uri::String)
    return with_dataset(gdal_path(uri)) do d
        gt = geotransform(d)
        (gt[3] == 0 && gt[5] == 0) || throw(ArgumentError(
            uri * " has a rotated geotransform; the optical geogrid assumes an axis-aligned grid"))
        nx, ny = rastersize(d)
        ImageFootprint((gt[1], gt[4]), (gt[2], gt[6]), (nx, ny), GFT.EPSG(projected_epsg(projection_wkt(d))))
    end
end

footprint_epsg(fp::Footprint) = Int(GFT.val(fp.crs))

"""
    warped_footprint(uri, reference, target) -> Footprint

The grid `ensure_same_projection` would warp the scene at `uri` onto — `target`'s projection at
`reference`'s spacing, `-tap`-aligned — read off a warped VRT, which GDAL builds from the source's header
alone.
"""
function warped_footprint(uri::String, reference::Footprint, target::Int)
    flags = ["-of", "VRT", "-t_srs", "EPSG:" * string(target),
             "-tr", string(abs(reference.spacing[1])), string(abs(reference.spacing[2])),
             "-r", "lanczos", "-tap"]
    return with_dataset(gdal_path(uri)) do src
        v = warp(src, flags, vsimem_name("warped_footprint") * ".vrt")
        try
            gt = geotransform(v)
            nx, ny = rastersize(v)
            ImageFootprint((gt[1], gt[4]), (gt[2], gt[6]), (nx, ny), GFT.EPSG(target))
        finally
            close(v)
        end
    end
end

"""
    pair_centroid(coord, epsg) -> (lon, lat)

The centre of the coregistered coordinate in WGS84 degrees, which is what the parameter-region lookup
takes: the region has to be the one covering the overlap.
"""
function pair_centroid(coord::ProjectedCoordinate{Float64}, epsg::Int)
    cx = coord.origin[1] + (coord.size[1] - 1) * coord.spacing[1] / 2
    cy = coord.origin[2] + (coord.size[2] - 1) * coord.spacing[2] / 2
    return with_transform(ct -> transform_point(ct, cx, cy), epsg, 4326)
end
