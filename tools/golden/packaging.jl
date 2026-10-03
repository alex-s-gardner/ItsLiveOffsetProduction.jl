# The `ItsLiveInput` metadata every packaging stage needs, regardless of which script calls it:
# `img_pair_info`, the georef, and the Sentinel-1 subswath-offset-bias correction. Shared by
# `julia_e2e.jl` (the byte-exact from-granule reproduction) and `e2e_run.jl` (which times packaging
# without reproducing the driver's byte quantization), so it lives here rather than in either.

using Dates
import Downloads
import JSON3
import GeoFormatTypes as GFT
import FastGeoProjections as FGP
using ItsLiveOffsetProduction
using OpticalDatasets

# The reference implementation being reproduced. Every golden product's `source` attribute names these
# two versions, and the container they come from is `ghcr.io/asfhyp3/hyp3-autorift:0.28.4`; neither is a
# property of the imagery, so neither can be derived from it.
const PLUGIN_NAME = "hyp3_autorift"
const PLUGIN_VERSION = "0.28.4"
const AUTORIFT_VERSION = "2.1.1"

# `testautoRIFT.py:1363-1365`, for an optical pair.
const DETECTION_METHOD = "feature"
const MOTION_COORDINATES = "map"

# `arImgDisp_*` negates `Dy0` as its first act and AutoRIFT.jl has no such internal step, so
# `MultichipResult.dy` runs opposite to the reference's `Dy` — which `testautoRIFT.py:800` passes
# straight into `netCDF_packaging` as `DY`, with no sign change of its own.
# Measured on the Landsat 8 case: this sign reaches 86.17% of `dy` exact against the reference's own,
# and the other 7.19%. `--compare-intermediate` prints both.
const DY_SIGN = -1

# ---------------------------------------------------------------------------
# Metadata, from the scenes' own STAC items
# ---------------------------------------------------------------------------

"""
    scene_identification(name) -> OpticalDatasets.Identification

The scene's STAC item, read and parsed into an `Identification`.

The acquisition *time of day* is the one field a Landsat product id does not carry, and the reference
takes it from the same place — `process.py:63-74` reads this item and parses `properties.datetime`. The
item is anonymous, unlike the imagery it points at.

`DateTime` is millisecond precision against the item's microseconds, so the time is truncated on the
way in; `ItsLiveOffsetProduction.ImagePairInfo`'s docstring records that the product's date attributes inherit that.
"""
function scene_identification(name::AbstractString)
    # A Sentinel-2 product name carries its own acquisition time to the second, so it needs no item at
    # all. The Landsat STAC call exists for the one field a Landsat id does not carry — the time of day
    # — so asking the Landsat collection for an `S2...` name is a 404, not a slower answer.
    startswith(name, "S2") && return sentinel2_identification(name)
    url = "https://landsatlook.usgs.gov/stac-server/collections/landsat-c2l1/items/$name"
    body = String(take!(Downloads.download(url, IOBuffer(); timeout = 60)))
    return open_optical(JSON3.read(body, Dict{String,Any}))
end

"""
    native_img_pair_info(c, s, roi_valid) -> ItsLiveOffsetProduction.ImagePairInfo

`IMG_INFO_DICT` for an optical pair, from the two product ids and their STAC items.

`path`, `row` and `collection_number` are written as `Float64` because the reference's own attributes
are doubles: it carries them through NumPy rather than as Python ints.
"""
function native_img_pair_info(c::GoldenCase, s::Setup, roi_valid::Real)
    early, late = acquisition_order(c)
    # Concurrently: these are two independent round trips to the same host, and run one after the
    # other they cost 0.47 s against 0.22 s together. A Sentinel-2 name needs no request at all, so
    # for an S2 pair both tasks are pure parsing and the spawn costs nothing. At `-t 1` this still
    # runs, just sequentially.
    t_early = Threads.@spawn scene_identification(early)
    t_late = Threads.@spawn scene_identification(late)
    id1, id2 = fetch(t_early), fetch(t_late)
    lon, lat = pair_centroid(s.pair.coordinate, s.epsg)
    extra = Dict{String,Any}(
        "id_img1" => id1.id, "id_img2" => id2.id,
        "sensor_img1" => id1.sensor, "sensor_img2" => id2.sensor,
        "correction_level_img1" => id1.correction_level,
        "correction_level_img2" => id2.correction_level)
    # Path, row, collection and processing date are Landsat's own. A Sentinel-2 product's
    # `img_pair_info` carries none of them — checked against the golden S2 products, whose variable has
    # exactly the six attributes above plus the ones `write_product` derives — and
    # `sentinel2_identification` leaves them `nothing`, so they are added only where they exist.
    for (n, id) in ((1, id1), (2, id2))
        id.path === nothing && continue
        extra["path_img$n"] = Float64(id.path)
        extra["row_img$n"] = Float64(id.row)
        extra["collection_number_img$n"] = Float64(id.collection_number)
        extra["collection_category_img$n"] = id.collection_category
        extra["processing_date_img$n"] = id.processing_date
    end
    return ItsLiveOffsetProduction.ImagePairInfo(id1.acquisition_time, id2.acquisition_time,
                                  id1.mission, id2.mission, id1.satellite, id2.satellite,
                                  Float64(roi_valid), round(lat; digits = 2), round(lon; digits = 2),
                                  extra)
end

# `netcdf_output.py:427,438`. The year is the wall clock's, as it is there.
#
# Sentinel-1 and Sentinel-2 share `mission_img1 == "S"` (`OpticalDatasets.sentinel2_identification`
# collapses both to the single-letter mission the rest of `ImagePairInfo` uses), so this doesn't
# distinguish them, but the clause doesn't need to: both are Copernicus Sentinel data, so both get
# it. The reference's own radar clause additionally names the ISCE3 version it ran — not reproduced
# here, since this package does not use ISCE3 and the string would be naming a tool that never ran.
# The attributed year is `acquisition_date_img1`'s, matching the Landsat clause's own choice of img1
# for the satellite it names — checked only on a pair where both acquisitions are the same year.
function product_source(info::ItsLiveOffsetProduction.ImagePairInfo)
    s = "NASA MEaSUREs ITS_LIVE project. Processed by ASF DAAC HyP3 $(Dates.year(Dates.now())) " *
        "using the $PLUGIN_NAME plugin version $PLUGIN_VERSION running autoRIFT version " *
        AUTORIFT_VERSION
    if startswith(info.mission_img1, "L")
        return s * ". Landsat-$(info.satellite_img1) images courtesy of the U.S. Geological Survey"
    end
    if info.mission_img1 == "S"
        return s * ". Contains modified Copernicus Sentinel data " *
               "$(Dates.year(info.acquisition_date_img1)), processed by ESA"
    end
    return s
end

"""
    native_georef(s::Setup) -> ItsLiveOffsetProduction.ItsLiveGeoref

The product grid, from the geogrid window's own geotransform.

`x`/`y` are **cell centres**, where a geotransform names the outer edge, so each gains half a pixel.
`pixel_size_x`/`pixel_size_y` are the *scene's* pixel size rather than the grid's — they scale a chip
size in pixels into metres, and a chip is cut in the imagery.
"""
function native_georef(s::Setup)
    gt = ImagePairGeometry.window_geotransform(s.grid, s.window)
    nx, ny = size(s.window)
    x = collect(gt[1] + gt[2] / 2 .+ (0:(nx - 1)) .* gt[2])
    y = collect(gt[4] + gt[6] / 2 .+ (0:(ny - 1)) .* gt[6])
    to_lonlat = FGP.Transformation(FGP.EPSG(Int(s.info.epsg)), FGP.EPSG(4326);
                                   always_xy = true)
    px = ImagePairGeometry.xsize(s.pair.coordinate)
    py = abs(s.pair.coordinate.spacing[2])
    return ItsLiveOffsetProduction.ItsLiveGeoref(x, y, ItsLiveOffsetProduction.cf_grid_mapping(GFT.EPSG(s.info.epsg)),
                                  (a, b) -> (to_lonlat(a, b)...,), px, py)
end

# ---------------------------------------------------------------------------
# Radar metadata
# ---------------------------------------------------------------------------

# Whether `c` reads its two scenes as a SAR acquisition (mosaic, RSLC or GSLC) rather than through
# `OpticalDatasets` — every platform `e2e_imagery` has a branch for.
radar_sensor(c::GoldenCase) = startswith(c.platform, "S1") || startswith(c.platform, "NISAR")

# Whether `c` writes the 16-variable `:radar` schema. NISAR-L2 is a radar sensor but a geocoded,
# map-grid product — `product.jl`'s `schema` calls it `:optical` for exactly that reason, and
# `setup` already routes it through the same `coregister`/`ProjectedCoordinate` path as Landsat and
# Sentinel-2 rather than through `_radar_setup`.
radar_pair(c::GoldenCase) = startswith(c.platform, "S1") || c.platform == "NISAR-L1"

# The pair's footprint midpoint in geodetic degrees. Mirrors `_radar_setup`'s own inline computation
# (`e2e.jl`): a radar footprint is solved for with `rdr2geo` rather than carried by a projection, so
# `pair_centroid` (which expects a `ProjectedCoordinate`) does not apply.
function radar_pair_centroid(s::Setup)
    b = footprint_bounds(IdentityTransform(), s.pair.coordinate)
    return ((b.X[1] + b.X[2]) / 2, (b.Y[1] + b.Y[2]) / 2)
end

"""
    radar_georef(s::Setup) -> ItsLiveOffsetProduction.ItsLiveGeoref

[`native_georef`](@ref) for a `RadarCoordinate` pair (S1 or NISAR-L1).

The grid itself (`x`/`y`, the CF mapping) is built identically to the optical path — `s.grid`/
`s.window` already come from the same `parameter_grid`/`grid_window` calls regardless of platform —
only `pixel_size_x` differs, since `s.pair.coordinate` carries no `.spacing` to read a scene pixel
size from.

`pixel_size_x` is the **ground**-range sample spacing, `coord.dr / sin(coord.incidence_angle)` —
slant range projected through the incidence angle at the scene centre, matching the reference's own
`rangePixelSize` to 6 significant figures on both a measured Sentinel-1 and a NISAR-L1 case (checked
against `capture_packaging`'s captured scalar, not assumed). `pixel_size_y` is unused for `:radar`:
`write_product` reads `input.georef.pixel_size_y` nowhere — the `azimuth_pixel_size` attribute and the
`vr`/`va` conversion both use a `pixel_size_y` recomputed from `offset2va`/`dt_seconds` instead (see
`ItsLiveGeoref`'s docstring) — so it is set equal to `pixel_size_x` rather than carrying a second,
unverified formula for a value nothing reads.
"""
function radar_georef(s::Setup)
    gt = ImagePairGeometry.window_geotransform(s.grid, s.window)
    nx, ny = size(s.window)
    x = collect(gt[1] + gt[2] / 2 .+ (0:(nx - 1)) .* gt[2])
    y = collect(gt[4] + gt[6] / 2 .+ (0:(ny - 1)) .* gt[6])
    to_lonlat = FGP.Transformation(FGP.EPSG(Int(s.info.epsg)), FGP.EPSG(4326);
                                   always_xy = true)
    coord = s.pair.coordinate
    px = coord.dr / sin(coord.incidence_angle)
    return ItsLiveOffsetProduction.ItsLiveGeoref(x, y, ItsLiveOffsetProduction.cf_grid_mapping(GFT.EPSG(s.info.epsg)),
                                  (a, b) -> (to_lonlat(a, b)...,), px, px)
end

"""
    s1_swath_bias(c, s, p) -> Union{ItsLiveOffsetProduction.SwathOffsetBias,Nothing}

The Sentinel-1 subswath-offset-bias correction's own inputs, or `nothing` where the correction has no
meaning — which is most of this manifest.

`nothing` for every platform but `S1-SLC`: a burst job opens one or two subswaths
([`burst_swaths`](@ref)), which `SLCDatasets.subswath_borders` refuses outright, matching
`ItsLiveOffsetProduction.SwathOffsetBias`'s own docstring ("or for a pair whose product was not opened
for all three subswaths"). For an `S1-SLC` pair it is `nothing` whenever the two acquisitions are the
same spacecraft — `subswath_borders`' own `same_platform` — since the correction removes a difference
between two *different* antennas' patterns and has no meaning within one. Of this manifest's 5
full-SLC cases, 4 are same-platform and only the fifth (`S1A` paired with `S1B`) gets a real one.

`p` is `AutoRIFT.params(s.geometry; ...)`, already computed by the caller; `grid_spacing_x` is its
`grid_spacing.X`, which the correction's coarse-grid smoothing step reads.
"""
function s1_swath_bias(c::GoldenCase, s::Setup, p)
    c.platform == "S1-SLC" || return nothing
    rp, sp = _s1_products(c, s.run)
    bd = SLCDatasets.subswath_borders(rp, sp)
    bd.same_platform && return nothing
    lx, ly = ItsLiveOffsetProduction.image_location(s.geometry)
    return ItsLiveOffsetProduction.SwathOffsetBias(lx, ly, Float64(p.grid_spacing.X), bd.ncols,
                                                   Float64(bd.border12), Float64(bd.border23),
                                                   bd.same_platform, bd.reference_platform)
end

# A `UtcTime`'s instant as a `DateTime`, truncated to the millisecond `ImagePairInfo.acquisition_date_img1`
# carries (see its docstring: the reference's own microsecond precision is lost the same way on the
# optical path).
_s1_datetime(t) = t.datetime + Dates.Millisecond(round(Int, t.seconds * 1000))

# Splits a Sentinel-1 product id (or `burst2safe`'s synthesized SAFE name, same layout —
# `S1A_IW_SLC__1SSH_<start>_<stop>_<orbit>_<datatake>_<productid>`, ten fields because `SLC__1SSH`
# contributes an empty one) into the fields `img_pair_info`'s `extra` carries verbatim. Measured
# against a real `capture_packaging` capture: field 1's tail is `satellite_img` ("1A" out of "S1A"),
# and fields 8/9/10 are `absolute_orbit_number`/`mission_data_take_ID`/`product_unique_ID` exactly.
function _s1_name_fields(name::AbstractString)
    parts = split(name, '_')
    length(parts) == 10 || throw(ArgumentError(
        "\"$name\" does not split into a Sentinel-1 product id's 10 underscore-delimited fields"))
    # `String`, not the bare `SubString` `split`/indexing return: a `SubString` that doesn't reach
    # the parent's own end carries no null terminator of its own, and the netCDF attribute write
    # reads straight through to the parent's — every field but the last (`product_id`, which does
    # reach the end) came out on disk with the rest of `name` appended until this was added.
    return (; satellite = String(parts[1][2:end]), orbit = String(parts[8]),
            datatake = String(parts[9]), product_id = String(parts[10]))
end

"""
    s1_img_pair_info(c, s, roi_valid) -> ItsLiveOffsetProduction.ImagePairInfo

`IMG_INFO_DICT` for a Sentinel-1 pair (full-SLC or burst), from the parsed products
[`_s1_products`](@ref) already builds for the geometry and mosaic stages.

`_s1_products` reads whichever SAFE or granule zip sits beside the run's outputs and returns the pair
reference-first (its own docstring: "earlier acquisition of the two") — annotation parsing needs no
unpacked raster, so this works even when the pixel data is not yet staged. Every `extra` field below
and `mission_img`/`satellite_img` were checked against a real `capture_packaging` capture of this
exact case; `sensor_img` ("C", Sentinel-1's C-band instrument) is the one constant not derivable from
either the id or the annotation.

**`acquisition_date_img1`/`img2` carry a measured, so far unexplained offset of order one second
against that same capture — `date_dt`, their difference, does not**: checked on this case, `s1_mosaic`'s
two `start`s differ from the reference's own `acquisition_date_img1`/`img2` by about +1.5 s each, but
the gap between them agrees with the reference's own `date_dt` to about 1.6 ms, well inside
`DateTime`'s millisecond rounding. So whatever `s1_mosaic`'s `start` is anchored to differs from the
reference's own epoch for the *absolute* instant by a nearly constant shift, while the physically
meaningful quantity — the interval between the two acquisitions, which is what feeds `dt_seconds`,
`vr`/`va` and the error-vector scaling — is unaffected. Worth checking against the reference's own
`IMG_INFO_DICT` population before trusting the absolute date attributes on a product this writes.
"""
function s1_img_pair_info(c::GoldenCase, s::Setup, roi_valid::Real)
    rp, sp = _s1_products(c, s.run)
    sws = burst_swaths(c)
    _, t_early = s1_mosaic(SafeSwaths(rp), sws)
    _, t_late = s1_mosaic(SafeSwaths(sp), sws)
    a1, a2 = annotation(rp, first(sws)), annotation(sp, first(sws))
    name(p) = replace(basename(p.path), r"\.(SAFE|zip)$" => "")
    f1, f2 = _s1_name_fields(name(rp)), _s1_name_fields(name(sp))
    lon, lat = radar_pair_centroid(s)
    extra = Dict{String,Any}(
        "id_img1" => name(rp), "id_img2" => name(sp),
        "sensor_img1" => "C", "sensor_img2" => "C",
        "flight_direction_img1" => lowercase(a1.pass_direction),
        "flight_direction_img2" => lowercase(a2.pass_direction),
        "absolute_orbit_number_img1" => f1.orbit, "absolute_orbit_number_img2" => f2.orbit,
        "mission_data_take_ID_img1" => f1.datatake, "mission_data_take_ID_img2" => f2.datatake,
        "product_unique_ID_img1" => f1.product_id, "product_unique_ID_img2" => f2.product_id)
    return ItsLiveOffsetProduction.ImagePairInfo(_s1_datetime(t_early), _s1_datetime(t_late), "S", "S",
                                  f1.satellite, f2.satellite, Float64(roi_valid),
                                  round(lat; digits = 2), round(lon; digits = 2), extra)
end

# A NISAR `zeroDopplerStartTime`/`zeroDopplerEndTime`-style HDF5 string
# ("2025-10-28T23:52:01.000000000", nanosecond precision) as a `DateTime`, truncated the same way
# `_s1_datetime` truncates a `UtcTime`.
function _nisar_datetime(s::AbstractString)
    whole, frac = split(s, '.'; limit = 2)
    dt = Dates.DateTime(whole, Dates.dateformat"yyyy-mm-ddTHH:MM:SS")
    return dt + Dates.Millisecond(round(Int, parse(Float64, "0." * frac) * 1000))
end

# `IMG_INFO_DICT` common to both NISAR schemas, from the two products' `Identification`.
# `mission_img`/`satellite_img` are the literal constants `"N"`/`1` the reference writes for every
# NISAR case regardless of which satellite flew it — checked against a real `capture_packaging`
# capture of both the L1 and the L2 case, not derived from `absolute_orbit` or any other field.
function _nisar_img_pair_info(early::AbstractString, late::AbstractString, id1, id2,
                              roi_valid::Real, lon::Real, lat::Real)
    extra = Dict{String,Any}(
        "id_img1" => early, "id_img2" => late,
        "sensor_img1" => "L", "sensor_img2" => "L",
        "flight_direction_img1" => lowercase(id1.pass_direction),
        "flight_direction_img2" => lowercase(id2.pass_direction),
        "absolute_orbit_number_img1" => id1.absolute_orbit,
        "absolute_orbit_number_img2" => id2.absolute_orbit)
    return ItsLiveOffsetProduction.ImagePairInfo(_nisar_datetime(id1.start_time), _nisar_datetime(id2.start_time),
                                  "N", "N", 1, 1, Float64(roi_valid),
                                  round(lat; digits = 2), round(lon; digits = 2), extra)
end

function nisar_l1_img_pair_info(c::GoldenCase, s::Setup, roi_valid::Real)
    run = rslc_run(c)
    run === nothing && error("$(basename(c.product)): no cached run holds both RSLC granules")
    early, late = acquisition_order(c)
    id1 = open_slc(joinpath(run, early * ".h5")).identification
    id2 = open_slc(joinpath(run, late * ".h5")).identification
    lon, lat = radar_pair_centroid(s)
    return _nisar_img_pair_info(early, late, id1, id2, roi_valid, lon, lat)
end

function nisar_l2_img_pair_info(c::GoldenCase, s::Setup, roi_valid::Real)
    run = gslc_run(c)
    run === nothing && error("$(basename(c.product)): no cached run holds both GSLC granules")
    early, late = acquisition_order(c)
    id1 = open_geocoded(joinpath(run, early * ".h5")).identification
    id2 = open_geocoded(joinpath(run, late * ".h5")).identification
    lon, lat = pair_centroid(s.pair.coordinate, s.epsg)
    return _nisar_img_pair_info(early, late, id1, id2, roi_valid, lon, lat)
end

"""
    radar_img_pair_info(c, s, roi_valid) -> ItsLiveOffsetProduction.ImagePairInfo

[`native_img_pair_info`](@ref) for a radar-sensor pair, dispatched to the platform's own metadata
source: `s1_img_pair_info` for Sentinel-1 (full-SLC or burst), [`nisar_l1_img_pair_info`](@ref) /
[`nisar_l2_img_pair_info`](@ref) for NISAR's RSLC and GSLC identification groups.
"""
function radar_img_pair_info(c::GoldenCase, s::Setup, roi_valid::Real)
    startswith(c.platform, "S1") && return s1_img_pair_info(c, s, roi_valid)
    c.platform == "NISAR-L2" && return nisar_l2_img_pair_info(c, s, roi_valid)
    return nisar_l1_img_pair_info(c, s, roi_valid)
end
