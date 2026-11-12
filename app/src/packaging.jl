# The `ItsLiveInput` metadata the product needs, and the product itself.


# `testautoRIFT.py:1363-1365`, for an optical pair.
const DETECTION_METHOD = "feature"
const MOTION_COORDINATES = "map"

"""
    scene_identification(scene) -> (OpticalDatasets.Identification, ItsLiveOffsetProduction.UtcMicros)

The scene's identity and its acquisition instant to the microsecond. A Sentinel-2 product name, like an
ASTER granule name, carries its acquisition time to the second; a Landsat product id does not carry the
time of day, which the reference takes from the scene's STAC item (`process.py:63-74`) and which the
scene's `_MTL.json` gives as the same instant at the same microsecond precision. The identification's
`acquisition_time` is that instant to the millisecond.
"""
function scene_identification(scene::SceneSource)
    id = startswith(scene.name, "S2") ? OpticalDatasets.sentinel2_identification(scene.name) :
         is_aster_id(scene.name) ? aster_identification(scene.name) : nothing
    id === nothing || return id, ItsLiveOffsetProduction.UtcMicros(id.acquisition_time)
    name, t = landsat_acquisition(scene.mtl)
    name == scene.name || throw(ArgumentError(scene.mtl * " describes " * name * ", not " * scene.name))
    return OpticalDatasets.landsat_identification(name, DateTime(t)), t
end

"""
    native_img_pair_info(s, roi_valid) -> ItsLiveOffsetProduction.ImagePairInfo

`IMG_INFO_DICT` for an optical pair, from the two scenes' identities, and each scene's pass
(`flight_direction_img1/2`, [`flight_direction`](@ref)) as a radar pair records its. `path`, `row` and
`collection_number` are written as `Float64` because the reference's own attributes are doubles. A
Sentinel-2 product carries none of path, row, collection and processing date. A NISAR L2 pair's comes from
its GSLCs' identification ([`nisar_img_pair_info`](@ref)).
"""
function native_img_pair_info(s::Setup, roi_valid::Float64)
    lon, lat = pair_centroid(s.pair.coordinate, s.epsg)
    eg, lg = s.early_gslc, s.late_gslc
    if eg !== nothing && lg !== nothing
        return nisar_img_pair_info(s.early.name, s.late.name, eg.product.identification, lg.product.identification,
                                   roi_valid, lon, lat)
    end
    (id1, t1), (id2, t2) = scene_identification(s.early), scene_identification(s.late)
    # Sized for its six entries, so the attributes iterate in the order a six-pair `Dict` constructor gives.
    extra = sizehint!(Dict{String,AttrValue}(), 6)
    a = Threads.@spawn flight_direction(s.early)
    b = Threads.@spawn flight_direction(s.late)
    extra["id_img1"] = id1.id
    extra["id_img2"] = id2.id
    extra["sensor_img1"] = id1.sensor
    extra["sensor_img2"] = id2.sensor
    extra["correction_level_img1"] = id1.correction_level
    extra["correction_level_img2"] = id2.correction_level
    _landsat_extra!(extra, 1, id1)
    _landsat_extra!(extra, 2, id2)
    extra["flight_direction_img1"] = _await(a)::String
    extra["flight_direction_img2"] = _await(b)::String
    return ItsLiveOffsetProduction.ImagePairInfo(t1, t2, id1.mission, id2.mission, id1.satellite, id2.satellite,
                                                 roi_valid, round(lat; digits = 2), round(lon; digits = 2), extra)
end

function _landsat_extra!(extra::Dict{String,AttrValue}, n::Int, id::OpticalDatasets.Identification)
    path, row, collection = id.path, id.row, id.collection_number
    category, processed = id.collection_category, id.processing_date
    (path === nothing || row === nothing || collection === nothing || category === nothing ||
     processed === nothing) && return extra
    k = string(n)
    extra["path_img" * k] = Float64(path)
    extra["row_img" * k] = Float64(row)
    extra["collection_number_img" * k] = Float64(collection)
    extra["collection_category_img" * k] = category
    extra["processing_date_img" * k] = processed
    return extra
end

"""
    product_source(info) -> String

The `source` attribute (`netcdf_output.py:427,438`). The year is the wall clock's, as it is there; the
attributed acquisition year is `acquisition_date_img1`'s.
"""
function product_source(info::ItsLiveOffsetProduction.ImagePairInfo)
    s = "NASA MEaSUREs ITS_LIVE project. Processed by ASF DAAC HyP3 " * string(Dates.year(Dates.now())) *
        " using autorift_itslive version " * app_version()
    startswith(info.mission_img1, "L") &&
        return s * ". Landsat-" * string(info.satellite_img1) * " images courtesy of the U.S. Geological Survey"
    info.mission_img1 == "S" && return s * ". Contains modified Copernicus Sentinel data " *
                                       string(Dates.year(info.acquisition_date_img1)) * ", processed by ESA"
    info.mission_img1 == "T" && return s * ". ASTER L1T data courtesy of NASA/METI/AIST/Japan Space Systems, " *
                                       "and U.S./Japan ASTER Science Team, distributed by NASA LP DAAC"
    return s
end

# Whether any point of `c` resolved.
resolved(c::Correlation) = any(!iszero, c.chip_size)

# `roi_valid_percentage` of `c`: 0 for a correlation that searched nothing, whose percentage has no denominator.
valid_percentage(c::Correlation) =
    any(!iszero, c.search) ? ItsLiveOffsetProduction.roi_valid_percentage(c.chip_size, c.search) : 0.0

"""
    write_named_product(output, early, late, roi, input) -> String

Write `input` to `output`, or, when `output` is empty, to the product's ITS_LIVE name for scenes `early`
and `late` and valid percentage `roi` in the working directory, and return the path written.
"""
function write_named_product(output::String, early::String, late::String, roi::Float64,
                             input::ItsLiveOffsetProduction.ItsLiveInput)
    path = isempty(output) ? ItsLiveOffsetProduction.product_filename(early, late, roi) : output
    isempty(dirname(path)) || mkpath(dirname(path))
    ItsLiveOffsetProduction.write_product(path, input)
    return path
end

"""
    package(s, c, output, lonlat, tide, tide_model)

Write the product of setup `s` and correlation `c` to `output`. `lonlat` is the grid CRS's transform to
longitude and latitude, whose type is fixed by the caller's branch.

The product grid is the geogrid window's own geotransform at cell centres, where a geotransform names the
outer edge. `pixel_size_x`/`pixel_size_y` are the *scene's* pixel size: they scale a chip size in pixels
into metres, and a chip is cut in the imagery.

A NISAR L2 pair's offsets over polar floating ice and ocean first have the tide and inverse-barometer motion
removed ([`tide_correct`](@ref)) when `tide` is set, with the region's own tide model or, when not empty,
`tide_model`. An ASTER pair's offsets have the L1T product's distortion and jitter removed
([`aster_correct`](@ref)), fitted over the stable surface; the constant shift is left for the product's own
stable-surface correction, and a pair the fit cannot be trusted on is refused ([`NoProduct`](@ref)).
`img_pair_info:geometric_correction` describes the fit and `geometric_correction_flag` is 1
(nasa-jpl/its_live_production#68). Any other optical pair is not corrected, nor is a correlation that resolved
nothing.
"""
function package(s::Setup, c::Correlation, output::String, lonlat::F, tide::Bool, tide_model::String,
                 no_velocity_reason::String = "") where {F}
    gt = ImagePairGeometry.window_geotransform(s.grid, s.window)
    nx, ny = size(s.window)
    x = collect(gt[1] + gt[2] / 2 .+ (0:(nx - 1)) .* gt[2])
    y = collect(gt[4] + gt[6] / 2 .+ (0:(ny - 1)) .* gt[6])
    px = Float64(ImagePairGeometry.xsize(s.pair.coordinate))
    py = abs(s.pair.coordinate.spacing[2])
    georef = ItsLiveOffsetProduction.ItsLiveGeoref(x, y, ItsLiveOffsetProduction.cf_grid_mapping(s.region.epsg),
                                                   lonlat, px, py)
    roi = valid_percentage(c)
    info = native_img_pair_info(s, roi)
    isempty(no_velocity_reason) || (info.extra["no_velocity_reason"] = no_velocity_reason)
    coeffs = ItsLiveOffsetProduction.coefficients(s.geometry)
    refv, ssm = ItsLiveOffsetProduction.reference_velocity(s.inputs)
    dx, dy = tide ? tide_correct(s, x, y, lonlat, c.dx, c.dy, info, tide_model) : (c.dx, c.dy)
    if is_aster(s) && resolved(c)
        # The offsets less the reference motion, so the fit sees only the geometry's error; added back after.
        dxr, dyr = ItsLiveOffsetProduction.reference_offsets(coeffs, refv)
        fit = aster_fit_mask(ssm, refv, dx)
        cx, cy = aster_correct(Float32.(dx .- dxr), Float32.(dy .- dyr), fit, x, y,
                               aster_track_frame(s.early.band, s.region.epsg))
        dx, dy = Float32.(cx .+ dxr), Float32.(cy .+ dyr)
        info.extra["geometric_correction"] = ASTER_GEOMETRIC_CORRECTION
        info.extra["geometric_correction_flag"] = Int64(1)
    end
    input = ItsLiveOffsetProduction.ItsLiveInput(
        :optical, DETECTION_METHOD, MOTION_COORDINATES, app_version(),
        replace(PARAMETER_SHAPEFILE, "/vsicurl/" => ""), product_source(info), nothing,
        georef, info, coeffs, refv, ssm, nothing, dx, dy, c.chip_size, c.scale, c.interpolated)
    return write_named_product(output, s.early.name, s.late.name, roi, input)
end

"""
    package(s, c, output; tide = true, tide_model = "")

[`package`](@ref) with the grid CRS's own transform to longitude and latitude — one concrete type per
branch.
"""
function package(s::Setup, c::Correlation, output::String; tide::Bool = true, tide_model::String = "",
                 no_velocity_reason::String = "")
    g = FGP.native_crs(s.region.epsg)
    g isa FGP.UTMCRS &&
        return package(s, c, output, FGP.Transformation(g, FGP.GeographicCRS()), tide, tide_model, no_velocity_reason)
    g isa FGP.PolarStereographicCRS &&
        return package(s, c, output, FGP.Transformation(g, FGP.GeographicCRS()), tide, tide_model, no_velocity_reason)
    throw(ArgumentError("a product grid is UTM or polar stereographic, not EPSG:" * string(s.region.epsg)))
end
