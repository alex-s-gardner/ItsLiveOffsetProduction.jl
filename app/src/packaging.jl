# The `ItsLiveInput` metadata the product needs, and the product itself.

# The reference implementation being reproduced. Every golden product's `source` attribute names these two
# versions, and the container they come from is `ghcr.io/asfhyp3/hyp3-autorift:0.28.4`.
const PLUGIN_NAME = "hyp3_autorift"
const PLUGIN_VERSION = "0.28.4"
const AUTORIFT_VERSION = "2.1.1"

# `testautoRIFT.py:1363-1365`, for an optical pair.
const DETECTION_METHOD = "feature"
const MOTION_COORDINATES = "map"

"""
    scene_identification(scene) -> OpticalDatasets.Identification

The scene's identity. A Sentinel-2 product name carries its acquisition time to the second; a Landsat
product id does not carry the time of day, which the reference takes from the scene's STAC item
(`process.py:63-74`) and which the scene's `_MTL.json` gives as the same instant at the same millisecond
precision.
"""
function scene_identification(scene::SceneSource)
    startswith(scene.name, "S2") && return OpticalDatasets.sentinel2_identification(scene.name)
    id, t = landsat_acquisition(scene.mtl)
    id == scene.name || throw(ArgumentError(scene.mtl * " describes " * id * ", not " * scene.name))
    return OpticalDatasets.landsat_identification(id, t)
end

"""
    native_img_pair_info(s, roi_valid) -> ItsLiveOffsetProduction.ImagePairInfo

`IMG_INFO_DICT` for an optical pair, from the two scenes' identities. `path`, `row` and
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
    id1, id2 = scene_identification(s.early), scene_identification(s.late)
    # Sized for its six entries, so the attributes iterate in the order a six-pair `Dict` constructor gives.
    extra = sizehint!(Dict{String,AttrValue}(), 6)
    extra["id_img1"] = id1.id
    extra["id_img2"] = id2.id
    extra["sensor_img1"] = id1.sensor
    extra["sensor_img2"] = id2.sensor
    extra["correction_level_img1"] = id1.correction_level
    extra["correction_level_img2"] = id2.correction_level
    _landsat_extra!(extra, 1, id1)
    _landsat_extra!(extra, 2, id2)
    return ItsLiveOffsetProduction.ImagePairInfo(id1.acquisition_time, id2.acquisition_time,
                                                 id1.mission, id2.mission, id1.satellite, id2.satellite,
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
        " using the " * PLUGIN_NAME * " plugin version " * PLUGIN_VERSION * " running autoRIFT version " *
        AUTORIFT_VERSION
    startswith(info.mission_img1, "L") &&
        return s * ". Landsat-" * string(info.satellite_img1) * " images courtesy of the U.S. Geological Survey"
    info.mission_img1 == "S" && return s * ". Contains modified Copernicus Sentinel data " *
                                       string(Dates.year(info.acquisition_date_img1)) * ", processed by ESA"
    return s
end

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
    package(s, c, output, lonlat)

Write the product of setup `s` and correlation `c` to `output`. `lonlat` is the grid CRS's transform to
longitude and latitude, whose type is fixed by the caller's branch.

The product grid is the geogrid window's own geotransform at cell centres, where a geotransform names the
outer edge. `pixel_size_x`/`pixel_size_y` are the *scene's* pixel size: they scale a chip size in pixels
into metres, and a chip is cut in the imagery.
"""
function package(s::Setup, c::Correlation, output::String, lonlat::F) where {F}
    gt = ImagePairGeometry.window_geotransform(s.grid, s.window)
    nx, ny = size(s.window)
    x = collect(gt[1] + gt[2] / 2 .+ (0:(nx - 1)) .* gt[2])
    y = collect(gt[4] + gt[6] / 2 .+ (0:(ny - 1)) .* gt[6])
    px = Float64(ImagePairGeometry.xsize(s.pair.coordinate))
    py = abs(s.pair.coordinate.spacing[2])
    georef = ItsLiveOffsetProduction.ItsLiveGeoref(x, y, ItsLiveOffsetProduction.cf_grid_mapping(s.region.epsg),
                                                   lonlat, px, py)
    roi = ItsLiveOffsetProduction.roi_valid_percentage(c.chip_size, c.search)
    info = native_img_pair_info(s, roi)
    coeffs = ItsLiveOffsetProduction.coefficients(s.geometry)
    refv, ssm = ItsLiveOffsetProduction.reference_velocity(s.inputs)
    input = ItsLiveOffsetProduction.ItsLiveInput(
        :optical, DETECTION_METHOD, MOTION_COORDINATES, AUTORIFT_VERSION,
        replace(PARAMETER_SHAPEFILE, "/vsicurl/" => ""), product_source(info), nothing,
        georef, info, coeffs, refv, ssm, nothing, c.dx, c.dy, c.chip_size, c.scale, c.interpolated)
    return write_named_product(output, s.early.name, s.late.name, roi, input)
end

"""
    package(s, c, output)

[`package`](@ref) with the grid CRS's own transform to longitude and latitude — one concrete type per
branch.
"""
function package(s::Setup, c::Correlation, output::String)
    g = FGP.native_crs(s.region.epsg)
    g isa FGP.UTMCRS && return package(s, c, output, FGP.Transformation(g, FGP.GeographicCRS()))
    g isa FGP.PolarStereographicCRS && return package(s, c, output, FGP.Transformation(g, FGP.GeographicCRS()))
    throw(ArgumentError("a product grid is UTM or polar stereographic, not EPSG:" * string(s.region.epsg)))
end
