# A job: the sensor, and the two scenes resolved from the inputs the caller gave — a path, an `s3://` or
# `https://` URL, a STAC item, or a scene ID (`sceneid.jl`).

"""
    SensorSpec

A sensor code the command line accepts: the platform name the chain dispatches on and the pattern the
reference scene's name follows. A Landsat pair may mix missions — an L8 scene against an L7 one — so the
secondary only has to be a Landsat scene.
"""
struct SensorSpec
    code::String
    platform::String
    reference::Regex
    secondary::Regex
end

const LANDSAT_SCENE = r"^L[TEOC]0[45789]_"
# A Sentinel-2 product's `.SAFE` directory in a path or URL, capturing the product name.
const S2_SAFE = r"(S2[ABC]_MSIL1C_[0-9A-Z_]+?)\.SAFE"

const SENSORS = (
    SensorSpec("L4", "L4", r"^LT04_", LANDSAT_SCENE),
    SensorSpec("L5", "L5", r"^LT05_", LANDSAT_SCENE),
    SensorSpec("L7", "L7", r"^LE07_", LANDSAT_SCENE),
    SensorSpec("L8", "L8", r"^LC08_", LANDSAT_SCENE),
    SensorSpec("L9", "L9", r"^LC09_", LANDSAT_SCENE),
    SensorSpec("S2", "S2", r"^S2[ABC]_MSIL1C_", r"^S2[ABC]_MSIL1C_"),
    SensorSpec("AST", "AST", r"^AST_L1T_", r"^AST_L1T_"),
    SensorSpec("S1", "S1-SLC", r"^S1[ABCD]_IW_SLC__", r"^S1[ABCD]_IW_SLC__"),
    SensorSpec("S1BURST", "S1-BURST", r"^S1[ABCD]_IW_SLC__", r"^S1[ABCD]_IW_SLC__"),
    SensorSpec("NL1", "NISAR-L1", r"^NISAR_L1_\w+_RSLC_", r"^NISAR_L1_\w+_RSLC_"),
    SensorSpec("NL2", "NISAR-L2", r"^NISAR_L2_\w+_GSLC_", r"^NISAR_L2_\w+_GSLC_"),
)

# The sensors the command line names that this build does not process.
const UNBUILT_SENSORS = String[]

# Whether `s` names a sensor the command line knows, built or not.
is_sensor_code(s::String) = (u = uppercase(s); any(x -> x.code == u, SENSORS) || u in UNBUILT_SENSORS)

function sensor_spec(code::String)
    for s in SENSORS
        s.code == code && return s
    end
    code in UNBUILT_SENSORS && throw(ArgumentError(
        "sensor " * code * " is not processed by this build; it processes " * join([s.code for s in SENSORS], ", ")))
    throw(ArgumentError("unknown sensor \"" * code * "\"; expected one of " *
                        join([s.code for s in SENSORS], ", ") *
                        (isempty(UNBUILT_SENSORS) ? "" : " or " * join(UNBUILT_SENSORS, ", "))))
end

"""
    SceneSource

One scene, resolved from its input: the granule `name` every per-granule rule reads (dates, filters,
metadata), where the correlated `band` is, and — for Landsat — where its `_MTL.json` and `_ANG.txt` are.
Each location is a local path, an `s3://` URI or an `https://` URL.
"""
struct SceneSource
    input::String
    name::String
    band::String
    mtl::String
    ang::String
end

"""
    Job

One pair to process: the sensor and its platform, the two scenes, the product name for messages, and the
OPERA frame ID a Sentinel-1 product records as `frame_img1`/`frame_img2` (empty for none).
"""
struct Job
    sensor::String
    platform::String
    reference::SceneSource
    secondary::SceneSource
    product::String
    frame_id::String
end

function Job(sensor::String, reference::String, secondary::String; subset::String = "", frame_id::String = "")
    spec = sensor_spec(sensor)
    !isempty(frame_id) && !startswith(spec.platform, "S1") && throw(ArgumentError(
        "--frame-id names a Sentinel-1 pair's OPERA frame; sensor " * sensor * " takes none"))
    given = (reference, secondary)
    spec.platform == "S1-BURST" && ((reference, secondary) = resolve_burst_ids(reference, secondary))
    if !isempty(subset)
        spec.platform == "S1-BURST" || throw(ArgumentError(
            "--subset selects the bursts of a Sentinel-1 burst pair (S1BURST); sensor " * sensor * " takes none"))
        (is_zip_list(reference) && is_zip_list(secondary)) || throw(ArgumentError(
            "--subset selects bursts out of SLC zips: give each scene as its zips, comma-separated"))
        reference, secondary = subset_inputs(reference, secondary, subset)
    elseif occursin(',', reference) && is_zip_list(reference)
        throw(ArgumentError("the reference names several SLC zips; a burst pair from zips needs --subset, " *
                            "the polygon its bursts are selected by"))
    end
    ref = scene_source(reference, spec.platform)
    sec = scene_source(secondary, spec.platform)
    occursin(spec.reference, ref.name) || throw(ArgumentError(
        "the reference input \"" * given[1] * "\" is scene " * ref.name * ", which is not a " * sensor * " scene"))
    occursin(spec.secondary, sec.name) || throw(ArgumentError(
        "the secondary input \"" * given[2] * "\" is scene " * sec.name * ", which cannot pair with a " *
        sensor * " scene"))
    return Job(sensor, spec.platform, ref, sec, ref.name * "_X_" * sec.name, frame_id)
end

"""
    scene_source(input, platform) -> SceneSource

The scene `input` names. A Landsat input is the band file the reference correlates — `_B8.TIF` for
Landsat 7, 8 and 9, `_B2.TIF` (green) for 4 and 5 — with its `_MTL.json` and `_ANG.txt` beside it, or a
landsatlook STAC item, whose assets name all three. A Sentinel-2 input is the `B08.jp2` band inside a
`.SAFE`, or an Earth Search STAC item, which is resolved through its `s2:product_uri` to the band in the
public `s3://sentinel-s2-l1c` bucket. An ASTER input is the granule's band-3N `_VNIR_B03N.tif`. A NISAR input is the granule's `.h5`, a Sentinel-1 input the product's `.SAFE`
directory or its `.zip`. A burst pair's input is the bursts' extractor URLs, or bursts of SLC zips
(`<zip>#IW<n>/<POL>/<index>`), comma-separated, whose scene is the product `burst2safe` would assemble from
them; or that product's `.SAFE`. Any but a burst pair's may also be the scene's ID (`sceneid.jl`): a Landsat
ID names its landsatlook STAC item, a Sentinel-2 ID its band in that bucket, a Sentinel-1 ID its zip at ASF
a NISAR ID the `.h5` CMR lists, and an ASTER ID the band-3N file CMR lists at LP DAAC.
"""
function scene_source(input::String, platform::String)
    if startswith(platform, "S1")
        is_s1_slc_id(input) && return SceneSource(input, s1_slc_name(input), s1_slc_url(s1_slc_name(input)), "", "")
        if platform == "S1-BURST" && is_burst_list(input)
            return SceneSource(input, SLCDatasets.asf_safe_name(burst_list(input, nothing)), input, "", "")
        end
        if platform == "S1-BURST" && is_zip_burst_list(input)
            return SceneSource(input, SLCDatasets.asf_safe_name(zip_burst_list(input, nothing)), input, "", "")
        end
        base = basename(rstrip(input, '/'))
        name = endswith(base, ".SAFE") ? base[1:(end - 5)] : endswith(base, ".zip") ? base[1:(end - 4)] :
               throw(ArgumentError("\"" * input * "\" is neither a Sentinel-1 product's .SAFE, its .zip nor its ID"))
        return SceneSource(input, String(name), input, "", "")
    end
    if startswith(platform, "NISAR")
        occursin(NISAR_ID, input) && return SceneSource(input, input, nisar_url(input), "", "")
        base = basename(input)
        endswith(base, ".h5") || throw(ArgumentError("\"" * input * "\" is neither a NISAR granule's .h5 nor its ID"))
        return SceneSource(input, base[1:(end - 3)], input, "", "")
    end
    platform == "AST" && return aster_scene_source(input)
    if platform == "S2"
        is_stac_item(input) && return _s2_from_item(input)
        occursin(S2_ID, input) && return SceneSource(input, input, s2_band_url(input), "", "")
        m = match(S2_SAFE, input)
        m === nothing && throw(ArgumentError("cannot tell which Sentinel-2 product \"" * input * "\" is in; " *
                                             "give the B08.jp2 inside its .SAFE, its STAC item or its ID"))
        endswith(lowercase(input), "b08.jp2") || throw(ArgumentError(
            "\"" * input * "\" is not the B08 band the Sentinel-2 chain correlates"))
        return SceneSource(input, String(something(m.captures[1])), input, "", "")
    end
    is_stac_item(input) && return _landsat_from_item(input)
    occursin(LANDSAT_ID, input) && return _landsat_from_item(LANDSAT_STAC_ITEMS * input)
    base = basename(input)
    m = match(LANDSAT_BAND_FILE, base)
    m === nothing && throw(ArgumentError("\"" * input * "\" is not a Landsat Collection 2 band file " *
                                         "(<scene>_B<n>.TIF), a STAC item or a scene ID"))
    name = String(something(m.captures[1]))
    band = "B" * String(something(m.captures[2]))
    want = landsat_band_file(_landsat_platform(name))
    band == want || throw(ArgumentError("\"" * input * "\" is band " * band * "; the chain correlates " *
                                        want * " for " * _landsat_platform(name)))
    stem = input[1:(end - length(base))] * name
    return SceneSource(input, name, input, stem * "_MTL.json", stem * "_ANG.txt")
end

function _landsat_from_item(input::String)
    item = read_stac_item(input)
    name = item.id
    occursin(LANDSAT_SCENE, name) || throw(ArgumentError(input * " is not a Landsat scene's STAC item"))
    key = scene_band(_landsat_platform(name)) === :green ? "green" : "pan"
    return SceneSource(input, name, asset_href(item, key), asset_href(item, "MTL.json"),
                       asset_href(item, "ANG.txt"))
end

function _s2_from_item(input::String)
    item = read_stac_item(input)
    uri = item.properties.product_uri
    isempty(uri) && throw(ArgumentError(input * " carries no s2:product_uri"))
    name = endswith(uri, ".SAFE") ? uri[1:(end - 5)] : uri
    return SceneSource(input, String(name), s2_band_url(String(name)), "", "")
end

const S2_BUCKET = "sentinel-s2-l1c"
# The bucket is public: objects are read over HTTPS with unsigned requests.
const S2_BUCKET_URL = "https://sentinel-s2-l1c.s3.eu-central-1.amazonaws.com/"

# `S2B_MSIL1C_20200612T150759_N0209_R025_T22WEB_20200612T184700`: the tile is field 6, `T<zone><band><square>`.
function _s2_tile(name::String)
    tile = String(split(name, '_')[6])
    (startswith(tile, "T") && length(tile) == 6) || throw(ArgumentError(
        "cannot read an MGRS tile out of \"" * name * "\"; expected field 6 like T22WEB"))
    return tile
end

# The fields of the bucket's `productInfo.json` and `tileInfo.json` the chain reads.
mutable struct S2Datastrip
    id::String
    path::String
    S2Datastrip() = new("", "")
end
mutable struct S2TileInfo
    path::String
    productName::String
    datastrip::S2Datastrip
    S2TileInfo() = new("", "", S2Datastrip())
end
mutable struct S2ProductInfo
    name::String
    tiles::Vector{S2TileInfo}
    S2ProductInfo() = new("", S2TileInfo[])
end
StructTypes.StructType(::Type{S2Datastrip}) = StructTypes.Mutable()
StructTypes.StructType(::Type{S2TileInfo}) = StructTypes.Mutable()
StructTypes.StructType(::Type{S2ProductInfo}) = StructTypes.Mutable()

# The bucket's object `key`, or `nothing` when there is none.
function _s2_object(key::String)
    status, body = ChunkManifests.http_request(S2_BUCKET_URL * key)
    status == 404 && return nothing
    status == 200 || throw(ErrorException("HTTP " * string(status) * " fetching s3://" * S2_BUCKET * "/" * key))
    return String(body)
end

"""
    s2_location(name) -> (tile, datastrip)

The keys in `s3://sentinel-s2-l1c` of Sentinel-2 product `name`'s tile directory, which holds its bands
and tile metadata (`metadata.xml`, the `.SAFE`'s `MTD_TL.xml`), and of its datastrip directory, which
holds the datastrip metadata (`metadata.xml`, the `.SAFE`'s `MTD_DS.xml`). Both trees are dated by the
product's sensing day, month and day unpadded: `tiles/<zone>/<band>/<square>/<y>/<m>/<d>/<n>/`, where
`<n>` counts the products of that tile and day, and `products/<y>/<m>/<d>/<product>/`.

The product's `productInfo.json` names its tile directory. A product from before December 2016 is filed
under its original name (`S2A_OPER_PRD_MSIL1C_PDMC_…`), not `name`; its tile is then the one whose
`tileInfo.json` names a product of the same spacecraft, relative orbit and processing baseline whose
production time, or whose datastrip's sensing start, is `name`'s last field. A tile directory several
products share with one datastrip holds the same granule, and the first is taken. A product the bucket
does not hold is an error.
"""
function s2_location(name::String)
    date, stem = _s2_keys(name)
    info = _s2_object("products/" * date * "/" * name * "/productInfo.json")
    info === nothing || return s2_product_tile(name, info)
    hits = S2TileInfo[]
    n = 0
    while true
        doc = _s2_object(stem * string(n) * "/tileInfo.json")
        doc === nothing && break
        x = JSON3.read(doc, S2TileInfo)
        s2_original_product(x, name) && !any(h -> h.datastrip.id == x.datastrip.id, hits) && push!(hits, x)
        n += 1
    end
    isempty(hits) && throw(ArgumentError(
        "Sentinel-2 product " * name * " is not in s3://" * S2_BUCKET * ": there is no products/" * date * "/" *
        name * "/productInfo.json, and none of the " * string(n) * " directories under " * stem *
        " holds it"))
    length(hits) == 1 || throw(ArgumentError(
        "Sentinel-2 product " * name * " matches " * string(length(hits)) * " datastrips under s3://" *
        S2_BUCKET * "/" * stem))
    return (tile = hits[1].path, datastrip = hits[1].datastrip.path)
end

# The bucket's date path of product `name`, `<y>/<m>/<d>`, and the prefix of its tile's directories that day.
function _s2_keys(name::String)
    tile = _s2_tile(name)
    t = split(name, '_')[3]
    date = string(parse(Int, t[1:4])) * "/" * string(parse(Int, t[5:6])) * "/" * string(parse(Int, t[7:8]))
    return date, "tiles/" * string(parse(Int, tile[2:3])) * "/" * tile[4:4] * "/" * tile[5:6] * "/" * date * "/"
end

# The tile and datastrip directories of product `name` that its `productInfo.json`, `info`, lists.
function s2_product_tile(name::String, info::AbstractString)
    stem = _s2_keys(name)[2]
    tiles = filter(x -> startswith(x.path, stem), JSON3.read(info, S2ProductInfo).tiles)
    length(tiles) == 1 || throw(ArgumentError(
        "expected one " * _s2_tile(name) * " tile in " * name * "'s productInfo.json, found " *
        string(length(tiles))))
    return (tile = tiles[1].path, datastrip = tiles[1].datastrip.path)
end

# Whether tile directory `x` (its `tileInfo.json`) holds product `name`, filed under its original name
# `S2A_OPER_PRD_MSIL1C_PDMC_<production>_R<orbit>_V<start>_<stop>`, from datastrip
# `S2A_OPER_MSI_L1C_DS_<centre>__<production>_S<sensing start>_N<baseline>`.
s2_original_product(x::AbstractString, name::String) = s2_original_product(JSON3.read(x, S2TileInfo), name)

function s2_original_product(x::S2TileInfo, name::String)
    f = split(name, '_')
    p = match(r"^(S2[A-D])_OPER_PRD_MSIL1C_PDMC_(\d{8}T\d{6})_(R\d{3})_V", x.productName)
    d = match(r"_S(\d{8}T\d{6})_N(\d{2})\.(\d{2})$", x.datastrip.id)
    (p === nothing || d === nothing) && return false
    _capture(p, 1) == f[1] && _capture(p, 3) == f[5] && "N" * _capture(d, 2) * _capture(d, 3) == f[4] ||
        return false
    return _capture(p, 2) == f[7] || _capture(d, 1) == f[7]
end

_capture(m::RegexMatch, i::Int) = something(m.captures[i])

"""
    s2_band_url(name) -> String

The `B08.jp2` of Sentinel-2 product `name` in the public `s3://sentinel-s2-l1c` bucket, as an HTTPS URL
([`s2_location`](@ref)).
"""
s2_band_url(name::String) = S2_BUCKET_URL * s2_location(name).tile * "/B08.jp2"
