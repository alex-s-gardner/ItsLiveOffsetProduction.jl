# A job: the sensor, and the two scenes resolved from the inputs the caller gave — a path, an `s3://` or
# `https://` URL, or a STAC item.

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

const SENSORS = (
    SensorSpec("L4", "L4", r"^LT04_", LANDSAT_SCENE),
    SensorSpec("L5", "L5", r"^LT05_", LANDSAT_SCENE),
    SensorSpec("L7", "L7", r"^LE07_", LANDSAT_SCENE),
    SensorSpec("L8", "L8", r"^LC08_", LANDSAT_SCENE),
    SensorSpec("L9", "L9", r"^LC09_", LANDSAT_SCENE),
    SensorSpec("S2", "S2", r"^S2[ABC]_MSIL1C_", r"^S2[ABC]_MSIL1C_"),
    SensorSpec("S1", "S1-SLC", r"^S1[ABCD]_IW_SLC__", r"^S1[ABCD]_IW_SLC__"),
    SensorSpec("S1BURST", "S1-BURST", r"^S1[ABCD]_IW_SLC__", r"^S1[ABCD]_IW_SLC__"),
    SensorSpec("NL1", "NISAR-L1", r"^NISAR_L1_\w+_RSLC_", r"^NISAR_L1_\w+_RSLC_"),
    SensorSpec("NL2", "NISAR-L2", r"^NISAR_L2_\w+_GSLC_", r"^NISAR_L2_\w+_GSLC_"),
)

# The sensors the command line names that this build does not process.
const UNBUILT_SENSORS = String[]

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

One pair to process: the sensor and its platform, the two scenes, and the product name for messages.
"""
struct Job
    sensor::String
    platform::String
    reference::SceneSource
    secondary::SceneSource
    product::String
end

function Job(sensor::String, reference::String, secondary::String)
    spec = sensor_spec(sensor)
    ref = scene_source(reference, spec.platform)
    sec = scene_source(secondary, spec.platform)
    occursin(spec.reference, ref.name) || throw(ArgumentError(
        "the reference input \"" * reference * "\" is scene " * ref.name * ", which is not a " * sensor * " scene"))
    occursin(spec.secondary, sec.name) || throw(ArgumentError(
        "the secondary input \"" * secondary * "\" is scene " * sec.name * ", which cannot pair with a " *
        sensor * " scene"))
    return Job(sensor, spec.platform, ref, sec, ref.name * "_X_" * sec.name)
end

"""
    scene_source(input, platform) -> SceneSource

The scene `input` names. A Landsat input is the band file the reference correlates — `_B8.TIF` for
Landsat 7, 8 and 9, `_B2.TIF` (green) for 4 and 5 — with its `_MTL.json` and `_ANG.txt` beside it, or a
landsatlook STAC item, whose assets name all three. A Sentinel-2 input is the `B08.jp2` band inside a
`.SAFE`, or an Earth Search STAC item, which is resolved through its `s2:product_uri` to the band on the
public Google Cloud mirror. A NISAR input is the granule's `.h5`, a Sentinel-1 input the product's `.SAFE`
directory or its `.zip`. A burst pair's input is the bursts' extractor URLs, comma-separated, whose scene
is the product `burst2safe` would assemble from them, or that product's `.SAFE`.
"""
function scene_source(input::String, platform::String)
    if startswith(platform, "S1")
        if platform == "S1-BURST" && is_burst_list(input)
            return SceneSource(input, SLCDatasets.asf_safe_name(burst_list(input, nothing)), input, "", "")
        end
        base = basename(rstrip(input, '/'))
        name = endswith(base, ".SAFE") ? base[1:(end - 5)] : endswith(base, ".zip") ? base[1:(end - 4)] :
               throw(ArgumentError("\"" * input * "\" is neither a Sentinel-1 product's .SAFE nor its .zip"))
        return SceneSource(input, String(name), input, "", "")
    end
    if startswith(platform, "NISAR")
        base = basename(input)
        endswith(base, ".h5") || throw(ArgumentError("\"" * input * "\" is not a NISAR granule's .h5"))
        return SceneSource(input, base[1:(end - 3)], input, "", "")
    end
    if platform == "S2"
        is_stac_item(input) && return _s2_from_item(input)
        m = match(r"(S2[ABC]_MSIL1C_[0-9A-Z_]+?)\.SAFE", input)
        m === nothing && throw(ArgumentError("cannot tell which Sentinel-2 product \"" * input * "\" is in; " *
                                             "give the B08.jp2 inside its .SAFE, or its STAC item"))
        endswith(lowercase(input), "b08.jp2") || throw(ArgumentError(
            "\"" * input * "\" is not the B08 band the Sentinel-2 chain correlates"))
        return SceneSource(input, String(something(m.captures[1])), input, "", "")
    end
    is_stac_item(input) && return _landsat_from_item(input)
    base = basename(input)
    m = match(r"^(L[TEOC]0[45789]_L1\w\w_\d{6}_\d{8}_\d{8}_\d{2}_\w{2})_B(\d+)\.TIF$"i, base)
    m === nothing && throw(ArgumentError("\"" * input * "\" is not a Landsat Collection 2 band file " *
                                         "(<scene>_B<n>.TIF) or a STAC item"))
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

const S2_MIRROR = "https://storage.googleapis.com/gcp-public-data-sentinel-2/tiles"

"""
    s2_band_url(name) -> String

The `B08.jp2` of Sentinel-2 product `name` on the public Google Cloud mirror. The granule directory inside
a `.SAFE` is not derivable from the product name, so the product's `manifest.safe` is read for it.
"""
function s2_band_url(name::String)
    # `S2B_MSIL1C_20200612T150759_N0209_R025_T22WEB_20200612T184700`: the tile is field 6, `T<zone><band>
    # <square>`, and the mirror splits it into three directory levels.
    tile = split(name, '_')[6]
    (startswith(tile, "T") && length(tile) == 6) || throw(ArgumentError(
        "cannot read an MGRS tile out of \"" * name * "\"; expected field 6 like T22WEB"))
    safe = S2_MIRROR * "/" * tile[2:3] * "/" * tile[4:4] * "/" * tile[5:6] * "/" * name * ".SAFE"
    manifest = read_text(safe * "/manifest.safe")
    # The band under IMG_DATA only: from processing baseline 04.00 the manifest also lists QI_DATA masks
    # whose names end in B08.jp2.
    hits = unique(String[m.match for m in eachmatch(r"GRANULE/[^\"<>]*/IMG_DATA/[^\"<>]*_B08\.jp2", manifest)])
    # A product from before December 2016 holds every tile of its datatake: its own is the one named.
    length(hits) > 1 && filter!(h -> occursin(tile * "_B08", h), hits)
    length(hits) == 1 || throw(ArgumentError("expected one B08 entry in " * name * "'s manifest.safe, found " *
                                             string(length(hits))))
    return safe * "/" * hits[1]
end
