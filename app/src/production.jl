# What a production run needs around the chain: the sensor read from the inputs when it is not named,
# Earthdata Login credentials from the environment, and the product, its browse image and thumbnail
# uploaded to S3 and published to the ITS_LIVE Open Data bucket.

"""
    infer_sensor(reference) -> String

The sensor code of the pair whose reference input is `reference`: `S1BURST` for bursts — burst IDs,
extractor URLs, bursts of SLC zips, or SLC zips given as a list, which only a burst pair reads — and
otherwise the code whose reference pattern ([`SENSORS`](@ref)) the scene name the input carries
matches. Refused: an input it cannot tell the sensor of, which needs the sensor named.
"""
function infer_sensor(reference::String)
    (is_burst_list(reference) || is_zip_burst_list(reference) ||
     occursin(S1_BURST_ID, strip(first(split(reference, ',')))) ||
     (occursin(',', reference) && is_zip_list(reference))) && return "S1BURST"
    name = _input_scene_name(reference)
    for s in SENSORS
        s.code == "S1BURST" || !occursin(s.reference, name) || return s.code
    end
    throw(ArgumentError("cannot tell the sensor from the reference \"" * reference * "\"; name it first, " *
                        "one of " * join([s.code for s in SENSORS], ", ")))
end

# The scene name an input carries: the ID itself, or the file or directory a path or URL ends in, read
# up to the extension. A Sentinel-2 band is named after its tile, so its `.SAFE` directory is read instead.
function _input_scene_name(input::String)
    m = match(S2_SAFE, input)
    m === nothing || return String(something(m.captures[1]))
    return String(first(split(basename(rstrip(input, '/')), '.')))
end

# ---------------------------------------------------------------------------
# Earthdata Login credentials
# ---------------------------------------------------------------------------

"""
    with_earthdata_netrc(f)

`f()`, with `EARTHDATA_USERNAME` and `EARTHDATA_PASSWORD`, when both are set, written as the
`urs.earthdata.nasa.gov` entry of a `.netrc` in a new private directory and `NETRC` pointed at it for
`f`'s length, so every Earthdata read — a token request and each login redirect libcurl follows — uses
them. The directory is removed when `f` returns or throws. With neither set, `f()` alone. Refused: one of
the two set without the other, and `NETRC` already set, whose file `f` would otherwise not read.
"""
function with_earthdata_netrc(f::F) where {F}
    user = get(ENV, "EARTHDATA_USERNAME", "")
    pass = get(ENV, "EARTHDATA_PASSWORD", "")
    isempty(user) && isempty(pass) && return f()
    (isempty(user) || isempty(pass)) && throw(ArgumentError(
        "EARTHDATA_USERNAME and EARTHDATA_PASSWORD are set together, but only one of them is"))
    haskey(ENV, "NETRC") && throw(ArgumentError(
        "NETRC names a .netrc and EARTHDATA_USERNAME and EARTHDATA_PASSWORD name credentials; set one or the other"))
    any(isspace, user * pass) && throw(ArgumentError(
        "EARTHDATA_USERNAME or EARTHDATA_PASSWORD holds white space, which a .netrc cannot carry"))
    # `mktempdir` makes the directory readable by its owner alone.
    dir = mktempdir(; cleanup = false)
    path = joinpath(dir, ".netrc")
    write(path, "machine " * SLCDatasets.EDL_HOST * " login " * user * " password " * pass * "\n")
    ENV["NETRC"] = path
    try
        return f()
    finally
        delete!(ENV, "NETRC")
        rm(dir; recursive = true, force = true)
    end
end

# ---------------------------------------------------------------------------
# Upload
# ---------------------------------------------------------------------------

"""
    Upload

A file to upload: its path, its content type and the HyP3 `file_type` tag it carries (`""` for none).
"""
struct Upload
    path::String
    type::String
    tag::String
end

"""
    upload_file(u, bucket, prefix) -> String

Copy file `u` to `s3://<bucket>/<prefix>/<file name>` (no prefix when `prefix` is empty), with its
content type and tag, and return that URI. The credentials are the AWS credential chain's, unless
[`publish_credentials!`](@ref) set the bucket's own.
"""
function upload_file(u::Upload, bucket::String, prefix::String)
    key = isempty(prefix) ? basename(u.path) : rstrip(prefix, '/') * "/" * basename(u.path)
    uri = "s3://" * bucket * "/" * key
    target = "/vsis3/" * bucket * "/" * key
    vsi_copy_file(u.path, target, ["Content-Type=" * u.type]) ||
        throw(_gdal_error("could not upload " * u.path * " to " * uri))
    isempty(u.tag) || vsi_set_tags(target, ["file_type=" * u.tag]) ||
        throw(_gdal_error("uploaded " * uri * " but could not tag it"))
    return uri
end

# ---------------------------------------------------------------------------
# Publishing to the Open Data bucket
# ---------------------------------------------------------------------------

# The Open Data directory of each platform's products (`hyp3_autorift.utils.PLATFORM_SHORTNAME_LONGNAME_MAPPING`).
function opendata_platform(platform::String)
    startswith(platform, "L") && return "landsatOLI"
    platform == "S2" && return "sentinel2"
    startswith(platform, "S1") && return "sentinel1"
    startswith(platform, "NISAR") && return "nisar"
    throw(ArgumentError("no Open Data directory is recorded for platform " * platform))
end

"""
    point_to_region(lat, lon) -> String

The ITS_LIVE region of a product centred at `(lat, lon)`, as `hyp3_autorift.utils.point_to_region` names
it: the hemispheres and the 10° cell, truncated toward zero, `N78W124` for `(78.6, -124.2)`. The pole
falls in the 80° cell and the antimeridian in the 170° one.
"""
function point_to_region(lat::Float64, lon::Float64)
    ns = lat >= 0.0 ? "N" : "S"
    ew = lon >= 0.0 ? "E" : "W"
    rlat = min(unsafe_trunc(Int, 10 * trunc(abs(lat / 10.0))), 80)
    rlon = min(unsafe_trunc(Int, 10 * trunc(abs(lon / 10.0))), 170)
    return ns * string(rlat; pad = 2) * ew * string(rlon; pad = 3)
end

"""
    opendata_prefix(platform, product) -> String

Where the product at `product` is published in the Open Data bucket,
`velocity_image_pair/<platform directory>/v02/<region>` (`hyp3_autorift.utils.get_opendata_prefix`), its
region from the centre `img_pair_info` records.
"""
function opendata_prefix(platform::String, product::String)
    lat, lon = ItsLiveOffsetProduction.read_center(product)
    return "velocity_image_pair/" * opendata_platform(platform) * "/v02/" * point_to_region(lat, lon)
end

"""
    publish_credentials!(bucket)

`PUBLISH_ACCESS_KEY_ID` and `PUBLISH_SECRET_ACCESS_KEY` set as the credentials of `/vsis3/<bucket>/`
alone, so publishing uses them and every other read and upload the credential chain's. Refused when
either is not set.
"""
function publish_credentials!(bucket::String)
    id = get(ENV, "PUBLISH_ACCESS_KEY_ID", "")
    secret = get(ENV, "PUBLISH_SECRET_ACCESS_KEY", "")
    (isempty(id) || isempty(secret)) && throw(ArgumentError(
        "--publish-bucket needs the bucket's credentials in PUBLISH_ACCESS_KEY_ID and PUBLISH_SECRET_ACCESS_KEY"))
    prefix = "/vsis3/" * bucket * "/"
    GDAL.vsisetpathspecificoption(prefix, "AWS_ACCESS_KEY_ID", id)
    GDAL.vsisetpathspecificoption(prefix, "AWS_SECRET_ACCESS_KEY", secret)
    # Long-lived keys: a session token from the credential chain belongs to other credentials.
    GDAL.vsisetpathspecificoption(prefix, "AWS_SESSION_TOKEN", "")
    return nothing
end

"""
    PublishInfo

Where a product was published, as `hyp3-autorift` records it in `publish_info.json` for the ITS_LIVE
metadata step to find.
"""
struct PublishInfo
    bucket::String
    prefix::String
    name::String
end
StructTypes.StructType(::Type{PublishInfo}) = StructTypes.Struct()

"""
    deliver(path, platform; bucket, bucket_prefix, publish_bucket) -> Vector{String}

The product at `path` with its browse image and thumbnail (drawn beside it when either upload is asked
for): published to `publish_bucket` under the product's Open Data prefix, and uploaded to `bucket` under
`bucket_prefix` — with `publish_info.json`, which records where it was published, when there is both. The
URIs written, in that order; none when both buckets are empty. The publish credentials and the prefix are
settled before anything is drawn.
"""
function deliver(path::String, platform::String; bucket::String, bucket_prefix::String, publish_bucket::String)
    isempty(bucket) && isempty(publish_bucket) && return String[]
    publishing = !isempty(publish_bucket)
    publishing && publish_credentials!(publish_bucket)
    prefix = publishing ? opendata_prefix(platform, path) : ""
    browse, thumb = write_browse(path)
    files = [Upload(path, "application/x-netcdf", "product"), Upload(browse, "image/png", "amp-browse"),
             Upload(thumb, "image/png", "amp-thumbnail")]
    uris = String[]
    if publishing
        for f in files
            push!(uris, upload_file(f, publish_bucket, prefix))
        end
    end
    isempty(bucket) && return uris
    if publishing
        info = joinpath(dirname(abspath(path)), "publish_info.json")
        write(info, JSON3.write(PublishInfo(publish_bucket, prefix, basename(path))))
        push!(files, Upload(info, "application/json", ""))
    end
    for f in files
        push!(uris, upload_file(f, bucket, bucket_prefix))
    end
    return uris
end
