# What a production run needs around the chain: the sensor read from the inputs when it is not named,
# Earthdata Login credentials from the environment, and the product uploaded to S3.

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
    upload_product(path, bucket, prefix) -> String

Copy the product at `path` to `s3://<bucket>/<prefix>/<file name>` (no prefix when `prefix` is empty) as
`application/x-netcdf`, with the credentials of the AWS credential chain, and return that URI. The
object is tagged `file_type=product`, as HyP3 tags its products.
"""
function upload_product(path::String, bucket::String, prefix::String)
    key = isempty(prefix) ? basename(path) : rstrip(prefix, '/') * "/" * basename(path)
    uri = "s3://" * bucket * "/" * key
    target = "/vsis3/" * bucket * "/" * key
    vsi_copy_file(path, target, ["Content-Type=application/x-netcdf"]) ||
        throw(_gdal_error("could not upload " * path * " to " * uri))
    vsi_set_tags(target, ["file_type=product"]) || throw(_gdal_error("uploaded " * uri * " but could not tag it"))
    return uri
end
