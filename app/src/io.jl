# Reading objects the chain names by URL or path: the transports, whole small objects (metadata, STAC
# items, an S2 manifest), and the typed forms of the JSON among them.

import JSON3, StructTypes

# The region of each bucket the scenes come from. S3 answers a request to any other region's endpoint
# with a redirect the transport does not follow.
const BUCKET_REGION = Dict("usgs-landsat" => "us-west-2", "its-live-data" => "us-west-2",
                           "sentinel-s2-l1c" => "eu-central-1")

# One S3 transport per bucket, kept for the run so its connections are reused.
const S3_TRANSPORTS = Dict{String,ChunkManifests.S3Transport}()
const S3_TRANSPORTS_LOCK = ReentrantLock()

function s3_transport(bucket::String)
    lock(S3_TRANSPORTS_LOCK)
    try
        t = get(S3_TRANSPORTS, bucket, nothing)
        t === nothing || return t
        region = get(BUCKET_REGION, bucket, "")
        isempty(region) && throw(ArgumentError("no region is recorded for s3://" * bucket))
        # Requester-pays: the account the credential chain names pays for the reads.
        t = ChunkManifests.S3Transport(bucket; region, requesterpays = true)
        S3_TRANSPORTS[bucket] = t
        return t
    finally
        unlock(S3_TRANSPORTS_LOCK)
    end
end

is_s3(uri::AbstractString) = startswith(uri, "s3://")
is_http(uri::AbstractString) = startswith(uri, "https://") || startswith(uri, "http://")

s3_bucket(uri::AbstractString) = String(first(split(uri[6:end], '/'; limit = 2)))

"""
    transport_for(uri) -> ChunkManifests.LeafTransport

The transport an `s3://` URI, an `http(s)://` URL or a local path is read through.
"""
function transport_for(uri::String)::ChunkManifests.LeafTransport
    is_s3(uri) && return s3_transport(s3_bucket(uri))
    is_http(uri) && return ChunkManifests.HTTPTransport()
    return ChunkManifests.LocalTransport()
end

"""
    gdal_path(uri) -> String

`uri` as GDAL opens it: `/vsis3/` for `s3://`, `/vsicurl/` for `http(s)://`, and a local path made
absolute.
"""
function gdal_path(uri::String)
    is_s3(uri) && return "/vsis3/" * uri[6:end]
    is_http(uri) && return "/vsicurl/" * uri
    return abspath(uri)
end

"""
    read_object(uri) -> Vector{UInt8}

The whole object at `uri`.
"""
function read_object(uri::String)
    is_s3(uri) || is_http(uri) || return read(uri)
    # One request without `Range`: these are small, and a STAC API answers a ranged request with the
    # whole generated document.
    return ChunkManifests.fetchall(transport_for(uri), uri)
end

read_text(uri::String) = String(read_object(uri))

# ---------------------------------------------------------------------------
# STAC items
# ---------------------------------------------------------------------------
#
# Only the fields the chain reads. `Mutable` structs: keys are matched by name in any order, and keys
# without a field are skipped.

mutable struct S3Alternate
    href::String
    S3Alternate() = new("")
end
mutable struct Alternates
    s3::S3Alternate
    Alternates() = new(S3Alternate())
end
mutable struct StacAsset
    href::String
    alternate::Alternates
    StacAsset() = new("", Alternates())
end
mutable struct StacProperties
    datetime::String
    product_uri::String          # Sentinel-2 only: `s2:product_uri`
    StacProperties() = new("", "")
end
mutable struct StacItem
    type::String
    stac_version::String
    id::String
    properties::StacProperties
    assets::Dict{String,StacAsset}
    StacItem() = new("", "", "", StacProperties(), Dict{String,StacAsset}())
end

StructTypes.StructType(::Type{S3Alternate}) = StructTypes.Mutable()
StructTypes.StructType(::Type{Alternates}) = StructTypes.Mutable()
StructTypes.StructType(::Type{StacAsset}) = StructTypes.Mutable()
StructTypes.StructType(::Type{StacProperties}) = StructTypes.Mutable()
StructTypes.StructType(::Type{StacItem}) = StructTypes.Mutable()
StructTypes.names(::Type{StacProperties}) = ((:product_uri, Symbol("s2:product_uri")),)

"""
    is_stac_item(input) -> Bool

Whether `input` names a STAC item rather than an image: a local `.json` file, a URL ending in `.json`, or
a STAC API item URL (`…/collections/<collection>/items/<id>`).
"""
is_stac_item(input::AbstractString) =
    endswith(lowercase(input), ".json") ||
    (is_http(input) && occursin("/collections/", input) && occursin("/items/", input))

function read_stac_item(uri::String)
    item = JSON3.read(read_text(uri), StacItem)
    item.type == "Feature" && !isempty(item.stac_version) || throw(ArgumentError(
        uri * " is not a STAC item: a STAC item is a GeoJSON Feature with a stac_version"))
    return item
end

# The S3 location of an item's asset where the provider gives one, else its href.
function asset_href(item::StacItem, key::String)
    a = get(item.assets, key, nothing)
    a === nothing && throw(ArgumentError("STAC item " * item.id * " has no " * key * " asset"))
    return isempty(a.alternate.s3.href) ? a.href : a.alternate.s3.href
end

# ---------------------------------------------------------------------------
# Landsat metadata
# ---------------------------------------------------------------------------

mutable struct MtlProductContents
    LANDSAT_PRODUCT_ID::String
    MtlProductContents() = new("")
end
mutable struct MtlImageAttributes
    DATE_ACQUIRED::String
    SCENE_CENTER_TIME::String
    MtlImageAttributes() = new("", "")
end
mutable struct MtlFile
    PRODUCT_CONTENTS::MtlProductContents
    IMAGE_ATTRIBUTES::MtlImageAttributes
    MtlFile() = new(MtlProductContents(), MtlImageAttributes())
end
mutable struct MtlDocument
    LANDSAT_METADATA_FILE::MtlFile
    MtlDocument() = new(MtlFile())
end

StructTypes.StructType(::Type{MtlProductContents}) = StructTypes.Mutable()
StructTypes.StructType(::Type{MtlImageAttributes}) = StructTypes.Mutable()
StructTypes.StructType(::Type{MtlFile}) = StructTypes.Mutable()
StructTypes.StructType(::Type{MtlDocument}) = StructTypes.Mutable()

"""
    landsat_acquisition(mtl_uri) -> (product_id, acquisition_time)

The product id and the acquisition instant from a Landsat `_MTL.json`: `DATE_ACQUIRED` and
`SCENE_CENTER_TIME` (`15:00:24.5318140Z`), truncated to the millisecond a `DateTime` holds — the same
instant, at the same precision, as the STAC item's `properties.datetime`.
"""
function landsat_acquisition(mtl_uri::String)
    doc = JSON3.read(read_text(mtl_uri), MtlDocument)
    m = doc.LANDSAT_METADATA_FILE
    id = m.PRODUCT_CONTENTS.LANDSAT_PRODUCT_ID
    isempty(id) && throw(ArgumentError(mtl_uri * " carries no LANDSAT_PRODUCT_ID"))
    return id, parse_utc(m.IMAGE_ATTRIBUTES.DATE_ACQUIRED * "T" * m.IMAGE_ATTRIBUTES.SCENE_CENTER_TIME)
end

"""
    parse_utc(s) -> DateTime

`YYYY-MM-DDTHH:MM:SS[.fff…][Z]`, truncated to the millisecond. Parsed field by field: a `DateFormat`
parse builds its error messages at run time, which a trimmed program cannot.
"""
function parse_utc(s::AbstractString)
    t = endswith(s, "Z") ? s[1:(end - 1)] : s
    date, time = split(t, 'T'; limit = 2)
    ymd = split(date, '-')
    hms = split(time, ':')
    (length(ymd) == 3 && length(hms) == 3) || throw(ArgumentError("not a UTC timestamp: " * String(s)))
    sec = split(hms[3], '.'; limit = 2)
    ms = 0
    if length(sec) == 2
        f = sec[2]
        ms = length(f) >= 3 ? parse(Int, f[1:3]) : parse(Int, f) * 10^(3 - length(f))
    end
    return DateTime(parse(Int, ymd[1]), parse(Int, ymd[2]), parse(Int, ymd[3]), parse(Int, hms[1]),
                    parse(Int, hms[2]), parse(Int, sec[1]), ms)
end

# `YYYYMMDD`, as a scene name or an object key carries a date: field by field, as `Dates.format` is not
# available to a trimmed program.
yyyymmdd(t::DateTime) = string(Dates.year(t); pad = 4) * string(Dates.month(t); pad = 2) * string(Dates.day(t); pad = 2)

function parse_yyyymmdd(s::AbstractString)
    length(s) == 8 || throw(ArgumentError("not a YYYYMMDD date: " * String(s)))
    return Date(parse(Int, s[1:4]), parse(Int, s[5:6]), parse(Int, s[7:8]))
end
