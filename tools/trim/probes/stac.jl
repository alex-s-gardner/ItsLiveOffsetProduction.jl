# Trim probe: typed parsing of real STAC items (a landsatlook Landsat item and an Earth Search
# Sentinel-2 item) with JSON3 + StructTypes, keeping only the fields the chain uses.
#
#   stac <landsat-item.json> <sentinel2-item.json>

import JSON3, StructTypes

say(s::String) = (write(Core.stdout, s, "\n"); nothing)

function message(err)
    m = err isa ErrorException ? err.msg : nothing
    return m isa String ? m : string(typeof(err))
end

mutable struct S3Alternate
    href::String
    S3Alternate() = new("")
end
mutable struct Alternates
    s3::S3Alternate
    Alternates() = new(S3Alternate())
end
mutable struct Asset
    href::String
    alternate::Alternates
    Asset() = new("", Alternates())
end
mutable struct Properties
    datetime::String
    product_uri::String          # Sentinel-2 only: `s2:product_uri`
    Properties() = new("", "")
end
mutable struct Item
    type::String
    stac_version::String
    id::String
    properties::Properties
    assets::Dict{String,Asset}
    Item() = new("", "", "", Properties(), Dict{String,Asset}())
end

# Mutable: keys are matched by name in any order, and keys without a field are skipped.
StructTypes.StructType(::Type{S3Alternate}) = StructTypes.Mutable()
StructTypes.StructType(::Type{Alternates}) = StructTypes.Mutable()
StructTypes.StructType(::Type{Asset}) = StructTypes.Mutable()
StructTypes.StructType(::Type{Properties}) = StructTypes.Mutable()
StructTypes.StructType(::Type{Item}) = StructTypes.Mutable()
StructTypes.names(::Type{Properties}) = ((:product_uri, Symbol("s2:product_uri")),)

read_item(path::String) = JSON3.read(read(path, String), Item)

function (@main)(args::Vector{String})::Cint
    try
        ls = read_item(args[1])
        say(string("landsat id=", ls.id, " type=", ls.type, " datetime=", ls.properties.datetime))
        say(string("  B8 s3=", ls.assets["pan"].alternate.s3.href))
        say(string("  MTL.json s3=", ls.assets["MTL.json"].alternate.s3.href))
        s2 = read_item(args[2])
        say(string("s2 id=", s2.id, " product_uri=", s2.properties.product_uri, " datetime=", s2.properties.datetime))
        say(string("  B08 href=", s2.assets["nir"].href))
        return Cint(0)
    catch err
        say(string("error: ", message(err)))
        return Cint(1)
    end
end
