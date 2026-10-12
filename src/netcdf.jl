# A netCDF-4 writer over the netCDF C library (NetCDF_jll): the calls `write_product` makes and no
# others. Every value written has a concrete Julia type, and so a fixed netCDF type, which keeps the
# writer usable from a `--trim`med program.
#
# Dimensions are declared and variables shaped in netCDF's own order, slowest first — `(time, y, x)`.
# A Julia matrix is column-major, so an `(x, y)` Julia array is the `(y, x)` netCDF variable's memory
# layout, and data are handed over in that layout.

import NetCDF_jll

const NC_NETCDF4 = Cint(0x1000)
const NC_CLOBBER = Cint(0)
const NC_GLOBAL = Cint(-1)
const NC_UNLIMITED = Csize_t(0)
const NC_CHUNKED = Cint(0)

const NC_BYTE = Cint(1)
const NC_CHAR = Cint(2)
const NC_SHORT = Cint(3)
const NC_INT = Cint(4)
const NC_FLOAT = Cint(5)
const NC_DOUBLE = Cint(6)
const NC_UBYTE = Cint(7)
const NC_USHORT = Cint(8)
const NC_UINT = Cint(9)
const NC_INT64 = Cint(10)
const NC_UINT64 = Cint(11)
const NC_STRING = Cint(12)

_nctype(::Type{Int8}) = NC_BYTE
_nctype(::Type{Int16}) = NC_SHORT
_nctype(::Type{Int32}) = NC_INT
_nctype(::Type{Int64}) = NC_INT64
_nctype(::Type{UInt8}) = NC_UBYTE
_nctype(::Type{UInt16}) = NC_USHORT
_nctype(::Type{UInt32}) = NC_UINT
_nctype(::Type{UInt64}) = NC_UINT64
_nctype(::Type{Float32}) = NC_FLOAT
_nctype(::Type{Float64}) = NC_DOUBLE

"""
    AttrValue

The types a netCDF attribute value may have here: a string, a number of one of netCDF's numeric types,
or a vector of `UInt8`. Each is written as its own netCDF type — a `String` as `NC_CHAR` or
`NC_STRING` (see [`NcFile`](@ref)), an `Int64` as `NC_INT64`, a `Float32` as `NC_FLOAT`.
"""
const AttrValue = Union{String,Float64,Float32,Int64,Int32,Int16,Int8,UInt64,UInt32,UInt16,UInt8,
                        Vector{UInt8}}

"""
    attrvalue(x) -> AttrValue

`x` as an attribute value: a string becomes a `String`; a value of an [`AttrValue`](@ref) type is
returned as it is; anything else is an error, since it has no netCDF type to be written as.
"""
attrvalue(x::AbstractString) = String(x)
attrvalue(x::AttrValue) = x
attrvalue(x) = throw(ArgumentError(
    "an attribute value must be a string, a netCDF numeric type or a Vector{UInt8}; got a " *
    string(nameof(typeof(x)))))

"""
    attrdict(d) -> Dict{String,AttrValue}

`d`'s entries as attribute values (see [`attrvalue`](@ref)), inserted in `d`'s iteration order.
"""
function attrdict(d::AbstractDict)
    out = Dict{String,AttrValue}()
    for (k, v) in d
        out[String(k)] = attrvalue(v)
    end
    return out
end
attrdict(d::Dict{String,AttrValue}) = d

function _nccheck(status::Cint, what::String)
    status == 0 && return nothing
    msg = unsafe_string(ccall((:nc_strerror, NetCDF_jll.libnetcdf), Cstring, (Cint,), status))
    throw(ErrorException("netCDF " * what * ": " * msg))
end

"""
    NcFile

A netCDF-4 file open for writing: its id, its dimensions by name, and whether its string attributes
are written as `NC_STRING` (`nc_string = true`) or `NC_CHAR`.
"""
struct NcFile
    ncid::Cint
    dims::Dict{String,Cint}
    nc_string::Bool
end

"""
    NcVar

A variable of an [`NcFile`](@ref), or the file's global attributes (`varid == NC_GLOBAL`).
"""
struct NcVar
    file::NcFile
    varid::Cint
end

_globals(f::NcFile) = NcVar(f, NC_GLOBAL)

function nc_create(path::String, nc_string::Bool)
    id = Ref{Cint}(0)
    _nccheck(ccall((:nc_create, NetCDF_jll.libnetcdf), Cint, (Cstring, Cint, Ptr{Cint}),
                   path, NC_NETCDF4 | NC_CLOBBER, id), "nc_create " * path)
    return NcFile(id[], Dict{String,Cint}(), nc_string)
end

nc_close(f::NcFile) =
    _nccheck(ccall((:nc_close, NetCDF_jll.libnetcdf), Cint, (Cint,), f.ncid), "nc_close")

"""
    with_netcdf(f, path, nc_string)

Create `path`, call `f(file)`, and close it whether or not `f` throws.
"""
function with_netcdf(f::F, path::String, nc_string::Bool) where {F}
    file = nc_create(path, nc_string)
    try
        f(file)
    finally
        nc_close(file)
    end
    return path
end

# `len = nothing` is the unlimited dimension.
function def_dim!(f::NcFile, name::String, len::Union{Int,Nothing})
    id = Ref{Cint}(0)
    _nccheck(ccall((:nc_def_dim, NetCDF_jll.libnetcdf), Cint, (Cint, Cstring, Csize_t, Ptr{Cint}),
                   f.ncid, name, len === nothing ? NC_UNLIMITED : Csize_t(len), id), "nc_def_dim " * name)
    f.dims[name] = id[]
    return id[]
end

"""
    def_var!(file, name, T, dims; fillvalue = nothing, chunks = nothing, deflate = 0, shuffle = false) -> NcVar

A variable of element type `T` over the named dimensions, in netCDF order (slowest first). `fillvalue`
sets `_FillValue`; `nothing` leaves the library's default fill and no attribute. `chunks` are chunk
sizes in the same order as `dims`; `deflate` a zlib level, `0` for none.
"""
function def_var!(f::NcFile, name::String, ::Type{T}, dims::Union{Tuple{Vararg{String}},Vector{String}};
                  fillvalue::Union{T,Nothing} = nothing,
                  chunks::Union{Nothing,Tuple{Vararg{Int}},Vector{Int}} = nothing,
                  deflate::Int = 0, shuffle::Bool = false) where {T}
    dimids = Cint[f.dims[d] for d in dims]
    id = Ref{Cint}(0)
    _nccheck(ccall((:nc_def_var, NetCDF_jll.libnetcdf), Cint, (Cint, Cstring, Cint, Cint, Ptr{Cint}, Ptr{Cint}),
                   f.ncid, name, _nctype(T), length(dimids), dimids, id), "nc_def_var " * name)
    varid = id[]
    if chunks !== nothing
        sizes = Csize_t[c for c in chunks]
        _nccheck(ccall((:nc_def_var_chunking, NetCDF_jll.libnetcdf), Cint, (Cint, Cint, Cint, Ptr{Csize_t}),
                       f.ncid, varid, NC_CHUNKED, sizes), "nc_def_var_chunking " * name)
    end
    if deflate > 0 || shuffle
        _nccheck(ccall((:nc_def_var_deflate, NetCDF_jll.libnetcdf), Cint, (Cint, Cint, Cint, Cint, Cint),
                       f.ncid, varid, Cint(shuffle), Cint(deflate > 0), Cint(deflate)),
                 "nc_def_var_deflate " * name)
    end
    if fillvalue !== nothing
        _nccheck(ccall((:nc_def_var_fill, NetCDF_jll.libnetcdf), Cint, (Cint, Cint, Cint, Ref{T}),
                       f.ncid, varid, Cint(0), fillvalue), "nc_def_var_fill " * name)
    end
    return NcVar(f, varid)
end

# Attributes, one method per value type.
function put_att!(v::NcVar, key::String, val::String)
    if v.file.nc_string
        p = Cstring[Base.unsafe_convert(Cstring, val)]
        GC.@preserve val p _nccheck(
            ccall((:nc_put_att_string, NetCDF_jll.libnetcdf), Cint,
                  (Cint, Cint, Cstring, Csize_t, Ptr{Cstring}), v.file.ncid, v.varid, key, 1, p),
            "nc_put_att_string " * key)
    else
        _nccheck(ccall((:nc_put_att_text, NetCDF_jll.libnetcdf), Cint, (Cint, Cint, Cstring, Csize_t, Cstring),
                       v.file.ncid, v.varid, key, sizeof(val), val), "nc_put_att_text " * key)
    end
    return v
end

function put_att!(v::NcVar, key::String, val::T) where {T<:Union{Int8,Int16,Int32,Int64,UInt8,UInt16,UInt32,UInt64,Float32,Float64}}
    _nccheck(ccall((:nc_put_att, NetCDF_jll.libnetcdf), Cint, (Cint, Cint, Cstring, Cint, Csize_t, Ref{T}),
                   v.file.ncid, v.varid, key, _nctype(T), 1, val), "nc_put_att " * key)
    return v
end

function put_att!(v::NcVar, key::String, val::Vector{UInt8})
    _nccheck(ccall((:nc_put_att, NetCDF_jll.libnetcdf), Cint, (Cint, Cint, Cstring, Cint, Csize_t, Ptr{UInt8}),
                   v.file.ncid, v.varid, key, NC_UBYTE, length(val), val), "nc_put_att " * key)
    return v
end

put_att!(v::NcVar, key::String, val::SubString{String}) = put_att!(v, key, String(val))

# A value of the `AttrValue` union, dispatched by an explicit test of each member so that every call
# is to a concrete method.
function put_attvalue!(v::NcVar, key::String, val::AttrValue)
    if val isa String
        put_att!(v, key, val)
    elseif val isa Float64
        put_att!(v, key, val)
    elseif val isa Float32
        put_att!(v, key, val)
    elseif val isa Int64
        put_att!(v, key, val)
    elseif val isa Int32
        put_att!(v, key, val)
    elseif val isa Int16
        put_att!(v, key, val)
    elseif val isa Int8
        put_att!(v, key, val)
    elseif val isa UInt64
        put_att!(v, key, val)
    elseif val isa UInt32
        put_att!(v, key, val)
    elseif val isa UInt16
        put_att!(v, key, val)
    elseif val isa UInt8
        put_att!(v, key, val)
    elseif val isa Vector{UInt8}
        put_att!(v, key, val)
    end
    return v
end

"""
    put_var!(v, data; start, count)

Write `data`, whose memory layout is the variable's netCDF layout, into the hyperslab at `start` (zero
based) of extent `count`, both in netCDF order.
"""
function put_var!(v::NcVar, data::Array{T}; start::Tuple{Vararg{Int}}, count::Tuple{Vararg{Int}}) where {T}
    prod(count) == length(data) || throw(DimensionMismatch(
        "a hyperslab of " * string(prod(count)) * " values cannot take " * string(length(data))))
    s = Csize_t[i for i in start]
    c = Csize_t[i for i in count]
    _nccheck(ccall((:nc_put_vara, NetCDF_jll.libnetcdf), Cint, (Cint, Cint, Ptr{Csize_t}, Ptr{Csize_t}, Ptr{T}),
                   v.file.ncid, v.varid, s, c, data), "nc_put_vara")
    return v
end
