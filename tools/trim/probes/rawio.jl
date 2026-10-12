# Arrays and scalars on disk without Serialization, so a trimmed program can read what a JIT session
# wrote. A directory holds one `<name>.bin` per array (raw column-major bytes) and an `index.txt` with
# one line per array, `<name> <eltype> <dim1> <dim2> ...`, plus a `meta.txt` of `<key> <value...>` lines.

# Explicit `open`/`close` rather than `open(f, path, mode)` blocks, and element type names from
# `_tname` rather than the printed type, so a trimmed program can write these too.
function _append(path::String, line::String)
    io = open(path, "a")
    try
        write(io, line)
        write(io, "\n")
    finally
        close(io)
    end
    return nothing
end

function writeraw(dir::String, name::String, A::AbstractArray)
    B = collect(A)
    io = open(joinpath(dir, name * ".bin"), "w")
    try
        write(io, B)
    finally
        close(io)
    end
    _append(joinpath(dir, "index.txt"), name * " " * _tname(eltype(B)) * " " * join(map(string, size(B)), " "))
    return nothing
end

# JIT-only: the values are printed with `string`.
writemeta(dir::AbstractString, key::AbstractString, values...) =
    _append(joinpath(dir, "meta.txt"), key * " " * join(map(string, values), " "))

# The whitespace-separated fields of the line in `file` whose first field is `key`.
function _fields(file::String, key::String)
    for line in eachline(file)
        parts = split(line)
        !isempty(parts) && parts[1] == key && return String[String(p) for p in parts[2:end]]
    end
    throw(ArgumentError(file * " has no entry " * key))
end

meta_strings(dir::String, key::String) = _fields(joinpath(dir, "meta.txt"), key)
meta_string(dir::String, key::String) = meta_strings(dir, key)[1]
meta_float(dir::String, key::String) = parse(Float64, meta_string(dir, key))
meta_int(dir::String, key::String) = parse(Int, meta_string(dir, key))
meta_floats(dir::String, key::String) = Float64[parse(Float64, s) for s in meta_strings(dir, key)]
meta_ints(dir::String, key::String) = Int[parse(Int, s) for s in meta_strings(dir, key)]

# Element type names as `index.txt` spells them; `string(::Type)` reaches the type printer, which a
# trimmed program does not have.
_tname(::Type{Float32}) = "Float32"
_tname(::Type{Float64}) = "Float64"
_tname(::Type{Int32}) = "Int32"
_tname(::Type{Int64}) = "Int64"
_tname(::Type{Bool}) = "Bool"
_tname(::Type{UInt16}) = "UInt16"

function readraw(dir::String, name::String, ::Type{T}) where {T}
    f = _fields(joinpath(dir, "index.txt"), name)
    f[1] == _tname(T) || throw(ArgumentError(name * " holds " * f[1] * ", not " * _tname(T)))
    length(f) == 3 || throw(ArgumentError(name * " is not a matrix"))
    A = Matrix{T}(undef, parse(Int, f[2]), parse(Int, f[3]))
    read!(joinpath(dir, name * ".bin"), A)
    return A
end
