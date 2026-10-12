# Trim probe: Zarr's multi-chunk read path — the `@async` producer, the `Channel` of compressed chunks
# and the `Threads.@spawn` decoders in `readblock!` — on an array whose type is fixed in code rather
# than parsed from its `.zarray` JSON.
#
#   zarr <directory-store> <i0> <i1> <j0> <j1>
#
# The store holds a 2-D Int16 zlib-compressed Zarr v2 array at its root, `.`-separated chunk keys.
# Prints the sum of `a[i0:i1, j0:j1]`.

import Zarr
using Zarr: ZArray, MetadataV2, ZlibCompressor, DirectoryStore, ChunkKeyEncoding

say(s::String) = (write(Core.stdout, s, "\n"); nothing)

function message(err)
    m = err isa ErrorException ? err.msg : nothing
    return m isa String ? m : string(typeof(err))
end

const Meta = MetadataV2{Int16,2,ZlibCompressor,Nothing}

function open_typed(dir::String, shape::NTuple{2,Int}, chunks::NTuple{2,Int})
    m = Meta(2, "array", shape, chunks, "<i2", ZlibCompressor(), Int16(0), 'C', nothing, ChunkKeyEncoding('.', false))
    return ZArray{Int16,2,DirectoryStore,Meta}(m, DirectoryStore(dir), "", Dict{String,Any}(), false)
end

function (@main)(args::Vector{String})::Cint
    try
        z = open_typed(args[1], (1000, 800), (128, 128))
        i0, i1, j0, j1 = (parse(Int, a) for a in args[2:5])
        w = z[i0:i1, j0:j1]
        say(string("window ", size(w, 1), "x", size(w, 2), " sum=", sum(Int, w), " threads=", Threads.nthreads()))
        return Cint(0)
    catch err
        say(string("error: ", message(err)))
        return Cint(1)
    end
end
