# A lazy image's values, computed once per tile and kept on disk.
#
# A blocked correlation sweeps every block once per pass — a coarse and a fine sweep per chip-size level —
# so a whole run reads the scene several times over. Where the values are derived rather than read, the
# multiplier lands on the derivation: a coregistered secondary is resampled four times over. A memory cache
# does not help, because the access order is a cyclic scan; so each tile is persisted the first time it is
# derived, and every later pass reads it back.
#
# The tiles live in a scratch file accessed by positioned reads and writes (`io = :pread`) or through a
# memory mapping (`:mmap`). Neither shares a file position, so tasks filling or reading different tiles
# never contend; `:pread` keeps the pages in the kernel's file cache without counting them toward the
# process's resident size, at the cost of one copy per read.
#
# **Two types, so the parent's type stops at the fill.** `TileCache{P}` derives tiles from its parent;
# [`fill_tiles!`](@ref) derives all of them, and what is left — every later read — is a [`TileFile`](@ref),
# whose type does not depend on how the tiles were made. Everything downstream of the fill is therefore
# compiled once, whatever chain of filters and warps produced the image.

"""
    TileFile <: AbstractMatrix{Float32}

The scratch file a [`TileCache`](@ref) keeps its tiles in, readable as the image once every tile it is
asked for is present. Deleted by `close`. The image starts `origin` rows and columns into the cached grid,
so a cache can align its tiles with a source's own chunks while the image keeps its own extent.
"""
struct TileFile <: AbstractMatrix{Float32}
    tile::Int
    # The cached grid's size, `origin` included.
    dims::Tuple{Int,Int}
    origin::Tuple{Int,Int}
    ntiles::Tuple{Int,Int}
    io::IOStream
    path::String
    # `:mmap` or `:pread`.
    mode::Symbol
    # With `:mmap`, the whole file, mapped once; every tile is a disjoint `Float32` range of it. Empty with
    # `:pread`.
    mapped::Vector{UInt8}
    # One lock per tile: two tasks must not derive the same absent tile at once, and the lock's release is
    # what makes a later reader see the write.
    locks::Vector{ReentrantLock}
    present::Vector{Bool}
end

"""
    TileCache(parent; tile = 512, dir, io = :pread, origin = (0, 0)) <: AbstractMatrix{Float32}

`parent`'s values, computed a tile at a time into a [`TileFile`](@ref) under `dir`. Reads the same values
`parent` would return; nothing about the tiling reaches the result. The filled file is `parent` less its
first `origin` rows and columns.
"""
struct TileCache{P<:AbstractMatrix} <: AbstractMatrix{Float32}
    parent::P
    file::TileFile
end

function TileCache(parent::AbstractMatrix; tile::Int = 512, dir::String, io::Symbol = :pread,
                   origin::Tuple{Int,Int} = (0, 0))
    tile > 0 || throw(ArgumentError("`tile` must be positive"))
    (0 <= origin[1] < size(parent, 1) && 0 <= origin[2] < size(parent, 2)) ||
        throw(ArgumentError("a tile cache's origin must lie inside its parent"))
    AutoRIFT._check_filter_cache_io(io)
    dims = size(parent)
    nt = AutoRIFT._tile_grid(dims, tile)
    nbytes = prod(nt) * tile * tile * sizeof(Float32)
    mkpath(dir)
    # A file of its own, created atomically: both images of a pair are cached in one directory and are
    # usually the same size, so a name derived from the dimensions would give them one file.
    path, stream = mktemp(dir; cleanup = false)
    truncate(stream, nbytes)
    mapped = io === :mmap ? Mmap.mmap(stream, Vector{UInt8}, nbytes) : UInt8[]
    file = TileFile(tile, dims, origin, nt, stream, path, io, mapped, [ReentrantLock() for _ in 1:prod(nt)],
                    zeros(Bool, prod(nt)))
    return TileCache(parent, file)
end

Base.size(f::TileFile) = (f.dims[1] - f.origin[1], f.dims[2] - f.origin[2])
Base.size(c::TileCache) = c.file.dims

# `Mmap` attaches its unmapping finalizer to the array's underlying `Memory` from Julia 1.11 on, so
# `finalize(arr)` alone does not unmap.
function _finalize_mapping!(mapped::Vector{UInt8})
    finalize(mapped)
    finalize(mapped.ref.mem)
    return nothing
end

"""
    close(f::TileFile)

Close the scratch file and delete it.
"""
function Base.close(f::TileFile)
    close(f.io)
    f.mode === :mmap && _finalize_mapping!(f.mapped)
    isfile(f.path) && rm(f.path)
    return nothing
end

_tile_bytes(f::TileFile) = f.tile * f.tile * sizeof(Float32)
_tile_offset(f::TileFile, ti::Int, tj::Int) = AutoRIFT._tile_grid_offset(f.ntiles, _tile_bytes(f), ti, tj)
_tile_span(f::TileFile, r::AbstractUnitRange) = AutoRIFT._tile_grid_span(f.tile, r)
_tile_extent(f::TileFile, ti::Int, tj::Int) = AutoRIFT._tile_grid_extent(f.dims, f.tile, ti, tj)
_tile_index(f::TileFile, ti::Int, tj::Int) = (tj - 1) * f.ntiles[1] + ti

# The `tile`-square plane backing one tile of a mapped file, as a view into the mapping.
function _tile_array(f::TileFile, ti::Int, tj::Int)
    off = _tile_offset(f, ti, tj)
    n = f.tile * f.tile
    return reshape(reinterpret(Float32, view(f.mapped, (off + 1):(off + 4n))), f.tile, f.tile)
end

# Make tile `(ti, tj)` present, deriving it from `c.parent` if it is not, under that tile's own lock.
function _ensure_tile!(c::TileCache, ti::Int, tj::Int)
    f = c.file
    k = _tile_index(f, ti, tj)
    lock(f.locks[k])
    try
        f.present[k] && return nothing
        rows, cols = _tile_extent(f, ti, tj)
        if f.mode === :mmap
            buf = _tile_array(f, ti, tj)
            fill!(buf, 0.0f0)
            copyto!(view(buf, 1:length(rows), 1:length(cols)), c.parent[rows, cols])
        else
            bytes = Vector{UInt8}(undef, _tile_bytes(f))
            buf = reshape(reinterpret(Float32, bytes), f.tile, f.tile)
            fill!(buf, 0.0f0)
            copyto!(view(buf, 1:length(rows), 1:length(cols)), c.parent[rows, cols])
            AutoRIFT._pwrite(f.io, bytes, _tile_offset(f, ti, tj))
        end
        f.present[k] = true
    finally
        unlock(f.locks[k])
    end
    return nothing
end

# A task's result, or the exception it failed with rather than the `TaskFailedException` wrapping it, so its
# message reaches the command line.
function _await(t::Task)
    try
        return fetch(t)
    catch e
        e isa TaskFailedException && throw(e.task.exception)
        rethrow()
    end
end

# Tiles kept in flight while a remote granule's chunks are read straight into them: a tile then costs
# little memory beyond its own plane, and with only one per thread the link idles while threads decode.
remote_tiles_in_flight() = 4 * Threads.nthreads()

"""
    fill_tiles!(c::TileCache; inflight = Threads.nthreads()) -> TileFile

Derive every tile of `c`, `inflight` at a time — a tile's derivation is a remote read and possibly a
filter and a warp — and return the file, which from here on is the image. Tiles are taken in order, so the
ones in flight are neighbours. A tile waiting on the network leaves its thread to the others, so a source
read remotely keeps more tiles in flight than there are threads; each costs its own buffers.
"""
function fill_tiles!(c::TileCache; inflight::Int = Threads.nthreads())
    nt = c.file.ntiles
    n = prod(nt)
    next = Threads.Atomic{Int}(1)
    tasks = [Threads.@spawn begin
                 while true
                     k = Threads.atomic_add!(next, 1)
                     k > n && break
                     _ensure_tile!(c, mod1(k, nt[1]), cld(k, nt[1]))
                 end
             end for _ in 1:min(n, max(1, inflight))]
    foreach(_await, tasks)
    return c.file
end

function Base.getindex(f::TileFile, irows::AbstractUnitRange, icols::AbstractUnitRange)
    checkbounds(f, irows, icols)
    # In the cached grid's indices from here on.
    rows, cols = irows .+ f.origin[1], icols .+ f.origin[2]
    out = Matrix{Float32}(undef, length(rows), length(cols))
    for tj in _tile_span(f, cols), ti in _tile_span(f, rows)
        f.present[_tile_index(f, ti, tj)] || throw(ArgumentError("a tile was read before it was derived"))
        trows, tcols = _tile_extent(f, ti, tj)
        wr = intersect(rows, trows)
        wc = intersect(cols, tcols)
        (isempty(wr) || isempty(wc)) && continue
        sr, sc = wr .- (first(trows) - 1), wc .- (first(tcols) - 1)
        if f.mode === :mmap
            copyto!(view(out, (wr .- first(rows) .+ 1), (wc .- first(cols) .+ 1)), view(_tile_array(f, ti, tj), sr, sc))
        else
            # The request's columns of the tile — contiguous in the file, so one positioned read — into a
            # buffer this task keeps.
            b = get!(() -> Matrix{Float32}(undef, f.tile, f.tile), task_local_storage(),
                     (:tilecache_read_buffer, f.tile))::Matrix{Float32}
            AutoRIFT._pread_into!(f.io, reinterpret(UInt8, vec(b)), 4 * f.tile * length(sc),
                                  _tile_offset(f, ti, tj) + 4 * f.tile * (first(sc) - 1))
            copyto!(view(out, (wr .- first(rows) .+ 1), (wc .- first(cols) .+ 1)), view(b, sr, 1:length(sc)))
        end
    end
    return out
end

Base.getindex(f::TileFile, i::Integer, j::Integer) = f[i:i, j:j][1, 1]
Base.getindex(f::TileFile, rows::AbstractUnitRange, j::Integer) = f[rows, j:j][:, 1]
Base.getindex(f::TileFile, i::Integer, cols::AbstractUnitRange) = f[i:i, cols][1, :]

# Reading one costs I/O, which is what it is for: a blocked run windows its reads.
AutoRIFT.ondisk(::TileFile) = true
