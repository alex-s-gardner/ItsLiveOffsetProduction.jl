# A lazy image's values, computed once per tile and kept on disk.
#
# **What this is for.** A blocked correlation sweeps every block once per pass — a coarse and a fine sweep
# per chip-size level — so a whole run reads 3.4-4.4x the scene where one pass's windows sum to 0.66-0.96x
# (`src/api.jl:262`). Where those values are *derived* rather than read, the multiplier lands on the
# derivation: a coregistered secondary is resampled four times over, and a NISAR granule's samples are
# decompressed four times over.
#
# **Why not a memory cache.** The access order is a cyclic scan, so a tile's next use is a whole sweep
# away and any cache smaller than the touched set has evicted it by then; `memory.md` replays the real
# request sequence through LRU and random replacement and finds no interior optimum. `cache_budget`
# therefore chooses between holding the *whole* pair and holding none of it — and a NISAR pair is 23 GB
# (L1) to 48 GB (L2) of `Float32`, which is not a choice.
#
# **What this does instead.** Persist each tile the first time it is asked for. The first pass pays the
# resample or the decompression; every later pass reads an uncompressed tile back. Resident memory stays
# proportional to the window being served rather than to the scene, because the tiles live in a scratch
# file and are read through a memory mapping rather than `read`/`write`, so the kernel's page cache holds
# them outside the process's own footprint, where it can reclaim them under pressure.
#
# **Tiles rather than windows, because the windows change.** Each pass reads a block grown by that level's
# halo, so no two passes ask for the same rectangle. A cache keyed on the request would miss every time;
# one keyed on a fixed tile grid hits whatever the request's shape.
#
# **Memory-mapped, not a shared `IOStream`.** An earlier version held one lock around every tile's
# `seek` plus transfer, on the reasoning that a tile's own derivation is the expensive part and the
# transfer is not. That is true of the transfer alone, but the lock does not serialize only the transfer
# — it serializes every task's *turn* at the one shared file position, and at the request rates a
# blocked NISAR run reaches, that queuing is the cost, not the few blocked microseconds a single `seek`
# takes. A memory mapping has no shared position: each tile owns disjoint bytes of it, so two tasks
# filling different tiles never contend, and a repeat read is the kernel's own page cache rather than a
# lock's queue.

using Mmap

"""
    TileCache(parent; tile = 512, dir) <: AbstractMatrix{Float32}

`parent`'s values, computed a tile at a time and cached in a memory-mapped scratch file.

Reads the same values `parent` would return — assert it on a window if you want that checked; nothing
about the tiling reaches the result, since a tile is filled by asking `parent` for exactly that tile.

`tile` is the side of the cache's grid in pixels. It wants to be at least the chunk or burst granularity
of whatever `parent` reads, so that filling one tile is one read of the source rather than a fraction of
one, and small enough that a request's tile-aligned footprint is not much larger than the request.

The file is deleted by [`close`](@ref). It is sparse: only the tiles a run touches are ever written, so a
correlation over part of a scene costs that part.
"""
struct TileCache{P<:AbstractMatrix} <: AbstractMatrix{Float32}
    parent::P
    tile::Int
    dims::Tuple{Int,Int}
    ntiles::Tuple{Int,Int}
    io::IOStream
    path::String
    # The whole file, mapped once. Every tile is a disjoint `Float32` range of this one array, so
    # filling it is an ordinary array write with no file-position state to contend over.
    mapped::Vector{UInt8}
    # One lock per tile: two tasks must not derive the same absent tile at once, and a `ReentrantLock`'s
    # release is what makes a later reader see the write rather than a stale or torn one. No lock is
    # taken for a read once a tile is present — see the header note above this struct.
    tiles::Vector{ReentrantLock}
    present::Vector{Bool}
    # Atomic because tiles now fill concurrently: these are diagnostics, but a lost count would
    # misreport the reuse the cache is here to deliver.
    filled::Threads.Atomic{Int}
    served::Threads.Atomic{Int}
end

function TileCache(parent::AbstractMatrix; tile::Integer = 512,
                   dir::AbstractString = mktempdir(; cleanup = false))
    tile > 0 || throw(ArgumentError("`tile` must be positive, got $tile"))
    dims = size(parent)
    nt = AutoRIFT._tile_grid(dims, tile)
    nbytes = prod(nt) * tile * tile * sizeof(Float32)
    mkpath(dir)
    # A file of its own, created atomically. Both images of a pair are cached in one directory and are
    # usually the same size, so a name derived from the dimensions would give them one file, and each
    # would read back the other's tiles.
    path, io = mktemp(dir; cleanup = false)
    truncate(io, nbytes)
    mapped = Mmap.mmap(io, Vector{UInt8}, nbytes)
    return TileCache(parent, Int(tile), dims, nt, io, path, mapped,
                     [ReentrantLock() for _ in 1:prod(nt)], zeros(Bool, prod(nt)),
                     Threads.Atomic{Int}(0), Threads.Atomic{Int}(0))
end

Base.size(c::TileCache) = c.dims

# See the identical helper on `FilterTileCache` in `src/tile.jl`: `Mmap.jl` attaches its
# `munmap`/`UnmapViewOfFile` finalizer to the array itself on Julia 1.10, and to that array's
# underlying `Memory` (`arr.ref.mem`) from 1.11 on — `finalize(arr)` alone is a silent no-op on the
# newer layout.
function _finalize_mapping!(mapped)
    finalize(mapped)
    hasfield(typeof(mapped), :ref) && finalize(mapped.ref.mem)
    return nothing
end

"""
    close(c::TileCache)

Close the scratch file and delete it.
"""
function Base.close(c::TileCache)
    close(c.io)
    # The mapping behind `c.mapped` outlives `close(c.io)` until GC finalizes it, and Windows
    # refuses to delete a file with an active mapping where POSIX allows it. Finalizing explicitly
    # here costs one `munmap` slightly earlier than GC would have, once per cache.
    _finalize_mapping!(c.mapped)
    isfile(c.path) && rm(c.path)
    return nothing
end

# Bytes per tile, and where a tile starts. Every tile is stored at full size even where the image ends
# inside it, so the offset is arithmetic rather than a lookup; the file is sparse, so the padding of the
# last row and column of tiles occupies nothing.
#
# The grid arithmetic itself (`_tile_span`/`_tile_extent`/the offset formula) depends only on `dims` and
# `tile`, never on what a tile holds, so it lives once in `AutoRIFT._tile_grid_span` etc. — shared with
# `FilterTileCache`, the same cache one layer down, for filtered results instead of raw reads.
_tile_bytes(c::TileCache) = c.tile * c.tile * sizeof(Float32)
_tile_offset(c::TileCache, ti::Integer, tj::Integer) =
    AutoRIFT._tile_grid_offset(c.ntiles, _tile_bytes(c), ti, tj)
_tile_span(c::TileCache, r::AbstractUnitRange) = AutoRIFT._tile_grid_span(c.tile, r)
_tile_extent(c::TileCache, ti::Integer, tj::Integer) = AutoRIFT._tile_grid_extent(c.dims, c.tile, ti, tj)

# The `c.tile`-square plane backing one tile, as a view into the mapping — reading or writing through
# it reads or writes the file directly, with no separate transfer step.
function _tile_array(c::TileCache, ti::Integer, tj::Integer)
    off = _tile_offset(c, ti, tj)
    n = c.tile * c.tile
    return reshape(reinterpret(Float32, view(c.mapped, (off + 1):(off + 4n))), c.tile, c.tile)
end

# Make tile `(ti, tj)` present, deriving it from `c.parent` if it is not already.
#
# **The derivation happens under this tile's own lock, and only this tile's.** Two tasks deriving
# different tiles never contend — they take different locks and write disjoint bytes of the mapping —
# which is what makes a blocked run's concurrent misses actually run concurrently.
function _ensure_tile!(c::TileCache, ti::Integer, tj::Integer)
    k = (tj - 1) * c.ntiles[1] + ti
    @lock c.tiles[k] begin
        if c.present[k]
            Threads.atomic_add!(c.served, 1)
            return nothing
        end
        buf = _tile_array(c, ti, tj)
        rows, cols = _tile_extent(c, ti, tj)
        fill!(buf, 0.0f0)
        copyto!(view(buf, 1:length(rows), 1:length(cols)), c.parent[rows, cols])
        c.present[k] = true
        Threads.atomic_add!(c.filled, 1)
    end
    return nothing
end

function Base.getindex(c::TileCache, rows::AbstractUnitRange, cols::AbstractUnitRange)
    checkbounds(c, rows, cols)
    out = Matrix{Float32}(undef, length(rows), length(cols))
    for tj in _tile_span(c, cols), ti in _tile_span(c, rows)
        _ensure_tile!(c, ti, tj)
        buf = _tile_array(c, ti, tj)
        trows, tcols = _tile_extent(c, ti, tj)
        # The part of this tile the request wants, in the tile's own indices and in the output's.
        wr = intersect(rows, trows)
        wc = intersect(cols, tcols)
        (isempty(wr) || isempty(wc)) && continue
        copyto!(view(out, (wr .- first(rows) .+ 1), (wc .- first(cols) .+ 1)),
                view(buf, (wr .- first(trows) .+ 1), (wc .- first(tcols) .+ 1)))
    end
    return out
end

Base.getindex(c::TileCache, i::Integer, j::Integer) = c[i:i, j:j][1, 1]
Base.getindex(c::TileCache, rows::AbstractUnitRange, j::Integer) = c[rows, j:j][:, 1]
Base.getindex(c::TileCache, i::Integer, cols::AbstractUnitRange) = c[i:i, cols][1, :]

# Reading one of these costs I/O, which is what it is for. A blocked run therefore windows its reads, and
# the automatic cache budget will decline to hold a scene this size in memory.
AutoRIFT.ondisk(::TileCache) = true

"""
    cache_report(c::TileCache) -> String

How many tiles were derived and how many were served from the mapping.

The ratio is what the cache bought: a run that swept the grid `n` times reads `n` tiles for every one it
derives, so `served / filled` approaching `2 * levels - 1` is the multiplicity removed.
"""
cache_report(c::TileCache) = AutoRIFT._tile_cache_report(c.filled[], c.served[], _tile_bytes(c))
