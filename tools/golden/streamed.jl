# An optical scene as the correlator reads it, streamed from the granule's own object a window at a time.
#
# Every type here is a lazy `AbstractMatrix` indexed `(row, col)` that answers a range `getindex` by
# reading only what that window needs from its parent. They compose in the order `process.py` applies
# the steps to the native scene: the granule, then `apply_landsat_filtering`'s native filter, then
# `ensure_same_projection`'s warp, then the crop to the overlap. Nothing holds a whole scene except
# [`ScratchMatrix`](@ref), and that on disk: `Destripe`'s band-reject is a transform of the entire array.

using Mmap
import AWS

# The object a GDAL virtual path names, as OpticalDatasets reads it: `scene_path` spells a streamed scene
# for GDAL, which the rest of the harness opens it with. A requester-pays bucket is read through a
# transport that says who pays, one per bucket and kept for the session so its connections are reused.
const _S3_TRANSPORTS = Dict{String,Any}()
const _S3_TRANSPORTS_LOCK = ReentrantLock()

# The region each bucket the scenes come from lives in. S3 answers a request to any other region's
# endpoint with a redirect that the transport does not follow.
const _BUCKET_REGION = Dict("usgs-landsat" => "us-west-2")

function _s3_transport(bucket::AbstractString)
    region = get(_BUCKET_REGION, bucket, nothing)
    region === nothing && error("no region recorded for s3://$bucket; add it to `_BUCKET_REGION`")
    return ChunkManifests.S3Transport(bucket; aws = AWS.AWSConfig(; region), requesterpays = true)
end

function _optical_source(path::AbstractString)
    if startswith(path, "/vsis3/")
        uri = "s3://" * path[8:end]
        bucket = first(split(path[8:end], '/'))
        t = @lock _S3_TRANSPORTS_LOCK get!(() -> _s3_transport(bucket), _S3_TRANSPORTS, bucket)
        return (uri,), (; transport = t)
    end
    startswith(path, "/vsicurl/") && return (path[10:end],), (;)
    return (path,), (;)
end

"""
    open_scene(path) -> OpticalDatasets.OpticalRaster

The band at `path` — a local file, or a `/vsis3/` or `/vsicurl/` path as `scene_path` spells it — opened
for chunked reads from the object itself. Only the header is read.
"""
function open_scene(path::AbstractString)
    args, kwargs = _optical_source(path)
    return OpticalDatasets.open_optical(args...; kwargs...)
end

# Scratch arrays to close once a run is done with them: their files are deleted by `close`.
const SCRATCH = Any[]
const SCRATCH_LOCK = ReentrantLock()
scratch!(x) = (@lock SCRATCH_LOCK push!(SCRATCH, x); x)
function close_scratch!()
    @lock SCRATCH_LOCK begin
        foreach(close, SCRATCH)
        empty!(SCRATCH)
    end
    return nothing
end

# Scalar and single-row/column reads, through the window method each type below defines.
macro _window_indexing(T)
    quote
        Base.getindex(a::$(esc(T)), i::Integer, j::Integer) = a[i:i, j:j][1, 1]
        Base.getindex(a::$(esc(T)), rows::AbstractUnitRange, j::Integer) = a[rows, j:j][:, 1]
        Base.getindex(a::$(esc(T)), i::Integer, cols::AbstractUnitRange) = a[i:i, cols][1, :]
        AutoRIFT.ondisk(::$(esc(T))) = true
    end
end

# `rows` grown by `reach` on both sides and clamped to `1:n`.
_grown(rows::AbstractUnitRange, reach::Integer, n::Integer) =
    max(first(rows) - reach, 1):min(last(rows) + reach, n)

# ---------------------------------------------------------------------------
# The native-scene filter
# ---------------------------------------------------------------------------

"""
    NativeFiltered(parent, method::AutoRIFT.WallisGapfill, name) <: AbstractMatrix{Float32}

`parent`, the scene `name`, with `method` applied as `native_scene` applies it to the whole scene — gaps
are the pixels reading zero — computed a window at a time.

A window is filtered from `parent` grown by `AutoRIFT.filter_reach(method)` on every side and clamped
to the scene, so each pixel sees the neighbourhood a whole-scene call gives it: its validity verdict is
the whole scene's. The fill noise is not. Each window draws from a generator seeded by the scene's name
and the window's position, so the noise is reproducible, differs between the two scenes of a pair and
between windows, and differs from a whole-scene call's; nothing else does. The reference draws from
NumPy's unseeded global generator, so its noise is independent between the scenes too.
"""
struct NativeFiltered{A<:AbstractMatrix} <: AbstractMatrix{Float32}
    parent::A
    method::AutoRIFT.WallisGapfill
    reach::Int
    seed::UInt64
end

NativeFiltered(parent::AbstractMatrix, m::AutoRIFT.WallisGapfill, name::AbstractString) =
    NativeFiltered{typeof(parent)}(parent, m, AutoRIFT.filter_reach(m), hash(String(name)))

Base.size(f::NativeFiltered) = size(f.parent)

function Base.getindex(f::NativeFiltered, rows::AbstractUnitRange, cols::AbstractUnitRange)
    checkbounds(f, rows, cols)
    gr, gc = _grown(rows, f.reach, size(f, 1)), _grown(cols, f.reach, size(f, 2))
    img = Float32.(f.parent[gr, gc])
    out, _ = AutoRIFT.preprocess(img, img .!= 0, f.method, hash((first(gr), first(gc)), f.seed))
    return out[rows .- (first(gr) - 1), cols .- (first(gc) - 1)]
end

@_window_indexing NativeFiltered

"""
    ScratchMatrix(A; dir) <: AbstractMatrix{Float32}

`A` copied into a memory-mapped scratch file under `dir`, which [`close`](@ref) deletes.

For a filter that has to see the whole scene at once, so its result is held on disk rather than in the
process's own memory; the kernel's page cache serves the reads.
"""
struct ScratchMatrix <: AbstractMatrix{Float32}
    mapped::Matrix{Float32}
    io::IOStream
    path::String
end

function ScratchMatrix(A::AbstractMatrix; dir::AbstractString)
    mkpath(dir)
    path, io = mktemp(dir; cleanup = false)
    mapped = Mmap.mmap(io, Matrix{Float32}, size(A))
    copyto!(mapped, A)
    return ScratchMatrix(mapped, io, path)
end

Base.size(m::ScratchMatrix) = size(m.mapped)
Base.getindex(m::ScratchMatrix, rows::AbstractUnitRange, cols::AbstractUnitRange) =
    m.mapped[rows, cols]
@_window_indexing ScratchMatrix

function Base.close(m::ScratchMatrix)
    close(m.io)
    _finalize_mapping!(m.mapped)
    isfile(m.path) && rm(m.path)
    return nothing
end

"""
    fft_filtered(path, name, fp; dir) -> ScratchMatrix

`native_scene`'s `:fft` filter — Wallis at width 5, then the `Destripe` band-reject — of the scene at
`path`, streamed whole and kept on disk. `fp` is the scene's native footprint, whose CRS and spacing
orient the band-reject along the scan.
"""
function fft_filtered(path::AbstractString, name::AbstractString, fp; dir::AbstractString)
    r = open_scene(path)
    img = Matrix{Float32}(undef, size(r))
    OpticalDatasets.read_window!(img, r, axes(r, 1), axes(r, 2))
    along, cross = orbit_scan_angles(scene_ang(name, joinpath(CACHE, "angcache")), footprint_epsg(fp);
                                     spacing = fp.spacing)
    return ScratchMatrix(fft_native_filter(img, along, cross); dir)
end

# ---------------------------------------------------------------------------
# The warp to the reference's projection
# ---------------------------------------------------------------------------

"""
    WarpedScene(parent, source, target; nodata = nothing) <: AbstractMatrix{Float32}

`parent`, a scene on the grid of footprint `source`, resampled onto the grid of footprint `target` with
GDAL's Lanczos kernel — `_warp_to`'s warp, done one window at a time.

A window is warped from the part of `parent` it maps into: the window's edges are transformed into the
source grid, bounded, grown by `margin` pixels for the kernel's support and clamped to the scene. GDAL
then warps that piece, carrying its own geotransform, onto exactly the window's extent. A target pixel
whose source location lies outside the scene is zero, as it is in `_warp_to`'s file.

The source is handed to GDAL in `parent`'s element type and with `nodata`, the source file's no-data
value or `nothing`: a granule's integer band warps to rounded integers with its no-data pixels left out
of the kernel, as `_warp_to` writes it, and a filtered `Float32` scene, which `process.py` writes with no
no-data value, stays `Float32`.
"""
struct WarpedScene{A<:AbstractMatrix} <: AbstractMatrix{Float32}
    parent::A
    source_gt::NTuple{6,Float64}
    source_epsg::Int
    target_gt::NTuple{6,Float64}
    target_epsg::Int
    dims::Tuple{Int,Int}
    nodata::Union{Nothing,Float64}
    margin::Int
end

function WarpedScene(parent::AbstractMatrix, source, target; nodata = nothing, margin::Integer = 16)
    reverse(source.size) == size(parent) || throw(DimensionMismatch(
        "the source footprint is $(source.size) (x, y) but the scene is $(size(parent)) (row, col)"))
    gt(fp) = (Float64(fp.origin[1]), Float64(fp.spacing[1]), 0.0,
              Float64(fp.origin[2]), 0.0, Float64(fp.spacing[2]))
    return WarpedScene{typeof(parent)}(parent, gt(source), footprint_epsg(source), gt(target),
                                       footprint_epsg(target), reverse(target.size),
                                       nodata === nothing ? nothing : Float64(nodata), margin)
end

Base.size(w::WarpedScene) = w.dims

# The source pixel window covering target rows `rows` and columns `cols`, or `nothing` if they map
# entirely off the scene.
function _source_window(w::WarpedScene, rows::AbstractUnitRange, cols::AbstractUnitRange)
    t = w.target_gt
    # Pixel *edges* of the window, sampled densely enough that the bound is the curved edge's.
    n = 33
    es = range(first(cols) - 1, last(cols); length = n)
    ns = range(first(rows) - 1, last(rows); length = n)
    pc = vcat(es, es, fill(first(cols) - 1, n), fill(last(cols), n))
    pr = vcat(fill(first(rows) - 1, n), fill(last(rows), n), ns, ns)
    xs = t[1] .+ pc .* t[2]
    ys = t[4] .+ pr .* t[6]
    zs = zeros(length(xs))
    ArchGDAL.createcoordtrans(ArchGDAL.importEPSG(w.target_epsg; order = :trad),
                              ArchGDAL.importEPSG(w.source_epsg; order = :trad)) do ct
        ArchGDAL.transform!(xs, ys, zs, ct) || error(
            "a window edge of the warped scene did not transform from EPSG:$(w.target_epsg) to " *
            "EPSG:$(w.source_epsg)")
    end
    s = w.source_gt
    sc = (xs .- s[1]) ./ s[2]
    sr = (ys .- s[4]) ./ s[6]
    nr, nc = size(w.parent)
    r0 = max(floor(Int, minimum(sr)) - w.margin + 1, 1)
    r1 = min(ceil(Int, maximum(sr)) + w.margin, nr)
    c0 = max(floor(Int, minimum(sc)) - w.margin + 1, 1)
    c1 = min(ceil(Int, maximum(sc)) + w.margin, nc)
    (r0 > r1 || c0 > c1) && return nothing
    return (r0:r1, c0:c1)
end

function Base.getindex(w::WarpedScene, rows::AbstractUnitRange, cols::AbstractUnitRange)
    checkbounds(w, rows, cols)
    win = _source_window(w, rows, cols)
    win === nothing && return zeros(Float32, length(rows), length(cols))
    sr, sc = win
    src = w.parent[sr, sc]
    s, t = w.source_gt, w.target_gt
    sgt = [s[1] + (first(sc) - 1) * s[2], s[2], 0.0, s[4] + (first(sr) - 1) * s[6], 0.0, s[6]]
    x0 = t[1] + (first(cols) - 1) * t[2]
    x1 = t[1] + last(cols) * t[2]
    y0 = t[4] + (first(rows) - 1) * t[6]
    y1 = t[4] + last(rows) * t[6]
    flags = ["-of", "MEM", "-t_srs", "EPSG:$(w.target_epsg)",
             "-te", string(min(x0, x1)), string(min(y0, y1)), string(max(x0, x1)), string(max(y0, y1)),
             "-ts", string(length(cols)), string(length(rows)), "-r", "lanczos"]
    srs = ArchGDAL.toWKT(ArchGDAL.importEPSG(w.source_epsg))
    out = ArchGDAL.create(ArchGDAL.getdriver("MEM"); width = size(src, 2), height = size(src, 1),
                          nbands = 1, dtype = eltype(src)) do ds
        ArchGDAL.setgeotransform!(ds, sgt)
        ArchGDAL.setproj!(ds, srs)
        ArchGDAL.write!(ds, permutedims(src), 1)
        w.nodata === nothing || ArchGDAL.setnodatavalue!(ArchGDAL.getband(ds, 1), w.nodata)
        ArchGDAL.gdalwarp([ds], flags; dest = "/vsimem/$(gensym())") do v
            ArchGDAL.read(v, 1)
        end
    end
    size(out) == (length(cols), length(rows)) || error(
        "GDAL warped a $(length(rows))x$(length(cols)) window to $(reverse(size(out))); the window's " *
        "extent and the target grid disagree")
    return Float32.(permutedims(out))
end

@_window_indexing WarpedScene

# ---------------------------------------------------------------------------
# The correlator's gap mask
# ---------------------------------------------------------------------------

"""
    NonzeroMask(img) <: AbstractMatrix{Bool}

`img .!= 0`, computed a window at a time: the gaps the reference's `wallis_fill` detects as zero
(`autoRIFT.py`'s `np.isclose(image, 0)`), stated as the validity mask `AutoRIFT.WallisGapfill` reads them
from.
"""
struct NonzeroMask{A<:AbstractMatrix} <: AbstractMatrix{Bool}
    parent::A
end

Base.size(m::NonzeroMask) = size(m.parent)
Base.getindex(m::NonzeroMask, rows::AbstractUnitRange, cols::AbstractUnitRange) =
    m.parent[rows, cols] .!= 0
Base.getindex(m::NonzeroMask, i::Integer, j::Integer) = m.parent[i, j] != 0
AutoRIFT.ondisk(m::NonzeroMask) = AutoRIFT.ondisk(m.parent)

# ---------------------------------------------------------------------------
# Filling a cache
# ---------------------------------------------------------------------------

"""
    fill_tiles!(c::TileCache) -> c

Derive every tile of `c`, concurrently. A tile's derivation is a remote read and possibly a filter and a
warp, and the no-data mask samples every tile, so doing it here up front is what keeps those reads
parallel rather than one tile after another inside the mask's row bands.
"""
function fill_tiles!(c::TileCache)
    nt = c.ntiles
    Threads.@threads :dynamic for k in 1:prod(nt)
        ti, tj = mod1(k, nt[1]), cld(k, nt[1])
        _ensure_tile!(c, ti, tj)
    end
    return c
end
