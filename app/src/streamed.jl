# An optical scene as the correlator reads it, streamed from the granule's own object a window at a time.
#
# Every type here is a lazy `AbstractMatrix` indexed `(row, col)` that answers a range `getindex` by
# reading only what that window needs from its parent. They compose in the order `process.py` applies the
# steps to the native scene: the granule, then `apply_landsat_filtering`'s native filter, then
# `ensure_same_projection`'s warp, then the crop to the overlap. Nothing holds a whole scene except
# [`ScratchMatrix`](@ref), and that on disk: `Destripe`'s band-reject is a transform of the entire array.

"""
    scene_eltype(uri) -> DataType

The sample type of the band at `uri`, read from its header.
"""
function scene_eltype(uri::String)
    t = transport_for(uri)
    lower = lowercase(uri)
    return endswith(lower, ".jp2") ? ChunkManifests.jp2_eltype(uri; transport = t) :
           ChunkManifests.tiff_eltype(uri; transport = t)
end

"""
    open_scene(uri, T) -> OpticalDatasets.OpticalRaster{T}

The band at `uri`, opened for chunked reads of `T` samples from the object itself. Only the header is
read.
"""
open_scene(uri::String, ::Type{T}) where {T} =
    OpticalDatasets.open_optical(uri, T; transport = transport_for(uri))

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
_grown(rows::AbstractUnitRange, reach::Int, n::Int) = max(first(rows) - reach, 1):min(last(rows) + reach, n)

"""
    NativeFiltered(parent, method::AutoRIFT.WallisGapfill, name) <: AbstractMatrix{Float32}

`parent`, the scene `name`, with `method` applied as `native_scene` applies it to the whole scene — gaps
are the pixels reading zero — computed a window at a time. A window is filtered from `parent` grown by the
filter's reach on every side, so each pixel's validity verdict is the whole scene's. The fill noise is
seeded by the scene's name and the window's position, so it is reproducible and differs between the two
scenes of a pair.
"""
struct NativeFiltered{A<:AbstractMatrix} <: AbstractMatrix{Float32}
    parent::A
    method::AutoRIFT.WallisGapfill
    reach::Int
    seed::UInt64
end

NativeFiltered(parent::AbstractMatrix, m::AutoRIFT.WallisGapfill, name::String) =
    NativeFiltered{typeof(parent)}(parent, m, AutoRIFT.filter_reach(m), hash(name))

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

`A` copied into a memory-mapped scratch file under `dir`, which `close` deletes. For a filter that has to
see the whole scene at once, so its result is held on disk rather than in the process's own memory.
"""
struct ScratchMatrix <: AbstractMatrix{Float32}
    mapped::Matrix{Float32}
    io::IOStream
    path::String
end

function ScratchMatrix(A::Matrix{Float32}; dir::String)
    mkpath(dir)
    path, io = mktemp(dir; cleanup = false)
    mapped = Mmap.mmap(io, Matrix{Float32}, size(A))
    copyto!(mapped, A)
    return ScratchMatrix(mapped, io, path)
end

Base.size(m::ScratchMatrix) = size(m.mapped)
Base.getindex(m::ScratchMatrix, rows::AbstractUnitRange, cols::AbstractUnitRange) = m.mapped[rows, cols]
@_window_indexing ScratchMatrix

function Base.close(m::ScratchMatrix)
    close(m.io)
    finalize(m.mapped)
    finalize(m.mapped.ref.mem)
    isfile(m.path) && rm(m.path)
    return nothing
end

"""
    fft_native_filter(img, along, cross) -> Matrix{Float32}

`process.py`'s `:fft` native filter of a whole scene: Wallis at width 5, then the `Destripe` band-reject
oriented by the scan's `along`/`cross` angles, with the pixels reading zero kept at zero through both.
"""
function fft_native_filter(img::Matrix{Float32}, along::Float64, cross::Float64)
    valid = img .!= 0
    w, _ = AutoRIFT.preprocess(img, valid, AutoRIFT.Wallis(5, 0.0))
    w[.!valid] .= 0.0f0
    d, _ = AutoRIFT.preprocess(w, valid, AutoRIFT.Destripe(; along_track = along, cross_track = cross))
    d[.!valid] .= 0.0f0
    return d
end

"""
    fft_filtered(raster, scene, fp; dir) -> ScratchMatrix

[`fft_native_filter`](@ref) of the whole of `raster`, the band of `scene`, kept on disk. `fp` is the
scene's native footprint, whose CRS and spacing orient the band-reject along the scan from the orbit in
the scene's `_ANG.txt`.
"""
function fft_filtered(raster::OpticalDatasets.OpticalRaster{T}, scene::SceneSource, fp::Footprint;
                      dir::String) where {T}
    img = Matrix{Float32}(undef, size(raster))
    OpticalDatasets.read_window!(img, raster, axes(raster, 1), axes(raster, 2))
    along, cross = orbit_scan_angles(read_text(scene.ang), footprint_epsg(fp); spacing = fp.spacing)
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
source grid, bounded, grown by `margin` pixels for the kernel's support and clamped to the scene. GDAL then
warps that piece, carrying its own geotransform, onto exactly the window's extent; a target pixel whose
source lies outside the scene is zero. The source keeps `parent`'s element type and `nodata`, the source
file's no-data value or `nothing`, so a granule's integer band warps to rounded integers with its no-data
pixels left out of the kernel, as `_warp_to` writes it.
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
    source_wkt::String
end

_footprint_gt(fp::Footprint) = (fp.origin[1], fp.spacing[1], 0.0, fp.origin[2], 0.0, fp.spacing[2])

function WarpedScene(parent::AbstractMatrix, source::Footprint, target::Footprint;
                     nodata::Union{Nothing,Float64} = nothing, margin::Int = 16)
    reverse(source.size) == size(parent) || throw(DimensionMismatch(
        "the source footprint and the scene it describes differ in size"))
    return WarpedScene{typeof(parent)}(parent, _footprint_gt(source), footprint_epsg(source),
                                       _footprint_gt(target), footprint_epsg(target), reverse(target.size),
                                       nodata, margin, epsg_wkt(footprint_epsg(source)))
end

Base.size(w::WarpedScene) = w.dims

# The source pixel window covering target rows `rows` and columns `cols`, or `nothing` if they map entirely
# off the scene.
function _source_window(w::WarpedScene, rows::AbstractUnitRange, cols::AbstractUnitRange)
    t = w.target_gt
    # Pixel *edges* of the window, sampled densely enough that the bound is the curved edge's.
    n = 33
    es = range(first(cols) - 1, last(cols); length = n)
    ns = range(first(rows) - 1, last(rows); length = n)
    pc = vcat(es, es, fill(Float64(first(cols) - 1), n), fill(Float64(last(cols)), n))
    pr = vcat(fill(Float64(first(rows) - 1), n), fill(Float64(last(rows)), n), ns, ns)
    xs = t[1] .+ pc .* t[2]
    ys = t[4] .+ pr .* t[6]
    zs = zeros(length(xs))
    ok = with_transform(ct -> transform_points!(ct, xs, ys, zs), w.target_epsg, w.source_epsg)
    ok || throw(ErrorException("a window edge of a warped scene did not transform from EPSG:" *
                               string(w.target_epsg) * " to EPSG:" * string(w.source_epsg)))
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
    return _warp_window(w, w.parent[sr, sc], sr, sc, rows, cols)
end

function _warp_window(w::WarpedScene, src::Matrix{T}, sr::UnitRange{Int}, sc::UnitRange{Int},
                      rows::AbstractUnitRange, cols::AbstractUnitRange) where {T}
    s, t = w.source_gt, w.target_gt
    sgt = (s[1] + (first(sc) - 1) * s[2], s[2], 0.0, s[4] + (first(sr) - 1) * s[6], 0.0, s[6])
    x0 = t[1] + (first(cols) - 1) * t[2]
    x1 = t[1] + last(cols) * t[2]
    y0 = t[4] + (first(rows) - 1) * t[6]
    y1 = t[4] + last(rows) * t[6]
    flags = ["-of", "MEM", "-t_srs", "EPSG:" * string(w.target_epsg),
             "-te", string(min(x0, x1)), string(min(y0, y1)), string(max(x0, x1)), string(max(y0, y1)),
             "-ts", string(length(cols)), string(length(rows)), "-r", "lanczos"]
    ds = mem_dataset(permutedims(src), sgt, w.source_wkt; nodata = w.nodata)
    out = try
        v = warp(ds, flags, vsimem_name("warped_window"))
        try
            nx, ny = rastersize(v)
            (nx, ny) == (length(cols), length(rows)) || throw(ErrorException(
                "GDAL warped a window to a different size; the window's extent and the target grid disagree"))
            read_window(v, T, 0, 0, nx, ny)
        finally
            close(v)
        end
    finally
        close(ds)
    end
    return Float32.(permutedims(out))
end

@_window_indexing WarpedScene

"""
    LazyWindow(parent, row0, col0, dims) <: AbstractMatrix{Float32}

The `dims` window of `parent` starting after row `row0` and column `col0`, read through `parent`'s own
range indexing — which a `view` of a lazy array would defer to one read per element.
"""
struct LazyWindow{A<:AbstractMatrix} <: AbstractMatrix{Float32}
    parent::A
    row0::Int
    col0::Int
    dims::Tuple{Int,Int}
end

Base.size(w::LazyWindow) = w.dims

function Base.getindex(w::LazyWindow, rows::AbstractUnitRange, cols::AbstractUnitRange)
    checkbounds(w, rows, cols)
    return Float32.(w.parent[(w.row0 .+ rows), (w.col0 .+ cols)])
end

@_window_indexing LazyWindow

"""
    NonzeroMask(img) <: AbstractMatrix{Bool}

`img .!= 0`, a window at a time: the gaps the reference's `wallis_fill` detects as zero
(`np.isclose(image, 0)`), as the validity mask `AutoRIFT.WallisGapfill` reads them from.
"""
struct NonzeroMask{A<:AbstractMatrix} <: AbstractMatrix{Bool}
    parent::A
end

Base.size(m::NonzeroMask) = size(m.parent)
# A `Matrix{Bool}` rather than the `BitMatrix` a broadcast comparison gives, so a window and a whole-mask
# read are one type.
Base.getindex(m::NonzeroMask, rows::AbstractUnitRange, cols::AbstractUnitRange) =
    map(!=(0), m.parent[rows, cols])::Matrix{Bool}
Base.getindex(m::NonzeroMask, i::Integer, j::Integer) = m.parent[i, j] != 0
AutoRIFT.ondisk(m::NonzeroMask) = AutoRIFT.ondisk(m.parent)
