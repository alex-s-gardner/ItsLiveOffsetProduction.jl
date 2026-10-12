# The correlator's two images for an optical pair, and the driver's no-data mask.

# Which raw-image pixel the no-data mask is read at, as an offset added to the geogrid's location band.
# The reference indexes `I1[yGrid - 1, xGrid - 1]` of a zero-based array over a band it read straight out of
# `window_location.tif` (`testautoRIFT.py:341,602`), so in one-based terms it samples the location value
# itself. Measured: with this offset the grown mask equals the reference's own `noDataMask` at every one of
# 5,503,691 points.
const NODATA_SAMPLE_OFFSET = 0

"""
    Imagery

The two images the correlator is handed, each a filled [`TileFile`](@ref) cropped to the overlap, in
acquisition order, and the driver's no-data mask over the grid.
"""
struct Imagery
    early::TileFile
    late::TileFile
    declined::Matrix{Bool}
end

Base.close(m::Imagery) = (close(m.early); close(m.late); nothing)

"""
    native_imagery(s::Setup; dir, cache_io) -> Imagery

The two images the correlator is handed, unfiltered by the correlator's own filter, and the driver's no-data
mask. Each image is the granule's object streamed a tile at a time through the steps `process.py` applies
to the native scene — the native filter, the warp to the earlier scene's projection, the crop to the
overlap — and kept in a tile cache under `dir`; no scene-sized array is held in memory except where the
`:fft` filter needs the whole scene. Not quantized to `UInt8`: the reference's `uniform_data_type` exists
for its own memory budget, and `ZNCC` removes each window's own mean and scale.

A point is declined where *either* scene reads zero at its location in the scene after the native filter and
the warp, which is what the reference tests (`testautoRIFT.py:337-349`), or where the geogrid placed it
outside the image. A NISAR L2 pair's images are its GSLCs' amplitudes ([`gslc_image`](@ref)).
"""
function native_imagery(s::Setup; dir::String, cache_io::Symbol)
    g = s.geometry
    lx, ly = permutedims(g.location_x), permutedims(g.location_y)
    sentinel = Int32(g.nodata.output)
    declined = _sentinel_declined(lx, ly, sentinel)
    want = (s.pair.coordinate.size[2], s.pair.coordinate.size[1])     # (rows, cols) of the overlap
    eg, lg = s.early_gslc, s.late_gslc
    early = eg === nothing ? scene_image(s, s.early, s.early_fp, s.pair.reference_offset, want; dir, cache_io) :
            gslc_image(eg, s.pair.reference_offset, want; dir, cache_io)
    _accumulate_zeros!(declined, early, lx, ly, sentinel, NODATA_SAMPLE_OFFSET)
    late = lg === nothing ? scene_image(s, s.late, s.late_fp, s.pair.secondary_offset, want; dir, cache_io) :
           gslc_image(lg, s.pair.secondary_offset, want; dir, cache_io)
    _accumulate_zeros!(declined, late, lx, ly, sentinel, NODATA_SAMPLE_OFFSET)
    return Imagery(early, late, declined)
end

"""
    scene_image(s, scene, fp, offset, want; dir, cache_io) -> TileFile

`scene` as the correlator reads it — the granule's band, its native filter, the warp when the pair is in
two projections — cropped at `offset` to `want` and filled into a tile cache. One concrete chain of types
per branch below; the tile file returned has the same type whichever branch made it.
"""
function scene_image(s::Setup, scene::SceneSource, fp::Footprint, offset::NTuple{2,Int}, want::Tuple{Int,Int};
                     dir::String, cache_io::Symbol)
    T = scene_eltype(scene.band)
    T === UInt16 && return scene_image(s, scene, fp, offset, want, open_scene(scene.band, UInt16); dir, cache_io)
    T === UInt8 && return scene_image(s, scene, fp, offset, want, open_scene(scene.band, UInt8); dir, cache_io)
    throw(ArgumentError(scene.band * " holds samples of a type the optical chain does not read"))
end

function scene_image(s::Setup, scene::SceneSource, fp::Footprint, offset::NTuple{2,Int}, want::Tuple{Int,Int},
                     src::OpticalDatasets.OpticalRaster{T}; dir::String, cache_io::Symbol) where {T}
    reverse(fp.size) == size(src) || throw(DimensionMismatch(
        scene.name * ": its geotransform and its band disagree on the scene's size"))
    nf = native_filter(s.job, scene.name)
    if nf === :none
        return _filled(s, scene, fp, offset, want, src, true; dir, cache_io)
    elseif nf === :wallis_fill
        return _filled(s, scene, fp, offset, want,
                       NativeFiltered(src, AutoRIFT.WallisGapfill(5, WALLIS_STD_CUTOFF), scene.name), false;
                       dir, cache_io)
    elseif nf === :fft
        m = fft_filtered(src, scene, fp; dir)
        try
            return _filled(s, scene, fp, offset, want, m, false; dir, cache_io)
        finally
            close(m)
        end
    end
    throw(ArgumentError("no route for native filter " * string(nf) * " on " * scene.name))
end

# The scene `img`, warped when the pair is in two projections, cropped and filled into a tile cache. A
# granule's own band warps with its no-data value; a filtered scene, which `process.py` writes with none,
# warps without (`raw`).
function _filled(s::Setup, scene::SceneSource, fp::Footprint, offset::NTuple{2,Int}, want::Tuple{Int,Int},
                 img::AbstractMatrix, raw::Bool; dir::String, cache_io::Symbol)
    if s.warp
        nodata = raw ? with_dataset(band_nodata, gdal_path(scene.band)) : nothing
        # Each scene onto the grid the warp writes it to: its own extent, at the earlier scene's spacing.
        target = warped_footprint(scene.band, s.early_fp, s.epsg)
        return _cropped(WarpedScene(img, fp, target; nodata), offset, want; dir, cache_io)
    end
    return _cropped(img, offset, want; dir, cache_io)
end

function _cropped(scene::AbstractMatrix, offset::NTuple{2,Int}, want::Tuple{Int,Int}; dir::String,
                  cache_io::Symbol)
    ox, oy = offset
    (size(scene, 1) >= oy + want[1] && size(scene, 2) >= ox + want[2]) || throw(DimensionMismatch(
        "the overlap does not fit in the scene; the crop and the geotransform disagree"))
    return fill_tiles!(TileCache(LazyWindow(scene, oy, ox, want); tile = 1024, dir, io = cache_io))
end

# A point outside the image is declined before any scene is read.
function _sentinel_declined(lx::Matrix{Int32}, ly::Matrix{Int32}, sentinel::Int32)
    out = Matrix{Bool}(undef, size(lx))
    Threads.@threads for j in axes(lx, 2)
        for i in axes(lx, 1)
            out[i, j] = lx[i, j] == sentinel || ly[i, j] == sentinel
        end
    end
    return out
end

"""
    _accumulate_zeros!(declined, scene, lx, ly, sentinel, offset; rows = 1024)

One scene's contribution to the no-data mask: a grid point is declined where the scene reads zero at its
location band plus `offset`, or where that sample falls outside the scene. Accumulated with `|=`, because
either scene reading zero is enough. The scene is read a band of `rows` rows at a time, each band once.
"""
function _accumulate_zeros!(out::Matrix{Bool}, scene::TileFile, lx::Matrix{Int32}, ly::Matrix{Int32},
                            sentinel::Int32, offset::Int; rows::Int = 1024)
    nr, nc = size(scene)
    for j in axes(lx, 2), i in axes(lx, 1)
        x, y = lx[i, j], ly[i, j]
        (out[i, j] || x == sentinel || y == sentinel) && continue
        r, c = Int(y) + offset, Int(x) + offset
        (1 <= r <= nr && 1 <= c <= nc) || (out[i, j] = true)
    end
    for r0 in 1:rows:nr
        band = r0:min(r0 + rows - 1, nr)
        b = scene[band, 1:nc]
        Threads.@threads for j in axes(lx, 2)
            for i in axes(lx, 1)
                x, y = lx[i, j], ly[i, j]
                (out[i, j] || x == sentinel || y == sentinel) && continue
                r = Int(y) + offset
                r in band || continue
                iszero(b[r - r0 + 1, Int(x) + offset]) && (out[i, j] = true)
            end
        end
    end
    return out
end
