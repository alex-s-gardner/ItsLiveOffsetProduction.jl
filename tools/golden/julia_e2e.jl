# From two granules to the ITS_LIVE product netCDF, with nothing taken from the reference.
#
#   AWS_PROFILE=itslive julia --project=tools/golden -t 12 tools/golden/julia_e2e.jl LC08_L1TP_009011_20200703
#   AWS_PROFILE=itslive julia --project=tools/golden -t 12 tools/golden/julia_e2e.jl \
#       LC08_L1TP_009011_20200703 --compare-intermediate
#
# **The optical granules are streamed, never staged.** Each scene is read from its own object — Landsat
# from `s3://usgs-landsat`, Sentinel-2 from the Google mirror over HTTPS — a chunk at a time, inside the
# timed region, as `process.py` reads them through GDAL's `/vsis3/` and `/vsicurl/`. `s3://usgs-landsat`
# is requester-pays, so `AWS_PROFILE` has to name credentials that can pay. What is kept locally is a
# run's scratch: the tiles of each image the correlator has already been handed, deleted when the run
# ends (`TileCache`, `streamed.jl`). A NISAR granule is read from ASF over HTTPS with an Earthdata Login
# token taken from `~/.netrc` (`nisar_url`). A Sentinel-1 pair is still read from local SAFE trees: a
# full-SLC granule is a `.zip` on ASF, expanded under the golden cache (`stage_safe`), and a burst pair's
# SAFE trees are in its run directory.
#
# **What was missing.** `e2e_run.jl` runs the chain from the granule to `dx`/`dy`, and
# `validate_itslive_write.jl` writes a product — but from the reference's own captured arrays. So no run
# has produced a product from the imagery alone, which is the claim this makes.
#
# Four things sit between those two scripts, all of them the ITS_LIVE driver's rather than the
# correlator's, and none reproduced by any rung of the ladder:
#
#   nodata      the imagery's zero mask sampled at each grid point, which zeroes the search limit and
#               both chip bounds before correlation (`testautoRIFT.py:337-403`).
#   chop        the grid truncated to a multiple of `max(ChipSizeMaxX) / ChipSize0X`
#               (`autoRIFT.py:767`), the result pasted back into an `origSize` array
#               (`testautoRIFT.py:792-805`). The dropped strip is nodata *in the product*, so it is not
#               a reporting detail.
#   metadata    the two STAC items, the parameter shapefile, and the `IMG_INFO_DICT` attribute family.
#   package     `write_product` of an `ItsLiveInput` built from all of the above.
#
# **`--compare-intermediate` is where a disagreement is localized.** `autoRIFT_intermediate.nc` holds the
# reference's `Dx`, `Dy`, `ChipSizeX`, `InterpMask` and `SearchLimitX` on the chopped grid and its
# `noDataMask` at `origSize` — every array this script derives, before any stable shift, on exactly the
# grids this produces them on. A product diff alone cannot say whether a difference entered at the mask,
# the chop or the correlator; that mode can.

include("e2e_run.jl")      # `setup`, `_chop_to`, `_cropped_scene`, the memory trace, `packaging.jl`
include("streamed.jl")     # the lazy, streamed optical scene the correlator reads

# NISAR granules are read from ASF a chunk at a time, as the optical ones are from their buckets.
NISAR_REMOTE[] = true
include("product.jl")
include("compare.jl")

using Dates, NCDatasets, OpticalDatasets
import AWSS3, ChunkManifests
# Loaded, not called: `cf_grid_mapping` comes from `ItsLiveOffsetProduction`, whose trigger is
# Rasters plus DimensionalData plus ArchGDAL plus DiskArrays, and a product cannot be written without
# the CF grid-mapping attributes it builds.
using DimensionalData, Rasters

# ---------------------------------------------------------------------------
# What is declared rather than derived
# ---------------------------------------------------------------------------

# Which raw-image pixel the no-data mask is read at, as an offset added to the geogrid's own location
# band. The reference indexes `I1[yGrid - 1, xGrid - 1]` with a zero-based NumPy array
# (`testautoRIFT.py:341`) over a band it read straight out of `window_location.tif` (`:602`), so in
# one-based terms it samples the location value itself. Measured: with this offset the grown mask equals
# the reference's own `noDataMask` at every one of 5,503,691 points, and `+1` reaches 99.9868%.
# `--compare-intermediate` scores both and warns if this one is not the better.
const NODATA_SAMPLE_OFFSET = 0

# `arImgDisp_*` negates `Dy0` as its first act and AutoRIFT.jl has no such internal step, so
# `MultichipResult.dy` runs opposite to the reference's `Dy` — which `testautoRIFT.py:800` passes
# straight into `netCDF_packaging` as `DY`, with no sign change of its own.
# Measured on the Landsat 8 case: this sign reaches 86.17% of `dy` exact against the reference's own,
# and the other 7.19%. `--compare-intermediate` prints both.
const DY_SIGN = -1

# Which run directory the native product is written into. A number no container run uses, so it never
# lands beside a reference product and is never picked up by `run_product`.
const RUN_FOR_OUTPUT = 900

# ---------------------------------------------------------------------------
# The imagery the correlator is handed
# ---------------------------------------------------------------------------

"""
    native_imagery(s::Setup; cache_io = :pread) -> (early, late, declined, declined_shifted, filter)

The two images the correlator is handed, unfiltered by the correlator's own filter, and the driver's
no-data mask.

**Each image is lazy, streamed and cached a tile at a time.** [`streamed_scene`](@ref) composes the
granule's object with the steps `process.py` applies to the native scene; the overlap is a crop of that
([`LazyWindow`](@ref)), and a [`TileCache`](@ref) keeps each tile once derived, because the correlator
reads every block once per pass. No scene-sized array is held in memory. The correlator applies
`filter` itself a block at a time — `native_run` passes it as `preprocess` and correlates blocked.

**Not quantized to `UInt8`.** The reference's `uniform_data_type` (`autoRIFT.py:359-384`) exists for its
own memory budget, not for anything the correlator needs — `AutoRIFT.jl`'s correlator takes any element
type, and `ZNCC` already removes each correlation window's own mean and scale, so a whole-image mean,
standard deviation and rescale onto 256 levels buys nothing here, at the cost `dev/CORRECTNESS.md` item
4b on `AutoRIFT.jl` quantifies — up to 7 px on one scene, a 31x tie-breaking disparity in the base
level's weakest correlation quartile. So this is a deliberate divergence from what the Python reference's
current `DataType = 0` default produces.

**The mask is sampled from the image the correlator is handed.** A point is declined where *either*
scene reads zero (`testautoRIFT.py:337-349`) in the scene after the native filter and the warp, which is
what the reference tests. Both sampling offsets are accumulated, so `--compare-intermediate` can score
the choice.
"""
function native_imagery(s::Setup; cache_io::Symbol = :pread)
    c = s.case
    m = correlator_filter(c)
    early, late = acquisition_order(c)
    # `setup` coregisters the scenes `scene_path` names, and it names them by acquisition — `:reference`
    # is the earlier — so `reference_offset` is the earlier scene's on a reversed job too.
    s.reference_path == scene_path(c, :reference) || error(
        "$(basename(c.product)): the setup's reference scene is $(s.reference_path), not the earlier " *
        "acquisition's $(scene_path(c, :reference)), so the coregistration offsets would be paired " *
        "with the wrong scenes")

    want = reverse(_coord_size(s.pair.coordinate))      # (ny, nx) of the overlap
    g = s.geometry
    lx, ly = permutedims(g.location_x), permutedims(g.location_y)
    sentinel = Int32(g.nodata.output)
    # Seeded from the geogrid's own out-of-image sentinel (`testautoRIFT.py:603`), before any imagery.
    declined = _sentinel_declined(lx, ly, sentinel)
    shifted = copy(declined)

    rfp, sfp = scene_footprint(s.reference_path), scene_footprint(s.secondary_path)
    # `ensure_same_projection` warps both scenes or neither: `-tap` can move even the reference's origin.
    warp = footprint_epsg(rfp) != footprint_epsg(sfp)
    images = map(((s.reference_path, early, rfp, s.pair.reference_offset),
                  (s.secondary_path, late, sfp, s.pair.secondary_offset))) do (path, name, fp, off)
        scene = streamed_scene(c, path, name, fp; warp_to = warp ? (rfp, s.epsg) : nothing)
        ox, oy = off
        size(scene, 1) >= oy + want[1] && size(scene, 2) >= ox + want[2] || throw(DimensionMismatch(
            "the overlap at offset $off does not fit in a $(size(scene)) scene; the crop and the " *
            "geotransform disagree"))
        img = scratch!(TileCache(LazyWindow(scene, oy, ox, want); tile = 1024,
                                 dir = joinpath(CACHE, "scratch"), io = cache_io))
        fill_tiles!(img)
        _accumulate_zeros_banded!((declined => NODATA_SAMPLE_OFFSET, shifted => NODATA_SAMPLE_OFFSET + 1),
                                  img, lx, ly, sentinel)
        img
    end
    return (images[1], images[2], declined, shifted, m)
end

"""
    streamed_scene(c, path, name, fp; warp_to = nothing) -> AbstractMatrix

The scene `name` at `path` as geogrid and the correlator see it: read from its object, with the native
filter `process.py` applies to it ([`native_filter`](@ref)), and — when `warp_to` is
`(reference_footprint, epsg)` — resampled onto the grid `_warp_to` would write. `fp` is the scene's own
footprint. Lazy throughout except for an `:fft` filter, which is computed whole and held on disk.
"""
function streamed_scene(c::GoldenCase, path::AbstractString, name::AbstractString, fp; warp_to = nothing)
    src = open_scene(path)
    reverse(fp.size) == size(src) || error(
        "$name: its geotransform describes $(fp.size) (x, y) but the band is $(size(src)) (row, col)")
    nf = native_filter(c, name)
    img = nf === nothing ? src :
          nf === :wallis_fill ? NativeFiltered(src, AutoRIFT.WallisGapfill(5, 0.25), name) :
          nf === :fft ? scratch!(fft_filtered(path, name, fp; dir = joinpath(CACHE, "scratch"))) :
          error("no route for native filter `$nf` on $name")
    warp_to === nothing && return img
    reference, epsg = warp_to
    # A filtered scene is written with no no-data value; the granule's own band carries one.
    nodata = nf === nothing ? scene_nodata(path) : nothing
    return WarpedScene(img, fp, warped_footprint(path, reference, epsg); nodata)
end

# The no-data value of the band at `path`, or `nothing`, read from its header.
function scene_nodata(path::AbstractString)
    gdal_network_setup()
    return ArchGDAL.read(ds -> ArchGDAL.getnodatavalue(ArchGDAL.getband(ds, 1)), path)
end

# A point outside the image is declined before any scene is read.
function _sentinel_declined(lx::AbstractMatrix{<:Integer}, ly::AbstractMatrix{<:Integer},
                            sentinel::Int32)
    out = Matrix{Bool}(undef, size(lx))
    Threads.@threads for j in axes(lx, 2)
        for i in axes(lx, 1)
            out[i, j] = lx[i, j] == sentinel || ly[i, j] == sentinel
        end
    end
    return out
end

"""
    _accumulate_zeros_banded!(outs, scene, lx, ly, sentinel; rows = 1024) -> outs

One scene's contribution to each no-data mask `out => offset` in `outs`: a grid point is declined where
the scene reads zero at its location band plus `offset`, or where that sample falls outside the scene.
Accumulated with `|=`, because either scene reading zero is enough. The scene is read a band of `rows`
rows at a time, each band once, so the largest array held is one band.
"""
function _accumulate_zeros_banded!(outs, scene::AbstractMatrix, lx::AbstractMatrix{<:Integer},
                                   ly::AbstractMatrix{<:Integer}, sentinel::Int32; rows::Integer = 1024)
    nr, nc = size(scene)
    # Off the scene is declined without a read.
    for (out, offset) in outs, j in axes(lx, 2), i in axes(lx, 1)
        x, y = lx[i, j], ly[i, j]
        (out[i, j] || x == sentinel || y == sentinel) && continue
        r, c = Int(y) + offset, Int(x) + offset
        (1 <= r <= nr && 1 <= c <= nc) || (out[i, j] = true)
    end
    for r0 in 1:rows:nr
        band = r0:min(r0 + rows - 1, nr)
        b = scene[band, 1:nc]
        for (out, offset) in outs
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
    end
    return outs
end

# ---------------------------------------------------------------------------
# The driver's no-data mask, chop, and origSize padding
# ---------------------------------------------------------------------------

"""
    dilated_nodata_mask(declined) -> Matrix{Bool}

`declined` grown by one point in every direction, which is what gates the *outputs*.

The driver masks twice with two different masks, and the difference is one dilation
(`testautoRIFT.py:519-521`, `cv2.dilate` with a 3x3 kernel and one iteration). The narrow mask zeroes
the search limit and the chip bounds *before* correlating; the grown one is applied to `DX`, `DY`,
`CHIPSIZEX` and `SEARCHLIMITX` *after* (`:806-818`), so it is what decides the product's coverage and
`ItsLiveOffsetProduction.roi_valid_percentage`'s denominator.

Measured on the Landsat 8 case: dilating closes the mask's agreement with the reference's own
`noDataMask` from 99.8399% to 100%, which is 8,811 points.

`_dilate` is `e2e.jl`'s separable sweep; at radius 1 it is the 3x3 square structuring element.
"""
dilated_nodata_mask(declined::AbstractMatrix{Bool}) = _dilate(declined, 1)

"""
    driver_chip_size_max(s::Setup, declined) -> Int

The largest x chip bound surviving the no-data mask, which is the maximum the driver hands the
correlator and so the top of its pyramid.

`AutoRIFT.params(g)` reads the maximum off the geometry alone, which is the larger number whenever the
imagery is gappy: a `PairGeometry` carries no imagery and cannot know which points the driver dropped.
"""
function driver_chip_size_max(s::Setup, declined::AbstractMatrix{Bool})
    g = s.geometry
    best = _max_outside(permutedims(g.chip_max_x), declined, Int32(g.nodata.output))
    best > 0 || error("every chip-size maximum is declined or missing, so no pyramid level could run")
    return best
end

function _max_outside(cmax::AbstractMatrix{<:Integer}, declined::AbstractMatrix{Bool},
                      sentinel::Int32)
    best = 0
    for i in eachindex(cmax, declined)
        v = cmax[i]
        (declined[i] || v == sentinel) && continue
        best = max(best, Int(v))
    end
    return best
end

# `round(x / 2) * 2`, the reference's way of keeping a chip extent even.
_even_chip(x::Real) = round(Int, x / 2) * 2

"""
    chopped_size(gridsize, factor) -> Tuple{Int,Int}

`gridsize` truncated to a multiple of `factor` on both axes (`autoRIFT.py:767`).
"""
chopped_size(gridsize::Tuple{Int,Int}, factor::Integer) =
    (fld(gridsize[1], factor) * factor, fld(gridsize[2], factor) * factor)

"""
    pad_to_orig(A, origsize, fill) -> Matrix

`A` pasted into the top-left corner of an `origsize` array of `fill` (`testautoRIFT.py:792-805`).

The reference allocates every output at the pre-chop size and writes the correlated block into it, so
the strip the chop dropped reaches the product as nodata. Padding it here is what makes the two
products the same shape, and what makes the dropped strip count the same way in
`ItsLiveOffsetProduction.roi_valid_percentage`.
"""
function pad_to_orig(A::AbstractMatrix, origsize::Tuple{Int,Int}, filler)
    out = fill(convert(eltype(A), filler), origsize)
    out[axes(A, 1), axes(A, 2)] .= A
    return out
end

# The search limit the driver gives the correlator, over the whole window: the geogrid band with the
# sentinel and the declined points zeroed, and the chopped-away strip zeroed too. This is
# `roi_valid_percentage`'s denominator.
function driver_search_limit(s::Setup, declined::AbstractMatrix{Bool}, chop::Tuple{Int,Int})
    g = s.geometry
    return _search_limit(permutedims(g.search_x), declined, Int32(g.nodata.output), chop)
end

function _search_limit(sx::AbstractMatrix{<:Integer}, declined::AbstractMatrix{Bool},
                       sentinel::Int32, chop::Tuple{Int,Int})
    out = zeros(Int32, size(sx))
    Threads.@threads for j in 1:chop[2]
        for i in 1:chop[1]
            v = sx[i, j]
            (declined[i, j] || v == sentinel || v <= 0) && continue
            out[i, j] = Int32(v)
        end
    end
    return out
end

# ---------------------------------------------------------------------------
# Radar imagery (metadata/georef/swath-bias live in packaging.jl, shared with e2e_run.jl)
# ---------------------------------------------------------------------------

"""
    radar_imagery(s::Setup; cache_io = :pread) -> (early, late, declined, declined_shifted, filter)

[`native_imagery`](@ref) for a radar pair: the correlator's two images, unfiltered by the correlator's
own filter, and the driver's no-data mask.

**The images are lazy and cached a tile at a time.** They are [`e2e_imagery`](@ref)'s — the reference
mosaic or RSLC amplitude and the secondary resampled onto it a window at a time (`ResampledMosaic`,
`ResampledRSLC`), or a GSLC pair's amplitude windows into the two granules — each behind a
[`TileCache`](@ref), so a tile is resampled or decoded once however many passes read it. No scene-sized
array is held in memory: a Sentinel-1 mosaic is up to 24043 x 67860 pixels, an RSLC 57760 x 50511 and a
GSLC overlap 54885 x 110085. The correlator applies `filter` a block at a time — `native_run` passes it
as `preprocess` and correlates blocked.

**No offset.** An optical pair's two scenes are the raw, uncropped granules, so `native_imagery` crops
each to the overlap at its own coregistration offset. A radar pair's secondary is instead resampled
onto the reference's own mosaic or RSLC grid, so both arrays share the geometry's grid and the location
bands apply to them with no shift.

**Not quantized to `UInt8`**, for the reason `native_imagery` gives. NISAR-L1 (RSLC) does not share L2's
`UInt8`-before-highpass defect either: `loadProduct` hard-casts NISAR to `UInt8` before filtering, which
saturates the high-pass response's negative half and leaves ~30% of an L2 capture's `in_I1` at exactly
128, but only 2.4% of an L1 capture's `in_I1` and 2.9% of `in_I2` sit there. So the plain `Float32`
highpass is the right one for every radar case.
"""
function radar_imagery(s::Setup; cache_io::Symbol = :pread)
    c = s.case
    m = correlator_filter(c)
    m isa AutoRIFT.Highpass || error(
        "no route for radar filter $m on $(basename(c.product)); every radar case reaches " *
        "AutoRIFT.Highpass through correlator_filter today")
    (ref, sec) = e2e_imagery(s; lazy = true)
    size(ref) == size(sec) || error("$(basename(c.product)): the reference image is $(size(ref)) and " *
                                    "the secondary $(size(sec)); they are not on the same grid")
    g = s.geometry
    lx, ly = permutedims(g.location_x), permutedims(g.location_y)
    sentinel = Int32(g.nodata.output)
    declined = _sentinel_declined(lx, ly, sentinel)
    shifted = copy(declined)
    images = map((ref, sec)) do img
        tc = scratch!(TileCache(img; tile = 1024, dir = joinpath(CACHE, "scratch"), io = cache_io))
        fill_tiles!(tc)
        _accumulate_zeros_banded!((declined => NODATA_SAMPLE_OFFSET, shifted => NODATA_SAMPLE_OFFSET + 1),
                                  tc, lx, ly, sentinel)
        tc
    end
    return (images[1], images[2], declined, shifted, m)
end

# ---------------------------------------------------------------------------
# The whole chain
# ---------------------------------------------------------------------------

struct NativeRun
    case::String
    stages::Vector{Pair{String,Float64}}
    # The highest sampled footprint inside each stage, bytes; empty when the run was not traced.
    stage_peaks::Vector{Pair{String,Int}}
    peak::Int
    floor::Int
    scene::Tuple{Int,Int}
    origsize::Tuple{Int,Int}
    chop::Tuple{Int,Int}
    chop_factor::Int
    chip_max::Int
    declined::Int
    searched::Int
    resolved::Int
    roi_valid::Float64
    product::String
end

total_seconds(r::NativeRun) = sum(last, r.stages)

"""
    native_run(c::GoldenCase; warm, trace, cache_io = :pread, filter_tile = 512, outname)
        -> (NativeRun, NamedTuple)

The chain on one case. The second return value carries the intermediate arrays, for
`--compare-intermediate`.

`cache_io` is how the scratch tile caches are accessed (`TileCache`'s `io`, `AutoRIFT`'s
`filter_cache_io`) and `filter_tile` the side of `AutoRIFT.FilterTileCache`'s tiles. `:pread` keeps the
cached pages out of the process's resident size, which a mapping counts in full: measured on a Sentinel-1
burst pair, 15 GiB resident against 32 GiB mapped. 512 px filtered tiles keep it as fast as the mapping,
since a block's read window then spans little more of each tile than it needs; at 1024 px the same run
correlates 14% slower.
"""
function native_run(c::GoldenCase; warm::Bool = false, trace::Bool = true, cache_io::Symbol = :pread,
                    filter_tile::Integer = 512, outname::AbstractString = "julia_native.nc")
    # Every scene is read from its source inside the timed region, as the reference reads it; an entry
    # left in `STAGED` by another script would re-route a scene to a local copy.
    empty!(STAGED)
    buf = zeros(UInt64, 64)
    stages = Pair{String,Float64}[]
    marks = UInt64[]            # the counter at each stage's end, for the per-stage peak
    box = Ref{Any}(nothing)

    work = function ()
        local s, p, r, b1, b2, filt, declined, shifted, search, grid, pts
        local cmax, scale, factor, orig, chop, chipx, roi, path
        empty!(marks)
        push!(marks, cpu_ticks())

        t = @elapsed s = setup(c; lazy_warp = true)
        push!(stages, "geometry" => t)
        push!(marks, cpu_ticks())

        # The no-data mask is accumulated inside the imagery pass, where the unfiltered values live.
        t = @elapsed (b1, b2, declined, shifted, filt) = radar_sensor(c) ? radar_imagery(s; cache_io) :
                                                         native_imagery(s; cache_io)
        push!(stages, "imagery" => t)
        push!(marks, cpu_ticks())

        t = @elapsed begin
            cmax = driver_chip_size_max(s, declined)
            scale = AutoRIFT.chip_size_scale(s.geometry)
            # `chip_size_max` overrides what the extension derives from the geometry alone: a caller's
            # keyword wins over the inner default, which is how the driver's lower, imagery-aware
            # maximum reaches `Params`.
            # Every pair's images arrive unfiltered (`native_imagery`, `radar_imagery`), so the
            # correlator filters them; `nothing` is a pair the reference hands the correlator with
            # only its native filter applied.
            kw = (; threaded = Threads.nthreads() > 1,
                  preprocess = filt === nothing ? :none : filt,
                  chip_size_max = (X = cmax, Y = _even_chip(cmax * scale)))
            p = AutoRIFT.params(s.geometry; kw...)
            factor = cmax ÷ p.chip_size_min.X
            # The subpixel ladder is keyed by chip size from the smallest up (`testautoRIFT.py:330-334`),
            # so a lower maximum drops its coarsest entries: one per level the driver's maximum leaves.
            nlevels = trailing_zeros(factor) + 1
            length(p.subpixel) > nlevels &&
                (p = AutoRIFT.params(s.geometry; kw..., subpixel = p.subpixel[1:nlevels]))
            pts = AutoRIFT.pointset(s.geometry;
                                    pixel_size = ImagePairGeometry.xsize(s.pair.coordinate))
            orig = size(pts.x)
            chop = chopped_size(orig, factor)
            # The grown mask gates the outputs; the narrow one gates the search. Both are the
            # driver's, and conflating them moves the product's coverage.
            grown = dilated_nodata_mask(declined)
            search = driver_search_limit(s, grown, chop)
            # Chop first, then decline: both are the driver's own, and a point outside the chop is not
            # correlated at all.
            grid = _chop_to(pts, chop[1], chop[2])
            keep = .!view(declined, 1:chop[1], 1:chop[2])
            grid = AutoRIFT.rebuild(grid; radius_x = grid.radius_x .* keep,
                                    radius_y = grid.radius_y .* keep)
        end
        push!(stages, "grid" => t)
        push!(marks, cpu_ticks())

        t = @elapsed r = begin
            # Blocked, so a block's read window is all a task holds, sized to the 512² chunks the
            # granules are stored in. `filter_tile` is `AutoRIFT.FilterTileCache`'s tile, which keeps
            # each filtered tile across the passes; that cache holds only a `Highpass` result, so any
            # other filter is applied afresh to every block.
            b = AutoRIFT.block_size_for(grid, p, size(b1); chunk = (512, 512))
            h = AutoRIFT.halo(grid, p, size(b1))
            # A block's read window is the block grown by the halo on every side.
            @printf("    blocks %d x %d (X x Y), halo %d x %d, read window %d x %d\n", b.X, b.Y, h.X, h.Y,
                    b.X + 2h.X, b.Y + 2h.Y)
            fct = filt isa AutoRIFT.Highpass ? filter_tile : nothing
            # `wallis_fill` takes its gaps from the pixels reading zero (`np.isclose(image, 0)`), which
            # AutoRIFT reads from the validity mask rather than from the values.
            pair = filt isa AutoRIFT.WallisGapfill ?
                   AutoRIFT.ImagePair(b2, b1; reference_valid = NonzeroMask(b2),
                                      secondary_valid = NonzeroMask(b1)) :
                   AutoRIFT.ImagePair(b2, b1)
            try
                AutoRIFT.autorift(pair, grid, p, (b.X, b.Y), 0, fct; filter_cache_io = cache_io)
            finally
                close_scratch!()
            end
        end
        push!(stages, "correlate" => t)
        push!(marks, cpu_ticks())

        t = @elapsed begin
            # Pasted into the pre-chop grid, then dropped wherever the grown mask or a zero search
            # limit says the driver reports nothing (`testautoRIFT.py:806-818`). The second condition
            # is not implied by the first: a point can carry no search range of its own.
            dx = pad_to_orig(r.dx, orig, NaN32)
            dy = dy_sign(c.platform) .* pad_to_orig(r.dy, orig, NaN32)
            chipx = pad_to_orig(r.chip_size, orig, 0)
            interp = pad_to_orig(Matrix{Bool}(r.interpolated), orig, false)
            drop = grown .| (search .== 0)
            dx[drop] .= NaN32
            dy[drop] .= NaN32
            chipx[drop] .= 0
            interp[drop] .= false

            roi = ItsLiveOffsetProduction.roi_valid_percentage(chipx, search)
            coeffs = ItsLiveOffsetProduction.coefficients(s.geometry)
            refv, ssm = ItsLiveOffsetProduction.reference_velocity(s.inputs)
            info = radar_sensor(c) ? radar_img_pair_info(c, s, roi) : native_img_pair_info(c, s, roi)
            # `radar_pair(c)` picks the 16-variable `:radar` schema; NISAR-L2 is a radar sensor
            # (`radar_sensor(c)`, above) but writes the plain 12-variable `:optical` one — see
            # `radar_pair`'s docstring.
            radar = radar_pair(c)
            pair_type = radar ? :radar : :optical
            coordinates = radar ? "radar, map" : MOTION_COORDINATES
            dt_seconds = radar ? s.pair.dt : nothing
            georef = radar ? radar_georef(s) : native_georef(s)
            swath_bias = radar ? s1_swath_bias(c, s, p) : nothing
            input = ItsLiveOffsetProduction.ItsLiveInput(
                pair_type, DETECTION_METHOD, coordinates, AUTORIFT_VERSION,
                replace(PARAMETER_SHAPEFILE, "/vsicurl/" => ""), product_source(info), dt_seconds,
                georef, info, coeffs, refv, ssm, swath_bias,
                dx, dy, chipx, scale, interp)
            path = joinpath(run_dir(c, RUN_FOR_OUTPUT), outname)
            mkpath(dirname(path))
            ItsLiveOffsetProduction.write_product(path, input)
        end
        push!(stages, "package" => t)
        push!(marks, cpu_ticks())

        box[] = (; s, p, r, declined, shifted, grown, search, orig, chop, factor, cmax, scale,
                 roi, chipx, scene = size(b1), path, filt)
        return nothing
    end

    if warm
        try
            work()
        finally
            close_scratch!()
        end
        empty!(stages)
        GC.gc(true); GC.gc(true)
    end
    floor_bytes = last(rusage!(buf))
    stage_peaks = Pair{String,Int}[]
    peak = try
        if trace
            _, tr, _ = with_trace(; interval = 0.05) do
                work()
            end
            for (k, (name, _)) in enumerate(stages)
                inside = (marks[k] .<= tr.tick .< marks[k + 1])
                push!(stage_peaks, name => any(inside) ? Int(maximum(tr.footprint[inside])) : 0)
            end
            maximum(tr.footprint)
        else
            work()
            last(rusage!(buf))
        end
    finally
        close_scratch!()
    end

    k = box[]
    r = NativeRun(c.product, stages, stage_peaks, Int(peak), Int(floor_bytes), k.scene, k.orig, k.chop,
                  k.factor, k.cmax, count(k.declined), count(!iszero, k.search),
                  count(!iszero, k.chipx), k.roi, k.path)
    return (r, k)
end

# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

function report(r::NativeRun)
    @printf("  %-56s\n", first(r.case, 56))
    peaks = Dict(r.stage_peaks)
    for (name, secs) in r.stages
        @printf("    %-12s %8.1f s", name, secs)
        haskey(peaks, name) && (peaks[name] > 0 ? @printf("   peak %6.2f GiB", peaks[name] / 2^30) :
                                print("   peak not sampled"))
        println()
    end
    @printf("    %-12s %8.1f s\n", "total", total_seconds(r))
    @printf("    scene %d x %d, grid %d x %d chopped to %d x %d (factor %d, chip max %d)\n",
            r.scene..., r.origsize..., r.chop..., r.chop_factor, r.chip_max)
    @printf("    declined %d, searched %d, resolved %d, ItsLiveOffsetProduction.roi_valid_percentage %.1f\n",
            r.declined, r.searched, r.resolved, r.roi_valid)
    @printf("    peak %.2f GiB (%.2f above the floor %.2f)\n", r.peak / 2^30,
            (r.peak - r.floor) / 2^30, r.floor / 2^30)
    @printf("    wrote %s\n", r.product)
    return flush(stdout)
end

# The reference's own pre-shift arrays, on the grids this script produces them on. `(y, x)` in the file
# is `(x, y)` out of NCDatasets, so every band transposes once here.
function read_intermediate(dir::AbstractString)
    path = joinpath(dir, "autoRIFT_intermediate.nc")
    isfile(path) || error("no autoRIFT_intermediate.nc in $dir; a container run writes it")
    return NCDataset(path) do ds
        (; dx = permutedims(Array(ds["Dx"].var[:, :])),
         dy = permutedims(Array(ds["Dy"].var[:, :])),
         chip = permutedims(Array(ds["ChipSizeX"].var[:, :])),
         interp = permutedims(Array(ds["InterpMask"].var[:, :])) .!= 0,
         search = permutedims(Array(ds["SearchLimitX"].var[:, :])),
         nodata = permutedims(Array(ds["noDataMask"].var[:, :])) .!= 0,
         scale = Float64(ds["ScaleChipSizeY"].var[]),
         orig = (Int(ds["origSizeY"].var[]), Int(ds["origSizeX"].var[])))
    end
end

function compare_intermediate(k, dir::AbstractString)
    ref = read_intermediate(dir)
    println("\n=== against the reference's own autoRIFT_intermediate.nc ===")
    @printf("  origSize      julia %s   reference %s   %s\n", k.orig, ref.orig,
            k.orig == ref.orig ? "equal" : "DIFFER")
    @printf("  chopped       julia %s   reference %s   %s\n", k.chop, size(ref.dx),
            k.chop == size(ref.dx) ? "equal" : "DIFFER")
    @printf("  ScaleChipSizeY julia %.6f  reference %.6f  %s\n", k.scale, ref.scale,
            k.scale ≈ ref.scale ? "equal" : "DIFFER")

    # The no-data mask, both sampling offsets, against the reference's own band. This is what decides
    # `NODATA_SAMPLE_OFFSET` rather than a reading of the indexing.
    # The reference's stored band is the *grown* mask, so that is what this scores. The narrow one is
    # reported beside it to show the dilation is what accounts for the difference.
    # Both offsets were accumulated during the imagery pass, so scoring them costs no imagery here.
    println("\n  no-data mask against the reference's own ($(count(ref.nodata)) set):")
    best, bestoff = -1.0, nothing
    for (off, narrow) in ((NODATA_SAMPLE_OFFSET, k.declined), (NODATA_SAMPLE_OFFSET + 1, k.shifted))
        for (label, m) in (("narrow", narrow), ("grown", dilated_nodata_mask(narrow)))
            agree = count(m .== ref.nodata) / length(m)
            @printf("    offset %+d %-7s set %8d  agree %.6f%%%s\n", off, label, count(m), 100agree,
                    (off == NODATA_SAMPLE_OFFSET && label == "grown") ? "   <= gates the outputs" : "")
            label == "grown" && agree > best && ((best, bestoff) = (agree, off))
        end
    end
    bestoff == NODATA_SAMPLE_OFFSET || @warn "a different sampling offset agrees better; \
        NODATA_SAMPLE_OFFSET should be $bestoff" in_use=NODATA_SAMPLE_OFFSET

    nr, nc = min.(k.chop, size(ref.dx))
    cut(A) = view(A, 1:nr, 1:nc)
    chip = cut(k.r.chip_size)
    rchip = cut(ref.chip)
    @printf("\n  ChipSizeX     %d of %d equal (%.4f%%)\n", count(chip .== rchip), length(rchip),
            100count(chip .== rchip) / length(rchip))
    interp = cut(Matrix{Bool}(k.r.interpolated))
    @printf("  InterpMask    %d of %d equal (%.4f%%)\n", count(interp .== cut(ref.interp)),
            length(rchip), 100count(interp .== cut(ref.interp)) / length(rchip))

    step = 1 / AutoRIFT.subpixel_at(k.p, 1).upsampling
    for (axis, mine, theirs) in ((:dx, cut(k.r.dx), cut(ref.dx)), (:dy, cut(k.r.dy), cut(ref.dy)))
        for sgn in (1, -1)
            st = endpoint_stage("  $axis (sign $(sgn > 0 ? '+' : '-'))", "autoRIFT_intermediate.nc",
                                mine, sgn .* theirs, step; gated = false)
            println("  ", st.name, "  ", st.detail)
        end
    end
    return nothing
end

function report_product_diff(k, c::GoldenCase, dir::AbstractString; golden::Bool = true)
    mine = read_product(k.path)
    println("\n=== the native product against the reference's own, from $(basename(dir)) ===")
    d = compare_products(mine, read_product(run_product(dir)))
    show(stdout, d)
    println("\nagrees_on_data: ", agrees_on_data(d))
    if golden && have_golden(c)
        println("\n=== the native product against ASF's golden file ===")
        dg = compare_products(mine, read_product(c))
        show(stdout, dg)
        println("\nagrees_on_data: ", agrees_on_data(dg))
    end
    return d
end

# The heap size the collector aims to stay under, GiB. Without one the collector lets the heap grow well
# past what the run holds live. Measured on NISAR L1, the largest golden case, read from ASF: the footprint
# peaks at 14.6 GiB at 12 and 16.4 GiB at 14, against a 16 GiB budget; from local granules, 16.4 GiB with
# no target. 12 costs about 3% in runtime (858 s against 827-849 s from local granules), because the
# collector runs more often. `0` leaves the collector's own default.
const HEAP_HINT_GIB = 12.0

"""
    set_heap_hint!(gib)

Set the collector's heap-size target to `gib` GiB, as `julia --heap-size-hint` does at startup — unless the
process was started with that flag, which takes precedence. `0` sets nothing.
"""
function set_heap_hint!(gib::Real)
    gib > 0 || return nothing
    Base.JLOptions().heap_size_hint == 0 || return nothing
    ccall(:jl_gc_set_max_memory, Cvoid, (UInt64,), round(UInt64, gib * 2^30))
    return nothing
end

function main(args)
    isempty(args) && error("usage: julia_e2e.jl <product-fragment> [--compare-intermediate] " *
                           "[--no-golden] [--warm] [--no-trace] [--cache-io pread|mmap] [--filter-tile N] " *
                           "[--heap-hint-gib N] [--lattice]")
    h = findfirst(==("--heap-hint-gib"), args)
    set_heap_hint!(isnothing(h) ? HEAP_HINT_GIB : parse(Float64, args[h + 1]))
    c = only(cases(first(args)))
    FORCE_LATTICE[] = "--lattice" in args
    reason = unsupported(c)
    isnothing(reason) || error("$(basename(c.product)): $reason")
    i = findfirst(==("--cache-io"), args)
    cache_io = isnothing(i) ? :pread : Symbol(args[i + 1])
    j = findfirst(==("--filter-tile"), args)
    filter_tile = isnothing(j) ? 512 : parse(Int, args[j + 1])
    r, k = native_run(c; warm = "--warm" in args, trace = !("--no-trace" in args), cache_io, filter_tile)
    println("  scratch caches: ", cache_io, ", filtered tiles ", filter_tile, " px")
    report(r)
    dir = resolve_run(c)
    "--compare-intermediate" in args && compare_intermediate(k, dir)
    report_product_diff(k, c, dir; golden = !("--no-golden" in args))
    return nothing
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS)
