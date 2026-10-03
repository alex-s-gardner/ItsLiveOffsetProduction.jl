# From two granules to the ITS_LIVE product netCDF, with nothing taken from the reference.
#
#   julia --project=tools/golden -t 12 tools/golden/julia_e2e.jl LC08_L1TP_009011_20200703
#   julia --project=tools/golden -t 12 tools/golden/julia_e2e.jl LC08_L1TP_009011_20200703 --compare-intermediate
#   AWS_PROFILE=itslive julia --project=tools/golden -t 12 tools/golden/julia_e2e.jl \
#       LC08_L1TP_009011_20200703 --stream
#
# **`--stream` is the arm to time against the container.** `process.py` hands GDAL `/vsis3/` paths and
# reads the granule over the network inside its own timed region; staging the scenes locally first, which
# is the default here, measures less work than the reference does. Reading over `/vsis3` needs
# `AWS_PROFILE` to name credentials that can pay, since `s3://usgs-landsat` is requester-pays.
#
# **What was missing.** `e2e_run.jl` runs the chain from the granule to `dx`/`dy`, and
# `validate_itslive_write.jl` writes a product — but from the reference's own captured arrays. So no run
# has produced a product from the imagery alone, which is the claim this makes.
#
# Five things sit between those two scripts, all of them the ITS_LIVE driver's rather than the
# correlator's, and none reproduced by any rung of the ladder:
#
#   bytes       the overlap filtered and quantized to 256 levels, which is what the reference's
#               correlator is handed (`uniform_data_type`). `e2e_run.jl` correlates the filtered
#               `Float32` field instead — a different C++ template with a different answer, which
#               `stages.jl`'s header measures at up to 7 px on one scene.
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
# the chop, the bytes or the correlator; that mode can.

include("e2e_run.jl")      # `setup`, `_chop_to`, `_cropped_scene`, the memory trace, the tile cache
include("product.jl")
include("compare.jl")

using Dates, NCDatasets, OpticalDatasets
using ItsLiveOffsetProduction
import Downloads
import GeoFormatTypes as GFT
import JSON3
# `FastGeoProjections` rather than PROJ for the lon/lat transform: the harness already builds
# `grid_transform` with it, and the two agree to 1e-14 degrees on every CRS these cases use.
# PROJ is not loaded at all on this path: `grid_transform` pulls it in only for its `--proj-only`
# cross-check.
import FastGeoProjections as FGP
# Loaded, not called: `cf_grid_mapping` comes from `ItsLiveOffsetProduction`, whose trigger is
# Rasters plus DimensionalData plus ArchGDAL plus DiskArrays, and a product cannot be written without
# the CF grid-mapping attributes it builds.
using DimensionalData, Rasters

# ---------------------------------------------------------------------------
# What is declared rather than derived
# ---------------------------------------------------------------------------

# The reference implementation being reproduced. Every golden product's `source` attribute names these
# two versions, and the container they come from is `ghcr.io/asfhyp3/hyp3-autorift:0.28.4`; neither is a
# property of the imagery, so neither can be derived from it.
const PLUGIN_NAME = "hyp3_autorift"
const PLUGIN_VERSION = "0.28.4"
const AUTORIFT_VERSION = "2.1.1"

# `testautoRIFT.py:1363-1365`, for an optical pair.
const DETECTION_METHOD = "feature"
const MOTION_COORDINATES = "map"

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
    native_imagery(s::Setup) -> (bytes_early, bytes_late, declined, declined_shifted, filter)

The two byte images the correlator is handed, and the driver's no-data mask, in one pass over each scene.

**One pass, because the scene is read once and nothing full-size is copied.** The window is read straight
into a `Float32` buffer, so the `UInt16` scene is never materialized and there is no separate conversion
pass over 291 million pixels; the filter writes into a second buffer with
[`AutoRIFT.highpass!`](@ref); and both buffers are reused for the second scene rather than reallocated.
Four full-scene arrays live at once where the allocating form held eight.

**The mask is accumulated here because the unfiltered values are only available here.** A point is
declined where *either* scene reads zero (`testautoRIFT.py:337-349`), and zero converts exactly, so the
`Float32` buffer answers the same question the raw integers would. Accumulating per scene is what lets the
buffer be reused: the alternative holds both scenes' pixels at once purely to compare them.

Both sampling offsets are accumulated, so `--compare-intermediate` can score the choice without the
imagery being retained for it.

Quantizing is not optional. `uniform_data_type` collapses the filtered field onto 256 levels before
`runAutorift` sees it, so production reaches the reference's `UInt8` correlator template; correlating the
`Float32` field reaches a different one.
"""
function native_imagery(s::Setup)
    c = s.case
    m = correlator_filter(c)
    m isa AutoRIFT.Highpass || error(
        "this script covers the optical `hps` path, whose scene is filtered inside the correlator; " *
        "$(basename(c.product)) uses $m, which `process.py` applies to the native scene first and " *
        "which rung 5.3 is what reproduces")
    early, late = acquisition_order(c)
    # The offsets come from `coregister` in the job's reference/secondary order, and the arrays are
    # named by acquisition. Equal on every case but the two reversed jobs, and a silent mismatch there
    # would pair each scene with the other's offset.
    first(c.reference) == early || error(
        "$(basename(c.product)) is a reversed job — its reference acquisition is $(first(c.reference)) " *
        "but the earlier acquisition is $early — so the coregistration offsets and the scene paths " *
        "would be paired the wrong way round here")

    want = reverse(_coord_size(s.pair.coordinate))      # (ny, nx) of the overlap
    g = s.geometry
    lx, ly = permutedims(g.location_x), permutedims(g.location_y)
    sentinel = Int32(g.nodata.output)
    # Seeded from the geogrid's own out-of-image sentinel (`testautoRIFT.py:603`), before any imagery.
    declined = _sentinel_declined(lx, ly, sentinel)
    shifted = copy(declined)

    field = Matrix{Float32}(undef, want)       # the read destination, and the mask's view of the scene
    filtered = Matrix{Float32}(undef, want)    # `highpass!` needs a destination distinct from its input
    eroded = Matrix{Bool}(undef, want)         # `_filtered!`'s mask buffer, reused across both scenes
    keep = trues(want)
    bytes = Vector{Matrix{UInt8}}()

    for (path, off) in ((s.reference_path, s.pair.reference_offset),
                        (s.secondary_path, s.pair.secondary_offset))
        _scene_window!(field, path, off, want)
        _accumulate_zeros!(declined, field, lx, ly, sentinel, NODATA_SAMPLE_OFFSET)
        _accumulate_zeros!(shifted, field, lx, ly, sentinel, NODATA_SAMPLE_OFFSET + 1)
        AutoRIFT.highpass!(filtered, field, keep, m.width)
        AutoRIFT._filtered!(filtered, eroded, keep, m.width)
        # `keep`, not `eroded`: the reference quantizes over the whole array, and `bytescale`'s mean and
        # standard deviation are what the mask restricts.
        push!(bytes, AutoRIFT.bytescale(filtered, keep))
    end
    return (bytes[1], bytes[2], declined, shifted, m)
end

"""
    _scene_window!(dest, path, off, want) -> dest

The overlap out of one scene, read concurrently into `dest`, in `(row, col)`.

`OpticalDatasets.read_window!` reads chunk-aligned slabs from several tasks at once and converts each as
it lands, where `_cropped_scene` reads the **whole band** in one blocking request and then copies it
twice. On a 338 MiB Landsat band over `/vsis3` that is 32.3 s against 54.9 s, and it writes one array
rather than allocating four.

`off` is `(x, y)`, which is `coregister`'s order; the raster is indexed `(row, col)`, so the two swap
here.
"""
function _scene_window!(dest::AbstractMatrix, path::AbstractString, off::NTuple{2,Int},
                        want::Tuple{Int,Int})
    nr, nc = want
    ox, oy = off
    r = OpticalDatasets.open_optical(path)
    try
        size(r, 1) >= oy + nr && size(r, 2) >= ox + nc || throw(DimensionMismatch(
            "the overlap at offset $off does not fit in a $(size(r)) scene; the crop and the " *
            "geotransform disagree"))
        return OpticalDatasets.read_window!(dest, r, (oy + 1):(oy + nr), (ox + 1):(ox + nc))
    finally
        close(r)
    end
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

# One scene's contribution: declined wherever it reads zero at a point, or the point samples outside the
# overlap. Accumulated with `|=` because either scene reading zero is enough.
function _accumulate_zeros!(out::AbstractMatrix{Bool}, scene::AbstractMatrix,
                            lx::AbstractMatrix{<:Integer}, ly::AbstractMatrix{<:Integer},
                            sentinel::Int32, offset::Int)
    nr, nc = size(scene)
    Threads.@threads for j in axes(lx, 2)
        for i in axes(lx, 1)
            out[i, j] && continue
            x, y = lx[i, j], ly[i, j]
            (x == sentinel || y == sentinel) && continue
            r, c = Int(y) + offset, Int(x) + offset
            if !(1 <= r <= nr && 1 <= c <= nc)
                out[i, j] = true
                continue
            end
            iszero(scene[r, c]) && (out[i, j] = true)
        end
    end
    return out
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
# Metadata, from the scenes' own STAC items
# ---------------------------------------------------------------------------

"""
    scene_identification(name) -> OpticalDatasets.Identification

The scene's STAC item, read and parsed into an `Identification`.

The acquisition *time of day* is the one field a Landsat product id does not carry, and the reference
takes it from the same place — `process.py:63-74` reads this item and parses `properties.datetime`. The
item is anonymous, unlike the imagery it points at.

`DateTime` is millisecond precision against the item's microseconds, so the time is truncated on the
way in; `ItsLiveOffsetProduction.ImagePairInfo`'s docstring records that the product's date attributes inherit that.
"""
function scene_identification(name::AbstractString)
    # A Sentinel-2 product name carries its own acquisition time to the second, so it needs no item at
    # all. The Landsat STAC call exists for the one field a Landsat id does not carry — the time of day
    # — so asking the Landsat collection for an `S2...` name is a 404, not a slower answer.
    startswith(name, "S2") && return sentinel2_identification(name)
    url = "https://landsatlook.usgs.gov/stac-server/collections/landsat-c2l1/items/$name"
    body = String(take!(Downloads.download(url, IOBuffer(); timeout = 60)))
    return open_optical(JSON3.read(body, Dict{String,Any}))
end

"""
    native_img_pair_info(c, s, roi_valid) -> ItsLiveOffsetProduction.ImagePairInfo

`IMG_INFO_DICT` for an optical pair, from the two product ids and their STAC items.

`path`, `row` and `collection_number` are written as `Float64` because the reference's own attributes
are doubles: it carries them through NumPy rather than as Python ints.
"""
function native_img_pair_info(c::GoldenCase, s::Setup, roi_valid::Real)
    early, late = acquisition_order(c)
    # Concurrently: these are two independent round trips to the same host, and run one after the
    # other they cost 0.47 s against 0.22 s together. A Sentinel-2 name needs no request at all, so
    # for an S2 pair both tasks are pure parsing and the spawn costs nothing. At `-t 1` this still
    # runs, just sequentially.
    t_early = Threads.@spawn scene_identification(early)
    t_late = Threads.@spawn scene_identification(late)
    id1, id2 = fetch(t_early), fetch(t_late)
    lon, lat = pair_centroid(s.pair.coordinate, s.epsg)
    extra = Dict{String,Any}(
        "id_img1" => id1.id, "id_img2" => id2.id,
        "sensor_img1" => id1.sensor, "sensor_img2" => id2.sensor,
        "correction_level_img1" => id1.correction_level,
        "correction_level_img2" => id2.correction_level)
    # Path, row, collection and processing date are Landsat's own. A Sentinel-2 product's
    # `img_pair_info` carries none of them — checked against the golden S2 products, whose variable has
    # exactly the six attributes above plus the ones `write_product` derives — and
    # `sentinel2_identification` leaves them `nothing`, so they are added only where they exist.
    for (n, id) in ((1, id1), (2, id2))
        id.path === nothing && continue
        extra["path_img$n"] = Float64(id.path)
        extra["row_img$n"] = Float64(id.row)
        extra["collection_number_img$n"] = Float64(id.collection_number)
        extra["collection_category_img$n"] = id.collection_category
        extra["processing_date_img$n"] = id.processing_date
    end
    return ItsLiveOffsetProduction.ImagePairInfo(id1.acquisition_time, id2.acquisition_time,
                                  id1.mission, id2.mission, id1.satellite, id2.satellite,
                                  Float64(roi_valid), round(lat; digits = 2), round(lon; digits = 2),
                                  extra)
end

# `netcdf_output.py:427,438`. The year is the wall clock's, as it is there.
function product_source(info::ItsLiveOffsetProduction.ImagePairInfo)
    s = "NASA MEaSUREs ITS_LIVE project. Processed by ASF DAAC HyP3 $(Dates.year(Dates.now())) " *
        "using the $PLUGIN_NAME plugin version $PLUGIN_VERSION running autoRIFT version " *
        AUTORIFT_VERSION
    startswith(info.mission_img1, "L") || return s
    return s * ". Landsat-$(info.satellite_img1) images courtesy of the U.S. Geological Survey"
end

"""
    native_georef(s::Setup) -> ItsLiveOffsetProduction.ItsLiveGeoref

The product grid, from the geogrid window's own geotransform.

`x`/`y` are **cell centres**, where a geotransform names the outer edge, so each gains half a pixel.
`pixel_size_x`/`pixel_size_y` are the *scene's* pixel size rather than the grid's — they scale a chip
size in pixels into metres, and a chip is cut in the imagery.
"""
function native_georef(s::Setup)
    gt = ImagePairGeometry.window_geotransform(s.grid, s.window)
    nx, ny = size(s.window)
    x = collect(gt[1] + gt[2] / 2 .+ (0:(nx - 1)) .* gt[2])
    y = collect(gt[4] + gt[6] / 2 .+ (0:(ny - 1)) .* gt[6])
    to_lonlat = FGP.Transformation(FGP.EPSG(Int(s.info.epsg)), FGP.EPSG(4326);
                                   always_xy = true)
    px = ImagePairGeometry.xsize(s.pair.coordinate)
    py = abs(s.pair.coordinate.spacing[2])
    return ItsLiveOffsetProduction.ItsLiveGeoref(x, y, ItsLiveOffsetProduction.cf_grid_mapping(GFT.EPSG(s.info.epsg)),
                                  (a, b) -> (to_lonlat(a, b)...,), px, py)
end

# ---------------------------------------------------------------------------
# Radar imagery and metadata
# ---------------------------------------------------------------------------

# Whether `c` reads its two scenes as a SAR acquisition (mosaic, RSLC or GSLC) rather than through
# `OpticalDatasets` — every platform `e2e_imagery` has a branch for.
radar_sensor(c::GoldenCase) = startswith(c.platform, "S1") || startswith(c.platform, "NISAR")

# Whether `c` writes the 16-variable `:radar` schema. NISAR-L2 is a radar sensor but a geocoded,
# map-grid product — `product.jl`'s `schema` calls it `:optical` for exactly that reason, and
# `setup` already routes it through the same `coregister`/`ProjectedCoordinate` path as Landsat and
# Sentinel-2 rather than through `_radar_setup`.
radar_pair(c::GoldenCase) = startswith(c.platform, "S1") || c.platform == "NISAR-L1"

"""
    radar_imagery(s::Setup) -> (bytes_early, bytes_late, declined, declined_shifted, filter)

[`native_imagery`](@ref) for a radar pair: the correlator's two byte images and the driver's no-data
mask, sourced from [`e2e_imagery`](@ref) (the same mosaic/RSLC/GSLC reader `e2e_run.jl` uses, already
validated stage by stage against the reference) rather than a GDAL scene file.

**No offset.** An optical pair's two scenes are the raw, uncropped granules, so `native_imagery` crops
each to the overlap at its own coregistration offset. A radar pair's secondary is instead resampled
onto the reference's own mosaic or RSLC grid (`ResampledMosaic`/`ResampledRSLC`), so both arrays
`e2e_imagery` returns already share the geometry's grid and the location bands apply to them with no
shift.

NISAR-L2 (GSLC) has no route here: its band is ~13 billion pixels, and reproducing the reference's
whole-array `uniform_data_type` quantization over that needs a memory-bounded, tile-at-a-time pass
nothing in this package implements yet. Erroring here is deliberate rather than attempting an
allocation sized for it.

**NISAR-L1 (RSLC) does not share L2's `UInt8`-before-highpass defect.** `rung_bytes`'s skip message
for NISAR-L2 names the symptom: `loadProduct` hard-casts NISAR to `UInt8` before filtering, so the
high-pass response's entire negative half saturates to zero and ~30% of the reference's own
`capture/in_I1` sits at exactly 128. Checked directly against a real NISAR-L1 capture rather than
assumed: only 2.4% of `in_I1` and 2.9% of `in_I2` sit at 128 — nowhere near that signature — so this
function's plain `Float32` highpass, the same one the optical and Sentinel-1 paths use, is the right
one for L1 and needs no `UInt8`-saturating variant.
"""
function radar_imagery(s::Setup)
    c = s.case
    c.platform == "NISAR-L2" && error(
        "$(basename(c.product)): no memory-bounded highpass+quantize pass exists yet for a GSLC's " *
        "~13-billion-pixel band; see radar_imagery's docstring")
    m = correlator_filter(c)
    m isa AutoRIFT.Highpass || error(
        "no route for radar filter $m on $(basename(c.product)); every radar case reaches " *
        "AutoRIFT.Highpass through correlator_filter today")
    (ref, sec) = e2e_imagery(s)
    want = size(ref)
    size(sec) == want || error("$(basename(c.product)): the reference image is $want and the " *
                               "secondary $(size(sec)); they are not on the same grid")

    g = s.geometry
    lx, ly = permutedims(g.location_x), permutedims(g.location_y)
    sentinel = Int32(g.nodata.output)
    declined = _sentinel_declined(lx, ly, sentinel)
    shifted = copy(declined)

    filtered = Matrix{Float32}(undef, want)
    eroded = Matrix{Bool}(undef, want)
    keep = trues(want)
    bytes = Vector{Matrix{UInt8}}()

    # Each side is read once through its own range-indexed `getindex` — a plain copy for the dense S1
    # mosaic, a block-at-a-time resample or HDF5 read for the lazy NISAR-L1/`ResampledMosaic` types —
    # rather than scalar-indexed, which for the lazy types would resample or decode one pixel at a time.
    for img in (ref, sec)
        field = img[axes(img, 1), axes(img, 2)]
        _accumulate_zeros!(declined, field, lx, ly, sentinel, NODATA_SAMPLE_OFFSET)
        _accumulate_zeros!(shifted, field, lx, ly, sentinel, NODATA_SAMPLE_OFFSET + 1)
        AutoRIFT.highpass!(filtered, field, keep, m.width)
        AutoRIFT._filtered!(filtered, eroded, keep, m.width)
        push!(bytes, AutoRIFT.bytescale(filtered, keep))
    end
    return (bytes[1], bytes[2], declined, shifted, m)
end

# The pair's footprint midpoint in geodetic degrees. Mirrors `_radar_setup`'s own inline computation
# (`e2e.jl`): a radar footprint is solved for with `rdr2geo` rather than carried by a projection, so
# `pair_centroid` (which expects a `ProjectedCoordinate`) does not apply.
function radar_pair_centroid(s::Setup)
    b = footprint_bounds(IdentityTransform(), s.pair.coordinate)
    return ((b.X[1] + b.X[2]) / 2, (b.Y[1] + b.Y[2]) / 2)
end

"""
    radar_georef(s::Setup) -> ItsLiveOffsetProduction.ItsLiveGeoref

[`native_georef`](@ref) for a `RadarCoordinate` pair (S1 or NISAR-L1).

The grid itself (`x`/`y`, the CF mapping) is built identically to the optical path — `s.grid`/
`s.window` already come from the same `parameter_grid`/`grid_window` calls regardless of platform —
only `pixel_size_x` differs, since `s.pair.coordinate` carries no `.spacing` to read a scene pixel
size from.

`pixel_size_x` is the **ground**-range sample spacing, `coord.dr / sin(coord.incidence_angle)` —
slant range projected through the incidence angle at the scene centre, matching the reference's own
`rangePixelSize` to 6 significant figures on both a measured Sentinel-1 and a NISAR-L1 case (checked
against `capture_packaging`'s captured scalar, not assumed). `pixel_size_y` is unused for `:radar`:
`write_product` reads `input.georef.pixel_size_y` nowhere — the `azimuth_pixel_size` attribute and the
`vr`/`va` conversion both use a `pixel_size_y` recomputed from `offset2va`/`dt_seconds` instead (see
`ItsLiveGeoref`'s docstring) — so it is set equal to `pixel_size_x` rather than carrying a second,
unverified formula for a value nothing reads.
"""
function radar_georef(s::Setup)
    gt = ImagePairGeometry.window_geotransform(s.grid, s.window)
    nx, ny = size(s.window)
    x = collect(gt[1] + gt[2] / 2 .+ (0:(nx - 1)) .* gt[2])
    y = collect(gt[4] + gt[6] / 2 .+ (0:(ny - 1)) .* gt[6])
    to_lonlat = FGP.Transformation(FGP.EPSG(Int(s.info.epsg)), FGP.EPSG(4326);
                                   always_xy = true)
    coord = s.pair.coordinate
    px = coord.dr / sin(coord.incidence_angle)
    return ItsLiveOffsetProduction.ItsLiveGeoref(x, y, ItsLiveOffsetProduction.cf_grid_mapping(GFT.EPSG(s.info.epsg)),
                                  (a, b) -> (to_lonlat(a, b)...,), px, px)
end

"""
    s1_swath_bias(c, s, p) -> Union{ItsLiveOffsetProduction.SwathOffsetBias,Nothing}

The Sentinel-1 subswath-offset-bias correction's own inputs, or `nothing` where the correction has no
meaning — which is most of this manifest.

`nothing` for every platform but `S1-SLC`: a burst job opens one or two subswaths
([`burst_swaths`](@ref)), which `SLCDatasets.subswath_borders` refuses outright, matching
`ItsLiveOffsetProduction.SwathOffsetBias`'s own docstring ("or for a pair whose product was not opened
for all three subswaths"). For an `S1-SLC` pair it is `nothing` whenever the two acquisitions are the
same spacecraft — `subswath_borders`' own `same_platform` — since the correction removes a difference
between two *different* antennas' patterns and has no meaning within one. Of this manifest's 5
full-SLC cases, 4 are same-platform and only the fifth (`S1A` paired with `S1B`) gets a real one.

`p` is `work()`'s own `AutoRIFT.params(s.geometry; ...)`, already computed by the time this is called;
`grid_spacing_x` is its `grid_spacing.X`, which the correction's coarse-grid smoothing step reads.
"""
function s1_swath_bias(c::GoldenCase, s::Setup, p)
    c.platform == "S1-SLC" || return nothing
    rp, sp = _s1_products(c, s.run)
    bd = SLCDatasets.subswath_borders(rp, sp)
    bd.same_platform && return nothing
    lx, ly = image_location(s.geometry)
    return ItsLiveOffsetProduction.SwathOffsetBias(lx, ly, Float64(p.grid_spacing.X), bd.ncols,
                                                   Float64(bd.border12), Float64(bd.border23),
                                                   bd.same_platform, bd.reference_platform)
end

# A `UtcTime`'s instant as a `DateTime`, truncated to the millisecond `ImagePairInfo.acquisition_date_img1`
# carries (see its docstring: the reference's own microsecond precision is lost the same way on the
# optical path).
_s1_datetime(t) = t.datetime + Dates.Millisecond(round(Int, t.seconds * 1000))

# Splits a Sentinel-1 product id (or `burst2safe`'s synthesized SAFE name, same layout —
# `S1A_IW_SLC__1SSH_<start>_<stop>_<orbit>_<datatake>_<productid>`, ten fields because `SLC__1SSH`
# contributes an empty one) into the fields `img_pair_info`'s `extra` carries verbatim. Measured
# against a real `capture_packaging` capture: field 1's tail is `satellite_img` ("1A" out of "S1A"),
# and fields 8/9/10 are `absolute_orbit_number`/`mission_data_take_ID`/`product_unique_ID` exactly.
function _s1_name_fields(name::AbstractString)
    parts = split(name, '_')
    length(parts) == 10 || throw(ArgumentError(
        "\"$name\" does not split into a Sentinel-1 product id's 10 underscore-delimited fields"))
    return (; satellite = parts[1][2:end], orbit = parts[8], datatake = parts[9], product_id = parts[10])
end

"""
    s1_img_pair_info(c, s, roi_valid) -> ItsLiveOffsetProduction.ImagePairInfo

`IMG_INFO_DICT` for a Sentinel-1 pair (full-SLC or burst), from the parsed products
[`_s1_products`](@ref) already builds for the geometry and mosaic stages.

`_s1_products` reads whichever SAFE or granule zip sits beside the run's outputs and returns the pair
reference-first (its own docstring: "earlier acquisition of the two") — annotation parsing needs no
unpacked raster, so this works even when the pixel data is not yet staged. Every `extra` field below
and `mission_img`/`satellite_img` were checked against a real `capture_packaging` capture of this
exact case; `sensor_img` ("C", Sentinel-1's C-band instrument) is the one constant not derivable from
either the id or the annotation.

**`acquisition_date_img1`/`img2` carry a measured, so far unexplained offset of order one second
against that same capture — `date_dt`, their difference, does not**: checked on this case, `s1_mosaic`'s
two `start`s differ from the reference's own `acquisition_date_img1`/`img2` by about +1.5 s each, but
the gap between them agrees with the reference's own `date_dt` to about 1.6 ms, well inside
`DateTime`'s millisecond rounding. So whatever `s1_mosaic`'s `start` is anchored to differs from the
reference's own epoch for the *absolute* instant by a nearly constant shift, while the physically
meaningful quantity — the interval between the two acquisitions, which is what feeds `dt_seconds`,
`vr`/`va` and the error-vector scaling — is unaffected. Worth checking against the reference's own
`IMG_INFO_DICT` population before trusting the absolute date attributes on a product this writes.
"""
function s1_img_pair_info(c::GoldenCase, s::Setup, roi_valid::Real)
    rp, sp = _s1_products(c, s.run)
    sws = burst_swaths(c)
    _, t_early = s1_mosaic(SafeSwaths(rp), sws)
    _, t_late = s1_mosaic(SafeSwaths(sp), sws)
    a1, a2 = annotation(rp, first(sws)), annotation(sp, first(sws))
    name(p) = replace(basename(p.path), r"\.(SAFE|zip)$" => "")
    f1, f2 = _s1_name_fields(name(rp)), _s1_name_fields(name(sp))
    lon, lat = radar_pair_centroid(s)
    extra = Dict{String,Any}(
        "id_img1" => name(rp), "id_img2" => name(sp),
        "sensor_img1" => "C", "sensor_img2" => "C",
        "flight_direction_img1" => lowercase(a1.pass_direction),
        "flight_direction_img2" => lowercase(a2.pass_direction),
        "absolute_orbit_number_img1" => f1.orbit, "absolute_orbit_number_img2" => f2.orbit,
        "mission_data_take_ID_img1" => f1.datatake, "mission_data_take_ID_img2" => f2.datatake,
        "product_unique_ID_img1" => f1.product_id, "product_unique_ID_img2" => f2.product_id)
    return ItsLiveOffsetProduction.ImagePairInfo(_s1_datetime(t_early), _s1_datetime(t_late), "S", "S",
                                  f1.satellite, f2.satellite, Float64(roi_valid),
                                  round(lat; digits = 2), round(lon; digits = 2), extra)
end

# A NISAR `zeroDopplerStartTime`/`zeroDopplerEndTime`-style HDF5 string
# ("2025-10-28T23:52:01.000000000", nanosecond precision) as a `DateTime`, truncated the same way
# `_s1_datetime` truncates a `UtcTime`.
function _nisar_datetime(s::AbstractString)
    whole, frac = split(s, '.'; limit = 2)
    dt = Dates.DateTime(whole, Dates.dateformat"yyyy-mm-ddTHH:MM:SS")
    return dt + Dates.Millisecond(round(Int, parse(Float64, "0." * frac) * 1000))
end

# `IMG_INFO_DICT` common to both NISAR schemas, from the two products' `Identification`.
# `mission_img`/`satellite_img` are the literal constants `"N"`/`1` the reference writes for every
# NISAR case regardless of which satellite flew it — checked against a real `capture_packaging`
# capture of both the L1 and the L2 case, not derived from `absolute_orbit` or any other field.
function _nisar_img_pair_info(early::AbstractString, late::AbstractString, id1, id2,
                              roi_valid::Real, lon::Real, lat::Real)
    extra = Dict{String,Any}(
        "id_img1" => early, "id_img2" => late,
        "sensor_img1" => "L", "sensor_img2" => "L",
        "flight_direction_img1" => lowercase(id1.pass_direction),
        "flight_direction_img2" => lowercase(id2.pass_direction),
        "absolute_orbit_number_img1" => id1.absolute_orbit,
        "absolute_orbit_number_img2" => id2.absolute_orbit)
    return ItsLiveOffsetProduction.ImagePairInfo(_nisar_datetime(id1.start_time), _nisar_datetime(id2.start_time),
                                  "N", "N", 1, 1, Float64(roi_valid),
                                  round(lat; digits = 2), round(lon; digits = 2), extra)
end

function nisar_l1_img_pair_info(c::GoldenCase, s::Setup, roi_valid::Real)
    run = rslc_run(c)
    run === nothing && error("$(basename(c.product)): no cached run holds both RSLC granules")
    early, late = acquisition_order(c)
    id1 = open_slc(joinpath(run, early * ".h5")).identification
    id2 = open_slc(joinpath(run, late * ".h5")).identification
    lon, lat = radar_pair_centroid(s)
    return _nisar_img_pair_info(early, late, id1, id2, roi_valid, lon, lat)
end

function nisar_l2_img_pair_info(c::GoldenCase, s::Setup, roi_valid::Real)
    run = gslc_run(c)
    run === nothing && error("$(basename(c.product)): no cached run holds both GSLC granules")
    early, late = acquisition_order(c)
    id1 = open_geocoded(joinpath(run, early * ".h5")).identification
    id2 = open_geocoded(joinpath(run, late * ".h5")).identification
    lon, lat = pair_centroid(s.pair.coordinate, s.epsg)
    return _nisar_img_pair_info(early, late, id1, id2, roi_valid, lon, lat)
end

"""
    radar_img_pair_info(c, s, roi_valid) -> ItsLiveOffsetProduction.ImagePairInfo

[`native_img_pair_info`](@ref) for a radar-sensor pair, dispatched to the platform's own metadata
source: `s1_img_pair_info` for Sentinel-1 (full-SLC or burst), [`nisar_l1_img_pair_info`](@ref) /
[`nisar_l2_img_pair_info`](@ref) for NISAR's RSLC and GSLC identification groups.
"""
function radar_img_pair_info(c::GoldenCase, s::Setup, roi_valid::Real)
    startswith(c.platform, "S1") && return s1_img_pair_info(c, s, roi_valid)
    c.platform == "NISAR-L2" && return nisar_l2_img_pair_info(c, s, roi_valid)
    return nisar_l1_img_pair_info(c, s, roi_valid)
end

# ---------------------------------------------------------------------------
# The whole chain
# ---------------------------------------------------------------------------

struct NativeRun
    case::String
    stages::Vector{Pair{String,Float64}}
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
    native_run(c::GoldenCase; threads, warm, trace, outname) -> (NativeRun, NamedTuple)

The chain on one case. The second return value carries the intermediate arrays, for
`--compare-intermediate`.
"""
function native_run(c::GoldenCase; warm::Bool = false, trace::Bool = true, stream::Bool = false,
                    outname::AbstractString = "julia_native.nc")
    # **`stream` is what makes a timing comparable to the container's.** `process.py` builds `/vsis3/`
    # paths and reads the granule over the network inside its own timed region, so a run that stages the
    # scenes locally first is measuring less work. Staged is the default because every *other* question
    # here is about computation; a benchmark against the reference wants this.
    if stream
        empty!(STAGED)                          # a prior call's entry would re-route to a local copy
    else
        stage(c)                                # the granules locally, before the clock starts
    end
    buf = zeros(UInt64, 64)
    stages = Pair{String,Float64}[]
    box = Ref{Any}(nothing)

    work = function ()
        local s, p, r, b1, b2, filt, declined, shifted, search, grid, pts
        local cmax, scale, factor, orig, chop, chipx, roi, path

        t = @elapsed s = setup(c)
        push!(stages, "geometry" => t)

        # The no-data mask is accumulated inside the imagery pass, where the unfiltered values live.
        t = @elapsed (b1, b2, declined, shifted, filt) = radar_sensor(c) ? radar_imagery(s) :
                                                         native_imagery(s)
        push!(stages, "imagery" => t)

        t = @elapsed begin
            cmax = driver_chip_size_max(s, declined)
            scale = AutoRIFT.chip_size_scale(s.geometry)
            # `chip_size_max` overrides what the extension derives from the geometry alone: a caller's
            # keyword wins over the inner default, which is how the driver's lower, imagery-aware
            # maximum reaches `Params`.
            p = AutoRIFT.params(s.geometry; threaded = Threads.nthreads() > 1, preprocess = :none,
                                chip_size_max = (X = cmax, Y = _even_chip(cmax * scale)))
            factor = cmax ÷ p.chip_size_min.X
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

        t = @elapsed r = AutoRIFT.autorift(b2, b1, grid, p)
        push!(stages, "correlate" => t)

        t = @elapsed begin
            # Pasted into the pre-chop grid, then dropped wherever the grown mask or a zero search
            # limit says the driver reports nothing (`testautoRIFT.py:806-818`). The second condition
            # is not implied by the first: a point can carry no search range of its own.
            dx = pad_to_orig(r.dx, orig, NaN32)
            dy = DY_SIGN .* pad_to_orig(r.dy, orig, NaN32)
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

        box[] = (; s, p, r, declined, shifted, grown, search, orig, chop, factor, cmax, scale,
                 roi, chipx, scene = size(b1), path, filt)
        return nothing
    end

    if warm
        work()
        empty!(stages)
        GC.gc(true); GC.gc(true)
    end
    floor_bytes = last(rusage!(buf))
    peak = if trace
        _, tr, _ = with_trace(; interval = 0.05) do
            work()
        end
        maximum(tr.footprint)
    else
        work()
        last(rusage!(buf))
    end

    k = box[]
    r = NativeRun(c.product, stages, Int(peak), Int(floor_bytes), k.scene, k.orig, k.chop,
                  k.factor, k.cmax, count(k.declined), count(!iszero, k.search),
                  count(!iszero, k.chipx), k.roi, k.path)
    return (r, k)
end

# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

function report(r::NativeRun)
    @printf("  %-56s\n", first(r.case, 56))
    for (name, secs) in r.stages
        @printf("    %-12s %8.1f s\n", name, secs)
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

function main(args)
    isempty(args) && error("usage: julia_e2e.jl <product-fragment> [--compare-intermediate] " *
                           "[--stream] [--no-golden] [--warm] [--no-trace]")
    c = only(cases(first(args)))
    reason = unsupported(c)
    isnothing(reason) || error("$(basename(c.product)): $reason")
    r, k = native_run(c; warm = "--warm" in args, trace = !("--no-trace" in args),
                      stream = "--stream" in args)
    report(r)
    dir = resolve_run(c)
    "--compare-intermediate" in args && compare_intermediate(k, dir)
    report_product_diff(k, c, dir; golden = !("--no-golden" in args))
    return nothing
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS)
