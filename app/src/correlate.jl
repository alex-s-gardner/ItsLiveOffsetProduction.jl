# The correlator's grid and the correlation, from the geometry, the imagery and the pair's filter.

# `round(x / 2) * 2`, the reference's way of keeping a chip extent even.
even_chip(x::Real) = round(Int, x / 2) * 2

"""
    driver_chip_size_max(g, declined) -> Int

The largest x chip bound surviving the no-data mask, which is the maximum the driver hands the correlator
and so the top of its pyramid. `AutoRIFT.params(g)` reads the maximum off the geometry alone, which is
larger whenever the imagery is gappy.
"""
function driver_chip_size_max(g::PairGeometry, declined::Matrix{Bool})
    cmax = permutedims(g.chip_max_x)
    sentinel = Int32(g.nodata.output)
    size(cmax) == size(declined) || throw(DimensionMismatch("the chip bounds and the no-data mask differ in size"))
    best = 0
    for i in eachindex(cmax)
        v = cmax[i]
        (declined[i] || v == sentinel) && continue
        best = max(best, Int(v))
    end
    best > 0 || throw(ArgumentError("every chip-size maximum is declined or missing, so no pyramid level could run"))
    return best
end

# `gridsize` truncated to a multiple of `factor` on both axes (`autoRIFT.py:767`).
chopped_size(gridsize::Tuple{Int,Int}, factor::Int) =
    (fld(gridsize[1], factor) * factor, fld(gridsize[2], factor) * factor)

# `pts` truncated to `nr x nc`, which `autoRIFT.py:809-819` does to every per-point array before correlating.
# Every field, because `PointSet` requires every field to share `x`'s axes.
function chop_to(pts, nr::Int, nc::Int)
    c(A) = A[1:nr, 1:nc]
    return AutoRIFT.rebuild(pts; x = c(pts.x), y = c(pts.y),
                            radius_x = c(pts.radius_x), radius_y = c(pts.radius_y),
                            dx_prior = c(pts.dx_prior), dy_prior = c(pts.dy_prior),
                            chip_size_x = c(pts.chip_size_x), chip_size_y = c(pts.chip_size_y),
                            chip_size_min_x = c(pts.chip_size_min_x),
                            chip_size_max_x = c(pts.chip_size_max_x),
                            positioned = c(pts.positioned))
end

# The search limit the driver gives the correlator, over the whole window: the geogrid band with the
# sentinel and the declined points zeroed, and the chopped-away strip zeroed too. This is
# `roi_valid_percentage`'s denominator.
function driver_search_limit(g::PairGeometry, declined::Matrix{Bool}, chop::Tuple{Int,Int})
    sx = permutedims(g.search_x)
    sentinel = Int32(g.nodata.output)
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

"""
    Correlation

What the correlation stage hands packaging: the correlator's result pasted into the pre-chop grid, the
masks that gate it, and the chip-size ratio.
"""
struct Correlation
    dx::Matrix{Float32}
    dy::Matrix{Float32}
    chip_size::Matrix{UInt16}
    interpolated::Matrix{Bool}
    search::Matrix{Int32}
    scale::Float64
    # The correlator's grid spacing in x, which the Sentinel-1 subswath bias correction smooths over.
    grid_spacing_x::Float64
end

# `A` pasted into the top-left corner of an `origsize` array of `filler` (`testautoRIFT.py:792-805`): the
# reference allocates every output at the pre-chop size, so the strip the chop dropped reaches the product
# as no-data.
function pad_to_orig(A::Matrix{T}, origsize::Tuple{Int,Int}, filler::T) where {T}
    out = fill(filler, origsize)
    out[axes(A, 1), axes(A, 2)] .= A
    return out
end

# `arImgDisp_*` negates `Dy0` as its first act and AutoRIFT.jl has no such internal step, so
# `MultichipResult.dy` runs opposite to the reference's `Dy`, which `testautoRIFT.py:800` passes straight
# into `netCDF_packaging` for an optical pair.
const DY_SIGN = -1.0f0

# The sign the correlator's `dy` is written with. A radar pair's `Dy` is negated once more before packaging
# (`testautoRIFT.py:785-786`, `if optical_flag == 0: Dy = -Dy`); a NISAR L2 pair runs with
# `optical_flag = 1`, as an optical one.
dy_sign(::Setup) = DY_SIGN

"""
    correlate(s, img, pre; filter_tile, cache_io, masked, scratch) -> Correlation

`run_job`'s grid and correlate stages for one correlator filter `pre`: the driver's chip-size maximum and
the parameters the geometry implies, the grid chopped to a multiple of the coarsest level and declined where
the no-data mask is set, and a blocked `autorift` over 512² read chunks. `filter_tile` is the side of
`AutoRIFT.FilterTileCache`'s tiles, which keeps each filtered tile across the passes in a file under
`scratch`; that cache holds only a `Highpass` result. `masked` hands `WallisGapfill` its gaps as the pixels
reading zero (`np.isclose(image, 0)`), which AutoRIFT reads from the validity masks rather than from the
values.

The outputs are pasted into the pre-chop grid and dropped wherever the grown mask or a zero search limit
says the driver reports nothing (`testautoRIFT.py:806-818`).
"""
function correlate(s, img::Imagery, pre::AutoRIFT.PreprocessMethod; filter_tile::Union{Nothing,Int},
                   cache_io::Symbol, masked::Bool, scratch::String)
    g = s.geometry
    declined = img.declined
    cmax = driver_chip_size_max(g, declined)
    scale = AutoRIFT.chip_size_scale(g)
    p = AutoRIFT.params(g; threaded = AutoRIFT.True(), preprocess = pre,
                        chip_size_max = (X = cmax, Y = even_chip(cmax * scale)))
    factor = cmax ÷ p.chip_size_min.X
    pts = AutoRIFT.pointset(g; pixel_size = Float64(ImagePairGeometry.xsize(s.pair.coordinate)))
    orig = size(pts.x)
    chop = chopped_size(orig, factor)
    # The grown mask gates the outputs; the narrow one gates the search. Both are the driver's, and
    # conflating them moves the product's coverage.
    grown = _dilate(declined, 1)
    search = driver_search_limit(g, grown, chop)
    grid = chop_to(pts, chop[1], chop[2])
    keep = .!view(declined, 1:chop[1], 1:chop[2])
    grid = AutoRIFT.rebuild(grid; radius_x = grid.radius_x .* keep, radius_y = grid.radius_y .* keep)
    b = AutoRIFT.block_size_for(grid, p, size(img.early); chunk = (512, 512))
    r = if masked
        pair = AutoRIFT.ImagePair(img.late, img.early; reference_valid = NonzeroMask(img.late),
                                  secondary_valid = NonzeroMask(img.early))
        AutoRIFT.autorift(pair, grid, p, (b.X, b.Y), 0, filter_tile; filter_cache_io = cache_io,
                          filter_cache_dir = scratch)
    else
        AutoRIFT.autorift(AutoRIFT.ImagePair(img.late, img.early), grid, p, (b.X, b.Y), 0, filter_tile;
                          filter_cache_io = cache_io, filter_cache_dir = scratch)
    end
    dx = pad_to_orig(Matrix{Float32}(r.dx), orig, NaN32)
    dy = dy_sign(s) .* pad_to_orig(Matrix{Float32}(r.dy), orig, NaN32)
    chipx = pad_to_orig(Matrix{UInt16}(r.chip_size), orig, UInt16(0))
    interp = pad_to_orig(Matrix{Bool}(r.interpolated), orig, false)
    for i in eachindex(dx)
        if grown[i] || search[i] == 0
            dx[i] = NaN32
            dy[i] = NaN32
            chipx[i] = 0
            interp[i] = false
        end
    end
    return Correlation(dx, dy, chipx, interp, search, scale, Float64(p.grid_spacing.X))
end
