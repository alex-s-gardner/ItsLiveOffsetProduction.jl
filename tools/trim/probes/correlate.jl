# The chain's geometry and correlation stages on inputs `dump_correlate.jl` wrote: coregistration, the
# geogrid over typed FastGeoProjections transforms, the correlator parameters the geometry implies, and
# `autorift` with the pair's filter. Writes the geogrid bands and the correlator's result to <outdir> in
# the dump's format, so that the JIT run, the trimmed run and the app's own result can be compared byte
# for byte.
#
#   correlate <dumpdir> <outdir>

using AutoRIFT, ImagePairGeometry
import FastGeoProjections as FGP

include(joinpath(@__DIR__, "rawio.jl"))

say(s::String) = (write(Core.stderr, s); write(Core.stderr, "\n"); nothing)

function footprint(d::String, k::String)
    o = meta_floats(d, k * "_origin")
    sp = meta_floats(d, k * "_spacing")
    sz = meta_ints(d, k * "_size")
    return ImageFootprint((o[1], o[2]), (sp[1], sp[2]), (sz[1], sz[2]), meta_int(d, k * "_epsg"))
end

# The twelve parameter windows, in `GeometryInputs`' own field order.
function geometry_inputs(d::String)
    r(k) = readraw(d, "in_" * k, Float64)
    return GeometryInputs(r("dem"), r("dhdx"), r("dhdy"), r("vx"), r("vy"), r("srx"), r("sry"),
                          r("csminx"), r("csminy"), r("csmaxx"), r("csmaxy"), r("ssm"))
end

# The geogrid over the window the transform's footprint bounds give. The only place the transform's type
# is known; the `PairGeometry` returned has a type independent of it.
function geogrid(grid::MapGrid, pair::CoregisteredPair, inputs::GeometryInputs, nd::Float64, tf::TransformPair)
    window = grid_window(grid, footprint_bounds(tf, pair.coordinate))
    return pairgeometry_blocked(grid, pair, InMemoryInputs(inputs, window); transform = () -> tf, window,
                                nodata = nodata_from(nd))
end

function geogrid(grid::MapGrid, pair::CoregisteredPair, inputs::GeometryInputs, nd::Float64,
                 g::FGP.NativeCRS, s::FGP.NativeCRS)
    f = FGP.Transformation(g, s)
    return geogrid(grid, pair, inputs, nd, TransformPair(f, inv(f)))
end

# Grid CRS to scene CRS, one concrete transform type per branch.
function geogrid(grid::MapGrid, pair::CoregisteredPair, inputs::GeometryInputs, nd::Float64,
                 gcode::Int, scode::Int)
    gcode == scode && return geogrid(grid, pair, inputs, nd, transform_pair(IdentityTransform()))
    g, s = FGP.native_crs(gcode), FGP.native_crs(scode)
    if g isa FGP.UTMCRS
        s isa FGP.UTMCRS && return geogrid(grid, pair, inputs, nd, g, s)
        s isa FGP.PolarStereographicCRS && return geogrid(grid, pair, inputs, nd, g, s)
    elseif g isa FGP.PolarStereographicCRS
        s isa FGP.UTMCRS && return geogrid(grid, pair, inputs, nd, g, s)
        s isa FGP.PolarStereographicCRS && return geogrid(grid, pair, inputs, nd, g, s)
    end
    throw(ArgumentError("no optical transform from EPSG:" * string(gcode) * " to EPSG:" * string(scode)))
end

struct NonzeroMask{A<:AbstractMatrix} <: AbstractMatrix{Bool}
    parent::A
end
Base.size(m::NonzeroMask) = size(m.parent)
Base.getindex(m::NonzeroMask, rows::AbstractUnitRange, cols::AbstractUnitRange) = m.parent[rows, cols] .!= 0
Base.getindex(m::NonzeroMask, i::Integer, j::Integer) = m.parent[i, j] != 0
AutoRIFT.ondisk(m::NonzeroMask) = AutoRIFT.ondisk(m.parent)

function max_outside(cmax::AbstractMatrix{<:Integer}, declined::AbstractMatrix{Bool}, sentinel::Int32)
    best = 0
    for i in eachindex(cmax, declined)
        v = cmax[i]
        (declined[i] || v == sentinel) && continue
        best = max(best, Int(v))
    end
    return best
end

even_chip(x::Real) = round(Int, x / 2) * 2

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

# `run_job`'s grid and correlate stages for one filter, written to `out`.
function correlate(out::String, g::PairGeometry, pixel_size::Float64, pre::AutoRIFT.PreprocessMethod,
                   early::Matrix{Float32}, late::Matrix{Float32}, declined::Matrix{Bool},
                   filter_tile::Union{Nothing,Int}, masked::Bool)
    cmax = max_outside(permutedims(g.chip_max_x), declined, Int32(g.nodata.output))
    cmax > 0 || throw(ArgumentError("every chip-size maximum is declined or missing"))
    scale = AutoRIFT.chip_size_scale(g)
    p = AutoRIFT.params(g; threaded = AutoRIFT.True(), preprocess = pre,
                        chip_size_max = (X = cmax, Y = even_chip(cmax * scale)))
    factor = cmax ÷ p.chip_size_min.X
    pts = AutoRIFT.pointset(g; pixel_size)
    orig = size(pts.x)
    chop = (fld(orig[1], factor) * factor, fld(orig[2], factor) * factor)
    grid = chop_to(pts, chop[1], chop[2])
    keep = .!view(declined, 1:chop[1], 1:chop[2])
    grid = AutoRIFT.rebuild(grid; radius_x = grid.radius_x .* keep, radius_y = grid.radius_y .* keep)
    b = AutoRIFT.block_size_for(grid, p, size(early); chunk = (512, 512))
    r = if masked
        pair = AutoRIFT.ImagePair(late, early; reference_valid = NonzeroMask(late),
                                  secondary_valid = NonzeroMask(early))
        AutoRIFT.autorift(pair, grid, p, (b.X, b.Y), 0, filter_tile; filter_cache_io = :pread)
    else
        AutoRIFT.autorift(AutoRIFT.ImagePair(late, early), grid, p, (b.X, b.Y), 0, filter_tile;
                          filter_cache_io = :pread)
    end
    writeraw(out, "r_dx", r.dx)
    writeraw(out, "r_dy", r.dy)
    writeraw(out, "r_chip_size", r.chip_size)
    writeraw(out, "r_interpolated", Matrix{Bool}(r.interpolated))
    say("correlated: cmax " * string(cmax) * ", levels " * string(trailing_zeros(factor) + 1) *
        ", block " * string(b.X) * "x" * string(b.Y))
    return nothing
end

function run(d::String, out::String)
    rfp, sfp = footprint(d, "reference"), footprint(d, "secondary")
    pair = coregister(rfp, sfp; dt = meta_float(d, "dt"))
    gt = meta_floats(d, "grid_geotransform")
    sz = meta_ints(d, "grid_size")
    grid = MapGrid((gt[1], gt[2], gt[3], gt[4], gt[5], gt[6]), (sz[1], sz[2]), meta_int(d, "grid_epsg"))
    inputs = geometry_inputs(d)
    g = geogrid(grid, pair, inputs, meta_float(d, "dem_nodata"), meta_int(d, "grid_epsg"),
                meta_int(d, "scene_epsg"))
    w = meta_ints(d, "window")
    (first(g.window)[1], first(g.window)[2], last(g.window)[1], last(g.window)[2]) == (w[1], w[2], w[3], w[4]) ||
        throw(ArgumentError("the geogrid window differs from the dumped one"))
    for (k, band) in (("location_x", g.location_x), ("location_y", g.location_y), ("offset_x", g.offset_x),
                      ("offset_y", g.offset_y), ("search_x", g.search_x), ("search_y", g.search_y),
                      ("chip_min_x", g.chip_min_x), ("chip_min_y", g.chip_min_y), ("chip_max_x", g.chip_max_x),
                      ("chip_max_y", g.chip_max_y), ("stable_surface", g.stable_surface))
        writeraw(out, "g_" * k, band)
    end
    for (k, band) in (("off2vx_dx", g.off2vx_dx), ("off2vx_dy", g.off2vx_dy), ("off2vy_dx", g.off2vy_dx),
                      ("off2vy_dy", g.off2vy_dy), ("scale_x", g.scale_x), ("scale_y", g.scale_y))
        writeraw(out, "g_" * k, band)
    end
    say("geogrid: " * string(nvalid(g)) * " of " * string(npoints(g)) * " points in the image")

    early = readraw(d, "image_early", Float32)
    late = readraw(d, "image_late", Float32)
    declined = readraw(d, "declined", Bool)
    px = Float64(ImagePairGeometry.xsize(pair.coordinate))
    f = meta_strings(d, "filter")
    if f[1] == "highpass"
        correlate(out, g, px, AutoRIFT.Highpass(parse(Int, f[2])), early, late, declined, 512, false)
    elseif f[1] == "wallis_fill"
        correlate(out, g, px, AutoRIFT.WallisGapfill(parse(Int, f[2]), parse(Float64, f[3])), early, late,
                  declined, nothing, true)
    elseif f[1] == "none"
        correlate(out, g, px, AutoRIFT.NoPreprocess(), early, late, declined, nothing, false)
    else
        throw(ArgumentError("unknown filter " * f[1]))
    end
    return nothing
end

function (@main)(args::Vector{String})::Cint
    length(args) == 2 || (say("usage: correlate <dumpdir> <outdir>"); return Cint(2))
    d, out = args[1], args[2]
    mkpath(out)
    isempty(readdir(out)) || (say(out * " is not empty"); return Cint(2))
    try
        run(d, out)
    catch err
        m = err isa ArgumentError ? err.msg : err isa ErrorException ? err.msg : nothing
        say("error: " * (m isa String ? m : string(nameof(typeof(err)))))
        return Cint(1)
    end
    return Cint(0)
end
