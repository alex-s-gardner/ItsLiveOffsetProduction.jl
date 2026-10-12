# Run one optical case through the JIT app as far as the correlation, and write everything the
# correlation stage reads — and what it produced — to a directory that `correlate.jl` reads back.
#
#   julia --project=app -t 6,1 tools/trim/probes/dump_correlate.jl <sensor> <reference> <secondary> <outdir>
#
# Inputs written: the two scene footprints, the parameter grid and window, the twelve parameter
# windows, the DEM no-data value, the correlator filter, the two images and the no-data mask.
# Outputs written: the geogrid bands and the correlator's `dx`, `dy`, `chip_size`, `interpolated`.
# For an L4/L5 case also the reference scene, its scan angles and its `:fft` native filter result.

using AutoriftItsLiveProduction, AutoRIFT, ImagePairGeometry, ArchGDAL
const A = AutoriftItsLiveProduction
include(joinpath(@__DIR__, "rawio.jl"))

sensor, ref, sec, out = ARGS
mkpath(out)
isempty(readdir(out)) || error("$out is not empty")

c = A.Job(sensor, [ref], [sec])
opts = A.RunOptions()
mkpath(opts.scratch)
A._CACHE_DIR[] = opts.scratch

try
    s = A.setup(c)
    rpath, spath, rfp, sfp = A.aligned_scenes(c)
    b1, b2, declined, _, filt = A.native_imagery(s; cache_io = :pread)

    for (k, fp) in (("reference", rfp), ("secondary", sfp))
        writemeta(out, k * "_origin", fp.origin...)
        writemeta(out, k * "_spacing", fp.spacing...)
        writemeta(out, k * "_size", fp.size...)
        writemeta(out, k * "_epsg", A.footprint_epsg(fp))
    end
    writemeta(out, "dt", s.pair.dt)
    writemeta(out, "scene_epsg", s.epsg)
    writemeta(out, "grid_epsg", s.info.epsg)
    writemeta(out, "grid_geotransform", s.grid.geotransform...)
    writemeta(out, "grid_size", s.grid.size...)
    writemeta(out, "window", first(s.window)[1], first(s.window)[2], last(s.window)[1], last(s.window)[2])
    nd = ArchGDAL.getnodatavalue(ArchGDAL.getband(ArchGDAL.read(s.info.paths.dem), 1))
    writemeta(out, "dem_nodata", nd === nothing ? 0.0 : Float64(nd))
    if filt isa AutoRIFT.Highpass
        writemeta(out, "filter", "highpass", filt.width)
    elseif filt isa AutoRIFT.WallisGapfill
        writemeta(out, "filter", "wallis_fill", filt.width, filt.min_std)
    elseif filt === nothing
        writemeta(out, "filter", "none")
    else
        error("no dump rule for filter $(typeof(filt))")
    end
    for k in fieldnames(typeof(s.inputs))
        v = getfield(s.inputs, k)
        v === nothing || writeraw(out, "in_" * string(k), v)
    end
    writeraw(out, "image_early", Matrix{Float32}(b1[axes(b1, 1), axes(b1, 2)]))
    writeraw(out, "image_late", Matrix{Float32}(b2[axes(b2, 1), axes(b2, 2)]))
    writeraw(out, "declined", declined)

    # `run_job`'s grid and correlate stages, verbatim.
    cmax = A.driver_chip_size_max(s, declined)
    scale = AutoRIFT.chip_size_scale(s.geometry)
    kw = (; threaded = Threads.nthreads() > 1,
          preprocess = filt === nothing ? :none : filt,
          chip_size_max = (X = cmax, Y = A._even_chip(cmax * scale)))
    p = AutoRIFT.params(s.geometry; kw...)
    factor = cmax ÷ p.chip_size_min.X
    nlevels = trailing_zeros(factor) + 1
    length(p.subpixel) > nlevels &&
        (p = AutoRIFT.params(s.geometry; kw..., subpixel = p.subpixel[1:nlevels]))
    pts = AutoRIFT.pointset(s.geometry; pixel_size = ImagePairGeometry.xsize(s.pair.coordinate))
    orig = size(pts.x)
    chop = A.chopped_size(orig, factor)
    grid = A._chop_to(pts, chop[1], chop[2])
    keep = .!view(declined, 1:chop[1], 1:chop[2])
    grid = AutoRIFT.rebuild(grid; radius_x = grid.radius_x .* keep, radius_y = grid.radius_y .* keep)
    b = AutoRIFT.block_size_for(grid, p, size(b1); chunk = (512, 512))
    fct = filt isa AutoRIFT.Highpass ? opts.filter_tile : nothing
    pair = filt isa AutoRIFT.WallisGapfill ?
           AutoRIFT.ImagePair(b2, b1; reference_valid = A.NonzeroMask(b2), secondary_valid = A.NonzeroMask(b1)) :
           AutoRIFT.ImagePair(b2, b1)
    r = AutoRIFT.autorift(pair, grid, p, (b.X, b.Y), 0, fct; filter_cache_io = :pread)

    g = s.geometry
    for k in (:location_x, :location_y, :offset_x, :offset_y, :search_x, :search_y, :chip_min_x,
              :chip_min_y, :chip_max_x, :chip_max_y, :stable_surface, :off2vx_dx, :off2vx_dy,
              :off2vy_dx, :off2vy_dy, :scale_x, :scale_y)
        writeraw(out, "g_" * string(k), getfield(g, k))
    end
    writemeta(out, "cmax", cmax)
    writemeta(out, "nlevels", nlevels)
    writemeta(out, "block", b.X, b.Y)
    writeraw(out, "r_dx", r.dx)
    writeraw(out, "r_dy", r.dy)
    writeraw(out, "r_chip_size", r.chip_size)
    writeraw(out, "r_interpolated", Matrix{Bool}(r.interpolated))

    if any(n -> A.native_filter(c, n) === :fft, vcat(c.reference, c.secondary))
        early, _ = A.acquisition_order(c)
        sc = A.open_scene(s.reference_path)
        img = Matrix{Float32}(undef, size(sc))
        A.OpticalDatasets.read_window!(img, sc, axes(sc, 1), axes(sc, 2))
        along, cross = A.orbit_scan_angles(A.scene_ang(early, joinpath(A.cache_dir(), "angcache")),
                                           A.footprint_epsg(rfp); spacing = rfp.spacing)
        writeraw(out, "fft_input", img)
        writemeta(out, "fft_angles", along, cross)
        writeraw(out, "fft_output", A.fft_native_filter(img, along, cross))
    end
finally
    A.close_scratch!()
end
println("dumped to $out")
