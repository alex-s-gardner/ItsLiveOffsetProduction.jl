# From the raw granule to `dx`/`dy`: the whole Julia chain on a golden case, timed and traced.
#
#   julia --project=tools/golden -t 12,1 tools/golden/e2e_run.jl S2B_MSIL1C_20200612
#   julia --project=tools/golden -t 12,1 tools/golden/e2e_run.jl --all --tsv e2e_julia.tsv
#
# **Every other harness here starts from a capture.** `e2e.jl` compares one stage at a time against the
# reference's own artifacts and then correlates the capture's imagery; `e2e_table.jl` joins two sweeps
# that both read captures. So the figure in `README.md` is a correlator comparison: the same bytes on
# both sides, neither arm deriving the grid. Nothing measured the chain that turns a granule into
# `dx`/`dy`.
#
# This does, on the pieces the ladder validates separately:
#
#   geometry    the scenes resolved, warped to one projection if they are not, the footprints
#               coregistered, the parameter region looked up and the geogrid solved — `setup`, which is
#               rungs 5.0 through 5.2 and 5.5 as one call
#   grid        the geogrid's own point set and parameters — `AutoRIFT.pointset`/`params` of it
#   imagery     the correlator's input built from the granule: the overlap read out of each optical
#               scene with the filter `process.py` applies to a native scene, or the reference mosaic
#               and a lazily resampled secondary for a radar burst pair
#   correlate   `autorift` at the block `block_size_for` picks
#
# Peak is the whole run's, sampled: the point of measuring a chain rather than its last stage is that
# the earlier stages are what a peak is usually made of.
#
# **The scenes are staged locally before the clock starts**, so what is timed is the computation rather
# than 708 MB of requester-pays egress. `--stream` leaves them where `scene_path` resolves them, which is
# what the reference container does — the two differ by the `imagery` stage's I/O and nothing else.
#
# **Argument order.** `arImgDisp_*` cuts its chip from its second argument and the driver calls it
# `arImgDisp(I2, I1)`, so `I1` is AutoRIFT.jl's *secondary* slot; `e2e.jl`'s rung 5.4 pins `in_I1` to the
# scene at the reference offset, so `I1` is the reference *acquisition*. The reference acquisition
# therefore goes in the secondary slot, which is what rung 5.7 does and what this repeats.

include("e2e.jl")
include("nisar.jl")      # the NISAR crop and the RSLC coregistration
include("tilecache.jl")  # a derived image kept on disk, so a later pass reads instead of re-deriving
include(joinpath(dirname(@__DIR__), "ab", "memtrace.jl"))

using Printf

argvalue(flag, default) = (i = findfirst(==(flag), ARGS);
                           isnothing(i) ? default : ARGS[i + 1])

# The flags that take a value, so the value is not mistaken for a case name.
const VALUED = ("--run", "--block", "--tsv", "--tile")

# The case fragments in an argument list: everything that is neither a flag nor a flag's value.
function positional(args)
    out = String[]
    skip = false
    for a in args
        if skip
            skip = false
        elseif startswith(a, "--")
            skip = a in VALUED
        else
            push!(out, a)
        end
    end
    return out
end

# ---------------------------------------------------------------------------
# What the chain reaches
# ---------------------------------------------------------------------------

"""
    stage(c::GoldenCase) -> Vector{String}

Put `c`'s scenes on local disk and have [`scene_path`](@ref) resolve to the copies.

Called before the clock starts, which is the point: the reference streams a Landsat band from a
requester-pays bucket and so does this chain, and 354 MB of egress inside a timed stage measures the
network. A full-SLC pair is 2.5 GB of ASF egress per acquisition and a `unzip` of it, so both trees are
expanded here too. A burst pair is already local — its SAFE trees are in the run directory — and a NISAR
pair is read from the products the driver downloaded, so neither has anything to stage.
"""
function stage(c::GoldenCase)
    early, late = acquisition_order(c)
    if c.platform == "S1-SLC"
        dirs = cached_runs(c, resolve_run(c))
        return [stage_safe(g; search = dirs) for g in (early, late)]
    end
    (startswith(c.platform, "L") || c.platform == "S2") || return String[]
    out = String[]
    for (which, name) in ((:reference, early), (:secondary, late))
        haskey(STAGED, name) && (push!(out, STAGED[name]); continue)
        STAGED[name] = stage_scene(scene_path(c, which))
        push!(out, STAGED[name])
    end
    return out
end

"""
    gslc_run(c::GoldenCase) -> Union{Nothing,String}

The cached run holding both of `c`'s GSLC granules, or `nothing`.

**Not necessarily the run the ladder compares against.** `resolve_run` prefers a run with a capture, and
the granules live wherever the driver downloaded them — 10.3 GiB each, so they are not copied per run.
"""
function gslc_run(c::GoldenCase)
    root = runs_dir(c)
    isdir(root) || return nothing
    want = [n * ".h5" for n in (first(c.reference), first(c.secondary))]
    for d in sort(readdir(root; join = true); rev = true)
        isdir(d) && all(f -> isfile(joinpath(d, f)), want) && return d
    end
    return nothing
end

"""
    s1_orbit_dir(c::GoldenCase) -> Union{Nothing,String}

The cached run holding an orbit file for each of `c`'s two acquisitions, or `nothing`.

A full-SLC pair's annotations come from the staged SAFE and its state vectors from a `POEORB` file
beside the reference's outputs, which is where the driver downloaded them.
"""
function s1_orbit_dir(c::GoldenCase)
    root = runs_dir(c)
    isdir(root) || return nothing
    for d in sort(readdir(root; join = true); rev = true)
        isdir(d) || continue
        ok = all((first(c.reference), first(c.secondary))) do g
            try
                !isempty(s1_orbit(d, g))
            catch
                false
            end
        end
        ok && return d
    end
    return nothing
end

"""
    rslc_run(c::GoldenCase) -> Union{Nothing,String}

The cached run holding both of `c`'s RSLC granules and the DEM the resample reads, or `nothing`.
"""
function rslc_run(c::GoldenCase)
    root = runs_dir(c)
    isdir(root) || return nothing
    want = [[n * ".h5" for n in (first(c.reference), first(c.secondary))]; "dem.tif"]
    for d in sort(readdir(root; join = true); rev = true)
        isdir(d) && all(f -> isfile(joinpath(d, f)), want) && return d
    end
    return nothing
end

"""
    unsupported(c::GoldenCase) -> Union{Nothing,String}

Why the chain from the granule does not exist for `c`, or `nothing` when it does.

Named per platform rather than discovered by failure, so a case this cannot measure says what is
missing instead of erroring somewhere inside a stage.
"""
function unsupported(c::GoldenCase)
    if c.platform == "NISAR-L1"
        run = rslc_run(c)
        isnothing(run) && return "no cached run holds both RSLC granules and a DEM; the pair is read \
                                  from the products themselves"
    end
    if c.platform == "NISAR-L2"
        isnothing(gslc_run(c)) &&
            return "no cached run holds both GSLC granules; a geocoded pair is read from the products \
                    themselves rather than from the driver's cropped copies"
    end
    if c.platform == "S1-SLC"
        isnothing(s1_orbit_dir(c)) &&
            return "no cached run holds the two orbit files a full-SLC pair's geometry needs"
    end
    return nothing
end

# ---------------------------------------------------------------------------
# The imagery, from the granule
# ---------------------------------------------------------------------------

"""
    native_scene(path, name, c) -> Matrix{Float32}

One optical scene, filtered as `process.py` filters it before geogrid sees it.

**The filter runs on the whole scene and the crop comes after**, which is the order `process.py` uses
and the order `Destripe` requires: its band-reject is a transform of the entire array, so filtering a
crop is a different operation rather than a cheaper one.

`:fft` is Wallis at width 5 followed by the band-reject, with the scan angles derived from the
granule's own `_ANG.txt` ephemeris — the same route gate `5.orbit` measures. `:wallis_fill` is the
gap-filling Wallis. A pair whose platforms name no native filter is read and returned.
"""
function native_scene(path::AbstractString, name::AbstractString, c::GoldenCase)
    ds = ArchGDAL.read(path)
    # ArchGDAL hands back `(x, y)`; every grid index here counts in `(row, col)`.
    img = Float32.(permutedims(ArchGDAL.read(ArchGDAL.getband(ds, 1))))
    m = native_filter(c, name)
    m === nothing && return img

    valid = img .!= 0
    if m === :wallis_fill
        f, _ = AutoRIFT.preprocess(img, valid, AutoRIFT.WallisGapfill(5, 0.25))
        return f
    end
    m === :fft || error("no route for native filter `$m` on $name")
    gt = ArchGDAL.getgeotransform(ds)
    along, cross = orbit_scan_angles(scene_ang(name, joinpath(CACHE, "angcache")), scene_epsg(ds);
                                     spacing = (gt[2], gt[6]))
    w, _ = AutoRIFT.preprocess(img, valid, AutoRIFT.Wallis(5, 0.0))
    w[.!valid] .= 0.0f0
    d, _ = AutoRIFT.preprocess(w, valid, AutoRIFT.Destripe(; along_track = along, cross_track = cross))
    d[.!valid] .= 0.0f0
    return d
end

# The same crop, kept lazy: a further index shift rather than a read.
function crop_lazy(w::LazyWindow, off::NTuple{2,Int}, want::Tuple{Int,Int})
    nr, nc = want
    ox, oy = off
    size(w, 1) >= oy + nr && size(w, 2) >= ox + nc || throw(DimensionMismatch(
        "the overlap at offset $off does not fit in a $(size(w)) grid; the crop and the geotransform " *
        "disagree"))
    return LazyWindow(w.parent, w.row0 + oy, w.col0 + ox, want)
end

# The overlap `coregister` found, out of a scene already on the correlation grid.
function crop_overlap(img::AbstractMatrix, off::NTuple{2,Int}, want::Tuple{Int,Int})
    nr, nc = want
    ox, oy = off
    size(img, 1) >= oy + nr && size(img, 2) >= ox + nc || throw(DimensionMismatch(
        "the overlap at offset $off does not fit in a $(size(img)) scene; the crop and the " *
        "geotransform disagree"))
    return img[(oy + 1):(oy + nr), (ox + 1):(ox + nc)]
end

"""
    e2e_imagery(s::Setup) -> (reference, secondary)

The correlator's two images, built from the granule, in acquisition order.

A radar burst pair returns the reference acquisition's mosaic and a [`ResampledMosaic`](@ref) — lazy, so
the secondary is resampled a block at a time by the correlator rather than materialized or written to
disk. An optical pair returns the two overlaps.
"""
function e2e_imagery(s::Setup)
    c = s.case
    if c.platform == "NISAR-L2"
        # **Geocoded, so there is nothing to coregister**: `nisar_isce3.process_gslc` crops both products
        # to their overlap, takes the amplitude and runs the projected geogrid over that. The crop is an
        # index shift into each granule's own grid rather than a copy — the overlap is billions of samples
        # — and the correlator reads its blocks straight from the HDF5 datasets.
        #
        # The *window* comes from the driver's `*_adjusted.tif` geotransform, six numbers naming where
        # its crop began; every sample comes from the granule. Deriving the window instead means
        # reproducing `crop_gslcs`' own intersection rule, which no rung checks yet.
        run = gslc_run(c)
        pair = map((:reference, :secondary)) do which
            name = which === :reference ? first(acquisition_order(c)) : last(acquisition_order(c))
            h5 = joinpath(run, name * ".h5")
            fp = gslc_footprint(h5)
            crop = which === :reference ? s.reference_path : s.secondary_path
            r0, c0 = gslc_window(fp, crop)
            ds = ArchGDAL.read(crop)
            LazyWindow(gslc_amplitude(h5), r0, c0,
                       (ArchGDAL.height(ds), ArchGDAL.width(ds)))
        end
        want = reverse(s.pair.coordinate.size)
        return (crop_lazy(pair[1], s.pair.reference_offset, want),
                crop_lazy(pair[2], s.pair.secondary_offset, want))
    end
    if c.platform == "NISAR-L1"
        # **The reference needs no resampling**: it defines the grid everything else is put on, so its
        # amplitude is the granule's own samples. The secondary is resampled onto it.
        run = rslc_run(c)
        early, late = acquisition_order(c)
        ref = SLCDatasets.amplitude(open_slc(joinpath(run, early * ".h5")))
        sec = ResampledRSLC(joinpath(run, early * ".h5"), joinpath(run, late * ".h5"),
                            dem_sampler(joinpath(run, "dem.tif")))
        size(ref) == size(sec) || error("the reference grid is $(size(ref)) and the resampled " *
                                        "secondary $(size(sec)); they are not the same grid")
        return (ref, sec)
    end
    if c.platform == "S1-SLC"
        # The same three-subswath mosaic the burst path builds, over the real granules: `radar_mosaic` of
        # the reference acquisition, and the secondary resampled onto it a block at a time. Its own dims
        # are checked against the geometry's, since a disagreement puts the point set on a grid the pixels
        # are not on.
        sws = collect(S1_SWATHS)
        dem = joinpath(s.run, "dem.tif")
        isfile(dem) || error("no dem.tif in $(s.run); the secondary's resample solves for terrain height")
        early, late = acquisition_order(c)
        # The `.SAFE` trees are wherever the driver downloaded them, which need not be the run the ladder
        # compares against — so every cached run of the case is searched before ASF is.
        safes = cached_runs(c, s.run)
        mk(g) = Sentinel1Product(stage_safe(g; search = safes); orbit = s1_orbit(s.run, g),
                                 polarization = lowercase(s1_polarization(g)), swaths = sws)
        rp, sp = mk(early), mk(late)
        ref = radar_mosaic(rp, sws)
        sec = ResampledMosaic(rp, sp, sws, dem_sampler(dem))
        size(ref) == size(sec) || error("the reference mosaic is $(size(ref)) and the resampled " *
                                        "secondary $(size(sec))")
        co = s.pair.coordinate
        # **Fail rather than correlate on a grid the point set is not on.** `s1_pair` takes the merged
        # shape from `cslc_grid` where the run kept a CSLC, and both acquisitions' CSLCs sit on the grid
        # COMPASS resampled onto — 1640 x 21458 against the annotation's 1504 x 21530 on IW1 — so the
        # reference mosaic is itself a resample there and not a merge of raw bursts. Reproducing it needs
        # the reference's own coregistration target, which no rung derives yet.
        (co.nlines, co.nsamples) == size(ref) || error(
            "the geogrid was built on a $(co.nlines) x $(co.nsamples) mosaic and the imagery is " *
            "$(size(ref)); the point set and the pixels are on different grids" *
            (isnothing(cslc_grid(s.run, sws)) ? "" :
             ". $(basename(s.run)) holds a CSLC, so the geogrid is on the resampled grid and this " *
             "mosaic is on the annotation's"))
        return (ref, sec)
    end
    if c.platform == "S1-BURST"
        rp, sp = _s1_products(c, s.run)
        sws = burst_swaths(c)
        dem = joinpath(s.run, "dem.tif")
        isfile(dem) || error("no dem.tif in $(s.run); the secondary's resample solves for terrain " *
                            "height, so the pair cannot be coregistered without one")
        return (radar_mosaic(rp, sws), ResampledMosaic(rp, sp, sws, dem_sampler(dem)))
    end
    early, late = acquisition_order(c)
    want = reverse(s.pair.coordinate.size)
    ref = crop_overlap(native_scene(s.reference_path, early, c), s.pair.reference_offset, want)
    sec = crop_overlap(native_scene(s.secondary_path, late, c), s.pair.secondary_offset, want)
    return (ref, sec)
end

# ---------------------------------------------------------------------------
# One case
# ---------------------------------------------------------------------------

struct E2EResult
    case::String
    platform::String
    stages::Vector{Pair{String,Float64}}
    cpu::Float64
    peak::Int
    floor::Int
    scene::Tuple{Int,Int}
    block::Tuple{Int,Int}
    npoints::Int
    measured::Int
    filter::String
end

total_seconds(r::E2EResult) = sum(last, r.stages)

"""
    run_case(c::GoldenCase; n = nothing, block = nothing, trace = true) -> E2EResult

The whole chain on one case, each stage timed and the run's peak sampled.

`block` overrides what [`AutoRIFT.block_size_for`](@ref) picks from the grid's own halo.
"""
function run_case(c::GoldenCase; n::Union{Integer,Nothing} = nothing,
                  block::Union{Integer,Nothing} = nothing, trace::Bool = true,
                  staged::Bool = true, warm::Bool = true, tile::Integer = 1024,
                  filter_cache::Bool = true)
    reason = unsupported(c)
    isnothing(reason) || error(reason)
    staged && stage(c)

    buf = zeros(UInt64, 64)
    hz = tick_rate()
    stages = Pair{String,Float64}[]
    caches = TileCache[]
    out = Ref{Any}(nothing)
    shape = Ref((0, 0))
    blk = Ref((0, 0))
    np = Ref(0)
    filt = Ref("")

    work = function ()
        t = @elapsed s = setup(c)
        push!(stages, "geometry" => t)

        t = @elapsed begin
            grid = AutoRIFT.pointset(s.geometry;
                                     pixel_size = ImagePairGeometry.xsize(s.pair.coordinate))
            m = correlator_filter(c)
            p = AutoRIFT.params(s.geometry; threaded = Threads.nthreads() > 1,
                                preprocess = isnothing(m) ? :none : m)
        end
        push!(stages, "grid" => t)
        np[] = length(grid.x)
        filt[] = string(isnothing(m) ? :none : m)

        # `filter_cache_tile` caches each tile's *filtered* result across the coarse and fine passes of
        # every chip-size level — see `AutoRIFT.FilterTileCache`. Only `Highpass` is covered, which is
        # what both NISAR cases use; every other platform's filter (`Wallis`, `WallisGapfill`) has no
        # such cache and still wants the `TileCache` wrap below, which caches the *raw* read instead.
        fct = (filter_cache && m isa AutoRIFT.Highpass) ? 1024 : nothing

        t = @elapsed begin
            (i1, i2) = e2e_imagery(s)
            # **Cache what is expensive to read, not what is large.** `AutoRIFT.ondisk` is the image's
            # own answer to whether a read costs I/O or a derivation; an optical overlap is a plain array
            # and is left alone, while a resampled mosaic or an HDF5 band is wrapped so the second pass
            # over it reads a tile back instead of deriving it again.
            #
            # Kept even when `fct` is set. `FilterTileCache`'s own tiles are padded by
            # `filter_reach(m)` on every side, so two adjacent filter tiles' padded reads overlap by a
            # few pixels — and each spans more than one of this cache's own tiles, since a NISAR HDF5
            # chunk is 512 px against this cache's 1024. Without `TileCache` here, that overlap and that
            # straddling both cost a fresh decode every time; measured removing it, CPU time rose 949 s
            # to 1457 s on the same 800-row crop. Both caches now memory-map rather than share one
            # locked `IOStream`, so stacking them costs no more contention than either alone would.
            (i1, i2) = map((i1, i2)) do img
                (tile > 0 && AutoRIFT.ondisk(img)) || return img
                # `tc`: assigning `c` here would rebind the case, which this closure captures.
                tc = TileCache(img; tile, dir = joinpath(CACHE, "scratch"))
                push!(caches, tc)
                return tc
            end
        end
        push!(stages, "imagery" => t)
        shape[] = size(i1)

        # The reference acquisition takes the secondary slot; see the note at the top of this file.
        #
        # NISAR RSLC bands are stored in 512x512 HDF5 chunks (`SLCDatasets/src/nisar.jl`), and
        # `block_size_for`'s `chunk` keyword exists to round the block up to a multiple of that rather
        # than pick a size that straddles chunk boundaries on every read.
        chunk = startswith(c.platform, "NISAR") ? (512, 512) : (1, 1)
        b = isnothing(block) ? AutoRIFT.block_size_for(grid, p, size(i1); chunk) :
            (; X = Int(block), Y = Int(block))
        blk[] = (b.X, b.Y)
        t = @elapsed out[] = AutoRIFT.autorift(i2, i1, grid, p, (b.X, b.Y), 0, fct)
        push!(stages, "correlate" => t)
        return nothing
    end

    # A cold process spends its first pass compiling — 9.1 s of a Sentinel-2 case's `geometry` stage and
    # 0.4 s of its `grid` — so the measured pass is the second. The peak is unaffected: both passes hold
    # the same arrays, and the collection below returns the first pass's before the second allocates.
    if warm
        work()
        empty!(stages)
        GC.gc(true); GC.gc(true)
    end
    floor_bytes = last(rusage!(buf))
    cpu0 = cpu_seconds!(buf, hz)

    peak = if trace
        _, tr, _ = with_trace(; interval = 0.05) do
            work()
        end
        maximum(tr.footprint)
    else
        work()
        last(rusage!(buf))
    end
    cpu = cpu_seconds!(buf, hz) - cpu0
    # `tc`, not `c`: the case is captured by the closure above, so a loop variable of the same name
    # writes into its box rather than shadowing it.
    for tc in caches
        @printf("    %-12s %s\n", "tile cache", cache_report(tc))
        close(tc)
    end

    return E2EResult(c.product, c.platform, stages, cpu, Int(peak), Int(floor_bytes),
                     shape[], blk[], np[], count(isfinite, out[].dx), filt[])
end

# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

function report(r::E2EResult)
    @printf("  %-56s %s\n", first(r.case, 56), r.platform)
    for (name, secs) in r.stages
        @printf("    %-12s %8.1f s\n", name, secs)
    end
    @printf("    %-12s %8.1f s   cpu %.1f\n", "total", total_seconds(r), r.cpu)
    @printf("    scene %d x %d, block %dx%d, %s, %d of %d points measured\n", r.scene...,
            r.block..., r.filter, r.measured, r.npoints)
    @printf("    peak %.2f GiB (%.2f above the floor %.2f)\n", r.peak / 2^30,
            (r.peak - r.floor) / 2^30, r.floor / 2^30)
    flush(stdout)
end

function tsv_line(r::E2EResult)
    st = Dict(r.stages)
    return join([r.case, r.platform,
                 @sprintf("%.2f", get(st, "geometry", NaN)),
                 @sprintf("%.2f", get(st, "grid", NaN)),
                 @sprintf("%.2f", get(st, "imagery", NaN)),
                 @sprintf("%.2f", get(st, "correlate", NaN)),
                 @sprintf("%.2f", total_seconds(r)), @sprintf("%.1f", r.cpu),
                 string(r.peak), string(r.floor),
                 "$(r.scene[1])x$(r.scene[2])", "$(r.block[1])x$(r.block[2])",
                 string(r.npoints), string(r.measured), r.filter], '\t')
end

const TSV_HEADER = join(["case", "platform", "geometry_s", "grid_s", "imagery_s", "correlate_s",
                         "total_s", "cpu_s", "peak_bytes", "floor_bytes", "scene", "block",
                         "npoints", "measured", "filter"], '\t')

function main(args)
    isempty(args) && error("usage: e2e_run.jl <product-name-fragment>... | --all " *
                           "[--run N] [--block N] [--no-trace] [--stream] [--cold] [--tsv FILE] " *
                           "[--no-filter-cache]")
    cs = "--all" in args ? cases() : vcat((cases(a) for a in positional(args))...)
    n = "--run" in args ? parse(Int, argvalue("--run", "0")) : nothing
    block = "--block" in args ? parse(Int, argvalue("--block", "0")) : nothing
    tsv = argvalue("--tsv", nothing)
    trace = !("--no-trace" in args)
    staged = !("--stream" in args)
    warm = !("--cold" in args)
    tile = "--no-cache" in args ? 0 : parse(Int, argvalue("--tile", "1024"))
    filter_cache = !("--no-filter-cache" in args)

    @printf("%d case(s), %d threads\n\n", length(cs), Threads.nthreads())
    rows = E2EResult[]
    for c in cs
        reason = unsupported(c)
        if !isnothing(reason)
            @printf("  %-56s %s\n    skipped: %s\n", first(c.product, 56), c.platform,
                    replace(reason, r"\s+" => " "))
            flush(stdout)
            continue
        end
        try
            r = run_case(c; n, block, trace, staged, warm, tile, filter_cache)
            push!(rows, r)
            report(r)
        catch err
            @printf("  %-56s %s\n    FAILED: %s\n", first(c.product, 56), c.platform,
                    first(sprint(showerror, err), 300))
            flush(stdout)
        end
    end

    if !isnothing(tsv) && !isempty(rows)
        open(tsv, "w") do io
            println(io, TSV_HEADER)
            foreach(r -> println(io, tsv_line(r)), rows)
        end
        @printf("\nwrote %d row(s) to %s\n", length(rows), tsv)
    end
    return nothing
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main(ARGS)
end
