# Timing harness for the Sentinel-1 / NISAR coregistration path — rung 5.2 of `e2e.jl` — run
# standalone against a cached golden run so each stage can be timed while trying speedups.
#
#   julia --project=tools/golden -t 8 tools/golden/bench_coregister.jl <product> [--run N] \
#       [--replay] [--repeat N]
#
# `<product>` is matched against a golden case the way `e2e.jl` matches it — enough of the name to
# be unambiguous. `--replay` uses the captured COMPASS offsets instead of solving them, which is
# what `rung_coregister` does when a run kept them; omitted, `secondary_mosaic` runs the production
# path: the offset-lattice solve (`coregistration_offset`/`geo2rdr`, interpolated) followed by the
# eight-tap sinc resample. `--repeat N` reruns the resample stage N times in this process, after one
# untimed warm-up call, and reports the minimum and median — the first call in a fresh process pays
# JIT compilation that the production run only pays once per process too, but which would otherwise
# swamp a single sub-second measurement.
#
# NISAR L1 has no mosaic to merge or resample — an RSLC is one acquisition on one radar grid — so
# only `nisar_l1_pair` is timed for it.

include("e2e.jl")

using Printf, Statistics

function _fmt(label::AbstractString, r::NamedTuple)
    gcpct = r.time > 0 ? 100 * r.gctime / r.time : 0.0
    @printf("%-42s %8.3f s   %6.2f%% gc   %7.3f GiB\n", label, r.time, gcpct, r.bytes / 2^30)
end

function bench_coregister(c::GoldenCase; n::Union{Integer,Nothing} = nothing, replay::Bool = false,
                          repeat::Integer = 1)
    run = n === nothing ? resolve_run(c) : run_dir(c, n)
    isdir(run) || error("no run at $run; run reference.jl or intermediate.jl first")

    if c.platform == "NISAR-L1"
        r = @timed nisar_l1_pair(c, run)
        _fmt("nisar_l1_pair (no mosaic, no resample)", r)
        return r
    end

    startswith(c.platform, "S1") ||
        error("$(c.platform) has no radar coregistration to time")

    r_pair = @timed c.platform == "S1-BURST" ? s1_burst_pair(c, run) : s1_pair(c, run)
    _fmt("pair geometry (s1_pair / s1_burst_pair)", r_pair)

    products = _source_products(run)
    products === nothing && error(
        "no source SAFE/zip beside $run; the mosaic and resample need pixels, which this run's \
         source product supplies")
    all(isdir, products) || error(
        "$(basename(first(products[.!isdir.(products)]))) is still zipped; unzip it beside $run \
         to read pixels a window at a time")

    rp, sp = _s1_products(c, run)
    sws = burst_swaths(c)

    refoff = replay ? _reference_offsets(run, sws) : nothing
    refgrid = replay ? _reference_grid(run, sws) : nothing
    r_ref = @timed radar_mosaic(rp, sws; offsets = refoff, grid = refgrid)
    _fmt("radar_mosaic (reference: merge, no resample)", r_ref)

    dem = joinpath(run, "dem.tif")
    isfile(dem) || error("no $dem; the coregistration solves for terrain height")
    r_dem = @timed dem_sampler(dem)
    _fmt("dem_sampler (load + parse DEM)", r_dem)
    sampler = r_dem.value

    secoff = replay ? _secondary_offsets(run, sws) : nothing
    secgrid = replay ? _secondary_grid(run, sws) : nothing
    label = replay && secoff !== nothing ? "secondary_mosaic (replaying captured offsets)" :
            "secondary_mosaic (solve + resample)"

    # One untimed warm-up, so a `--repeat` series reports steady-state cost rather than paying JIT
    # compilation on every one of its entries.
    r_sec = @timed secondary_mosaic(rp, sp, sws, sampler; offsets = secoff, grid = secgrid)
    _fmt("$label [warm-up]", r_sec)

    times = Float64[r_sec.time]
    for i in 2:repeat
        r = @timed secondary_mosaic(rp, sp, sws, sampler; offsets = secoff, grid = secgrid)
        _fmt("$label [run $i]", r)
        push!(times, r.time)
    end
    if repeat > 1
        @printf("%-42s %8.3f s (min)   %8.3f s (median)\n", label, minimum(times), median(times))
    end

    total = r_pair.time + r_ref.time + r_dem.time + minimum(times)
    @printf("%-42s %8.3f s\n", "total (pair + reference + dem + best resample)", total)
    return (; pair = r_pair.time, reference = r_ref.time, dem = r_dem.time, resample = times)
end

function main(args)
    isempty(args) && error("usage: bench_coregister.jl <product> [--run N] [--replay] [--repeat N]")
    c = only(cases(args[1]))
    n = nothing
    i = findfirst(==("--run"), args)
    i === nothing || (n = parse(Int, args[i + 1]))
    repeat = 1
    j = findfirst(==("--repeat"), args)
    j === nothing || (repeat = parse(Int, args[j + 1]))
    replay = "--replay" in args
    @info "benchmarking" product=c.product platform=c.platform threads=Threads.nthreads() replay repeat
    bench_coregister(c; n, replay, repeat)
    return nothing
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS)
