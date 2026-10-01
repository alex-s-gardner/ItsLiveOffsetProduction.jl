# The golden comparison, one case or one phase at a time.
#
#   julia --project=tools/golden tools/golden/run.jl --selftest        # the harness's own gate
#   julia --project=tools/golden tools/golden/run.jl --reproducibility S2B_MSIL1C_20200612
#   julia --project=tools/golden tools/golden/run.jl --status
#
# `--reproducibility` is the measurement every tolerance rests on: run the reference container twice
# on one granule and diff the two products against each other and against golden. Whatever differs
# between two runs of unchanged code is irreducible, and no tolerance below that floor is meaningful
# for any implementation.
#
# Results are written to `results/<product>.<kind>.json` with the machine and versions that produced
# them, so a regression is a diff rather than a memory.

include("manifest.jl")
include("product.jl")
include("compare.jl")
include("reference.jl")

using JSON3, Dates, Printf

results_dir() = joinpath(@__DIR__, "results")

"""
    provenance() -> Dict

What produced a result: machine, OS, Julia, and the AutoRIFT.jl commit. Recorded beside every
comparison, since a number without the code that made it cannot be checked later.
"""
function provenance()
    commit = try
        strip(read(`git -C $(@__DIR__) rev-parse --short HEAD`, String))
    catch
        "unknown"
    end
    dirty = try
        !isempty(strip(read(`git -C $(@__DIR__) status --porcelain`, String)))
    catch
        false
    end
    return Dict("machine" => Sys.MACHINE, "cpu" => Sys.CPU_NAME, "ncores" => Sys.CPU_THREADS,
                "julia" => string(VERSION), "autorift_commit" => commit * (dirty ? "-dirty" : ""),
                "image" => IMAGE, "recorded" => string(now()))
end

"""
    record(name, kind, payload)

Write one comparison to `results/`. `kind` names what was compared — `"reproducibility"` for two
container runs against each other, and later `"julia"` for AutoRIFT.jl against golden.
"""
function record(name::AbstractString, kind::AbstractString, payload::Dict)
    mkpath(results_dir())
    path = joinpath(results_dir(), "$name.$kind.json")
    open(path, "w") do io
        JSON3.pretty(io, Dict("provenance" => provenance(), "comparison" => payload))
    end
    return path
end

# A `ProductDiff` as plain data, so a result file is readable without this code.
function as_dict(d::ProductDiff)
    return Dict(
        "a" => d.a, "b" => d.b,
        "variables" => [Dict("name" => v.name, "n_both" => v.n_both,
                             "only_a" => v.only_a, "only_b" => v.only_b,
                             "exact_fraction" => exact_fraction(v),
                             "max_abs" => v.max_abs, "p99" => v.p99, "bias" => v.bias)
                        for v in d.vars],
        "x_max" => d.x_max, "y_max" => d.y_max,
        "time_delta_s" => d.time_delta,
        "attrib_diffs" => Dict(k => [string(v[1]), string(v[2])] for (k, v) in d.attrib_diffs),
        "missing_vars" => d.missing_vars,
        "identical" => identical(d),
        # The one to read: `identical` also requires the time coordinate to match, which no two runs
        # of anything can manage.
        "agrees_on_data" => agrees_on_data(d),
    )
end

"""
    reproducibility(c::GoldenCase; threads = 8) -> ProductDiff

Run the reference twice on `c` and compare the two products.

What differs here is the floor: it is produced by unchanged code on identical inputs, so nothing
below it can be attributed to an implementation. Also compares run 1 against golden, which says
whether this machine reproduces ASF's output at all — an arm64 container against whatever built the
golden data.
"""
function reproducibility(c::GoldenCase; threads::Integer = 8)
    d1 = run_reference(c; n = 1, threads)
    d2 = run_reference(c; n = 2, threads)

    p1 = read_product(run_product(d1))
    p2 = read_product(run_product(d2))
    floor_diff = compare_products(p1, p2)

    println("\n=== two container runs against each other (the reproducibility floor) ===")
    println(floor_diff)
    payload = Dict("floor" => as_dict(floor_diff))

    if have_golden(c)
        g = read_product(c)
        vs_golden = compare_products(p1, g)
        println("\n=== container run 1 against golden ===")
        println(vs_golden)
        payload["vs_golden"] = as_dict(vs_golden)
    else
        @warn "golden product not cached; run fetch.jl" product=c.product
    end

    path = record(c.product, "reproducibility", payload)
    println("\nrecorded: ", path)
    return floor_diff
end

"""
    status()

What is cached and what is reachable, per case. The first thing to run in a fresh checkout.
"""
function status()
    @printf("%-42s %-9s %-6s %-8s %-8s\n", "product", "platform", "phase", "golden", "runs")
    for c in cases()
        nruns = isdir(runs_dir(c)) ? length(readdir(runs_dir(c))) : 0
        @printf("%-42s %-9s %-6d %-8s %-8d\n",
                first(c.product, 41), c.platform, c.phase,
                have_golden(c) ? "yes" : "no", nruns)
    end
    println("\ncache: ", CACHE)
    println("image: ", IMAGE, image_present() ? " (pulled)" : " (NOT pulled)")
    return nothing
end

function main(args)
    if "--selftest" in args
        include(joinpath(@__DIR__, "selftest.jl"))
        return nothing
    end
    if "--status" in args || isempty(args)
        return status()
    end
    i = findfirst(==("--reproducibility"), args)
    if i !== nothing
        t = 8
        j = findfirst(==("--threads"), args); j === nothing || (t = parse(Int, args[j + 1]))
        reproducibility(only(cases(args[i + 1])); threads = t)
        return nothing
    end
    error("usage: run.jl [--status | --selftest | --reproducibility <name> [--threads N]]")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS)
