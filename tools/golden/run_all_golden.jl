# Sweep every golden case: capture the reference's real `netCDF_packaging` call, write the Julia
# equivalent, and diff against the reference's own freshly-produced product.
#
#   julia --project=tools/golden tools/golden/run_all_golden.jl
#   julia --project=tools/golden tools/golden/run_all_golden.jl --skip S2B_MSIL1C_20200612,S1A_IW_SLC__1SSH_20150828
#
# One case's failure does not stop the sweep; failures are recorded and reported at the end.

include(joinpath(@__DIR__, "validate_itslive_write.jl"))

function all_case_fragments()
    return [c.product for c in cases()]
end

function run_all(; skip = String[])
    fragments = all_case_fragments()
    results = NamedTuple[]
    for frag in fragments
        any(s -> occursin(s, frag), skip) && continue
        println("\n" * "="^100)
        println("case: ", frag)
        println("="^100)
        t0 = time()
        try
            r = run_case(frag)
            elapsed = time() - t0
            push!(results, (; frag, ok = true, r, elapsed, err = nothing))
            println("physics_ok=", r.physics_ok, "  agrees_on_data=", r.agrees,
                    "  schema=", r.schema_mine, "/", r.schema_theirs,
                    "  cropped=", r.cropped_mine, "/", r.cropped_theirs,
                    "  (", round(elapsed / 60; digits = 1), " min)")
            if !isempty(r.diff.attrib_diffs) || !isempty(r.diff.missing_vars)
                show(stdout, r.diff)
                println()
            end
        catch e
            elapsed = time() - t0
            push!(results, (; frag, ok = false, r = nothing, elapsed, err = e))
            println("FAILED after ", round(elapsed / 60; digits = 1), " min: ", sprint(showerror, e))
        end
    end

    println("\n" * "="^100)
    println("SUMMARY")
    println("="^100)
    for res in results
        if !res.ok
            println(rpad(res.frag[1:min(60, end)], 62), "  ERROR: ", sprint(showerror, res.err)[1:min(80, end)])
            continue
        end
        r = res.r
        bad_vars = filter(v -> !agrees(v), r.diff.vars)
        status = r.physics_ok && isempty(bad_vars) ? "OK" : "CHECK"
        println(rpad(res.frag[1:min(60, end)], 62), "  ", status,
                "  physics=", r.physics_ok, "  bad_vars=", [v.name for v in bad_vars],
                "  missing=", r.diff.missing_vars)
    end
    return results
end

function main(args)
    skip = String[]
    i = findfirst(==("--skip"), args)
    i !== nothing && (skip = split(args[i + 1], ","))
    run_all(; skip)
    return nothing
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS)
