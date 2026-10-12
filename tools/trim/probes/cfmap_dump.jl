# `cf_grid_mapping` for a fixed set of EPSG codes, one `<code>.tsv` per code (`<key>\t<type>\t<value>`),
# so two versions of the function can be compared.
#
#   julia --project=<env> cfmap_dump.jl <outdir>
using ItsLiveOffsetProduction
out = ARGS[1]
mkpath(out)
for code in (3413, 3031, 32607, 32622, 32645, 32719, 32601, 32760)
    d = ItsLiveOffsetProduction.cf_grid_mapping(code)
    open(joinpath(out, string(code) * ".tsv"), "w") do io
        for k in sort!(collect(keys(d)))
            println(io, k, '\t', typeof(d[k]), '\t', d[k])
        end
    end
end
println("wrote ", out)
