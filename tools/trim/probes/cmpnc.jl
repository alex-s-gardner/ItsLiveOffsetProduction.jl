# Compare two netCDF products: the header as `ncdump -h -s` prints it, line by line in any order (so
# attribute order is free but every type and value must match), and every variable's data exactly.
# `date_created` and `_NCProperties` are excluded from the header, and `v_error` is compared only where
# `v` is nonzero: elsewhere it is a Monte Carlo draw.
#
#   julia --project=tools/golden cmpnc.jl <a.nc> <b.nc>
using NCDatasets
a, b = ARGS
header(f) = sort!(filter(l -> !occursin(r":date_created =|:_NCProperties =|^netcdf ", l),
                         readlines(`ncdump -h -s $f`)))
ha, hb = header(a), header(b)
ok = true
if ha != hb
    global ok = false
    println("HEADER DIFFERS")
    for l in setdiff(ha, hb); println("  only in a: ", l); end
    for l in setdiff(hb, ha); println("  only in b: ", l); end
else
    println("header identical (", length(ha), " lines)")
end
NCDataset(a) do da
    NCDataset(b) do db
        for k in keys(da)
            x, y = Array(da[k].var[:]), Array(db[k].var[:])
            if k == "v_error"
                v = Array(da["v"].var[:])
                m = v .!= 0
                same = x[m] == y[m]
                println(rpad(k, 18), same ? "identical where v != 0" : "DIFFERS where v != 0",
                        " (", count(x .!= y), " of ", length(x), " differ overall)")
                same || (global ok = false)
            else
                same = isequal(x, y)
                println(rpad(k, 18), same ? "identical" : "DIFFERS: $(count(.!isequal.(x, y))) of $(length(x))")
                same || (global ok = false)
            end
        end
    end
end
println(ok ? "PRODUCTS MATCH" : "PRODUCTS DIFFER")
