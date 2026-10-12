# How closely two products agree where they are not expected to be identical: per variable, each one's
# valid-pixel count, the count valid in both, and over those the median and 95th percentile of |a - b| and
# the share that are exactly equal.
#
#   julia --project=tools/golden agree.jl <a.nc> <b.nc> [variables...]
using NCDatasets, Statistics
a, b = ARGS[1], ARGS[2]
vars = length(ARGS) > 2 ? ARGS[3:end] : ["vx", "vy", "v", "chip_size_width"]
NCDataset(a) do da
    NCDataset(b) do db
        size(da["vx"]) == size(db["vx"]) || error("the products are on different grids: ", size(da["vx"]),
                                                  " vs ", size(db["vx"]))
        for k in vars
            x = Float64.(coalesce.(Array(da[k]), NaN))
            y = Float64.(coalesce.(Array(db[k]), NaN))
            va, vb = isfinite.(x), isfinite.(y)
            both = va .& vb
            d = abs.(x[both] .- y[both])
            q(p) = isempty(d) ? NaN : round(quantile(d, p); sigdigits = 3)
            println(rpad(k, 16), "valid a ", count(va), "  b ", count(vb), "  both ", count(both),
                    "  |a-b| median ", q(0.5), " p95 ", q(0.95), "  equal ",
                    isempty(d) ? "-" : string(round(100 * count(iszero, d) / length(d); digits = 2)) * "%")
        end
    end
end
