# Compare every array a probe wrote with the array of the same name in another directory, byte for
# byte; report each that differs with its count of differing elements.
#
#   julia cmpraw.jl <dir_a> <dir_b>
a, b = ARGS
ok = true
for line in eachline(joinpath(a, "index.txt"))
    name = split(line)[1]
    fb = joinpath(b, name * ".bin")
    if !isfile(fb)
        println(rpad(name, 20), " missing in ", b); global ok = false; continue
    end
    x, y = read(joinpath(a, name * ".bin")), read(fb)
    if x == y
        println(rpad(name, 20), " identical")
    else
        global ok = false
        T = Dict("Float32" => Float32, "Float64" => Float64, "Int32" => Int32, "UInt16" => UInt16,
                 "Bool" => Bool)[split(line)[2]]
        u, v = reinterpret(T, x), reinterpret(T, y)
        n = length(u) == length(v) ? count(i -> !isequal(u[i], v[i]), eachindex(u)) : -1
        println(rpad(name, 20), " DIFFERS: ", n, " of ", length(u), " elements")
    end
end
println(ok ? "ALL IDENTICAL" : "DIFFERENCES")
