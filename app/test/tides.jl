# Tide heights against pyTMD's TMD3 prediction (`tools/tidal/pytmd_reference.py synthetic`): random
# complex constants for the ten CATS2008 constituents at six points and thirteen times from 2014 to 2026.

using Dates: DateTime
using AutoriftItsLiveProduction: tide_days, tide_height, TidePrediction

@testset "tide heights match pyTMD" begin
    lines = readlines(joinpath(@__DIR__, "data", "tide_pytmd_synthetic.csv"))
    @test split(lines[1], ',')[5:6] == ["m2_re", "m2_im"]
    for l in lines[2:end]
        f = split(l, ',')
        t = DateTime(f[2][1:23])
        @test tide_days(t) ≈ parse(Float64, f[3]) atol = 1.0e-9
        z = ntuple(i -> complex(parse(Float64, f[3 + 2i]), parse(Float64, f[4 + 2i])), 10)
        @test tide_height(z, t) ≈ parse(Float64, f[4]) atol = 1.0e-8
    end
    # One `TidePrediction` serves every point at its instant.
    z = ntuple(i -> complex(0.1i, -0.05i), 10)
    t = DateTime(2020, 1, 13, 3, 12)
    @test tide_height(z, TidePrediction(t)) == tide_height(z, t)
end
