using Test

@testset "ItsLiveOffsetProduction" begin
    # The netCDF schemas, the metadata arithmetic, and the radar corrections — everything that runs
    # without a network. The granule-to-product pipeline is `tools/golden/julia_e2e.jl`, which needs
    # AWS credentials and is run by hand.
    include("itslive.jl")
end
