using Test

# `Pkg.test(; test_args = ["tidenetwork"])` runs the checks of the tide correction's published inputs,
# which read S3 and ERA5, in place of the offline suite.
if "tidenetwork" in ARGS
    include("tidenetwork.jl")
else
    @testset "AutoriftItsLiveProduction" begin
        # The parts of the chain that run without a network: the command line, how an input names its scene,
        # and geometry rules a run applies before it reads a scene. The end-to-end chain is
        # `tools/release/smoke.sh` and `tools/golden`, which read real granules.
        include("production.jl")
        include("browse.jl")
        include("catalog.jl")
        include("sceneid.jl")
        include("antimeridian.jl")
        include("tides.jl")
        include("tidecorrection.jl")
        include("tidegeometry.jl")
        include("pairquality.jl")
        include("flightdirection.jl")
    end
end
