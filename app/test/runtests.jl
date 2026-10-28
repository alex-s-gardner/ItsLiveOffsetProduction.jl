using Test

@testset "AutoriftItsLiveProduction" begin
    # The parts of the chain that run without a network: the command line, how an input names its scene,
    # and geometry rules a run applies before it reads a scene. The end-to-end chain is
    # `tools/release/smoke.sh` and `tools/golden`, which read real granules.
    include("production.jl")
    include("sceneid.jl")
    include("antimeridian.jl")
end
