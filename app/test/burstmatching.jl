# Which secondary burst a reference burst is resampled from: the one at the same time after the ascending
# node, whatever its index. The times are those of a real S1A × S1C pair (2025-07-11 × 2025-07-17, track
# 141), whose IW1 starts a burst earlier in the S1C product.
using AutoriftItsLiveProduction: matching_burst, S1_ANX_MATCH

anx(t0, n; cycle = 2.758) = (; burst_anx_time = [t0 + (i - 1) * cycle for i in 1:n])

@testset "burst matching" begin
    @testset "a secondary that starts a burst earlier is read one burst on" begin
        ra, sa = anx(1817.31, 9), anx(1814.55, 9)
        @test [matching_burst(ra, sa, i) for i in 1:9] == [2:9; nothing]
    end

    @testset "a secondary that starts a burst later is read one burst back" begin
        ra, sa = anx(1814.55, 9), anx(1817.31, 9)
        @test [matching_burst(ra, sa, i) for i in 1:9] == [nothing; 1:8]
    end

    @testset "same framing pairs burst for burst" begin
        ra, sa = anx(1815.0, 9), anx(1815.0 + 0.04, 9)
        @test [matching_burst(ra, sa, i) for i in 1:9] == collect(1:9)
    end

    @testset "a burst more than S1_ANX_MATCH away is no match" begin
        ra, sa = anx(1815.0, 1), anx(1815.0 + 1.5 * S1_ANX_MATCH, 1)
        @test matching_burst(ra, sa, 1) === nothing
    end
end
