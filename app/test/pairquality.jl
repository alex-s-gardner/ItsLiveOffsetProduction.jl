# The pair-quality diagnostics and the no-result error: an image without data or without structure is
# warned about, a brightness change alone is not, and a pair that resolves nothing ends the run with its
# diagnostics in the message.

using Random: Xoshiro
using AutoriftItsLiveProduction: Imagery, TileCache, fill_tiles!, pair_quality, quality_warnings,
                                 quality_report, check_resolved, Correlation, PairQuality, ImageSample

# A 2048² scene with structure at every scale above the 32-px block, and speckle on top.
function scene(rng; n = 2048)
    base = zeros(Float32, n, n)
    for k in 1:40
        cx, cy, r = rand(rng, 1:n), rand(rng, 1:n), rand(rng, 40:300)
        a = rand(rng, Float32) * 2
        for j in max(1, cy - r):min(n, cy + r), i in max(1, cx - r):min(n, cx + r)
            base[i, j] += a
        end
    end
    return base .+ 1
end
speckle(rng, a) = a .* Float32.(-log.(rand(rng, Float32, size(a)...) .+ eps(Float32)))
tiles(a, dir) = fill_tiles!(TileCache(a; tile = 1024, dir, io = :pread))

@testset "pair quality: same scene, dark scene, empty and flat images" begin
    rng = Xoshiro(38)
    dir = mktempdir()
    s = scene(rng)
    e = tiles(speckle(rng, s), dir)
    good = Imagery(e, tiles(speckle(rng, s), dir), falses(1, 1))
    q = pair_quality(good)
    @test q.early.valid > 0.99 && q.late.valid > 0.99
    @test abs(q.db) < 0.5
    @test isempty(quality_warnings(q))

    # The same scene 34 dB darker still correlates: brightness is reported, not warned about.
    dark = Imagery(e, tiles(speckle(rng, s) .* Float32(10^(-34 / 20)), dir), falses(1, 1))
    q = pair_quality(dark)
    @test q.db ≈ -34 atol = 1
    @test isempty(quality_warnings(q))

    # An image of no data, and one of a single value, leave a chip nothing to match.
    empty = Imagery(e, tiles(zeros(Float32, size(s)), dir), falses(1, 1))
    q = pair_quality(empty)
    @test q.late.valid == 0 && isnan(q.db)
    @test quality_warnings(q) == ["the late image has no data in any sample window"]
    flat = Imagery(tiles(fill(3.0f0, size(s)), dir), e, falses(1, 1))
    @test quality_warnings(pair_quality(flat)) == ["the early image has one value throughout the sample windows"]

    # The report names both scenes and their numbers.
    r = quality_report("SCENE_A", "SCENE_B", pair_quality(dark))
    @test occursin("early SCENE_A median", r) && occursin("late SCENE_B median", r) && occursin("dB", r)
    close(good); close(dark); close(empty); close(flat)
end

@testset "a pair that resolves nothing ends the run" begin
    none = Correlation(fill(NaN32, 4, 4), fill(NaN32, 4, 4), zeros(UInt16, 4, 4), falses(4, 4), ones(Int32, 4, 4), 1.0, 1.0)
    chip = zeros(UInt16, 4, 4)
    chip[1, 1] = 32                      # one resolved point is a result
    some = Correlation(fill(NaN32, 4, 4), fill(NaN32, 4, 4), chip, falses(4, 4), ones(Int32, 4, 4), 1.0, 1.0)
    q = PairQuality(ImageSample(0.00306, 1.0, false), ImageSample(0.166, 1.0, false), 34.7)
    @test check_resolved(some, "A", "B", q) === nothing
    @test_throws "no point of the pair resolved, so no product is written. early A median 0.00306, data in 100.0% of samples; late B median 0.166, data in 100.0% of samples; brightness difference 34.7 dB" check_resolved(none, "A", "B", q)
    q = PairQuality(ImageSample(NaN, 0.0, false), ImageSample(0.166, 1.0, false), NaN)
    @test_throws "no point of the pair resolved, so no product is written; the early image has no data in any sample window. early A median n/a" check_resolved(none, "A", "B", q)
end
