# The correlator's per-sensor settings: the displacement filter and the image filter.

using AutoriftItsLiveProduction: dispfilt, DISPFILT_STRICT, DISPFILT_OPTICAL, correlator_filter, Job,
                                 SceneSource, WALLIS, HIGHPASS, WALLIS_FILL, NO_FILTER
using AutoriftItsLiveProduction.AutoRIFT: GardnerFilter
import AutoriftItsLiveProduction

job(platform, ref, sec) = Job(platform, platform, SceneSource(ref, ref, "", "", ""),
                              SceneSource(sec, sec, "", "", ""), ref * "_X_" * sec, "")

@testset "each sensor's displacement filter" begin
    for p in ("S2", "S1-SLC", "S1-BURST")
        @test dispfilt(p) == DISPFILT_STRICT
    end
    for p in ("L4", "L5", "L7", "L8", "L9", "AST", "NISAR-L1", "NISAR-L2")
        @test dispfilt(p) == DISPFILT_OPTICAL
    end
    @test (DISPFILT_STRICT.mad_scale, DISPFILT_STRICT.agree_tolerance) == (3.0, 0.15)
    @test (DISPFILT_OPTICAL.mad_scale, DISPFILT_OPTICAL.agree_tolerance) == (4.0, 0.15)
    # A sensor without its own entry gets the optical filter.
    @test dispfilt("RADARSAT-1") == DISPFILT_OPTICAL
    # The coefficients not set per sensor are AutoRIFT's defaults.
    d = GardnerFilter()
    for f in (DISPFILT_STRICT, DISPFILT_OPTICAL)
        @test (f.window, f.iterations, f.min_agree_fraction) == (d.window, d.iterations, d.min_agree_fraction)
    end
end

@testset "the image filter each pair is correlated with" begin
    s2 = "S2A_MSIL1C_20200626T204021_N0209_R014_T07VEG_20200626T235742"
    @test correlator_filter(job("S2", s2, s2)) === WALLIS
    l8 = "LC08_L1TP_009011_20200703_20200913_02_T1"
    @test correlator_filter(job("L8", l8, l8)) === HIGHPASS
    @test correlator_filter(job("L7", l8, "LE07_L1TP_061018_20120428_20200909_02_T1")) === WALLIS_FILL
    @test correlator_filter(job("L5", "LT05_L1TP_060018_19851028_20200918_02_T1", l8)) === NO_FILTER
end

@testset "a product names the version that made it" begin
    v = AutoriftItsLiveProduction.app_version()
    # Unbuilt, the checkout's `git describe`: a tag, the commits since it, and the commit.
    @test occursin(r"^v\d+\.\d+\.\d+(-\d+-g[0-9a-f]+)?(-dirty)?$", v)
end
