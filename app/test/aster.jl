# ASTER L1T (`aster.jl`): the scene ID and its band file, the CMR record's band link, the identity the product
# records, the saturation mask the correlator reads, and the stable-surface correction on a synthetic scene.

using AutoriftItsLiveProduction: cmr_granules, aster_acquisition, aster_date, granule_link,
                                 aster_scene_source, aster_identification, UnsaturatedMask, aster_correct,
                                 aster_fit_mask, infer_sensor, opendata_platform
import ItsLiveOffsetProduction
import AutoriftItsLiveProduction
using Dates: DateTime
using Statistics: median

const AST_UR = "AST_L1T_00407222024142636_20251005235025"
const AST_DIR = "https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/AST_L1T.004/" * AST_UR * "/" * AST_UR
const AST_CMR = """
{"hits": 1, "items": [{"umm": {"GranuleUR": "$AST_UR",
  "RelatedUrls": [{"URL": "$(AST_DIR)_VNIR_B01.tif", "Type": "GET DATA"},
                  {"URL": "$(AST_DIR)_VNIR_B03N.tif", "Type": "GET DATA"},
                  {"URL": "s3://lp-prod-protected/AST_L1T.004/$AST_UR/$(AST_UR)_VNIR_B03N.tif", "Type": "GET DATA VIA DIRECT ACCESS"},
                  {"URL": "$(AST_DIR).cmr.json", "Type": "EXTENDED METADATA"}]}}]}
"""

@testset "an ASTER granule's name, band and identity" begin
    @test aster_acquisition(AST_UR) == DateTime(2024, 7, 22, 14, 26, 36)
    @test aster_acquisition("AST_L1T_00310032002055404_20150425130109") == DateTime(2002, 10, 3, 5, 54, 4)
    @test aster_date(AST_UR) == "20240722"
    @test_throws "is not an ASTER L1T granule name" aster_acquisition("AST_L1T_0040722")
    g = only(cmr_granules(AST_CMR, [AST_UR]; catalog = "LP DAAC's catalog (CMR)"))
    @test granule_link(g, "_VNIR_B03N.tif") == AST_DIR * "_VNIR_B03N.tif"
    @test_throws "lists 0 https:// _VNIR_B02.tif links, not one" granule_link(g, "_VNIR_B02.tif")
    @test_throws "LP DAAC's catalog (CMR) lists no granule AST_L1T_00407222024142636_20251005235026" cmr_granules(
        AST_CMR, ["AST_L1T_00407222024142636_20251005235026"]; catalog = "LP DAAC's catalog (CMR)")
    # A band file names its scene; any other ASTER file is refused.
    s = aster_scene_source("/data/" * AST_UR * "_VNIR_B03N.tif")
    @test (s.name, s.band, s.mtl, s.ang) == (AST_UR, "/data/" * AST_UR * "_VNIR_B03N.tif", "", "")
    @test_throws "the band the chain correlates" aster_scene_source("/data/" * AST_UR * "_VNIR_B02.tif")
    @test infer_sensor(AST_UR) == "AST"
    @test infer_sensor(AST_DIR * "_VNIR_B03N.tif") == "AST"
    @test opendata_platform("AST") == "aster"
    id = aster_identification(AST_UR)
    @test (id.mission, id.satellite, id.sensor, id.correction_level) == ("T", "1", "AST", "L1T")
    @test id.acquisition_time.datetime == DateTime(2024, 7, 22, 14, 26, 36)
end

@testset "saturated samples and their neighbours are invalid to the correlator" begin
    a = fill(UInt8(100), 12, 10)
    a[6, 5] = 0xff
    a[1, 1] = 0x00
    m = UnsaturatedMask(a)
    @test size(m) == (12, 10)
    full = m[1:12, 1:10]
    @test !full[1, 1]                                   # no data
    @test all(!full[i, j] for i in 4:8, j in 3:7)       # within 2 of the saturated sample
    @test full[3, 5] && full[9, 5] && full[6, 2] && full[6, 8]
    @test count(!, full) == 1 + 25
    # A window reads its margin, so its values agree with the whole-mask read.
    @test m[5:9, 6:10] == full[5:9, 6:10]
    @test m[6, 7] == full[6, 7]
end

@testset "the ASTER correction removes a smooth distortion and along-track jitter, keeping the constant" begin
    nr, nc = 300, 280
    x = collect(500_000.0 .+ 120.0 .* (0:(nc - 1)))
    y = collect(7_600_000.0 .- 120.0 .* (0:(nr - 1)))
    frame = ((0.0, 1.0), (1.0, 0.0))                     # the track runs north, detectors east
    S = [y[i] for i in 1:nr, j in 1:nc]
    X = [x[j] for i in 1:nr, j in 1:nc]
    jitter = 0.4 .* sin.(2π .* S ./ 4700.0)
    tilt = 0.3 .* (X .- x[1]) ./ 1e4
    dx = Float32.(0.25 .+ jitter .+ tilt)                 # cross-track offset (east): constant + jitter + ramp
    dy = Float32.(fill(-0.1, nr, nc))
    ssm = trues(nr, nc)
    cx, cy = aster_correct(dx, dy, ssm, x, y, frame)
    # Jitter and ramp removed; the stable cells' median level is kept for the product's stable shift.
    @test maximum(abs.(cx .- median(dx))) < 0.03
    @test maximum(abs.(cy .+ 0.1f0)) < 0.01
    # The fit's cells: the stable mask when it holds a measurement, else the slowest quarter of the reference.
    refv = ItsLiveOffsetProduction.ReferenceVelocity(Float32[j for i in 1:nr, j in 1:nc], zeros(Float32, nr, nc))
    @test aster_fit_mask(ssm, refv, dx) == ssm
    slow = aster_fit_mask(falses(nr, nc), refv, dx)
    @test all(slow[:, 1:70]) && !any(slow[:, 72:end])
end

@testset "the ASTER correction is refused where it cannot be trusted" begin
    nr, nc = 300, 280
    x = collect(500_000.0 .+ 120.0 .* (0:(nc - 1)))
    y = collect(7_600_000.0 .- 120.0 .* (0:(nr - 1)))
    frame = ((0.0, 1.0), (1.0, 0.0))
    dx = Float32.([0.25 + 0.4 * sin(2π * y[i] / 4700.0) for i in 1:nr, j in 1:nc])
    dy = fill(-0.1f0, nr, nc)
    # A quadrant without stable ground: the fit would extrapolate its cubics there.
    empty_corner = trues(nr, nc); empty_corner[151:end, 141:end] .= false
    @test_throws AutoriftItsLiveProduction.NoProduct aster_correct(dx, dy, empty_corner, x, y, frame)
    @test_throws "a quadrant of the scene holds 0 stable cells, fewer than the 1000" aster_correct(
        dx, dy, empty_corner, x, y, frame)
    # The floor is on each quadrant: one cell above it in the emptiest is accepted, one below refused.
    function corner(n)
        m = trues(nr, nc); m[151:end, 141:end] .= false
        m[CartesianIndices((151:300, 141:280))[1:n]] .= true
        return m
    end
    @test aster_correct(dx, dy, corner(1000), x, y, frame) isa Tuple
    @test_throws "holds 999 stable cells" aster_correct(dx, dy, corner(999), x, y, frame)
    # A residual the fit's terms cannot express (a 1.6 px step across the scene): refused on held-out ground.
    stepped = dx .+ Float32.([isodd(div(i - 1, 40)) ? 1.6 : 0.0 for i in 1:nr, j in 1:nc])
    @test_throws "px on held-out stable ground, more than 1.2" aster_correct(stepped, dy, trues(nr, nc), x, y, frame)
    # Stable ground too sparse for any 1 km block to hold 10 cells: the residual cannot be measured.
    sparse = falses(nr, nc); sparse[1:3:end, 1:3:end] .= true
    @test_throws "cannot be measured" aster_correct(dx, dy, sparse, x, y, frame)
    # No measured offsets at all.
    @test_throws "no measured offsets" aster_correct(fill(NaN32, nr, nc), dy, trues(nr, nc), x, y, frame)
end

