# The geogrid repackaging: an `ImagePairGeometry.PairGeometry` becomes the product structs.
#
# Each of these is a rename and a transpose rather than a computation, so every assertion checks that a
# band arrived under the right name, in `[row, col]` order, with the geometry's sentinel turned into
# `NaN`.

using ItsLiveOffsetProduction
using ImagePairGeometry
using ImagePairGeometry: nodata_from
using StaticArrays: SVector
using Test

"""A geometry over a grid the pair covers, with every input band present."""
function ipg_case(; csminy = 360.0, ssm = 1.0)
    grid = MapGrid(geotransform = (295000.0, 120.0, 0.0, 7805000.0, 0.0, -120.0),
                   size = (200, 200), crs = 32624)
    # **Deliberately non-square.** A square footprint gives a square grid window, and a square window
    # hides a transposed band: the shapes agree, the correlator runs, and the answer is wrong by two
    # thirds of its agreement. Every assertion below that compares against a band therefore has a shape
    # to disagree about.
    fp = ImageFootprint(origin = (300000.0, 7800000.0), spacing = (30.0, -30.0),
                        size = (400, 560))
    pair = coregister(fp, fp; dt = 91 * 86400.0)
    win = grid_window(grid, footprint_bounds(IdentityTransform(), pair.coordinate))
    n = size(win)
    inputs = GeometryInputs(dem = fill(500.0, n), dhdx = fill(0.02, n), dhdy = fill(-0.01, n),
                            vx = fill(120.0, n), vy = fill(-80.0, n),
                            srx = fill(400.0, n), sry = fill(300.0, n),
                            csminx = fill(240.0, n), csminy = fill(csminy, n),
                            csmaxx = fill(480.0, n), csmaxy = fill(720.0, n),
                            ssm = fill(ssm, n))
    r = pairgeometry(grid, pair, inputs; window = win, nodata = nodata_from(-32767.0))
    return r, grid, pair, win
end

# One geometry at the default inputs, shared by the testsets that do not vary them: building it runs
# the whole geometry kernel over the window.
const IPG_R, IPG_GRID, IPG_PAIR, IPG_WIN = ipg_case()
const SENTINEL = Int32(-32767)

@testset "coefficients and reference_velocity repackage the geogrid result" begin
    c = ItsLiveOffsetProduction.coefficients(IPG_R)
    @test c isa GeogridCoefficients
    # `[x,y]` to `[row,col]`, the same swap `pointset` makes.
    @test size(c.offset2vx_1) == reverse(size(IPG_R))
    @test size(c.scale_factor_1) == reverse(size(IPG_R))
    # Projected: no radar-only band. See "coefficients on a radar PairGeometry" below for the
    # radar-only bands themselves.
    @test c.offset2vr === nothing
    @test c.offset2va === nothing

    valid = findall(!=(SENTINEL), permutedims(IPG_R.location_x))
    invalid = findall(==(SENTINEL), permutedims(IPG_R.location_x))
    @test all(k -> c.offset2vx_1[k] == permutedims(IPG_R.off2vx_dx)[k], valid)
    @test all(k -> c.scale_factor_1[k] == permutedims(IPG_R.scale_x)[k], valid)
    # A geometry-invalid point becomes `NaN`, not the raw `-32767.0` sentinel: `ItsLiveAutoRIFT`'s
    # own helpers (`_finite_median`, every velocity-conversion formula) test `isnan`, not the value.
    @test all(k -> isnan(c.offset2vx_1[k]), invalid)
    @test all(k -> isnan(c.scale_factor_1[k]), invalid)

    n = size(IPG_WIN)
    inputs = GeometryInputs(dem = fill(500.0, n), dhdx = fill(0.02, n), dhdy = fill(-0.01, n),
                            vx = fill(120.0, n), vy = fill(-80.0, n),
                            srx = fill(400.0, n), sry = fill(300.0, n),
                            csminx = fill(240.0, n), csminy = fill(360.0, n),
                            csmaxx = fill(480.0, n), csmaxy = fill(720.0, n), ssm = fill(1.0, n))
    refv, mask = ItsLiveOffsetProduction.reference_velocity(inputs)
    @test refv isa ReferenceVelocity
    @test size(refv.vx) == reverse(size(IPG_R)) == size(mask)
    @test eltype(refv.vx) == Float32
    # Repackaged, not masked: the raw raster value passes straight through, transposed only — unlike
    # `coefficients`, since `inputs`' own sentinel is whatever its file declares, not `PairGeometry`'s
    # uniform one.
    @test all(==(Float32(120.0)), refv.vx)
    @test all(==(Float32(-80.0)), refv.vy)
    @test all(mask)   # `ssm = 1.0` everywhere in this fixture

    # No `ssm` at all: every point is stable-mask-false rather than an error, since a caller who
    # never supplied a stable-surface mask has that mask empty, not missing.
    inputs_no_ssm = GeometryInputs(dem = fill(500.0, n), dhdx = fill(0.02, n), dhdy = fill(-0.01, n),
                                   vx = fill(120.0, n), vy = fill(-80.0, n))
    _, mask_no_ssm = ItsLiveOffsetProduction.reference_velocity(inputs_no_ssm)
    @test size(mask_no_ssm) == reverse(size(IPG_R))
    @test !any(mask_no_ssm)

    @test_throws "vx" ItsLiveOffsetProduction.reference_velocity(GeometryInputs(dem = fill(500.0, n)))
end

"""A minimal radar `PairGeometry`, built directly rather than through `pairgeometry`.

`coefficients` only reads `off2vx_dr`/`off2vy_dr` and the coordinate's type — it runs no radar
solve — so a self-consistent orbit and a converging `geo2rdr` are not needed to test it, only a
`RadarCoordinate` that type-checks. Building one this way avoids the reference-fixture machinery
`radar_geogrid.jl` in ImagePairGeometry.jl's own test suite needs to test the solve itself, which is
not what this file is testing.
"""
function radar_pair_geometry()
    ts = [0.0, 10.0, 20.0, 30.0]
    orb = Orbit(; time = ts, position = [SVector(1.0 + i, 2.0 + i, 3.0 + i) for i in 0:3],
               velocity = [SVector(0.1, 0.1, 0.1) for _ in 0:3])
    coord = RadarCoordinate(; starting_range = 800_000.0, dr = 2.3, sensing_start = 0.0,
                            prf = 1717.13, nsamples = 10, nlines = 10, look_side = LookRight,
                            wavelength = 0.055, orbit = orb, incidence_angle = deg2rad(35))
    n = (4, 5)   # (x, y), deliberately non-square for the same reason `ipg_case` is
    sentinel32, sentinel64 = Int32(-32767), -32767.0
    loc = fill(sentinel32, n)
    off = fill(sentinel32, n)
    ss = fill(sentinel32, n)
    f = fill(sentinel64, n)
    dr = fill(sentinel64, n)
    dr[1, 1], dr[2, 2] = 2.5, 3.5   # two valid points, so both "present" and "absent" are exercised
    gt = (0.0, 120.0, 0.0, 0.0, 0.0, -120.0)
    return PairGeometry(loc, loc, off, off, off, off, off, off, off, off, ss,
                        f, f, f, f, dr, dr, f, f, gt, nothing, CartesianIndices((1:4, 1:5)),
                        nodata_from(sentinel64), coord)
end

@testset "coefficients on a radar PairGeometry" begin
    g = radar_pair_geometry()
    @test g.coordinate isa RadarCoordinate
    c = ItsLiveOffsetProduction.coefficients(g)
    @test size(c.offset2vr) == size(c.offset2va) == reverse(size(g))

    valid = permutedims(g.off2vx_dr) .!= -32767.0
    @test count(valid) == 2
    @test c.offset2vr[valid] == permutedims(g.off2vx_dr)[valid]
    @test c.offset2va[valid] == permutedims(g.off2vy_dr)[valid]
    @test all(isnan, c.offset2vr[.!valid])
    @test all(isnan, c.offset2va[.!valid])
end

@testset "image_location is a plain transposed repackaging" begin
    loc_x, loc_y = ItsLiveOffsetProduction.image_location(IPG_R)
    @test size(loc_x) == size(loc_y) == reverse(size(IPG_R))
    valid = findall(!=(SENTINEL), permutedims(IPG_R.location_x))
    invalid = findall(==(SENTINEL), permutedims(IPG_R.location_x))
    @test all(k -> loc_x[k] == permutedims(IPG_R.location_x)[k], valid)
    @test all(k -> loc_y[k] == permutedims(IPG_R.location_y)[k], valid)
    @test all(k -> isnan(loc_x[k]), invalid)
    @test all(k -> isnan(loc_y[k]), invalid)
end

