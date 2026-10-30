# The floating-ice correction's geometry, checked against first principles rather than against itself.
#
# The pair is a real Sentinel-1 IW2 pair, the orbits of S1A 2015-11-20 and 2015-12-14 over the Antarctic
# Peninsula (`data/s1_orbit/rdr2rdr.json`, copied from ImagePairGeometry.jl's isce3-checked fixture), so
# the geometry is a real side-looking one at Antarctic latitude rather than a toy. Every check builds a
# ground point, moves it the way the ocean moves floating ice, solves the acquisitions with the
# zero-Doppler solver production uses, and asserts that the correction removes exactly the change the
# vertical motion caused — in range for a radar-geometry pair, and as geocoded map position for a NISAR
# L2 one — independently of the formulas in `range_offset_correction` and `geocoded_shift`.

using JSON3
using Random: Xoshiro
using Statistics: median
import AutoRIFT
using AutoriftItsLiveProduction: range_offset_correction, cos_incidence, _interp_linear, grounding_taper,
                                 grounding_weights, distance_to, GROUNDING_TAPER_LENGTH
import ImagePairGeometry as IPG
import AutoriftItsLiveProduction as A

const S1FX = JSON3.read(read(joinpath(@__DIR__, "data", "s1_orbit", "rdr2rdr.json"), String))

function s1_coordinate(o)
    orbit = IPG.Orbit(; time = collect(Float64.(o.time)), position = [IPG.SVector{3,Float64}(p...) for p in o.pos],
                      velocity = [IPG.SVector{3,Float64}(v...) for v in o.vel])
    kw = (; orbit, starting_range = Float64(o.r0), dr = Float64(o.dr), sensing_start = Float64(o.t0),
          prf = Float64(o.prf), nsamples = 24845, nlines = 1500, look_side = IPG.LookRight,
          wavelength = Float64(o.wavelength))
    return IPG.RadarCoordinate(; kw..., incidence_angle = IPG.incidence_angle(; kw...))
end

const S1REF = s1_coordinate(S1FX.ref)
const S1SEC = s1_coordinate(S1FX.sec)
const EL = IPG.Ellipsoid()

# Slant range and azimuth time at which `c` images the ECEF point `xyz`.
function solve(c, xyz)
    t0 = IPG.midtime(c)
    pos, vel = IPG.interpolate(c.orbit, t0)
    p = IPG.geo2rdr(c.orbit, xyz, t0, t0, pos, vel, 60)
    return p.range, p.aztime
end

# The ground point the reference images at range sample `s` of its middle line, at ellipsoidal height `h`.
function ground(c, s, h)
    az = c.sensing_start + (c.nlines / 2) / c.prf
    return IPG.rdr2geo(c.orbit, EL, az, c.starting_range + s * c.dr; height = h, wavelength = c.wavelength,
                       side = c.look_side)
end

xyz(llh) = IPG.lonlat_to_xyz(EL, llh)

# The early acquisition's range offset a correlator measures for a ground point that moves from `a` to `b`:
# where the late acquisition sees `b`, minus where the early one sees `a`, in early-grid samples, with the
# orbit difference removed by the same point's offset when it does not move (what the geometry's prior
# accounts for).
function measured_dx(a, b)
    r_early, _ = solve(S1REF, xyz(a))
    r_late, _ = solve(S1SEC, xyz(b))
    r_late0, _ = solve(S1SEC, xyz(a))
    return ((r_late - r_early) - (r_late0 - r_early)) / S1REF.dr
end

@testset "the real pair is a Sentinel-1 IW pair over Antarctic ice" begin
    llh = ground(S1REF, 12000, 0.0)
    @test -72 < rad2deg(llh[2]) < -64           # Antarctic Peninsula latitude
    @test S1REF.dr ≈ 2.3295621147 atol = 1.0e-9   # IW slant-range sample spacing
    θ = acos(cos_incidence(S1REF)[2][end ÷ 2])
    @test 30 < rad2deg(θ) < 46
end

@testset "a purely vertical displacement is removed exactly" begin
    for s in (500, 6000, 12000, 18000, 24000), h in (0.0, 50.0), dz in (-1.5, -0.3, 0.4, 2.0)
        a = ground(S1REF, s, h)
        b = IPG.SVector(a[1], a[2], a[3] + dz)
        dx = measured_dx(a, b)
        # The sign alone: raising the surface brings it nearer, a negative range offset.
        @test sign(dx) == -sign(dz)
        samples, cosθ = cos_incidence(S1REF)
        corr = range_offset_correction(dz, _interp_linear(samples, cosθ, Float64(s)), S1REF.dr)
        # What is left is cos θ interpolated between 65 swath samples rather than solved at this one:
        # 79 µm of range per metre of dz at worst on this swath, 1e-4 of the correction.
        @test abs(dx + corr) < 2.0e-4 * abs(dz) / S1REF.dr
    end
end

@testset "horizontal motion survives the correction unchanged" begin
    # A point moving 20 m east and 10 m north while rising 1 m: after the correction the offset equals
    # that of the same horizontal motion with no rise, so the ice's velocity is not altered.
    for s in (2000, 12000, 22000)
        a = ground(S1REF, s, 0.0)
        e = IPG.SVector(-sin(a[1]), cos(a[1]), 0.0)
        n = IPG.SVector(-sin(a[2]) * cos(a[1]), -sin(a[2]) * sin(a[1]), cos(a[2]))
        moved = xyz(a) + 20.0 * e + 10.0 * n
        flat = IPG.xyz_to_lonlat(EL, moved)
        flat = IPG.SVector(flat[1], flat[2], a[3])
        risen = IPG.SVector(flat[1], flat[2], a[3] + 1.0)
        samples, cosθ = cos_incidence(S1REF)
        corr = range_offset_correction(1.0, _interp_linear(samples, cosθ, Float64(s)), S1REF.dr)
        @test measured_dx(a, risen) + corr ≈ measured_dx(a, flat) atol = 1.0e-4
    end
end

@testset "cos θ is the ellipsoid normal's, across the whole swath" begin
    samples, cosθ = cos_incidence(S1REF)
    @test issorted(cosθ; rev = true)                  # incidence grows with range
    # The incidence angle independently: the angle between the line of sight and the geodetic normal at
    # the ground point, solved at 9 range positions and compared with the interpolated table.
    for s in range(0.0, S1REF.nsamples - 1; length = 9)
        llh = ground(S1REF, s, 0.0)
        p = xyz(llh)
        t0 = IPG.midtime(S1REF)
        pos0, vel0 = IPG.interpolate(S1REF.orbit, t0)
        sat = IPG.geo2rdr(S1REF.orbit, p, t0, t0, pos0, vel0, 60).position
        los = (sat - p) / sqrt(sum(abs2, sat - p))
        nrm = IPG.SVector(cos(llh[2]) * cos(llh[1]), cos(llh[2]) * sin(llh[1]), sin(llh[2]))
        @test _interp_linear(samples, cosθ, s) ≈ sum(los .* nrm) atol = 2.0e-5
    end
end

@testset "the correlator's dx is late minus early, in range samples" begin
    # `measured_dx` above assumes the convention `correlate` relies on: with the images passed late
    # first, `dx` is a feature's column in the late image minus its column in the early one. A surface
    # that rose comes back nearer, at a lower range sample, so it must read negative here.
    rng = Xoshiro(20261029)
    early = Float32.(randn(rng, 256, 256))
    for _ in 1:2
        early .= (early .+ circshift(early, (1, 0)) .+ circshift(early, (0, 1))) ./ 3
    end
    late = circshift(early, (0, -2))                       # every feature two range samples nearer
    # Late first, as `correlate` passes them.
    r = AutoRIFT.autorift(late, early; chip_size = 32, chip_size_max = 32, grid_spacing = 16, search_radius = 8)
    @test median(filter(!isnan, r.dx)) ≈ -2 atol = 0.05
    @test median(filter(!isnan, r.dy)) ≈ 0 atol = 0.05
end

# Where a geocoded (GSLC) image places a ground point that sits at height `h_true` when imaged: the point on
# the fixed geocoding surface `h_dem` with the same range and zero-Doppler time.
function geocoded_position(c, llh_true, h_dem)
    t0 = IPG.midtime(c)
    pos0, vel0 = IPG.interpolate(c.orbit, t0)
    q = IPG.geo2rdr(c.orbit, xyz(llh_true), t0, t0, pos0, vel0, 60)
    return IPG.rdr2geo(c.orbit, EL, q.aztime, q.range; height = h_dem, wavelength = c.wavelength,
                       side = c.look_side)
end

@testset "a geocoded image's apparent motion from a rise is removed exactly" begin
    # NISAR L2 lies on a map grid made by geocoding each acquisition onto a DEM that does not move with
    # the tide, so a rise of dz shows as horizontal motion toward the sensor. The test geocodes the early
    # point and the late (moved and risen) point onto the fixed surface, as the GSLC processor does, and
    # checks that `geocoded_shift` predicts the late point's displacement exactly, in polar stereographic
    # metres as the product grid measures it.
    to_ps = A.FGP.Transformation(A.FGP.GeographicCRS(), A.FGP.native_crs(3031))
    ps(llh) = to_ps(rad2deg(llh[1]), rad2deg(llh[2]))
    for s in (1000, 12000, 23000), dz in (-1.5, 0.6, 2.0)
        P = ground(S1REF, s, 50.0)
        # 30 m east, 20 m north between acquisitions, as ice moving about 1 km/yr over 12 days.
        e, n = A._east(P[1]), A._north(P[1], P[2])
        moved = IPG.xyz_to_lonlat(EL, xyz(P) .+ 30.0 .* e .+ 20.0 .* n)
        early = geocoded_position(S1REF, P, 50.0)
        late = geocoded_position(S1REF, IPG.SVector(moved[1], moved[2], 50.0 + dz), 50.0)
        flat = geocoded_position(S1REF, IPG.SVector(moved[1], moved[2], 50.0), 50.0)
        t0 = IPG.midtime(S1REF)
        pos0, vel0 = IPG.interpolate(S1REF.orbit, t0)
        r = IPG.geo2rdr(S1REF.orbit, xyz(P), t0, t0, pos0, vel0, 60)
        δ = A.geocoded_shift(dz, -r.look / sqrt(sum(abs2, r.look)), r.velocity, e, n)
        # The correction moves the late point back by δ; what remains is the horizontal motion alone.
        corrected = IPG.xyz_to_lonlat(EL, xyz(late) .- δ)
        x_corr, y_corr = ps(corrected)
        x_flat, y_flat = ps(flat)
        # 0.05 mm per metre of dz at worst on this swath: the first-order shift against the exact solve.
        @test hypot(x_corr - x_flat, y_corr - y_flat) < 1.0e-4 * abs(dz)
        # And it was worth correcting: the uncorrected late point is off by about dz·cot θ.
        x_late, y_late = ps(late)
        @test hypot(x_late - x_flat, y_late - y_flat) > 0.9 * abs(dz)
        # The shift points toward the sensor for a rise: along the horizontal line of sight.
        los = -r.look / sqrt(sum(abs2, r.look))
        @test sign(sum(δ .* los)) == sign(dz)
    end
end

@testset "a map shift becomes offsets in the images' pixels, not the product grid's" begin
    # The map path correlates the GSLCs on their own grid (NISAR: 2.5 m × 5 m), and `dx`/`dy` are in those
    # pixels; the 120 m product grid's spacing would make the correction 48 and 24 times too small.
    nisar = (2.5, -5.0)
    @test A.map_offset_pixels(1.3, -0.65, nisar) == (0.52, -0.13)
    # Through the velocity operator this is the velocity the shift implies over 12 days, whatever the grid.
    for sp in (nisar, (10.0, -10.0), (15.0, -15.0))
        px, py = A.map_offset_pixels(1.3, -0.65, sp)
        @test px * sp[1] ≈ 1.3
        @test py * abs(sp[2]) ≈ -0.65
    end
end

@testset "the correlator's map dx and DY_SIGN·dy are east and north" begin
    # The map path writes `dx` and `DY_SIGN·dy` (`Correlation.dy`) as motion along +x and +y of a north-up
    # grid; `tide_correct(::Setup, …)` subtracts the geocoded shift in that frame.
    rng = Xoshiro(3)
    early = Float32.(randn(rng, 256, 256))
    for _ in 1:2
        early .= (early .+ circshift(early, (1, 0)) .+ circshift(early, (0, 1))) ./ 3
    end
    late = circshift(early, (-3, 2))                     # 2 columns east, 3 rows north
    r = AutoRIFT.autorift(late, early; chip_size = 32, chip_size_max = 32, grid_spacing = 16, search_radius = 8)
    @test median(filter(!isnan, r.dx)) ≈ 2 atol = 0.05
    @test A.DY_SIGN * median(filter(!isnan, r.dy)) ≈ 3 atol = 0.05
end

@testset "the taper is zero at the grounding line and one past it" begin
    @test grounding_taper(0.0) == 0
    @test grounding_taper(GROUNDING_TAPER_LENGTH) == 1
    @test grounding_taper(10 * GROUNDING_TAPER_LENGTH) == 1
    @test grounding_taper(GROUNDING_TAPER_LENGTH / 2) == 0.5
    d = range(0, GROUNDING_TAPER_LENGTH; length = 101)
    @test issorted(grounding_taper.(d))
    # Zero slope at both ends: no kink where the corrected and uncorrected regions meet.
    h = 1.0e-3
    @test grounding_taper(h) / h < 1.0e-5
    @test (1 - grounding_taper(GROUNDING_TAPER_LENGTH - h)) / h < 1.0e-5
end

@testset "distance to grounded ice is Euclidean" begin
    f = falses(41, 31)
    f[11, 16] = true
    d = distance_to(f, 120.0)
    for i in axes(f, 1), j in axes(f, 2)
        @test d[i, j] ≈ 120 * hypot(i - 11, j - 16)
    end
    @test all(>(1.0e10), distance_to(falses(5, 5), 120.0))
end

@testset "neither the ice front nor the ocean beyond it is tapered" begin
    # Columns 1–5 grounded, 6–60 floating, 61–80 ocean; 120 m pixels, one-pixel margin on each side.
    nx, ny = 82, 12
    water = falses(nx, ny)
    water[7:end, :] .= true
    w = grounding_weights(water, trues(nx, ny), 1, 120.0)
    @test all(iszero, w[1:5, :])                                   # grounded: no correction
    @test all(w[6, :] .≈ grounding_taper(120.0))                   # first floating pixel, 120 m out
    @test all(isone, w[6 + ceil(Int, GROUNDING_TAPER_LENGTH / 120):60, :])
    @test all(isone, w[60:end, :])                                 # the ice front and the open water
end

@testset "ocean beside a coast is tapered like floating ice" begin
    # Columns 1–40 land, 41–80 ocean: the correction rises over 2 km from the coastline.
    water = falses(82, 3)
    water[42:end, :] .= true
    w = grounding_weights(water, trues(82, 3), 1, 120.0)
    @test all(iszero, w[1:40, 1])
    @test w[41, 1] ≈ grounding_taper(120.0)
    @test all(isone, w[41 + ceil(Int, GROUNDING_TAPER_LENGTH / 120):80, 1])
end

@testset "the edge of the masks' data is not a grounding line" begin
    # Columns 1–40 ocean, 41–80 beyond the masks: the ocean is not tapered toward the unknown side, and the
    # unknown side is not corrected.
    water = falses(82, 3)
    water[1:41, :] .= true
    known = falses(82, 3)
    known[1:41, :] .= true
    w = grounding_weights(water, known, 1, 120.0)
    @test all(isone, w[1:40, 1])
    @test all(iszero, w[41:80, 1])
end
