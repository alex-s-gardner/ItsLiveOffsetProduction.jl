# The pieces of the floating-ice correction that run without a network: the model's constants sampled as
# pyTMD samples them, ERA5's hours and grid, and the range change a raised surface makes.

using Dates: DateTime, Minute
using JSON3
using AutoriftItsLiveProduction: TideConstants, tide_point, tide_height, TIDE_BANDS, CATS_PROJECTION,
                                 _era5_hours, _era5_at, _interp_linear, _hour_string
import ImagePairGeometry as IPG

@testset "tide model constants are bilinear and flexure-weighted" begin
    # A 3×3 model window, 2 km cells, flexure 1 except one grounded corner cell.
    bands = zeros(Float32, 3, 3, length(TIDE_BANDS))
    bands[:, :, 1] .= 1        # hRe_m2
    bands[:, :, 2] .= 0.5      # hIm_m2
    bands[:, :, 21] .= 1
    bands[1, 1, 21] = 0
    k = TideConstants(bands, (0.0, 0.0), (2000.0, -2000.0))
    # At a cell centre the value is that cell's, flexure applied.
    z, f = tide_point(k, 2000.0, -2000.0)
    @test z[1] == complex(1.0, -0.5) && f == 1
    # Midway between the grounded cell and its three neighbours: a quarter of the weight is grounded.
    z, f = tide_point(k, 1000.0, -1000.0)
    @test f ≈ 0.75 && z[1] ≈ complex(0.75, -0.375) && all(iszero, z[2:end])
    @test_throws "outside the tide model window" tide_point(k, 5000.0, -1000.0)
end

@testset "heights from the converted model match pyTMD on the original" begin
    # 4×4-cell windows of `CATS2008_v2023_z.tif` around twelve ice-shelf points, as `tide_constants`
    # reads them, against pyTMD's prediction from the model's own netCDF with its flexure applied
    # (`tools/tidal/tide_pytmd_cats.csv`): the conversion, the band order, the sign of `hIm`, the
    # projection and the interpolation together.
    for p in JSON3.read(read(joinpath(@__DIR__, "data", "tide_pytmd_cats.json"), String))
        bands = reshape(Float32.(p.bands), Tuple(p.size))
        k = TideConstants(bands, Tuple(Float64.(p.origin)), Tuple(Float64.(p.spacing)))
        x, y = CATS_PROJECTION(Float64(p.lon), Float64(p.lat))
        z, f = tide_point(k, x, y)
        @test f ≈ p.flexure atol = 1.0e-9
        for (t, h) in zip(p.times, p.height)
            @test tide_height(z, DateTime(t[1:23])) ≈ h atol = 1.0e-6
        end
    end
end

@testset "CATS2008_v2023 is on its own projection" begin
    # pyproj's `+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=-70 +datum=WGS84`, in metres, for two of the
    # reference points in `tools/tidal/tide_pytmd_cats.csv`.
    x, y = CATS_PROJECTION(-160.0, -81.0)
    @test x ≈ -979_806.32689977188 atol = 1.0e-6
    @test abs(y) < 1.0e-6
    x, y = CATS_PROJECTION(71.0, -70.5)
    @test x ≈ 1_345_856.3267148638 atol = 1.0e-6
    @test y ≈ -1_661_994.1509648934 atol = 1.0e-6
end

@testset "ERA5 hours and grid" begin
    h, f = _era5_hours(DateTime(1900, 1, 1, 5, 30))
    @test (h, f) == (5, 0.5)
    @test _era5_hours(DateTime(2020, 1, 13, 3))[1] == 1_052_187
    # Rows from 60°S, every longitude: pressure = 1000 lon + lat index, so bilinear is exact.
    p = Float32[1000 * i + j for i in 0:1439, j in 0:9]
    row0 = 600
    @test _era5_at(p, row0, 10.0, -60.5) ≈ 1000 * 40 + 2
    # West longitudes wrap to the grid's 0–360°, and the last column interpolates with the first.
    @test _era5_at(p, row0, -0.125, -60.0) ≈ (1000 * 1439 + 0) / 2
    @test_throws "outside the ERA5 rows" _era5_at(p, row0, 0.0, -75.0)
    @test _hour_string(DateTime(2026, 10, 19, 7) + Minute(59)) == "2026-10-19T07"
end

@testset "raising the surface shortens the range by cos θ" begin
    # A circular orbit, as ImagePairGeometry's tests build one.
    radius, incl = 7.0e6, deg2rad(98.0)
    omega = sqrt(3.986004418e14 / radius^3)
    ts = collect(-2000.0 .+ 10.0 .* (0:400))
    pos = [IPG.SVector(radius * cos(omega * t), radius * sin(omega * t) * cos(incl), radius * sin(omega * t) * sin(incl))
           for t in ts]
    vel = [IPG.SVector(-radius * omega * sin(omega * t), radius * omega * cos(omega * t) * cos(incl),
                       radius * omega * cos(omega * t) * sin(incl)) for t in ts]
    orb = IPG.Orbit(; time = ts, position = pos, velocity = vel)
    el = IPG.Ellipsoid()
    az, rg = 60.0 + 200 / 486.4863103, 8.0e5 + 250 * 2.3295621147
    llh = IPG.rdr2geo(orb, el, az, rg; height = 0.0, wavelength = 0.05546576, side = IPG.LookRight)
    sat, v = IPG.interpolate(orb, az)
    r(h) = IPG.geo2rdr(orb, IPG.lonlat_to_xyz(el, IPG.SVector(llh[1], llh[2], h)), az, az, sat, v).range
    θ = IPG.incidence_angle(orb, el, az, rg; wavelength = 0.05546576, side = IPG.LookRight, zrange = (0.0,))
    @test r(1.0) - r(0.0) ≈ -cos(θ) atol = 1.0e-5
    @test _interp_linear([0.0, 10.0, 20.0], [1.0, 3.0, 2.0], 15.0) == 2.5
end
