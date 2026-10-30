# The pieces of the floating-ice correction that run without a network: the model's constants sampled as
# pyTMD samples them, ERA5's hours and grid, and the range change a raised surface makes.

using Dates: DateTime, Minute
using JSON3
using AutoriftItsLiveProduction: TideConstants, tide_point, tide_height, TIDE_BANDS, CATS_PROJECTION,
                                 GR1KM_PROJECTION, CATS_MODEL, GR1KM_MODEL, BEDMACHINE_MASK, Region,
                                 tide_setting, tide_constants, covers, mask_window,
                                 epsg_wkt as A_epsg_wkt,
                                 _era5_hours, _era5_at, _interp_linear, _hour_string
import ImagePairGeometry as IPG
import GDAL

# A GeoTIFF of `bands[x, y, b]` on EPSG `epsg` with each band described by `names[b]`, through GDAL itself.
function write_bands(path::String, bands::Array{Float32,3}, gt::NTuple{6,Float64}, names::Vector{String};
                     epsg::Int = 3413)
    nx, ny, nb = size(bands)
    d = GDAL.gdalcreate(GDAL.gdalgetdriverbyname("GTiff"), path, nx, ny, nb, GDAL.GDT_Float32, C_NULL)
    GDAL.gdalsetgeotransform(d, collect(gt))
    GDAL.gdalsetprojection(d, A_epsg_wkt(epsg))
    for b in 1:nb
        band = GDAL.gdalgetrasterband(d, b)
        a = bands[:, :, b]
        GDAL.gdalrasterio(band, GDAL.GF_Write, 0, 0, nx, ny, a, nx, ny, GDAL.GDT_Float32, 0, 0)
        GDAL.gdalsetdescription(band, names[b])
    end
    GDAL.gdalclose(d)
    return path
end

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

@testset "heights from converted Gr1kmTM match pyTMD on the original" begin
    # As for CATS2008_v2023, from `Gr1kmTM_v1_z.tif` and pyTMD on the model's OTIS files, at points on six
    # Greenland ice tongues (`tools/tidal/pytmd_reference.py gr1km`): the OTIS conversion, the zero
    # long-period bands, the sign of `hIm`, the northern projection and the OTIS nodal corrections.
    for p in JSON3.read(read(joinpath(@__DIR__, "data", "tide_pytmd_gr1km.json"), String))
        bands = reshape(Float32.(p.bands), Tuple(p.size))
        k = TideConstants(bands, Tuple(Float64.(p.origin)), Tuple(Float64.(p.spacing)))
        x, y = GR1KM_PROJECTION(Float64(p.lon), Float64(p.lat))
        z, f = tide_point(k, x, y)
        @test f == 1
        @test all(iszero, z[9:10])                                  # mf, mm
        for (t, h) in zip(p.times, p.height)
            @test tide_height(z, DateTime(t[1:23])) ≈ h atol = 1.0e-6
        end
    end
end

@testset "each polar region has its own tide model and water masks" begin
    region(name, epsg) = Region(name, epsg, ntuple(_ -> "", Val(12)), "floating.tif", "ocean.tif")
    s = tide_setting(region("SPS", 3031), "")
    @test s.model == CATS_MODEL && s.projection === CATS_PROJECTION
    @test [(m.path, m.classes) for m in s.masks] == [("floating.tif", [0x01]), ("ocean.tif", [0x01])]
    g = tide_setting(region("NPS", 3413), "")
    @test g.model == GR1KM_MODEL && g.projection === GR1KM_PROJECTION
    # Greenland's floating ice and ocean come from BedMachine v7, not the region's own masks.
    @test [(m.path, m.classes) for m in g.masks] == [(BEDMACHINE_MASK, [0x00, 0x03])]
    @test tide_setting(region("NPS", 3413), "/vsimem/other.tif").model == "/vsimem/other.tif"
    @test tide_setting(region("N07", 32607), "") === nothing
    @test tide_setting(region("S19", 32719), "") === nothing
end

@testset "points beyond the model's extent are found, not interpolated" begin
    # A 4×3-cell model, 1 km cells, origin cell centre at (500, -500): a window over points partly off it
    # is clipped, and only points with four cell centres around them inside it are covered.
    path = "/vsimem/tide_extent.tif"
    write_bands(path, ones(Float32, 4, 3, length(TIDE_BANDS)), (0.0, 1000.0, 0.0, 0.0, 0.0, -1000.0),
                collect(TIDE_BANDS))
    k = tide_constants(path, [1500.0, 9000.0], [-1500.0, -1500.0])
    @test k !== nothing && size(k.bands, 1) == 4
    @test covers(k, 1500.0, -1500.0)
    @test !covers(k, 3600.0, -1500.0)                               # past the last cell centre
    @test tide_constants(path, [50_000.0], [-50_000.0]) === nothing # wholly outside
end

@testset "a class raster is read by map position, wherever its origin" begin
    # Product grid: 120 m pixels from (0, 0). Class raster: a 6×4 block starting 3 columns and 2 rows in,
    # value 10x + y at its pixel (x, y).
    grid = IPG.MapGrid((0.0, 120.0, 0.0, 0.0, 0.0, -120.0), (40, 40), 3413)
    path = "/vsimem/classes.tif"
    classes = Float32[10i + j for i in 1:6, j in 1:4]
    write_bands(path, reshape(classes, 6, 4, 1), (360.0, 120.0, 0.0, -240.0, 0.0, -120.0), ["classes"])
    # Window columns 5–6, rows 3–4 (one-based), grown by one pixel: grid columns 4–7, rows 2–5.
    m, inside = mask_window(path, grid, CartesianIndices((5:6, 3:4)), 1)
    @test size(m) == size(inside) == (4, 4)
    @test !inside[1, 1] && m[1, 1] == 0        # grid (4, 2) is above the raster's first row
    @test inside[1, 2] && m[1, 2] == 0x0b      # grid (4, 3) is raster pixel (1, 1)
    @test inside[4, 4] && m[4, 4] == 0x2b      # grid (7, 5) is raster pixel (4, 3)
    # A window wholly beside the raster is not covered.
    _, inside = mask_window(path, grid, CartesianIndices((30:31, 30:31)), 1)
    @test !any(inside)
    # A raster off the lattice, or on another projection, is refused.
    write_bands("/vsimem/off.tif", ones(Float32, 2, 2, 1), (50.0, 120.0, 0.0, 0.0, 0.0, -120.0), ["c"])
    @test_throws "not on the product grid's lattice" mask_window("/vsimem/off.tif", grid,
                                                                  CartesianIndices((1:2, 1:2)), 0)
    write_bands("/vsimem/south.tif", ones(Float32, 2, 2, 1), (0.0, 120.0, 0.0, 0.0, 0.0, -120.0), ["c"]; epsg = 3031)
    @test_throws "is on EPSG:3031; the product grid is on EPSG:3413" mask_window("/vsimem/south.tif", grid,
                                                                                  CartesianIndices((1:2, 1:2)), 0)
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
