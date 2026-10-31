# The tide correction's published inputs, read where production reads them: the tide models and the
# BedMachine mask in s3://its-live-data/tidal/ and ERA5 in Google's ARCO-ERA5 store. Not part of
# `runtests.jl`, which needs no network; run by `.github/workflows/TideData.yml`, or by hand:
#
#   julia --project=app -e 'import Pkg; Pkg.test(; test_args = ["tidenetwork"])'
#
# Each check compares against values computed independently: the heights against pyTMD's on the models'
# original files (the fixtures `runtests.jl` also checks against windows of the same COGs), the pressures
# against xarray's interpolation of the same ERA5 store.

using Test
using Dates: DateTime, Day, Hour, UTC, now
using JSON3
using AutoriftItsLiveProduction: CATS_MODEL, GR1KM_MODEL, BEDMACHINE_MASK, CATS_PROJECTION, GR1KM_PROJECTION,
                                 TIDE_BANDS, tide_constants, tide_point, tide_height, era5_pressure,
                                 ifs_pressure, pair_pressure, ERA5, IFS,
                                 mask_window, gdal_network_setup, with_dataset, band_count, band_description
import ImagePairGeometry as IPG

gdal_network_setup()

@testset "the published tide models predict pyTMD's heights" begin
    for (model, projection, fixture) in ((CATS_MODEL, CATS_PROJECTION, "tide_pytmd_cats.json"),
                                         (GR1KM_MODEL, GR1KM_PROJECTION, "tide_pytmd_gr1km.json"))
        with_dataset(model) do d
            @test band_count(d) == length(TIDE_BANDS)
            @test all(band_description(d, b) == TIDE_BANDS[b] for b in eachindex(TIDE_BANDS))
        end
        for p in JSON3.read(read(joinpath(@__DIR__, "data", fixture), String))
            x, y = projection(Float64(p.lon), Float64(p.lat))
            k = tide_constants(model, [x], [y])
            z, f = tide_point(k, x, y)
            @test f ≈ p.flexure atol = 1.0e-6
            for (t, h) in zip(p.times, p.height)
                @test tide_height(z, DateTime(t[1:23])) ≈ h atol = 1.0e-6
            end
        end
    end
end

@testset "the published BedMachine mask has its classes where expected" begin
    # NPS 120 m lattice; one product pixel on 79N's floating tongue and one on Greenland's ice sheet.
    grid = IPG.MapGrid((-4108807.5, 120.0, 0.0, 4108807.5, 0.0, -120.0), (68480, 68480), 3413)
    at(x, y) = (i = round(Int, (x - grid.geotransform[1]) / 120 + 0.5);
                j = round(Int, (y - grid.geotransform[4]) / -120 + 0.5);
                first(mask_window(BEDMACHINE_MASK, grid, CartesianIndices((i:i, j:j)), 0))[1, 1])
    @test at(GR1KM_PROJECTION(-20.9, 79.53)...) == 0x03       # floating ice
    @test at(GR1KM_PROJECTION(-40.0, 72.5)...) == 0x02        # grounded ice
    @test at(GR1KM_PROJECTION(-10.0, 75.0)...) == 0x00        # Greenland Sea
end

@testset "ERA5 mean sea level pressure matches xarray's" begin
    # xarray, `msl.interp(time, longitude, latitude)` on the same store.
    for (lon, lat, t, p) in ((-20.9, 79.53, DateTime(2019, 2, 17, 8, 35, 8), 103115.15932590279),
                             (-160.0, -81.0, DateTime(2015, 11, 20, 8, 2, 15), 96400.59306640625))
        @test only(era5_pressure(t, [lon], [lat])) ≈ p rtol = 1.0e-6
    end
    # An hour ERA5 does not yet reach reads as missing, not as an error.
    @test era5_pressure(DateTime(2049, 1, 1), [0.0], [-80.0]) === nothing
end

# ERA5's preliminary record (ERA5T) runs about six days behind real time; a pair whose later scene is newer
# than its end is corrected for the tide alone. A record that stops advancing would leave every new pair
# without the inverse barometer, so it is caught here rather than in the products.
const ERA5_MAX_LAG_DAYS = 10

@testset "ERA5 reaches within $(ERA5_MAX_LAG_DAYS) days of now" begin
    @test era5_pressure(now(UTC) - Day(ERA5_MAX_LAG_DAYS), [0.0], [-80.0]) !== nothing
end

@testset "IFS pressure matches GDAL's bilinear read of the same runs" begin
    # `gdallocationinfo -wgs84 -r bilinear` on the runs' surface msl records, 2026-10-24 00Z and 06Z;
    # 01:30 is a quarter of the way from the first to the second.
    lons, lats = [-20.9, -60.9], [79.53, 80.82]
    p00, p06 = [101304.3505, 101108.1905], [100707.3115, 101191.9835]
    @test ifs_pressure(DateTime(2026, 10, 24), lons, lats) ≈ p00 rtol = 1.0e-7
    @test ifs_pressure(DateTime(2026, 10, 24, 1, 30), lons, lats) ≈ 0.75 .* p00 .+ 0.25 .* p06 rtol = 1.0e-7
    @test ifs_pressure(DateTime(2020, 1, 1), lons, lats) === nothing      # before the archive
    @test ifs_pressure(DateTime(2099, 1, 1), lons, lats) === nothing      # not yet published
end

# IFS is the inverse barometer of every pair newer than ERA5: a run more than a day old not yet published,
# or a pair older than ERA5's end not taking ERA5, is caught here.
@testset "a pair takes ERA5 where ERA5 reaches it and IFS up to the last day" begin
    lons, lats = [-20.9], [79.53]
    @test last(pair_pressure(DateTime(2025, 1, 1), DateTime(2025, 1, 13), lons, lats)) == ERA5
    t = floor(now(UTC), Hour(6)) - Day(1)
    @test last(pair_pressure(t - Day(12), t, lons, lats)) == IFS
end
