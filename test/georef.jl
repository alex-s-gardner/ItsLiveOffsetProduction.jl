# CF grid-mapping attributes, and the `ItsLiveGeoref` a `Rasters.Raster` implies.

using ItsLiveOffsetProduction
using Rasters
using Rasters: EPSG, Projected, X, Y
using DimensionalData
using DimensionalData.Lookups: ForwardOrdered, Intervals, Regular, ReverseOrdered, Start
using DimensionalData: lookup
import GeoFormatTypes as GFT
using Test

# A north-up projected grid. Only the dims and the CRS matter here, so the values are zeros.
function _grid_raster(n; res = 120.0, epsg = 3413)
    x = X(Projected(0.0:res:(res * (n - 1)); order = ForwardOrdered(), span = Regular(res),
                    sampling = Intervals(Start()), crs = EPSG(epsg)))
    y = Y(Projected((res * (n - 1)):(-res):0.0; order = ReverseOrdered(), span = Regular(-res),
                    sampling = Intervals(Start()), crs = EPSG(epsg)))
    return Raster(zeros(Float32, n, n), (y, x))
end

@testset "cf_grid_mapping matches the reference's own CF attributes" begin
    # Every value here is the real `mapping` variable attributes of a captured golden ITS_LIVE
    # product (a polar-stereographic, EPSG:3413 case) — confirmed key-for-key and value-for-value
    # against `ncdump -h` on that file, not asserted from documentation. `OSRExportToCF1`, GDAL's own
    # generic CF exporter, was tried first and rejected: it reflects a *current* GDAL's naming
    # (`standard_parallel`, a `long_name`, no `scale_factor_at_projection_origin`), not the older one
    # the reference container's GDAL uses — see `cf_grid_mapping`'s docstring.
    d = cf_grid_mapping(3413)
    @test d["grid_mapping_name"] == "polar_stereographic"
    @test d["latitude_of_origin"] == 70.0
    @test d["latitude_of_projection_origin"] == 90.0
    @test d["straight_vertical_longitude_from_pole"] == -45.0
    @test d["false_easting"] == 0.0
    @test d["false_northing"] == 0.0
    @test d["semi_major_axis"] == 6378137.0
    @test d["inverse_flattening"] == 298.257223563
    @test d["scale_factor_at_projection_origin"] === Int64(1)
    @test d["spatial_epsg"] === Int64(3413)
    @test d["proj4text"] ==
          "+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
    @test startswith(d["spatial_ref"], "PROJCS[\"WGS 84 / NSIDC Sea Ice Polar Stereographic North\"")
    @test d["crs_wkt"] == d["spatial_ref"]
    # `write_product` always recomputes this from the grid actually written, never from the caller's dict.
    @test !haskey(d, "GeoTransform")

    # The southern grid: same variant, opposite-signed pole.
    ds = cf_grid_mapping(3031)
    @test ds["grid_mapping_name"] == "polar_stereographic"
    @test ds["latitude_of_projection_origin"] == -90.0

    # A CRS with no authority code at all — a custom or stripped WKT, not a hypothetical case.
    # `ArchGDAL.toEPSG` *throws* rather than returning a sentinel for this (confirmed by calling it),
    # so `spatial_epsg` being merely absent, instead of the whole function erroring, is the behavior
    # under test here.
    nocode_wkt = """PROJCS["custom polar",GEOGCS["WGS 84",DATUM["WGS_1984",\
    SPHEROID["WGS 84",6378137,298.257223563]],PRIMEM["Greenwich",0],\
    UNIT["degree",0.0174532925199433]],PROJECTION["Polar_Stereographic"],\
    PARAMETER["latitude_of_origin",70],PARAMETER["central_meridian",-45],\
    PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1]]"""
    GFT = Rasters.GeoFormatTypes
    nocode = cf_grid_mapping(GFT.WellKnownText(GFT.CRS(), nocode_wkt))
    @test !haskey(nocode, "spatial_epsg")
    @test nocode["grid_mapping_name"] == "polar_stereographic"

    # An EPSG code, an `Integer`, or a `GeoFormatTypes.GeoFormat` all work.
    @test cf_grid_mapping(EPSG(3413)) == d
    @test_throws ArgumentError cf_grid_mapping("not a crs")

    # A projection outside ITS_LIVE's two grid types is refused rather than silently mis-described.
    # `transverse_mercator` itself is exercised (key set only — not confirmed against a captured UTM
    # case, unlike `polar_stereographic` above; see the docstring).
    @test_throws "Polar_Stereographic and Transverse_Mercator" cf_grid_mapping(4326)
    utm = cf_grid_mapping(32624)
    @test utm["grid_mapping_name"] == "transverse_mercator"
    @test haskey(utm, "scale_factor_at_central_meridian")
    @test haskey(utm, "longitude_of_central_meridian")
    @test haskey(utm, "latitude_of_projection_origin")
end

@testset "ItsLiveGeoref built from a Raster" begin
    ref = _grid_raster(20; res = 120.0, epsg = 3413)
    g = ItsLiveGeoref(ref)
    @test g isa ItsLiveGeoref
    @test g.x == collect(lookup(ref, X))
    @test g.y == collect(lookup(ref, Y))
    @test g.pixel_size_x == 120.0
    @test g.pixel_size_y == 120.0
    @test g.mapping_attrs["grid_mapping_name"] == "polar_stereographic"

    lon, lat = g.lonlat(g.x[1], g.y[1])
    @test -180 <= lon <= 180
    @test -90 <= lat <= 90

    # A radar scene's own ground pixel size overrides the grid's — see `ItsLiveGeoref`'s docstring on
    # why `write_product` needs that rather than the projected grid's own spacing.
    g2 = ItsLiveGeoref(ref; pixel_size_x = 2.3, pixel_size_y = 13.9)
    @test (g2.pixel_size_x, g2.pixel_size_y) == (2.3, 13.9)
end

