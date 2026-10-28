# Longitudes across the antimeridian: a burst subset's polygon and footprints (`s1subset.jl`), and the
# geodetic windows of a radar pair (`radar.jl`).

using AutoriftItsLiveProduction: unwrap, subset_area, with_polygons, ring_overlap, wrap_lon, NearLon,
                                 geoid_sampler, mem_dataset, epsg_wkt

ring(lon0, lon1, lat0, lat1) = [(lon0, lat0), (lon1, lat0), (lon1, lat1), (lon0, lat1), (lon0, lat0)]

@testset "a footprint across ±180 is made continuous" begin
    r = ring(179.8, -179.9, -77.5, -77.0)
    @test unwrap(r, 180.0) == ring(179.8, 180.1, -77.5, -77.0)
    @test unwrap(r, -179.0) == ring(-180.2, -179.9, -77.5, -77.0)
    @test unwrap(ring(10.0, 11.0, 0.0, 1.0), 10.5) == ring(10.0, 11.0, 0.0, 1.0)
    # Moved by whole turns to sit next to the centre, whichever side it was given on.
    @test unwrap(ring(-170.0, -169.0, 0.0, 1.0), 175.0) == ring(190.0, 191.0, 0.0, 1.0)
end

@testset "a subset polygon split at ±180 is rejoined" begin
    split = "MULTIPOLYGON (((170 -70, 180 -70, 180 -60, 170 -60, 170 -70)), " *
            "((-180 -70, -170 -70, -170 -60, -180 -60, -180 -70)))"
    polys, center = subset_area(split)
    @test abs(center) ≈ 180
    @test extrema(x for p in polys for (x, _) in first(p)) == (center > 0 ? (170.0, 190.0) : (-190.0, -170.0))
    with_polygons(polys) do area
        # Across the join, wholly inside; half over the polygon's east edge, half covered.
        @test ring_overlap(area, unwrap(ring(179.0, -179.0, -65.0, -64.0), center)) ≈ 1
        @test ring_overlap(area, unwrap(ring(-175.0, -165.0, -65.0, -64.0), center)) ≈ 0.5
        @test ring_overlap(area, unwrap(ring(0.0, 1.0, -65.0, -64.0), center)) == 0
    end
    # A polygon away from the antimeridian keeps its longitudes.
    polys, center = subset_area("POLYGON ((-54.5 68.1, -47.5 69.2, -48.6 67.7, -54.5 68.1))")
    @test first(first(polys)) == [(-54.5, 68.1), (-47.5, 69.2), (-48.6, 67.7), (-54.5, 68.1)]
    @test -55 < center < -47
end

@testset "a subset polygon that cannot be read one way is refused" begin
    @test_throws "split one crossing the antimeridian" subset_area("POLYGON ((170 -70, -170 -70, -170 -60, 170 -60, 170 -70))")
    @test_throws "a subset spans at most 180°" subset_area(
        "MULTIPOLYGON (((0 0, 10 0, 10 1, 0 0)), ((170 0, 179 0, 179 1, 170 0)), ((-100 0, -90 0, -90 1, -100 0)))")
    @test_throws "is not a POLYGON or MULTIPOLYGON" subset_area("LINESTRING (0 0, 1 1)")
    @test_throws "encloses no area" subset_area("POLYGON ((0 0, 1 1, 0 0))")
    @test_throws "is not a WKT geometry" subset_area("POLYGON ((0 0")
end

@testset "longitudes move by whole turns" begin
    @test wrap_lon(187.5) == -172.5
    @test wrap_lon(-180.0) == -180.0
    @test wrap_lon(179.0) == 179.0
    @test NearLon(183.0)(-176.0, 1.0) == (184.0, 1.0)
    @test NearLon(0.0)(10.0, 1.0) == (10.0, 1.0)
end

@testset "a geoid window across ±180 joins the grid's two ends" begin
    # A global 1° grid of 100 sin(longitude), its cell centres on whole degrees from -180 as EGM2008's
    # are. The field is continuous across ±180 and odd about it, so a sample drawn from the wrong end of the
    # grid has the wrong sign.
    lons = -180.0:1.0:179.0
    lats = 90.0:-1.0:-90.0
    z = Float64[100 * sind(lon) for lon in lons, _ in lats]
    d = mem_dataset(z, (-180.5, 1.0, 0.0, 90.5, 0.0, -1.0), epsg_wkt(4326))
    try
        across = geoid_sampler(d, (178.2, 183.6), (-78.0, -76.0); margin_px = 2)
        for lon in (177.0, 179.0, 179.5, 180.0, -179.5, -179.0, -177.0, 180.5, 183.0)
            @test across(lon, -77.0) ≈ 100 * sind(lon) atol = 0.01
        end
        inside = geoid_sampler(d, (10.0, 12.0), (-1.0, 1.0); margin_px = 2)
        @test inside(11.0, 0.0) ≈ 100 * sind(11.0)
        @test inside(10.5, 0.5) ≈ (100 * sind(10.0) + 100 * sind(11.0)) / 2
    finally
        close(d)
    end
end
