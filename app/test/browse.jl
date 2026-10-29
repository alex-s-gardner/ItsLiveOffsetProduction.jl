# The browse image and thumbnail, against the reference's own. `data/browse_v.npy` is a 230 × 157 grid of
# speeds holding every value from 0 to 700 m/yr once and the fill a fifth of the time; `browse_ref.png`
# and `browse_ref_thumb.png` are what `hyp3_autorift.image.make_browse` and
# `hyp3lib.image.create_thumbnail` (hyp3-autorift 0.28.4: matplotlib 3.10.3, SciPy 1.16.0, Pillow 11.3.0)
# draw of it.

using AutoriftItsLiveProduction: browse_rgba, thumbnail_rgba, thumbnail_size, write_png, point_to_region,
                                 opendata_platform, PublishInfo
# The app's own dependencies, which the test environment reaches through it.
using AutoriftItsLiveProduction: GDAL, JSON3

const BROWSE_DATA = joinpath(@__DIR__, "data")

# A `.npy` of little-endian `int16` in C order, as an array indexed `[x, y]`.
function read_npy_i16(path, nx, ny)
    raw = read(path)
    header = Int(reinterpret(UInt16, raw[9:10])[1])
    occursin("'<i2'", String(raw[11:(10 + header)])) || error(path * " is not little-endian int16")
    return collect(reshape(reinterpret(Int16, raw[(11 + header):end]), nx, ny))
end

function read_png(path)
    h = GDAL.gdalopen(path, GDAL.GA_ReadOnly)
    nx, ny = Int(GDAL.gdalgetrasterxsize(h)), Int(GDAL.gdalgetrasterysize(h))
    out = Array{UInt8,3}(undef, 4, nx, ny)
    GDAL.gdaldatasetrasterio(h, GDAL.GF_Read, 0, 0, nx, ny, out, nx, ny, GDAL.GDT_Byte, 4, C_NULL, 4, 4nx, 1)
    GDAL.gdalclose(h)
    return out
end

@testset "the browse image and thumbnail are the reference's, pixel for pixel" begin
    v = read_npy_i16(joinpath(BROWSE_DATA, "browse_v.npy"), 230, 157)
    browse = browse_rgba(v)
    @test browse == read_png(joinpath(BROWSE_DATA, "browse_ref.png"))
    thumb = thumbnail_rgba(browse)
    @test thumb == read_png(joinpath(BROWSE_DATA, "browse_ref_thumb.png"))
    mktempdir() do dir
        @test read_png(write_png(joinpath(dir, "b.png"), browse)) == browse
    end
end

@testset "a thumbnail keeps the aspect ratio Pillow keeps" begin
    @test thumbnail_size(2048, 1536) == (100, 75)
    @test thumbnail_size(230, 157) == (100, 68)
    @test thumbnail_size(157, 230) == (68, 100)
    @test thumbnail_size(3000, 10) == (100, 1)
    @test thumbnail_size(80, 60) == (80, 60)
end

@testset "a product is published under its platform and region" begin
    @test point_to_region(81.79, -45.0) == "N80W040"
    @test point_to_region(78.6, -124.2) == "N70W120"
    @test point_to_region(-77.2, 180.0) == "S70E170"
    @test point_to_region(90.0, 0.0) == "N80E000"
    @test point_to_region(-0.5, -0.5) == "S00W000"
    @test [opendata_platform(p) for p in ("L5", "L9", "S2", "S1-SLC", "S1-BURST", "NISAR-L1", "NISAR-L2")] ==
          ["landsatOLI", "landsatOLI", "sentinel2", "sentinel1", "sentinel1", "nisar", "nisar"]
    @test JSON3.write(PublishInfo("its-live-data", "velocity_image_pair/landsatOLI/v02/N80W040", "a_X_b_G0120V02_P020.nc")) ==
          "{\"bucket\":\"its-live-data\",\"prefix\":\"velocity_image_pair/landsatOLI/v02/N80W040\"," *
          "\"name\":\"a_X_b_G0120V02_P020.nc\"}"
end
