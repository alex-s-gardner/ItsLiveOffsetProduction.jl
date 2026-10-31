# The command line a production run is given: the sensor read from the inputs when not named, the
# output path optional, the upload options, and Earthdata credentials from the environment.

using AutoriftItsLiveProduction: parse_args, infer_sensor, with_earthdata_netrc

const L8A = "LC08_L1TP_021248_20210507_20210517_02_T1"
const L8B = "LC08_L1TP_021248_20220611_20220617_02_T1"

@testset "the sensor and output path are optional" begin
    a = parse_args([L8A, L8B])
    @test (a.sensor, a.reference, a.secondary, a.output, a.bucket, a.bucket_prefix, a.publish_bucket) ==
          ("", L8A, L8B, "", "", "", "")
    a = parse_args(["l8", L8A, L8B, "out.nc", "--bucket", "b", "--bucket-prefix", "p/q", "--publish-bucket", "its-live-data"])
    @test (a.sensor, a.output, a.bucket, a.bucket_prefix, a.publish_bucket) == ("L8", "out.nc", "b", "p/q", "its-live-data")
    # HyP3 passes a job parameter it has no value for as `None` (#47); an empty value means the same.
    @test parse_args([L8A, L8B, "--publish-bucket", "None"]).publish_bucket == ""
    @test parse_args([L8A, L8B, "out.nc"]).output == "out.nc"
    @test_throws "expected [<sensor>] <reference> <secondary> [<output.nc>], got 1" parse_args([L8A])
    @test_throws "got 5 positional arguments" parse_args(["L8", L8A, L8B, "a.nc", "b.nc"])
    @test_throws "--bucket-prefix needs --bucket" parse_args([L8A, L8B, "--bucket-prefix", "p"])
    @test_throws "give one of them" parse_args([L8A, L8B, "--no-tide-correction", "--tide-model", "/vsimem/m.tif"])
    a = parse_args([L8A, L8B, "--tide-model", "/vsimem/m.tif"])
    @test (a.options.tide, a.options.tide_model) == (true, "/vsimem/m.tif")
    a = parse_args([L8A, L8B, "--no-tide-correction"])
    @test (a.options.tide, a.options.tide_model) == (false, "")
end

@testset "the sensor is read from the reference" begin
    zip = "https://datapool.asf.alaska.edu/SLC/SB/S1B_IW_SLC__1SDH_20180630T095955_20180630T100022_011601_015537_A96A.zip"
    for (input, sensor) in (
            (L8A, "L8"), ("LT05_L1TP_146037_19940407_20200913_02_T1", "L5"),
            ("https://landsatlook.usgs.gov/stac-server/collections/landsat-c2l1/items/LE07_L1GT_060115_20080309_20200913_02_T2", "L7"),
            ("s3://usgs-landsat/collection02/LC09_L1TP_021248_20220611_20220617_02_T1_B8.TIF", "L9"),
            ("https://storage.googleapis.com/gcp-public-data-sentinel-2/tiles/19/X/DJ/S2A_MSIL1C_20200624T181931_N0209_R127_T19XDJ_20200624T201221.SAFE/GRANULE/L1C_T19XDJ_A026150_20200624T181926/IMG_DATA/T19XDJ_20200624T181931_B08.jp2", "S2"),
            ("S2A_MSIL1C_20200624T181931_N0209_R127_T19XDJ_20200624T201221", "S2"),
            (zip, "S1"), ("/data/S1A_IW_SLC__1SDV_20211118T161558_20211118T161625_040628_04D1E5_D6DE.SAFE/", "S1"),
            (zip * "," * zip, "S1BURST"), (zip * "#IW2/HH/3", "S1BURST"),
            ("S1_342239_IW3_20211202T161351_VV_3876-BURST,S1_342240_IW3_20211202T161354_VV_3876-BURST", "S1BURST"),
            ("https://sentinel1-burst.asf.alaska.edu/S1A_IW_SLC__1SDV_20211202T161342_20211202T161411_040832_04D8FA_3876/IW3/VV/3.tiff", "S1BURST"),
            ("NISAR_L1_PR_RSLC_008_155_D_055_7700_SHNA_A_20251226T225459_20251226T225536_X05009_N_P_J_001", "NL1"),
            ("https://x/NISAR_L2_PR_GSLC_003_170_D_053_7700_SHNA_A_20251028T235201_20251028T235238_X05009_N_P_J_001.h5", "NL2"))
        @test infer_sensor(input) == sensor
    end
    @test_throws "cannot tell the sensor from the reference \"item.json\"" infer_sensor("item.json")
end

@testset "a Sentinel-1 SLC with --subset is a burst pair (#38)" begin
    poly = "POLYGON ((-55.13492 66.810363, -54.473561 68.324585, -48.482901 67.948708, -49.497015 66.462691, -55.13492 66.810363))"
    r1, r2 = "S1B_IW_SLC__1SDH_20180630T095955_20180630T100022_011601_015537_A96A",
             "S1B_IW_SLC__1SDH_20180630T095929_20180630T095957_011601_015537_164D"
    zip = "https://datapool.asf.alaska.edu/SLC/SB/" * r1 * ".zip"
    # SLC scene IDs as -r lists them, one ID, and a zip: each is S1 alone and S1BURST with a polygon.
    for input in (r1 * "," * r2, r1, zip)
        @test infer_sensor(input; subset = poly) == "S1BURST"
    end
    @test infer_sensor(r1) == "S1"
    # The command of #38, as parsed: no sensor given, so the polygon decides it.
    a = parse_args(["out.nc", "-r", r1, r2, "-s", "S1B_IW_SLC__1SDH_20180712T095955_20180712T100022_011776_015AA2_1C61",
                    "S1B_IW_SLC__1SDH_20180712T095930_20180712T095958_011776_015AA2_0514", "--subset", poly])
    @test isempty(a.sensor)
    @test infer_sensor(a.reference; a.subset) == "S1BURST"
    # A polygon does not change any other sensor; Job refuses it for them, as for an explicit S1.
    @test infer_sensor(L8A; subset = poly) == "L8"
end

@testset "an option given None, empty, or an empty geometry is unset (#47)" begin
    plain = parse_args([L8A, L8B])
    for unset in ("None", "", " ")
        a = parse_args([L8A, L8B, "--dem", unset, "--orbit-dir", unset, "--subset", unset, "--bucket", unset,
                        "--bucket-prefix", unset, "--publish-bucket", unset, "--tide-model", unset])
        @test (a.subset, a.bucket, a.bucket_prefix, a.publish_bucket) == ("", "", "", "")
        @test (a.options.dem, a.options.orbit_dir, a.options.tide, a.options.tide_model) ==
              (plain.options.dem, plain.options.orbit_dir, true, "")
    end
    # So does every other option that takes a value, keeping its default.
    a = parse_args([L8A, L8B, "--scratch", "None", "--cache-io", "None", "--filter-tile", "None", "--heap-hint-gib", "None"])
    @test (a.options.scratch, a.options.cache_io, a.options.filter_tile, a.heap_hint_gib) ==
          (plain.options.scratch, plain.options.cache_io, plain.options.filter_tile, plain.heap_hint_gib)
    for empty in ("POLYGON EMPTY", "MULTIPOLYGON EMPTY", "polygon empty")
        @test parse_args([L8A, L8B, "--subset", empty]).subset == ""
    end
    # A Sentinel-1 SLC pair with an unset --subset stays an SLC pair, as HyP3 runs one.
    r1 = "S1B_IW_SLC__1SDH_20180630T095955_20180630T100022_011601_015537_A96A"
    a = parse_args(["-r", r1, "-s", "S1B_IW_SLC__1SDH_20180712T095955_20180712T100022_011776_015AA2_1C61",
                    "--subset", "None", "--publish-bucket", "None"])
    @test infer_sensor(a.reference; a.subset) == "S1"
    # An unset --tide-model with --no-tide-correction is no conflict.
    @test parse_args([L8A, L8B, "--tide-model", "None", "--no-tide-correction"]).options.tide == false
    # Anything else is still a value.
    @test parse_args([L8A, L8B, "--bucket", "my-bucket"]).bucket == "my-bucket"
end

@testset "Earthdata credentials from the environment" begin
    withenv("EARTHDATA_USERNAME" => nothing, "EARTHDATA_PASSWORD" => nothing, "NETRC" => nothing) do
        @test with_earthdata_netrc(() -> get(ENV, "NETRC", "")) == ""
    end
    withenv("EARTHDATA_USERNAME" => "someone", "EARTHDATA_PASSWORD" => "s3cret", "NETRC" => nothing) do
        path, text, mode = with_earthdata_netrc(() -> (ENV["NETRC"], read(ENV["NETRC"], String), filemode(dirname(ENV["NETRC"]))))
        @test text == "machine urs.earthdata.nasa.gov login someone password s3cret\n"
        @test mode & 0o077 == 0
        # Removed when the run ends, and NETRC with it.
        @test !ispath(dirname(path)) && !haskey(ENV, "NETRC")
        @test_throws "boom" with_earthdata_netrc(() -> error("boom"))
        @test !haskey(ENV, "NETRC")
    end
    withenv("EARTHDATA_USERNAME" => "someone", "EARTHDATA_PASSWORD" => nothing, "NETRC" => nothing) do
        @test_throws "only one of them is" with_earthdata_netrc(() -> nothing)
    end
    withenv("EARTHDATA_USERNAME" => "someone", "EARTHDATA_PASSWORD" => "s3cret", "NETRC" => "/x/.netrc") do
        @test_throws "set one or the other" with_earthdata_netrc(() -> nothing)
    end
end

@testset "the scenes may be given as --reference and --secondary lists" begin
    r1, r2 = "S1_342239_IW3_20211202T161351_VV_3876-BURST", "S1_342240_IW3_20211202T161354_VV_3876-BURST"
    s1, s2 = "S1_342239_IW3_20211214T161351_VV_91C4-BURST", "S1_342240_IW3_20211214T161353_VV_91C4-BURST"
    a = parse_args(["--reference", r1, r2, "--secondary", s1, s2, "--bucket", "b"])
    @test (a.sensor, a.reference, a.secondary, a.output, a.bucket) == ("", r1 * "," * r2, s1 * "," * s2, "", "b")
    # A list runs to the next option, as argparse reads `nargs = '+'`: the output comes before it.
    a = parse_args(["S1BURST", "pair.nc", "-r", r1 * " " * r2, "-s", s1 * "," * s2, "--quiet"])
    @test (a.sensor, a.reference, a.secondary, a.output) == ("S1BURST", r1 * "," * r2, s1 * "," * s2, "pair.nc")
    @test !a.options.verbose
    # HyP3 renders a list it was not given as `None`.
    @test parse_args(["-r", L8A, "None", "-s", L8B]).reference == L8A
    @test_throws "--reference needs at least one scene" parse_args(["--reference", "--secondary", L8B])
    @test_throws "only one of them is" parse_args(["--reference", L8A])
    @test_throws "with --reference and --secondary, expected [<sensor>] [<output.nc>], got 2" parse_args(
        ["x.nc", "y.nc", "-r", L8A, "-s", L8B])
end
