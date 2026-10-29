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
    # HyP3 passes a job parameter it has no value for as `None`.
    @test parse_args([L8A, L8B, "--publish-bucket", "None"]).publish_bucket == ""
    @test parse_args([L8A, L8B, "out.nc"]).output == "out.nc"
    @test_throws "expected [<sensor>] <reference> <secondary> [<output.nc>], got 1" parse_args([L8A])
    @test_throws "got 5 positional arguments" parse_args(["L8", L8A, L8B, "a.nc", "b.nc"])
    @test_throws "--bucket-prefix needs --bucket" parse_args([L8A, L8B, "--bucket-prefix", "p"])
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
