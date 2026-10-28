# Scenes named by scene ID (`sceneid.jl`): IDs whose input is derived from the ID alone, and CMR records
# read out of a search response, here a real one cut to the fields read.

using AutoriftItsLiveProduction: scene_source, resolve_burst_ids, cmr_granules, burst_spec, granule_h5

const BURST_UR = "S1_342239_IW3_20211202T161351_VV_3876-BURST"
const NISAR_UR = "NISAR_L1_PR_RSLC_008_155_D_055_7700_SHNA_A_20251226T225459_20251226T225536_X05009_N_P_J_001"
const NISAR_DIR = "https://nisar.asf.earthdatacloud.nasa.gov/NISAR/NISAR_L1_RSLC_BETA_V1/" * NISAR_UR * "/" * NISAR_UR

attr(name, value) = """{"Name": "$name", "Values": ["$value"]}"""
const CMR_RESPONSE = """
{"hits": 2, "took": 30, "items": [
  {"meta": {"concept-id": "G2923495330-ASF"},
   "umm": {"GranuleUR": "$BURST_UR",
           "AdditionalAttributes": [$(attr("BURST_INDEX", 3)), $(attr("BYTE_OFFSET", 438958075)),
                                    $(attr("POLARIZATION", "VV")), $(attr("SUBSWATH_NAME", "IW3"))],
           "RelatedUrls": [{"URL": "https://sentinel1-burst.asf.alaska.edu/S1A_IW_SLC__1SDV_20211202T161342_20211202T161411_040832_04D8FA_3876/IW3/VV/3.tiff",
                            "Type": "USE SERVICE API", "Format": "GeoTIFF"}],
           "InputGranules": ["S1A_IW_SLC__1SDV_20211202T161342_20211202T161411_040832_04D8FA_3876-SLC"]}},
  {"meta": {"concept-id": "G0-ASF"},
   "umm": {"GranuleUR": "$NISAR_UR",
           "RelatedUrls": [{"URL": "$NISAR_DIR.h5", "Type": "GET DATA"},
                           {"URL": "s3://sds-n-cumulus-prod-nisar-products/$NISAR_UR.h5", "Type": "GET DATA VIA DIRECT ACCESS"},
                           {"URL": "$NISAR_DIR.h5.iso.xml", "Type": "EXTENDED METADATA"}]}}]}
"""

@testset "a scene's ID names its input where the layout fixes it" begin
    slc = "S1B_IW_SLC__1SDH_20180630T095955_20180630T100022_011601_015537_A96A"
    url = "https://datapool.asf.alaska.edu/SLC/SB/" * slc * ".zip"
    for id in (slc, slc * "-SLC")
        s = scene_source(id, "S1-SLC")
        @test (s.input, s.name, s.band) == (id, slc, url)
    end
    # A burst pair's SLC IDs become their zips, which `--subset` selects bursts out of; other inputs pass.
    @test resolve_burst_ids(slc * ", " * slc, "/data/x.SAFE") == (url * "," * url, "/data/x.SAFE")
    @test resolve_burst_ids("https://a/b.zip,https://a/c.zip", "x.SAFE") == ("https://a/b.zip,https://a/c.zip", "x.SAFE")
    @test_throws "mixes burst IDs, SLC IDs, paths or URLs" resolve_burst_ids(BURST_UR * "," * url, url)
    @test_throws "mixes burst IDs, SLC IDs, paths or URLs" resolve_burst_ids(url, BURST_UR * "," * slc)
    # An input that is no ID is read as before, and refused as before.
    @test_throws "a STAC item or a scene ID" scene_source("LC08_L1TP_021248", "L8")
    @test_throws "nor its ID" scene_source("S1A_IW_SLC__1SDV_20211118T161558", "S1-SLC")
    # An ID of another sensor's scene is no input of this one.
    @test_throws "cannot tell which Sentinel-2 product" scene_source("LC08_L1TP_021248_20210507_20210517_02_T1", "S2")
end

@testset "a CMR record names a burst's SLC zip and a granule's file" begin
    burst, nisar = cmr_granules(CMR_RESPONSE, [BURST_UR, NISAR_UR])
    @test burst_spec(burst) == "https://datapool.asf.alaska.edu/SLC/SA/" *
                               "S1A_IW_SLC__1SDV_20211202T161342_20211202T161411_040832_04D8FA_3876.zip#IW3/VV/3"
    @test granule_h5(nisar) == NISAR_DIR * ".h5"
    bare = """{"items": [{"umm": {"GranuleUR": "$BURST_UR", "InputGranules": ["$(only(burst.inputs))"]}}]}"""
    @test_throws "CMR's record of $BURST_UR has no single SUBSWATH_NAME" burst_spec(only(cmr_granules(bare, [BURST_UR])))
    @test_throws "lists 0 https:// .h5 links, not one" granule_h5(burst)
    @test_throws "lists no granule S1_000000_IW1_20211202T161351_VV_3876-BURST" cmr_granules(
        CMR_RESPONSE, [BURST_UR, "S1_000000_IW1_20211202T161351_VV_3876-BURST"])
end
