# Where a Sentinel-2 product is in `s3://sentinel-s2-l1c`, from the bucket's `productInfo.json` and
# `tileInfo.json`; here real ones cut to the fields read.

using AutoriftItsLiveProduction: s2_product_tile, s2_original_product

tileinfo(path, product, datastrip) = """
{"path": "$path", "productName": "$product", "datastrip": {"id": "$datastrip", "path": "products/x/datastrip/0"},
 "tileGeometry": {"type": "Polygon"}}"""

@testset "a product's tile directory is the one its productInfo.json lists" begin
    name = "S2A_MSIL1C_20200624T181931_N0209_R127_T19XDJ_20200624T201221"
    info = """
    {"name": "$name", "path": "products/2020/6/24/$name",
     "tiles": [{"path": "tiles/19/X/DJ/2020/6/24/0", "utmZone": 19, "latitudeBand": "X", "gridSquare": "DJ",
                "datastrip": {"id": "S2A_OPER_MSI_L1C_DS_SGS__20200624T201221_S20200624T181926_N02.09",
                              "path": "products/2020/6/24/$name/datastrip/0"}}]}"""
    @test s2_product_tile(name, info) ==
          (tile = "tiles/19/X/DJ/2020/6/24/0", datastrip = "products/2020/6/24/$name/datastrip/0")
    # Month and day are unpadded, and so is a zone below 10.
    veg = "S2A_MSIL1C_20200626T204021_N0209_R014_T07VEG_20200626T235742"
    @test s2_product_tile(veg, """{"tiles": [{"path": "tiles/7/V/EG/2020/6/26/1", "datastrip": {"path": "d"}}]}""").tile ==
          "tiles/7/V/EG/2020/6/26/1"
    @test_throws "expected one T19XDJ tile in $name's productInfo.json, found 0" s2_product_tile(name,
        """{"tiles": [{"path": "tiles/19/X/DK/2020/6/24/0", "datastrip": {"path": "d"}}]}""")
end

@testset "a product from before December 2016 is found by its original name" begin
    # Filed as S2A_OPER_PRD_MSIL1C_PDMC_<production>_R<orbit>_V<start>_<stop>; the compact name's last field
    # is the production time, or, for products ESA reprocessed, the datastrip's sensing start.
    wea = "S2A_MSIL1C_20160330T150931_N0201_R025_T22WEA_20160401T035851"
    @test s2_original_product(tileinfo("tiles/22/W/EA/2016/3/30/0",
                                       "S2A_OPER_PRD_MSIL1C_PDMC_20160401T035851_R025_V20160330T150931_20160330T150931",
                                       "S2A_OPER_MSI_L1C_DS_SGS__20160330T205642_S20160330T150931_N02.01"), wea)
    wwp = "S2A_MSIL1C_20160403T131312_N0201_R081_T27WWP_20160403T131310"
    @test s2_original_product(tileinfo("tiles/27/W/WP/2016/4/3/0",
                                       "S2A_OPER_PRD_MSIL1C_PDMC_20160404T093534_R081_V20160403T131310_20160403T131310",
                                       "S2A_OPER_MSI_L1C_DS_MTI__20160403T195930_S20160403T131310_N02.01"), wwp)
    # Another orbit's product on the same tile and day, another baseline, and a product filed under its ID.
    xmj = "S2A_MSIL1C_20160408T171929_N0201_R012_T23XMJ_20160409T015112"
    @test !s2_original_product(tileinfo("tiles/23/X/MJ/2016/4/8/0",
                                        "S2A_OPER_PRD_MSIL1C_PDMC_20160409T004000_R014_V20160408T203925_20160408T203925",
                                        "S2A_OPER_MSI_L1C_DS_SGS__20160408T223829_S20160408T203925_N02.01"), xmj)
    @test !s2_original_product(tileinfo("tiles/22/W/EA/2016/3/30/0",
                                        "S2A_OPER_PRD_MSIL1C_PDMC_20160401T035851_R025_V20160330T150931_20160330T150931",
                                        "S2A_OPER_MSI_L1C_DS_SGS__20160330T205642_S20160330T150931_N02.04"), wea)
    @test !s2_original_product(tileinfo("tiles/22/W/EA/2016/3/30/0", wea,
                                        "S2A_OPER_MSI_L1C_DS_SGS__20160330T205642_S20160330T150931_N02.01"), wea)
end
