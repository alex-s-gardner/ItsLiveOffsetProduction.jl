# A product's catalog records against the ITS_LIVE metadata step's own. `data/catalog/` holds four golden
# products — Landsat 8, Sentinel-1, Sentinel-2 and NISAR GSLC — cut to the grid and `img_pair_info` the
# records are made from, and the STAC item, premet and spatial file `itslive-metadata`'s
# `cryoforge/generate.py` (commit 442ed01) writes of each, published under its product name (below) at
# `s3://its-live-data/velocity_image_pair/test/v02/N00E000/`. The files are named for their sensor: a
# product name runs to 196 characters, past what a Windows checkout can write.

using AutoriftItsLiveProduction: write_catalog_records, parse_pair_date, parse_created_date, iso_utc,
                                 _np_round, _py_round, ItsLiveOffsetProduction

const CATALOG_DATA = joinpath(@__DIR__, "data", "catalog")
const CATALOG_PREFIX = "s3://its-live-data/velocity_image_pair/test/v02/N00E000/"
const CATALOG_PRODUCTS = (
    landsat8 = "LC08_L1TP_009011_20200703_20200913_02_T1_X_LC08_L1TP_009011_20200820_20200905_02_T1_G0120V02_P078",
    sentinel1 = "S1A_IW_SLC__1SSH_20150828T162412_20150828T162431_007465_00A4AF_DC3E_X_" *
                "S1A_IW_SLC__1SSH_20150909T162413_20150909T162443_007640_00A97D_922B_G0120V02_P086",
    sentinel2 = "S2A_MSIL1C_20200626T204021_N0209_R014_T07VEG_20200626T235742_X_" *
                "S2B_MSIL1C_20200701T204019_N0209_R014_T07VEG_20200701T223901_G0120V02_P085",
    nisar_gslc = "NISAR_L2_PR_GSLC_003_170_D_053_7700_SHNA_A_20251028T235201_20251028T235238_X05009_N_P_J_001_X_" *
                 "NISAR_L2_PR_GSLC_007_170_D_053_7700_SHNA_A_20251215T235203_20251215T235240_X05009_N_P_J_001_G0120V02_P093")

@testset "a product's catalog records are the metadata step's" begin
    @test sort(filter(endswith(".nc"), readdir(CATALOG_DATA))) == sort([string(k) * ".nc" for k in keys(CATALOG_PRODUCTS)])
    for (stem, name) in pairs(CATALOG_PRODUCTS)
        stem = string(stem)
        mktempdir() do dir
            p = joinpath(dir, stem * ".nc")
            cp(joinpath(CATALOG_DATA, stem * ".nc"), p)
            item, pm, sp = write_catalog_records(p, CATALOG_PREFIX * name * ".nc")
            @test basename(item) == name * ".stac.json"
            # The reference's item, compacted as `json.dumps(…, separators = (",", ":"))` writes it.
            @test read(item, String) == read(joinpath(CATALOG_DATA, stem * ".stac.compact.json"), String)
            @test read(pm, String) == read(joinpath(CATALOG_DATA, stem * ".nc.premet"), String)
            @test read(sp, String) == read(joinpath(CATALOG_DATA, stem * ".nc.spatial"), String)
        end
    end
end

@testset "a file that is no product is refused, naming what it lacks" begin
    IL = ItsLiveOffsetProduction
    mktempdir() do dir
        p = joinpath(dir, "x.nc")
        IL.with_netcdf(p, true) do f
            IL.def_dim!(f, "x", 2)
            IL.put_var!(IL.def_var!(f, "x", Float64, ("x",)), [0.0, 120.0]; start = (0,), count = (2,))
            IL.def_var!(f, "y", Float64, ("x",))
            IL.put_att!(IL.def_var!(f, "mapping", Int8, String[]), "spatial_epsg", Int64(3413))
        end
        @test_throws "has no mapping:GeoTransform" IL.read_summary(p)
    end
end

@testset "dates and numbers are printed as the reference prints them" begin
    @test iso_utc(parse_pair_date("20200703T15:00:24.531814")) == "2020-07-03T15:00:24.531814Z"
    @test iso_utc(parse_pair_date("20200820T15:00:40.28767")) == "2020-08-20T15:00:40.287670Z"
    @test iso_utc(parse_pair_date("20200626T20:40:21.")) == "2020-06-26T20:40:21Z"
    @test iso_utc(parse_created_date("01-Sep-2026 23:51:18")) == "2026-09-01T23:51:18Z"
    @test_throws "is not a YYYYMMDDTHH:MM:SS.f date" parse_pair_date("2020-07-03")
    # numpy's round scales and rounds half to even; Python's rounds the exact binary value.
    @test _np_round(-45.00015, 4) == -45.0002
    @test _py_round(-45.00015, 4) == -45.0001
    @test _py_round(2.675, 2) == 2.67
    @test _py_round(0.125, 2) == 0.12
end
