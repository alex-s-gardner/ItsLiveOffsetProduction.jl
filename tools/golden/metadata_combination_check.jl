# Does `SLCDatasets.Identification`/`OpticalDatasets.Identification` actually combine into an
# `ImagePairInfo` the way a driver is expected to? Neither package depends on the other or
# on ItsLiveOffsetProduction, and `AutoRIFT.jl` ships no adapter for either (see `ImagePairInfo`'s docstring) — the
# combination is caller code, demonstrated and checked here against real captured `img_pair_info`
# values rather than assumed to fit.
#
#   julia --project=tools/golden tools/golden/metadata_combination_check.jl

import Pkg
Pkg.activate(@__DIR__)

using AutoRIFT
using ItsLiveOffsetProduction
using OpticalDatasets
using SLCDatasets
using Dates
using Test

# --- Optical: OpticalDatasets carries everything img_pair_info needs directly. ---

@testset "Landsat identification combines into ItsLiveOffsetProduction.ImagePairInfo" begin
    # Real captured `img_pair_info` for `LC08_L1TP_009011_20200703_..._X_..._20200820_...` — see
    # `tools/golden/capture_packaging.py`. `datetime` is the STAC `properties.datetime` for the same
    # scene (`landsatlook.usgs.gov`'s STAC item), not invented — see `OpticalDatasets`' own tests for
    # where it came from.
    id1 = landsat_identification("LC08_L1TP_009011_20200703_20200913_02_T1",
                                 DateTime("2020-07-03T15:00:24.531"))
    id2 = landsat_identification("LC08_L1TP_009011_20200820_20200905_02_T1",
                                 DateTime("2020-08-20T15:00:40.288"))

    extra = Dict{String,Any}(
        "id_img1" => id1.id, "id_img2" => id2.id,
        "sensor_img1" => id1.sensor, "sensor_img2" => id2.sensor,
        "correction_level_img1" => id1.correction_level,
        "correction_level_img2" => id2.correction_level,
        "path_img1" => id1.path, "path_img2" => id2.path,
        "row_img1" => id1.row, "row_img2" => id2.row,
        "collection_number_img1" => id1.collection_number,
        "collection_number_img2" => id2.collection_number,
        "collection_category_img1" => id1.collection_category,
        "collection_category_img2" => id2.collection_category,
        "processing_date_img1" => id1.processing_date, "processing_date_img2" => id2.processing_date)
    info = ItsLiveOffsetProduction.ImagePairInfo(id1.acquisition_time, id2.acquisition_time, id1.mission, id2.mission,
                                  id1.satellite, id2.satellite, 78.4, 69.57, -49.22, extra)

    # Every value below is the real captured `img_pair_info` for this pair.
    @test info.mission_img1 == "L" && info.mission_img2 == "L"
    @test info.satellite_img1 == "8" && info.satellite_img2 == "8"
    @test info.extra["sensor_img1"] == "C"
    @test info.extra["path_img1"] == 9 && info.extra["row_img1"] == 11
    @test info.extra["collection_number_img1"] == 2
    @test info.extra["collection_category_img1"] == "T1"
    @test info.acquisition_date_img1 == DateTime("2020-07-03T15:00:24.531")
end

@testset "Sentinel-2 identification combines into ItsLiveOffsetProduction.ImagePairInfo" begin
    # Real captured `img_pair_info` for the golden `S2B_..._X_S2A_...` pair.
    id1 = sentinel2_identification("S2B_MSIL1C_20200612T150759_N0209_R025_T22WEB_20200612T184700")
    id2 = sentinel2_identification("S2A_MSIL1C_20200627T150921_N0209_R025_T22WEB_20200627T170912")

    extra = Dict{String,Any}("id_img1" => id1.id, "id_img2" => id2.id,
                             "correction_level_img1" => id1.correction_level,
                             "correction_level_img2" => id2.correction_level,
                             "sensor_img1" => id1.sensor, "sensor_img2" => id2.sensor)
    info = ItsLiveOffsetProduction.ImagePairInfo(id1.acquisition_time, id2.acquisition_time, id1.mission, id2.mission,
                                  id1.satellite, id2.satellite, 78.5, 68.9, -49.6, extra)

    @test info.mission_img1 == "S" && info.satellite_img1 == "2B"
    @test info.mission_img2 == "S" && info.satellite_img2 == "2A"
    @test info.extra["sensor_img1"] == "MSI"
    @test info.extra["correction_level_img1"] == "L1C"
end

# --- Radar: SLCDatasets.Identification alone is not enough. ---

@testset "SLCDatasets.Identification splits into mission_img1/satellite_img1, but not the rest" begin
    # Real captured `img_pair_info` for a Sentinel-1 golden pair (`S1A_IW_SLC..._007465..._X_..._007640...`).
    # `SLCDatasets.Identification.mission` is the platform id ("S1A"), which `ImagePairInfo` splits in
    # two — this is the one piece of the combination that is genuinely non-trivial (not a field
    # rename) and worth pinning down rather than assuming.
    id = SLCDatasets.Identification("S1A", "SLC", 7465, "ascending", "right",
                                    "2015-08-28T16:24:14.425753", "2015-08-28T16:24:31.000000", "")
    mission, satellite = string(id.mission[1]), id.mission[2:end]
    @test mission == "S"           # real captured `mission_img1`
    @test satellite == "1A"        # real captured `satellite_img1`
    @test lpad(string(id.absolute_orbit), 6, '0') == "007465"   # real captured `absolute_orbit_number_img1`

    # `sensor_img1` ("C", for C-band — coincidentally the same letter Landsat's OLI/TIRS uses, for an
    # unrelated reason) is not one of `Identification`'s 8 fields at all; it is a fixed constant for
    # every Sentinel-1 product, which a driver would hardcode rather than read from here.
    #
    # `mission_data_take_ID_img1` ("00A4AF") and `product_unique_ID_img1` ("DC3E") are not
    # `Identification` fields either — they are the trailing hex segments of the granule id string
    # itself (`S1A_IW_SLC__1SSH_..._007465_00A4AF_DC3E`), which `SLCDatasets` does not currently parse
    # out. A driver wanting them has to split the id string by hand; `SLCDatasets` exposes the parsed
    # burst/orbit geometry that string implies, not the string's own trailing fields. Documented here
    # as a real, disclosed gap — not something this check papers over by fabricating a parser.
end
