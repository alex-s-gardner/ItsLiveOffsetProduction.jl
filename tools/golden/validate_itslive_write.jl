# Validate `write_product` against a real captured `netCDF_packaging` call, from the reference
# container itself. See `capture_packaging.py`/`.jl`.
#
#   julia --project=tools/golden tools/golden/validate_itslive_write.jl <case-fragment>
#
# `run_case` is also `include`d by `run_all_golden.jl` to sweep every golden case.

import Pkg
Pkg.activate(@__DIR__)

using AutoRIFT, ItsLiveOffsetProduction, NCDatasets, Dates, JSON3
import FastGeoProjections as FGP

include(joinpath(@__DIR__, "manifest.jl"))
include(joinpath(@__DIR__, "reference.jl"))
include(joinpath(@__DIR__, "capture_packaging.jl"))
include(joinpath(@__DIR__, "product.jl"))
include(joinpath(@__DIR__, "compare.jl"))

# The private helpers below are reached through the module rather than imported: none is exported.
const ITSLIVE = ItsLiveOffsetProduction

# Keep millisecond precision where the reference string has it — DateTime cannot hold more, but
# there is no reason to throw away what it can hold.
function _parse_ref_date(str)
    length(str) <= 17 && return DateTime(str[1:17], dateformat"yyyymmddTHH:MM:SS")
    frac = rpad(str[19:min(end, 21)], 3, '0')
    return DateTime(str[1:17] * "." * frac, dateformat"yyyymmddTHH:MM:SS.sss")
end

const HANDLED_IMG_PAIR_INFO_KEYS = (
    "acquisition_date_img1", "acquisition_date_img2", "mission_img1", "mission_img2",
    "satellite_img1", "satellite_img2", "time_standard_img1", "time_standard_img2",
    "date_center", "date_dt", "latitude", "longitude", "roi_valid_percentage",
    "autoRIFT_software_version",
)

"""
    run_case(fragment; n = 101, threads = 8, force = false) -> NamedTuple

Capture `netCDF_packaging`'s real arguments for the golden case matching `fragment`, build the
equivalent `ItsLiveOffsetProduction.ItsLiveInput`, write it, and diff against the reference's own freshly-produced product
(not the possibly-older cached golden file, so `source`/`date_created`-adjacent fields are apples to
apples). Returns a result NamedTuple; never throws for a comparison mismatch, only for a genuine
failure to produce a comparable pair (capture, write, or read failure).
"""
function run_case(fragment::AbstractString; n::Integer = 101, threads::Integer = 8, force = false)
    c = only(cases(fragment))
    dir = run_dir(c, n)
    cap = capture_packaging_call(c; n, threads, force)
    k = read_packaging_capture(cap)
    a, s, info = k.arrays, k.scalars, k.img_pair_info

    is_radar = haskey(a, "offset2vr")
    pair_type = is_radar ? :radar : :optical

    dx_mean_shift = Float64(s["dx_mean_shift"])
    dy_mean_shift = Float64(s["dy_mean_shift"])
    dx_mean_shift1 = Float64(s["dx_mean_shift1"])
    dy_mean_shift1 = Float64(s["dy_mean_shift1"])
    applied = Int(s["stable_shift_applied"])
    shift_x = applied == 1 ? dx_mean_shift : applied == 2 ? dx_mean_shift1 : 0.0
    shift_y = applied == 1 ? dy_mean_shift : applied == 2 ? dy_mean_shift1 : 0.0
    dx_raw = Float32.(a["DX"] .+ shift_x)
    dy_raw = Float32.(a["DY"] .+ shift_y)

    coeffs = ItsLiveOffsetProduction.GeogridCoefficients(
        Float64.(a["offset2vx_1"]), Float64.(a["offset2vx_2"]),
        Float64.(a["offset2vy_1"]), Float64.(a["offset2vy_2"]),
        Float64.(a["scale_factor_1"]), Float64.(a["scale_factor_2"]),
        is_radar ? Float64.(a["offset2vr"]) : nothing,
        is_radar ? Float64.(a["offset2va"]) : nothing)
    refv = ItsLiveOffsetProduction.ReferenceVelocity(Float32.(a["VXref"]), Float32.(a["VYref"]))
    ssm = a["SSM"] .!= 0

    shift = ITSLIVE._stable_shift_correct(dx_raw, dy_raw, coeffs, refv, ssm)

    vx_ref, vy_ref = Float32.(a["VX"]), Float32.(a["VY"])
    both = .!isnan.(shift.vx) .& .!isnan.(vx_ref)
    vx_diff = shift.vx[both] .- vx_ref[both]
    vy_diff = shift.vy[both] .- vy_ref[both]
    physics_ok = isempty(vx_diff) || (maximum(abs, vx_diff) == 0 && maximum(abs, vy_diff) == 0)

    fresh_product = run_product(dir)

    mapping_attrs = Dict{String,Any}()
    img_pair_info_full = Dict{String,Any}()
    NCDataset(fresh_product) do ds
        for (kk, vv) in ds["mapping"].attrib
            kk in ("_FillValue", "GeoTransform") && continue
            mapping_attrs[String(kk)] = vv
        end
        for (kk, vv) in ds["img_pair_info"].attrib
            kk == "_FillValue" && continue
            img_pair_info_full[String(kk)] = vv
        end
        global source_attr = ds.attrib["source"]
    end

    tran = Float64.(s["tran"])
    x0, xres = tran[1] + tran[2] / 2, tran[2]
    y0, yres = tran[4] + tran[6] / 2, tran[6]
    ny, nx = size(a["VX"])
    x = collect(x0 .+ (0:(nx - 1)) .* xres)
    y = collect(y0 .+ (0:(ny - 1)) .* yres)

    to_lonlat = FGP.Transformation(FGP.EPSG(Int(s["epsg"])), FGP.EPSG(4326);
                                   always_xy = true)
    lonlat(px, py) = (to_lonlat(px, py)...,)

    georef = ItsLiveOffsetProduction.ItsLiveGeoref(x, y, mapping_attrs, lonlat,
                                     Float64(s["rangePixelSize"]), Float64(s["azimuthPixelSize"]))

    # `extra` is read from the *final produced product*, not the pre-patch capture — this is what
    # makes it pick up a sensor-specific post-`netCDF_packaging` patch (S1's `frame_img1`/`frame_img2`)
    # automatically, without hardcoding which sensors have one.
    extra = Dict{String,Any}(kk => vv for (kk, vv) in img_pair_info_full
                              if !(kk in HANDLED_IMG_PAIR_INFO_KEYS) && kk != "standard_name")

    img_pair_info = ItsLiveOffsetProduction.ImagePairInfo(
        _parse_ref_date(info["acquisition_date_img1"]), _parse_ref_date(info["acquisition_date_img2"]),
        info["mission_img1"], info["mission_img2"], info["satellite_img1"], info["satellite_img2"],
        Float64(info["roi_valid_percentage"]), Float64(info["latitude"]), Float64(info["longitude"]),
        extra)

    chip_size_x = round.(UInt16, a["CHIPSIZEX"])
    interp_mask = a["INTERPMASK"] .!= 0
    scale_chip_size_y = NCDataset(joinpath(dir, "autoRIFT_intermediate.nc")) do ds
        Float64(only(ds["ScaleChipSizeY"][:]))
    end

    input = ItsLiveOffsetProduction.ItsLiveInput(
        pair_type, String(s["detection_method"]), String(s["coordinates"]),
        String(info["autoRIFT_software_version"]), String(s["parameter_file"]), String(source_attr),
        is_radar ? Float64(s["dt"]) : nothing,
        georef, img_pair_info, coeffs, refv, ssm, nothing,
        dx_raw, dy_raw, chip_size_x, scale_chip_size_y, interp_mask)

    mine_path = joinpath(dir, "julia_output.nc")
    ItsLiveOffsetProduction.write_product(mine_path, input)

    mine = read_product(mine_path)
    theirs = read_product(fresh_product)
    d = compare_products(mine, theirs)
    type_diffs = attrib_type_diffs(mine_path, fresh_product)

    return (; product = c.product, platform = c.platform, pair_type, physics_ok,
            schema_mine = schema(mine), schema_theirs = schema(theirs),
            cropped_mine = cropped(mine), cropped_theirs = cropped(theirs),
            diff = d, agrees = agrees_on_data(d), type_diffs)
end

function main(args)
    isempty(args) && error("usage: validate_itslive_write.jl <case-fragment> [--force]")
    r = run_case(args[1]; force = "--force" in args)
    println("product:  ", r.product)
    println("platform: ", r.platform, "   pair_type: ", r.pair_type)
    println("physics (isolated VX/VY check): ", r.physics_ok ? "exact" : "MISMATCH")
    println("schema: mine=", r.schema_mine, "  reference=", r.schema_theirs)
    println("cropped: mine=", r.cropped_mine, "  reference=", r.cropped_theirs)
    show(stdout, r.diff)
    println()
    println("agrees_on_data: ", r.agrees)
    if isempty(r.type_diffs)
        println("attribute storage types: all shared attributes match")
    else
        println("attribute storage types differing (", length(r.type_diffs), "):")
        for k in sort(collect(keys(r.type_diffs)))
            ta, tb = r.type_diffs[k]
            println("  ", rpad(k, 34), " nc_type ", ta, "  vs  ", tb)
        end
    end
    return nothing
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS)
