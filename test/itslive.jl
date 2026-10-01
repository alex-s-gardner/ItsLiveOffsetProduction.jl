# The Sentinel-1 subswath-offset-bias correction, and `write_product`'s own wiring of it.
#
# `_swath_offset_bias_correct`/`_rotate_vel2radar` are reached through the module rather than imported:
# neither is exported and neither needs to be, since nothing outside `write_product` calls them.

using ItsLiveOffsetProduction
using NCDatasets
using Statistics
using Dates
using Random
using Rasters
using ArchGDAL
using DimensionalData
using Test

const ITSLIVE = ItsLiveOffsetProduction

# A flat displacement-to-velocity operator (`vx = dx`, `vy = dy`) over a wide, thin "image" — 3000
# range columns, so a ±500-column border band is a fraction of the width rather than the whole array.
# A real Sentinel-1 mosaic is tens of thousands of columns wide; this is the smallest width the
# algorithm's own margins still make sense on.
function _swath_case(; ny = 5, nx = 3000)
    loc_x = repeat(reshape(Float64.(0:(nx - 1)), 1, nx), ny, 1)
    loc_y = repeat(reshape(Float64.(0:(ny - 1)), ny, 1), 1, nx)
    coeffs = GeogridCoefficients(ones(ny, nx), zeros(ny, nx), zeros(ny, nx), ones(ny, nx),
                                          ones(ny, nx), ones(ny, nx), fill(0.5, ny, nx),
                                          fill(0.5, ny, nx))
    return loc_x, loc_y, coeffs
end

@testset "the bias subtraction is unconditional and cumulative past the second border" begin
    # `netcdf_output.py:1412-1447` appends the same `output_ref` entries whichever branch its
    # `nanstd` check takes — only whether the coarse-grid smoothing runs afterward depends on it, not
    # whether the bias is subtracted at all. High-variance data keeps that std check false at both
    # borders here, isolating the subtraction alone.
    loc_x, loc_y, coeffs = _swath_case()
    Random.seed!(1)
    dx = Float32.(50 .* randn(size(loc_x)))
    dy = Float32.(50 .* randn(size(loc_x)))
    border12, border23, ncols = 1000.0, 2000.0, size(loc_x, 2)

    sw = SwathOffsetBias(loc_x, loc_y, 1.0, ncols, border12, border23, false, "S1A")
    dxc, dyc = ITSLIVE._swath_offset_bias_correct(dx, dy, coeffs, sw, 1.0)
    out1, out2, out3, out4 = ITSLIVE.SWATH_OFFSET_BIAS_REF

    # One column well inside each of the three subswaths: IW1 untouched, IW2 shifted once, IW3
    # shifted by both borders' entries — the "mask2 is a superset of mask1" cumulative effect that is
    # this function's one easiest-to-get-backwards property (`netcdf_output.py:1449-1455`'s two
    # sequential, overlapping-region subtractions).
    iw1, iw2, iw3 = 501, 1501, 2501
    @test dxc[1, iw1] == dx[1, iw1] && dyc[1, iw1] == dy[1, iw1]
    @test dxc[1, iw2] - dx[1, iw2] ≈ -out1
    @test dyc[1, iw2] - dy[1, iw2] ≈ -out2
    @test dxc[1, iw3] - dx[1, iw3] ≈ -out1 - out3
    @test dyc[1, iw3] - dy[1, iw3] ≈ -out2 - out4

    # Sentinel-1B: the reference's own table sign-flipped (`netcdf_output.py:1373-1374`).
    sw_b = SwathOffsetBias(loc_x, loc_y, 1.0, ncols, border12, border23, false, "S1B")
    dxc_b, = ITSLIVE._swath_offset_bias_correct(dx, dy, coeffs, sw_b, 1.0)
    @test dxc_b[1, iw2] - dx[1, iw2] ≈ out1

    # A same-platform pair is not corrected at all, and the arrays it returns are the inputs
    # themselves — not equal copies, since a caller checking that no work was done should see it.
    sw_same = SwathOffsetBias(loc_x, loc_y, 1.0, ncols, border12, border23, true, "S1A")
    dxc_s, dyc_s = ITSLIVE._swath_offset_bias_correct(dx, dy, coeffs, sw_same, 1.0)
    @test dxc_s === dx && dyc_s === dy
end

@testset "low-variance borders trigger the coarse-grid smoothing" begin
    # Uniform dx/dy keeps `nanstd` at zero at both borders, so both trigger `_rotate_vel2radar`.
    loc_x, loc_y, coeffs = _swath_case()
    dx = fill(1.0f0, size(loc_x))
    dy = fill(2.0f0, size(loc_x))
    ncols = size(loc_x, 2)
    border12, border23 = 1000.0, 2000.0

    sw = SwathOffsetBias(loc_x, loc_y, 10.0, ncols, border12, border23, false, "S1A")
    dxc, dyc = ITSLIVE._swath_offset_bias_correct(dx, dy, coeffs, sw, 1.0)
    out1, out2, out3, out4 = ITSLIVE.SWATH_OFFSET_BIAS_REF

    # The coarse-grid roundtrip touches every point (`rotate_vel2radar` rebuilds `output_vel_x1` from
    # the smoothed coarse grid unconditionally), so nothing here should ever end up `NaN` starting
    # from finite input.
    @test !any(isnan, dxc) && !any(isnan, dyc)
    # Far from either border the step-function value survives exactly; near one, it is a smooth ramp
    # strictly between the two neighboring step levels rather than equal to either.
    @test dxc[1, 100] == 1.0f0
    @test dxc[1, end] ≈ 1.0f0 - out1 - out3
    band = dxc[1, 990:1010]
    @test minimum(band) > 1.0f0 && maximum(band) < 1.0f0 - out1
end

@testset "_rotate_vel2radar preserves the original NaN mask" begin
    # `netcdf_output.py:1279-1280`: a point that was never measured must not be filled in by the
    # coarse-grid roundtrip, even though the roundtrip has no other way to represent "no data" than
    # to leave whatever the smoothed coarse cell holds. The pre-existing NaN band (columns 1:5) is far
    # from the border's own ±500-column smoothing band (centred at 1500) so the two do not interact —
    # a border band this wide relative to the array's own width would otherwise wipe the entire row,
    # which is exactly the scale trap this file's header note on `_swath_case` exists to avoid.
    loc_x, loc_y, _ = _swath_case()
    vx = fill(1.0f0, size(loc_x))
    vy = fill(2.0f0, size(loc_x))
    vx[:, 1:5] .= NaN32
    vy[:, 1:5] .= NaN32
    out_x, out_y = ITSLIVE._rotate_vel2radar(loc_x, loc_y, vx, vy, 1500.0, 1.0, 1.0)
    @test all(isnan, out_x[:, 1:5]) && all(isnan, out_y[:, 1:5])
    @test !any(isnan, out_x[:, 6:end])
end

@testset "a coarse-cell collision resolves in the reference's own row-major order" begin
    # Four output points that all round to the *same* coarse cell (`grid_spacing_x = 100.0` against
    # an image only 56 pixels wide, so `x_grid`/`y_grid` have a single bin) — a degenerate case, but
    # one that isolates the property cleanly: the winner must be whichever point Python's own
    # `for irow: for icol:` visits *last*, i.e. `(row = 2, col = 2)`, not whichever value happens to
    # be largest. Julia's native `eachindex` on a `Matrix` visits row fastest, column slowest — the
    # transpose of that order — so this is exactly the case that order fix in `_rotate_vel2radar`
    # matters for.
    loc_x = Float64[0.0 50.0; 5.0 55.0]
    loc_y = Float64[0.0 0.0; 3.0 3.0]
    vel_x = Float32[100.0 20.0; 30.0 40.0]   # the largest value sits at (1,1), visited *first*
    vel_y = Float32[9.0 2.0; 3.0 4.0]

    out_x, out_y = ITSLIVE._rotate_vel2radar(loc_x, loc_y, vel_x, vel_y, 1000.0, 100.0, 1.0)
    @test all(==(vel_x[2, 2]), out_x)
    @test all(==(vel_y[2, 2]), out_y)
    @test vel_x[2, 2] != maximum(vel_x)   # confirms this pins traversal order, not magnitude
end

@testset "write() applies the correction before stable-shift correction, radar schema" begin
    x = collect(0.0:120.0:(120.0 * 49))
    y = collect(reverse(0.0:120.0:(120.0 * 29))) .+ 1.0e6
    ny, nx = length(y), length(x)
    georef = ItsLiveGeoref(x, y, Dict{String,Any}("grid_mapping_name" => "polar_stereographic"),
                                    (lon, lat) -> (0.0, 0.0), 2.3, 13.9)
    info = ImagePairInfo(DateTime(2020, 1, 1), DateTime(2020, 1, 13), "S", "S", "1A", "1B",
                                  80.0, 65.0, -45.0, Dict{String,Any}())
    coeffs = GeogridCoefficients(ones(ny, nx), zeros(ny, nx), zeros(ny, nx), ones(ny, nx),
                                          ones(ny, nx), ones(ny, nx), fill(0.3, ny, nx),
                                          fill(0.3, ny, nx))
    refv = ReferenceVelocity(zeros(Float32, ny, nx), zeros(Float32, ny, nx))
    loc_x = repeat(reshape(Float64.(0:(nx - 1)), 1, nx), ny, 1)
    loc_y = repeat(reshape(Float64.(0:(ny - 1)), ny, 1), 1, nx)
    sw = SwathOffsetBias(loc_x, loc_y, 1.0, nx, 20.0, 40.0, false, "S1A")

    input = ItsLiveInput(:radar, "feature", "radar", "1.0.0", "params.shp", "unit-test",
                                  12.0 * 86400, georef, info, coeffs, refv, trues(ny, nx), sw,
                                  fill(0.5f0, ny, nx), fill(0.2f0, ny, nx), fill(UInt16(32), ny, nx),
                                  1.0, falses(ny, nx))

    path = tempname() * ".nc"
    try
        write_product(path, input)
        @test isfile(path)
        NCDataset(path) do ds
            @test haskey(ds, "vx")
            @test haskey(ds, "vr")   # the radar-only schema addition
        end

        # `swath_bias = nothing` must still work — every other test in this repo that builds an
        # `ItsLiveInput` passes it, and this is the path they all take.
        input_nothing = ItsLiveInput(:radar, "feature", "radar", "1.0.0", "params.shp",
                                              "unit-test", 12.0 * 86400, georef, info, coeffs, refv,
                                              trues(ny, nx), nothing, fill(0.5f0, ny, nx),
                                              fill(0.2f0, ny, nx), fill(UInt16(32), ny, nx), 1.0,
                                              falses(ny, nx))
        path2 = tempname() * ".nc"
        write_product(path2, input_nothing)
        @test isfile(path2)
        rm(path2)
    finally
        isfile(path) && rm(path)
    end
end

@testset "ItsLiveGeoref(::Raster) survives write() unchanged, on the uncropped path" begin
    # `cf_grid_mapping` was checked in isolation (`extensions.jl`) against real captured `mapping`
    # attributes, and `ItsLiveGeoref(::Raster)` was checked in isolation against a synthetic raster
    # (`extensions.jl` again) — neither test runs the two together through `write`. This does, using
    # the *uncropped* path deliberately: `write` copies `georef.x`/`.y` straight through there rather
    # than cropping and re-aligning them (see `ImagePairInfo`'s docstring on `roi_valid_percentage`),
    # which is what makes "does the output still hold what the Raster gave it" answerable without
    # also re-deriving `write`'s crop/pad arithmetic here.
    #
    # `x[1]`/`y[1]` are a real captured golden case's own values (a polar-stereographic, EPSG:3413
    # case) rather than arbitrary ones, so this exercises the same coordinates `cf_grid_mapping`'s own
    # test does.
    xs = range(-3_440_587.5; step = 120.0, length = 10)
    ys = range(245_707.5; step = -120.0, length = 10)
    n = length(xs)
    r = Raster(zeros(Float32, n, n), (X(xs), Y(ys)); crs = Rasters.EPSG(3413))
    georef = ItsLiveGeoref(r)

    info = ImagePairInfo(DateTime(2020, 1, 1), DateTime(2020, 1, 17), "L", "L", "8", "8",
                                  0.9, 65.0, -45.0, Dict{String,Any}())   # floors to 0: uncropped
    coeffs = GeogridCoefficients(ones(n, n), zeros(n, n), zeros(n, n), ones(n, n),
                                          ones(n, n), ones(n, n), nothing, nothing)
    refv = ReferenceVelocity(zeros(Float32, n, n), zeros(Float32, n, n))
    input = ItsLiveInput(:optical, "feature", "map", "1.0.0", "params.shp", "unit-test",
                                  nothing, georef, info, coeffs, refv, trues(n, n), nothing,
                                  fill(0.5f0, n, n), fill(0.2f0, n, n), fill(UInt16(32), n, n), 1.0,
                                  falses(n, n))
    path = tempname() * ".nc"
    try
        write_product(path, input)
        NCDataset(path) do ds
            @test collect(ds["x"][:]) == georef.x
            @test collect(ds["y"][:]) == georef.y
            expected = cf_grid_mapping(3413)
            for (k, v) in expected
                got = ds["mapping"].attrib[k]
                @test (v isa AbstractString ? v == got : Float64(v) == Float64(got))
            end
        end
    finally
        isfile(path) && rm(path)
    end
end

@testset "write() refuses a swath_bias on an optical pair" begin
    # The correction has no meaning without a range/azimuth image axis, which an optical pair's
    # `dx`/`dy` are not measured against the same way a radar pair's are — the same reasoning
    # `write` already applies to `offset2vr`/`offset2va`, extended to this new field.
    x = collect(0.0:120.0:(120.0 * 9))
    y = collect(reverse(0.0:120.0:(120.0 * 9))) .+ 1.0e6
    n = length(y)
    georef = ItsLiveGeoref(x, y, Dict{String,Any}("grid_mapping_name" => "polar_stereographic"),
                                    (lon, lat) -> (0.0, 0.0), 120.0, 120.0)
    info = ImagePairInfo(DateTime(2020, 1, 1), DateTime(2020, 1, 17), "L", "L", "8", "8",
                                  80.0, 65.0, -45.0, Dict{String,Any}())
    coeffs = GeogridCoefficients(ones(n, n), zeros(n, n), zeros(n, n), ones(n, n),
                                          ones(n, n), ones(n, n), nothing, nothing)
    refv = ReferenceVelocity(zeros(Float32, n, n), zeros(Float32, n, n))
    loc = repeat(reshape(Float64.(0:(n - 1)), 1, n), n, 1)
    sw = SwathOffsetBias(loc, loc, 1.0, n, 3.0, 6.0, false, "S1A")

    input = ItsLiveInput(:optical, "feature", "map", "1.0.0", "params.shp", "unit-test",
                                  nothing, georef, info, coeffs, refv, trues(n, n), sw,
                                  fill(0.5f0, n, n), fill(0.2f0, n, n), fill(UInt16(32), n, n), 1.0,
                                  falses(n, n))
    @test_throws "requires swath_bias to be nothing" write_product(tempname() * ".nc", input)
end

@testset "roi_valid_percentage quantizes the fraction, not the percentage" begin
    # The reference quantizes to three decimals as a *fraction* and only then scales by 100
    # (`testautoRIFT.py:1366`), so the result lands on a whole tenth of a percent. 8 of 9 is
    # 0.8888…, which quantizes to 0.889 and scales to 88.9; quantizing the percentage instead would
    # keep 88.889.
    search = UInt8[1 1 1; 1 1 1; 1 1 1]
    chip = UInt16[16 16 16; 16 16 16; 16 16 0]
    @test roi_valid_percentage(chip, search) == 88.9

    # A point outside the ROI is in neither count, so zeroing it in both arrays raises the fraction.
    search2 = UInt8[1 1 1; 1 1 1; 1 1 0]
    @test roi_valid_percentage(chip, search2) == 100.0

    # Every point resolved, and none.
    @test roi_valid_percentage(fill(UInt16(16), 4, 4), fill(UInt8(1), 4, 4)) == 100.0
    @test roi_valid_percentage(zeros(UInt16, 4, 4), fill(UInt8(1), 4, 4)) == 0.0

    # An empty ROI has no denominator, and a mismatched grid is a caller error rather than a ratio.
    @test_throws "no point was searched" roi_valid_percentage(zeros(UInt16, 2, 2),
                                                                      zeros(UInt8, 2, 2))
    @test_throws DimensionMismatch roi_valid_percentage(zeros(UInt16, 2, 2),
                                                                 ones(UInt8, 3, 3))
end
