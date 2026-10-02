# Velocity, stable-shift correction, error estimates, and the ITS_LIVE product netCDF.
#
# Computes velocity, stable-shift correction, error estimates, and the radar conversion matrix from raw
# pixel offsets, and packages the result as an ITS_LIVE product netCDF — everything the reference
# implementation (`hyp3-autorift`'s `testautoRIFT.py`, `netcdf_output.py`, and `crop.py`) does across
# three files and two writes, done here from one [`ItsLiveInput`](@ref) and one write.
#
# Every formula below is transcribed from the reference, not re-derived: `testautoRIFT.py:892-1049` for
# stable-shift correction and velocity conversion, `netcdf_output.py:224-1211` for the error/stable-shift
# attribute family and the product schema, `crop.py:61-267` for the crop/align/chunk step. Two things
# the reference does are deliberately not reproduced, and are not silent gaps:
#
#   - **Sentinel-1 subswath-offset-bias correction** (`cal_swath_offset_bias`) needs burst/swath
#     geometry this package cannot obtain, so S1 products are not corrected for it.
#   - **The range-projected/observed velocity fusion** (`netcdf_output.py:321-411,1003-1081`, the
#     `VXP`/`VYP`/`alpha_sp`/`alpha_ref` block) computes `stable_shift_applied_p` and related values
#     that the reference never writes to any netCDF variable or attribute — verified by reading every
#     `setncattr`/`var[:]` call in that file. It is dead code there, so it is not ported here.

# ---------------------------------------------------------------------------
# Constants transcribed verbatim from the reference.
# ---------------------------------------------------------------------------

const NODATA_I16 = Int16(-32767)
const NODATA_F32 = Float32(-32767)
const CHUNK_SIZE = 512       # pixels, post-crop chunk edge (crop.py CHUNK_SIZE)
const PIXEL_SIZE = 120.0     # metres, crop.py's alignment grid — applied to every schema, not just optical
# `crop.py`'s own constant reads `...00:00:00Z`, but its comment says to use "exactly what xarray
# encodes" rather than the constant, and `xr.coding.times.encode_cf_datetime` normalizes the "Z" to
# an explicit offset when it re-serializes — confirmed against a real captured product, not assumed.
const TIME_UNITS = "seconds since 1980-01-06T00:00:00+00:00"
const CALENDAR = "proleptic_gregorian"
const GPS_EPOCH = DateTime(1980, 1, 6)

const TITLE = "autoRIFT surface velocities"
const AUTHOR = "Alex S. Gardner, JPL/NASA; Yang Lei, GPS/Caltech"
const INSTITUTION = "NASA Jet Propulsion Laboratory (JPL), California Institute of Technology"
const REFERENCES = join([
    "When using this data, please acknowledge the source (see global source attribute) and cite:",
    "* Gardner, A. S., Moholdt, G., Scambos, T., Fahnestock, M., Ligtenberg, S., van den Broeke, M.,",
    "  and Nilsson, J., 2018. Increased West Antarctic and unchanged East Antarctic ice discharge over",
    "  the last 7 years. The Cryosphere, 12, p.521. https://doi.org/10.5194/tc-12-521-2018",
    "* Lei, Y., Gardner, A. and Agram, P., 2021. Autonomous Repeat Image Feature Tracking (autoRIFT)",
    "  and Its Application for Tracking Ice Displacement. Remote Sensing, 13(4), p.749.",
    "  https://doi.org/10.3390/rs13040749",
    "",
    "Additionally, DOI's are provided for the software used to generate this data:",
    "* autoRIFT: https://doi.org/10.5281/zenodo.4025445",
    "* HyP3 autoRIFT plugin: https://doi.org/10.5281/zenodo.4037016",
    "* HyP3 processing environment: https://doi.org/10.5281/zenodo.3962581",
], "\n")

# The reference hardcodes this rather than fitting it per run (`testautoRIFT.py:1133-1138,1221,1229-1234`).
const ERROR_VECTOR_OPTICAL = (vx = 25.5, vy = 25.5)
const ERROR_VECTOR_RADAR = (
    vx = (slope = 0.0356, intercept = 0.5194),
    vy = (slope = 0.0501, intercept = 1.1638),
    vr = (slope = 0.0266, intercept = 0.3319),
    va = (slope = 0.0622, intercept = 1.3701),
)

const CHIP_SIZE_COORDINATES = (
    radar = "radar geometry: width = range, height = azimuth",
    optical = "image projection geometry: width = x, height = y",
)

const ERROR_STATIONARY_DESCRIPTION = "RMSE over stable surfaces, stationary or slow-flowing " *
    "surfaces with velocity < 15 meter/year identified from an external mask"
const ERROR_MODELED_DESCRIPTION = "1-sigma error calculated using a modeled error-dt relationship"
const ERROR_SLOW_DESCRIPTION = "RMSE over slowest 25% of retrieved velocities"
const STABLE_SHIFT_FLAG_DESCRIPTION = "flag for applying velocity bias correction: 0 = no " *
    "correction; 1 = correction from overlapping stable surface mask (stationary or slow-flowing " *
    "surfaces with velocity < 15 meter/year)(top priority); 2 = correction from slowest 25% of " *
    "overlapping velocities (second priority)"
const DR_TO_VR_FACTOR_DESCRIPTION = "multiplicative factor that converts slant range pixel " *
    "displacement dr to slant range velocity vr"

# ---------------------------------------------------------------------------
# Small numeric helpers.
# ---------------------------------------------------------------------------

function _finite_median(x)
    v = filter(!isnan, vec(x))
    return isempty(v) ? NaN : median(v)
end

# Only ever read when `mask` is known (by the caller) to intersect a finite region — the reference
# computes the empty-input case unconditionally too, but never writes it (`netcdf_output.py`'s
# `if stable_count != 0`/`else` branches around every use), so returning 0.0 here is never observed.
function _median_over_mask(diff::AbstractArray, mask::AbstractArray{Bool})
    v = diff[mask]
    v = filter(!isnan, v)
    return isempty(v) ? 0.0 : median(v)
end

function _std_over_mask(diff::AbstractArray, mask::AbstractArray{Bool})
    v = diff[mask]
    v = filter(!isnan, v)
    # `netcdf_output.py`'s `np.std` is uncorrected (divides by n, not n - 1); `corrected = false`
    # matches it. The two conventions agree to O(1/n) and the mask here is sometimes only tens to
    # hundreds of pixels, where that order is a visible fraction of an error of a few tens of m/year.
    return isempty(v) ? NaN : std(v; corrected = false)
end

# ---------------------------------------------------------------------------
# Sentinel-1 subswath-offset-bias correction (`netcdf_output.py:1214-1282,1345-1472`), applied to
# `dx`/`dy` before anything else — including stable-shift correction — sees them.
#
# The reference's own empirical table, only ever called with these four values
# (`testautoRIFT.py:912`); sign-flipped when the reference scene is Sentinel-1B
# (`netcdf_output.py:1373-1374`). Units match `dx`/`dy` — pixels, not velocity — despite the name
# `swath_offset_bias_ref` suggesting otherwise; the reference subtracts it from `DX`/`DY` directly.
# ---------------------------------------------------------------------------

const SWATH_OFFSET_BIAS_REF = (-0.01, 0.019, -0.0068, 0.006)

# Linear interpolation of the interior NaN runs of one row, leaving a run at either end unfilled — no
# valid neighbor to interpolate from there, matching `pandas.DataFrame.interpolate`'s own default (no
# extrapolation). A border band is hundreds of pixels from either edge of a real Sentinel-1 mosaic, so
# this case is never reached in practice; it is handled rather than assumed away only because leaving a
# real edge case silently wrong would violate more than it would save.
function _interpolate_row!(row::AbstractVector{Float32})
    n = length(row)
    i = 1
    while i <= n
        if isnan(row[i])
            j = i
            while j <= n && isnan(row[j])
                j += 1
            end
            if i > 1 && j <= n
                left, right, span = row[i - 1], row[j], Float32(j - (i - 1))
                for k in i:(j - 1)
                    row[k] = left + (right - left) * (k - (i - 1)) / span
                end
            end
            i = j
        else
            i += 1
        end
    end
    return row
end

# One border's worth of smoothing: NaN a ±500-column band around `border` on a coarse grid built at
# the output grid's own spacing (in original-image pixels), interpolate across it per row, then read
# each output point back from its nearest coarse cell — `rotate_vel2radar` in the reference. Multiple
# output points that map to the same coarse cell end up with the same smoothed value, and a point
# whose `loc_x`/`loc_y` collide with another's in scatter is overwritten by whichever the reference's
# own row-major iteration visits last; both are properties of the algorithm being reproduced, not
# defects in this port of it.
function _rotate_vel2radar(loc_x::AbstractMatrix{Float64}, loc_y::AbstractMatrix{Float64},
                            vel_x::AbstractMatrix{Float32}, vel_y::AbstractMatrix{Float32},
                            border::Float64, grid_spacing_x::Float64, scale_chip_size_y::Float64)
    ncols = round(Int, maximum(v -> isnan(v) ? -Inf : v, loc_x)) + 1
    nrows = round(Int, maximum(v -> isnan(v) ? -Inf : v, loc_y)) + 1
    skip_x = grid_spacing_x
    skip_y = round(Int, grid_spacing_x * scale_chip_size_y)
    x_grid = collect(0:skip_x:(ncols - 1))
    y_grid = collect(0:skip_y:(nrows - 1))

    nearest(grid, v) = argmin(k -> abs(grid[k] - v), eachindex(grid))

    coarse_x = fill(NaN32, length(y_grid), length(x_grid))
    coarse_y = fill(NaN32, length(y_grid), length(x_grid))
    # Row, then column within it — `netcdf_output.py`'s own `for irow: for icol:` order, not Julia's
    # column-major `eachindex`. The two orders agree on every point that lands in its own coarse
    # cell, and disagree only on which of two *colliding* points wins — rare, but a real difference
    # from the reference if left to Julia's native order rather than matched to Python's.
    for r in axes(loc_x, 1), c in axes(loc_x, 2)
        lx, ly = loc_x[r, c], loc_y[r, c]
        (isnan(lx) || isnan(ly)) && continue
        rc, cc = nearest(y_grid, ly), nearest(x_grid, lx)
        coarse_x[rc, cc] = vel_x[r, c]
        coarse_y[rc, cc] = vel_y[r, c]
    end

    shift = 500.0
    band = findall(xg -> border - shift < xg < border + shift, x_grid)
    coarse_x[:, band] .= NaN32
    coarse_y[:, band] .= NaN32
    foreach(_interpolate_row!, eachrow(coarse_x))
    foreach(_interpolate_row!, eachrow(coarse_y))

    out_x, out_y = copy(vel_x), copy(vel_y)
    for idx in eachindex(loc_x, loc_y, vel_x, vel_y)
        lx, ly = loc_x[idx], loc_y[idx]
        (isnan(lx) || isnan(ly)) && continue
        r, c = nearest(y_grid, ly), nearest(x_grid, lx)
        out_x[idx], out_y[idx] = coarse_x[r, c], coarse_y[r, c]
    end
    # A point that was already unmeasured stays unmeasured: the coarse-grid roundtrip can otherwise
    # fill it from a neighbor's cell (`netcdf_output.py:1279-1280`).
    out_x[isnan.(vel_x)] .= NaN32
    out_y[isnan.(vel_y)] .= NaN32
    return out_x, out_y
end

"""
    _swath_offset_bias_correct(dx, dy, coeffs, sw) -> (dx, dy)

`dx`/`dy` with the Sentinel-1 subswath-offset-bias correction applied — a no-op, returning `dx`/`dy`
unchanged, for a same-platform pair (`sw.same_platform`).

The correction itself is unconditional on the local velocity spread at each border: only whether the
coarse-grid smoothing runs afterward depends on it (`netcdf_output.py:1412-1447` appends the same
`output_ref` entries whichever branch the `nanstd` check takes — confirmed by reading both branches,
not assumed from the check existing).
"""
function _swath_offset_bias_correct(dx::AbstractMatrix{Float32}, dy::AbstractMatrix{Float32},
                                     coeffs::GeogridCoefficients, sw::SwathOffsetBias,
                                     scale_chip_size_y::Float64)
    sw.same_platform && return (dx, dy)

    sign = endswith(sw.reference_platform, "B") ? -1.0 : 1.0
    out1, out2, out3, out4 = sign .* SWATH_OFFSET_BIAS_REF

    loc_x = sw.location_x
    vx = Float32.(coeffs.offset2vx_1 .* (dx .* coeffs.scale_factor_1) .+
                  coeffs.offset2vx_2 .* (dy .* coeffs.scale_factor_2))
    vy = Float32.(coeffs.offset2vy_1 .* (dx .* coeffs.scale_factor_1) .+
                  coeffs.offset2vy_2 .* (dy .* coeffs.scale_factor_2))

    shift = 500.0
    _std(v) = (f = filter(!isnan, v); isempty(f) ? NaN : std(f))
    function stds_below_25(border)
        near = (loc_x .> border - shift) .& (loc_x .< border)
        far = (loc_x .> border) .& (loc_x .< border + shift)
        return _std(vx[near]) < 25 && _std(vx[far]) < 25 &&
               _std(vy[near]) < 25 && _std(vy[far]) < 25
    end
    flag12 = stds_below_25(sw.border12)
    flag23 = stds_below_25(sw.border23)

    dxc, dyc = copy(dx), copy(dy)
    m1 = (loc_x .> sw.border23) .& (loc_x .< sw.ncols)
    dxc[m1] .-= out3
    dyc[m1] .-= out4
    m2 = (loc_x .> sw.border12) .& (loc_x .< sw.ncols)
    dxc[m2] .-= out1
    dyc[m2] .-= out2

    if flag12
        dxc, dyc = _rotate_vel2radar(loc_x, sw.location_y, dxc, dyc, sw.border12, sw.grid_spacing_x,
                                     scale_chip_size_y)
    end
    if flag23
        dxc, dyc = _rotate_vel2radar(loc_x, sw.location_y, dxc, dyc, sw.border23, sw.grid_spacing_x,
                                     scale_chip_size_y)
    end
    return dxc, dyc
end

# ---------------------------------------------------------------------------
# Stable-shift correction and displacement-to-velocity conversion.
#
# `testautoRIFT.py:989-1049`: `DXref`/`DYref` are the inverse of the `offset2v*` system applied to the
# external reference velocity; the "slowest 25%" mask is derived from that reference, not supplied;
# the stable-surface mask is external and supplied. `dx`/`dy` are corrected by a single scalar median
# shift (over whichever mask has any stable coverage), then converted to velocity with the corrected
# values — the same conversion the correction step used to decide the shift in the first place.
# ---------------------------------------------------------------------------

function _stable_shift_correct(dx::AbstractMatrix{Float32}, dy::AbstractMatrix{Float32},
                                coeffs::GeogridCoefficients, refv::ReferenceVelocity,
                                ssm::AbstractMatrix{Bool})
    finite_dx = .!isnan.(dx)
    stable_count = count(ssm .& finite_dx)

    vtemp = sqrt.(refv.vx .^ 2 .+ refv.vy .^ 2)
    finite_v = filter(!isnan, vec(vtemp))
    ssm1 = if isempty(finite_v)
        falses(size(vtemp))
    else
        threshold = quantile(finite_v, 0.25)
        vtemp .<= threshold
    end
    stable_count1 = count(ssm1 .& finite_dx)

    denom = coeffs.offset2vx_1 .* coeffs.offset2vy_2 .- coeffs.offset2vx_2 .* coeffs.offset2vy_1
    dxref = (coeffs.offset2vy_2 ./ denom .* refv.vx .- coeffs.offset2vx_2 ./ denom .* refv.vy) ./
            coeffs.scale_factor_1
    dyref = (coeffs.offset2vx_1 ./ denom .* refv.vy .- coeffs.offset2vy_1 ./ denom .* refv.vx) ./
            coeffs.scale_factor_2

    dx_mean_shift = stable_count != 0 ? _median_over_mask(dx .- dxref, ssm) : 0.0
    dy_mean_shift = stable_count != 0 ? _median_over_mask(dy .- dyref, ssm) : 0.0
    dx_mean_shift1 = stable_count1 != 0 ? _median_over_mask(dx .- dxref, ssm1) : 0.0
    dy_mean_shift1 = stable_count1 != 0 ? _median_over_mask(dy .- dyref, ssm1) : 0.0

    stable_shift_applied = stable_count != 0 ? 1 : (stable_count1 != 0 ? 2 : 0)
    dxc, dyc = if stable_shift_applied == 1
        dx .- dx_mean_shift, dy .- dy_mean_shift
    elseif stable_shift_applied == 2
        dx .- dx_mean_shift1, dy .- dy_mean_shift1
    else
        dx, dy
    end

    vx = Float32.(coeffs.offset2vx_1 .* (dxc .* coeffs.scale_factor_1) .+
                  coeffs.offset2vx_2 .* (dyc .* coeffs.scale_factor_2))
    vy = Float32.(coeffs.offset2vy_1 .* (dxc .* coeffs.scale_factor_1) .+
                  coeffs.offset2vy_2 .* (dyc .* coeffs.scale_factor_2))

    return (; vx, vy, dxc, dyc, dxref, dyref, ssm1,
            stable_count, stable_count1, stable_shift_applied,
            dx_mean_shift, dy_mean_shift, dx_mean_shift1, dy_mean_shift1)
end

# ---------------------------------------------------------------------------
# Error / stable-shift attribute family, identical in form for vx, vy, vr, va
# (`netcdf_output.py:569-654,678-764,829-904,924-998`). `error_model` is either a bare number
# (optical) or a `(slope, intercept)` named tuple (radar) — the reference's own asymmetry, not a
# simplification made here.
# ---------------------------------------------------------------------------

function _error_family(comp, compref, ssm, ssm1, stable_count, stable_count1, stable_shift_applied,
                        mean_shift, mean_shift1, error_model, date_dt_days)
    error_stationary = stable_count != 0 ? _std_over_mask(comp .- compref, ssm) : NaN
    error_slow = stable_count1 != 0 ? _std_over_mask(comp .- compref, ssm1) : NaN
    error_modeled = if error_model isa NamedTuple
        (error_model.slope * date_dt_days + error_model.intercept) / date_dt_days * 365
    else
        error_model / date_dt_days * 365
    end
    error = stable_shift_applied == 1 ? error_stationary :
            stable_shift_applied == 2 ? error_slow : error_modeled
    # The reference's `else` branch is the bare Python literal `0` (`netcdf_output.py:633` etc.), not
    # a rounded float, so this is an `Int64` to match the `NC_INT64` it writes — `_round1` dispatches
    # on that to leave it alone, where `mean_shift`/`mean_shift1` are always `Float64`.
    stable_shift = stable_shift_applied == 2 ? mean_shift1 :
                   stable_shift_applied == 1 ? mean_shift : 0
    stable_shift_stationary = stable_count != 0 ? mean_shift : NaN
    stable_shift_slow = stable_count1 != 0 ? mean_shift1 : NaN
    return (; error, error_stationary, error_modeled, error_slow,
            stable_shift, stable_shift_stationary, stable_shift_slow)
end

# The per-component "mean shift as a velocity", `netcdf_output.py:270-288,309-319`: for vx/vy this is
# the *field* `offset2v_1*dx_mean_shift + offset2v_2*dy_mean_shift` (no `scale_factor` — unlike the
# velocity conversion itself), median over the stable mask; for vr/va it is `dx_mean_shift*offset2vr`
# (or `dy_mean_shift*offset2va`), median over *all* finite values, with no stable-mask restriction at
# all. Both asymmetries are the reference's, transcribed rather than reconciled.
function _velocity_mean_shift(offset2v_1, offset2v_2, dx_mean_shift, dy_mean_shift, mask, gate)
    gate || return 0.0
    # `netcdf_output.py` casts this field to float32 before masking and taking the median
    # (`temp = vx_mean_shift.astype(np.float32)`) even though `offset2v_1`/`offset2v_2` are float64 —
    # confirmed against a real capture, where skipping this cast changes the on-disk attribute type.
    field = Float32.(offset2v_1 .* dx_mean_shift .+ offset2v_2 .* dy_mean_shift)
    return _median_over_mask(field, mask)
end

function _radar_mean_shift(offset2v, mean_shift_scalar)
    return _finite_median(mean_shift_scalar .* offset2v)
end

# ---------------------------------------------------------------------------
# Chunk-alignment crop/pad (`crop.py:61-124`).
# ---------------------------------------------------------------------------

function _aligned_min(val, grid_spacing)
    nearest = floor(val / grid_spacing) * grid_spacing
    difference = val - nearest
    pixel_misalignment = mod(difference, PIXEL_SIZE)
    padding = difference - pixel_misalignment
    return val - padding, round(Int, padding / PIXEL_SIZE)
end

function _aligned_max(val, grid_spacing)
    nearest = ceil(val / grid_spacing) * grid_spacing
    difference = nearest - val
    pixel_misalignment = mod(difference, PIXEL_SIZE)
    padding = difference - pixel_misalignment
    return val + padding, round(Int, padding / PIXEL_SIZE)
end

function _valid_bbox(v::AbstractMatrix)
    valid = .!isnan.(v)
    any(valid) || return nothing
    rows = findall(vec(any(valid; dims = 2)))
    cols = findall(vec(any(valid; dims = 1)))
    return first(rows), last(rows), first(cols), last(cols)
end

function _crop_and_align(x::Vector{Float64}, y::Vector{Float64}, v::AbstractMatrix)
    bbox = _valid_bbox(v)
    bbox === nothing && return nothing
    iy1, iy2, ix1, ix2 = bbox
    xc, yc = x[ix1:ix2], y[iy1:iy2]         # x ascending, y descending
    grid_x_min, grid_x_max = xc[1], xc[end]
    grid_y_max, grid_y_min = yc[1], yc[end]

    grid_spacing = CHUNK_SIZE * PIXEL_SIZE
    ax_min, left_pad = _aligned_min(grid_x_min, grid_spacing)
    ax_max, right_pad = _aligned_max(grid_x_max, grid_spacing)
    ay_min, bottom_pad = _aligned_min(grid_y_min, grid_spacing)
    ay_max, top_pad = _aligned_max(grid_y_max, grid_spacing)

    nx = round(Int, (ax_max - ax_min) / PIXEL_SIZE) + 1
    ny = round(Int, (ay_max - ay_min) / PIXEL_SIZE) + 1
    xa = ax_min .+ (0:(nx - 1)) .* PIXEL_SIZE
    ya = ay_max .- (0:(ny - 1)) .* PIXEL_SIZE

    return (; iy1, iy2, ix1, ix2, xa, ya, top_pad, bottom_pad, left_pad, right_pad, nx, ny)
end

function _pad(data::AbstractMatrix{T}, top, bottom, left, right, fillvalue) where {T}
    ny, nx = size(data)
    out = fill(T(fillvalue), ny + top + bottom, nx + left + right)
    out[(top + 1):(top + ny), (left + 1):(left + nx)] .= data
    return out
end

function _crop_pad(data::AbstractMatrix, align, fillvalue)
    cropped = data[align.iy1:align.iy2, align.ix1:align.ix2]
    return _pad(cropped, align.top_pad, align.bottom_pad, align.left_pad, align.right_pad, fillvalue)
end

# ---------------------------------------------------------------------------
# Quantization to on-disk dtype (`netcdf_output.py`'s `np.round(np.clip(...)).astype(...)` calls).
# ---------------------------------------------------------------------------

_quantize_i16(data::AbstractArray{<:AbstractFloat}, nodata) =
    round.(Int16, clamp.(ifelse.(nodata, Float64(NODATA_I16), data), -32768.0, 32767.0))
_quantize_u16(data::AbstractArray{<:Real}) = round.(UInt16, clamp.(Float64.(data), 0.0, 65535.0))
_fill_f32(data::AbstractArray{<:AbstractFloat}, nodata) = ifelse.(nodata, NODATA_F32, Float32.(data))

# ---------------------------------------------------------------------------
# netCDF plumbing.
# ---------------------------------------------------------------------------

# Whether `_setattr!` forces a string attribute to `NC_STRING`, scoped around one `_write_file!`
# call by `_with_netcdf_create` so concurrent `write_product` calls on different tasks cannot race.
#
# The reference writes a cropped product's string attributes as `NC_STRING` (confirmed against a
# real product's `ncdump -h`: `string x:standard_name = "..."`) but the uncropped (`P000`) product's
# as `NC_CHAR` — confirmed against both a real cropped and a real uncropped product. `netCDF_packaging`
# writes both with the plain Python `netCDF4` library, which defaults to `NC_CHAR`; only the cropped
# path's extra reopen-and-patch in `crop.py` goes through `xarray`, which is where the upgrade to
# `NC_STRING` happens. So this package's own on-disk type has to track the same cropped/uncropped
# split, not pick one unconditionally.
const _FORCE_NC_STRING = ScopedValue(false)

# NCDatasets matches `NC_STRING` only when the value is wrapped in a one-element `Vector{String}`
# (its own `defAttrib` docstring) — a bare `String` writes `NC_CHAR` instead. `Int64` goes through
# `nc_put_att` directly with `NC_INT64`, since the plain `x.attrib[key] = val` path narrows a bare
# `Int64` to `Int32`.
function _setattr!(x, key, val)
    if val isa Int64 && hasproperty(x, :var)
        raw = x.var
        NCDatasets.nc_put_att(raw.ds.ncid, raw.varid, key, NCDatasets.NC_INT64, Int64[val])
    elseif val isa AbstractString && _FORCE_NC_STRING[]
        x.attrib[key] = [val]
    else
        x.attrib[key] = val
    end
    return x
end

# Opens `path` for writing with `_FORCE_NC_STRING` scoped to `cropped`, so every `_setattr!` call
# `f` makes sees the on-disk string type the reference itself would use for this product's schema.
function _with_netcdf_create(f, path::AbstractString, cropped::Bool)
    with(_FORCE_NC_STRING => cropped) do
        NCDatasets.NCDataset(path, "c") do ds
            f(ds)
        end
    end
end

# `data` is (y, x) — AutoRIFT's own convention. NCDatasets keeps a netCDF `(y,x)` variable's Julia-side
# dims as `(x,y)` (memory layout, not index order; verified against a round-trip `ncdump -h`), so
# every write here transposes.
function _put!(v, data::AbstractMatrix)
    arr = permutedims(data)
    if ndims(v) == 3
        v[:, :, 1] = arr
    else
        v[:, :] = arr
    end
    return v
end

function _defvar_scalar!(ds, name, fillvalue, attrs, dimnames = ())
    v = NCDatasets.defVar(ds, name, typeof(fillvalue), dimnames; fillvalue)
    for (k, val) in attrs
        _setattr!(v, k, val)
    end
    return v
end

function _defvar_coord!(ds, name, data, standard_name, description; fillvalue = nothing)
    v = NCDatasets.defVar(ds, name, Float64, (name,); fillvalue)
    _setattr!(v, "standard_name", standard_name)
    _setattr!(v, "description", description)
    _setattr!(v, "units", "m")
    v[:] = data
    return v
end

function _defvar_2d!(ds, name, ::Type{T}, fillvalue, dimnames, chunksizes) where {T}
    return NCDatasets.defVar(ds, name, T, dimnames; fillvalue = T(fillvalue),
                              deflatelevel = 2, shuffle = true, chunksizes)
end

_round1(x::Integer) = x
_round1(x::AbstractFloat) = isnan(x) ? x : round(x; digits = 1)

# Python's `str(float)` never uses scientific notation in the range these coordinates occupy, and
# every one of them is an exact multiple of half the 120 m pixel size, so one decimal digit is exact
# — but Julia's default `string`/interpolation switches to scientific notation past ~1e5, which
# `"$(x)"` in the `GeoTransform` string would otherwise silently reproduce as a value Python would
# never write. Built from `abs(x)` and reattached sign, rather than formatting the signed value
# directly: `"$(whole).$(tenth)"` for a negative `whole` reads as `whole - tenth/10`, not
# `whole + tenth/10`, which is off by a whole unit exactly at the `.5` values these coordinates are.
function _geo_num(x::Float64)
    s = x < 0 ? "-" : ""
    tenths = round(Int, abs(x) * 10)
    whole, tenth = divrem(tenths, 10)
    return string(s, whole, ".", tenth)
end

function _write_velocity_error_attribs!(v, stats, description, stable_shift_applied,
                                         stable_count, stable_count1)
    _setattr!(v, "error", _round1(stats.error))
    _setattr!(v, "error_description", description)
    _setattr!(v, "error_stationary", _round1(stats.error_stationary))
    _setattr!(v, "error_stationary_description", ERROR_STATIONARY_DESCRIPTION)
    _setattr!(v, "error_modeled", _round1(stats.error_modeled))
    _setattr!(v, "error_modeled_description", ERROR_MODELED_DESCRIPTION)
    _setattr!(v, "error_slow", _round1(stats.error_slow))
    _setattr!(v, "error_slow_description", ERROR_SLOW_DESCRIPTION)
    _setattr!(v, "stable_shift", _round1(stats.stable_shift))
    _setattr!(v, "stable_shift_flag", stable_shift_applied)
    _setattr!(v, "stable_shift_flag_description", STABLE_SHIFT_FLAG_DESCRIPTION)
    _setattr!(v, "stable_shift_stationary", _round1(stats.stable_shift_stationary))
    _setattr!(v, "stable_count_stationary", stable_count)
    _setattr!(v, "stable_shift_slow", _round1(stats.stable_shift_slow))
    _setattr!(v, "stable_count_slow", stable_count1)
    return v
end

# ---------------------------------------------------------------------------
# Metadata assembly (`testautoRIFT.py`'s `date_dt`/`date_center`/`satellite` computation).
# ---------------------------------------------------------------------------

const MISSION_NAMES = Dict("L" => "Landsat ", "S" => "Sentinel-", "N" => "NISAR")

# `testautoRIFT.py:1366`. The ratio is quantized to three decimals as a *fraction* and only then scaled
# to a percentage, so `0.784` becomes `78.4` — quantizing the percentage instead would keep a digit the
# reference discards. `round` is half-to-even on both sides.
function roi_valid_percentage(chip_size_x::AbstractArray,
                                        search_limit_x::AbstractArray)
    axes(chip_size_x) == axes(search_limit_x) || throw(DimensionMismatch(
        "chip_size_x and search_limit_x must cover the same grid: $(axes(chip_size_x)) vs " *
        "$(axes(search_limit_x))"))
    searched = count(!iszero, search_limit_x)
    searched == 0 && throw(ArgumentError(
        "every search limit is zero, so no point was searched and the valid fraction has no " *
        "denominator"))
    resolved = count(!iszero, chip_size_x)
    return round(Int, resolved / searched * 1000) / 1000 * 100
end

function _satellite_attribute(info::ImagePairInfo)
    haskey(MISSION_NAMES, info.mission_img1) ||
        throw(ArgumentError("unrecognized mission code $(info.mission_img1)"))
    haskey(MISSION_NAMES, info.mission_img2) ||
        throw(ArgumentError("unrecognized mission code $(info.mission_img2)"))
    # `string`, not `*`: `satellite_img1`/`satellite_img2` carry whatever type the source metadata
    # gives (a bare integer for NISAR), matching Python's f-string, which stringifies either way.
    s1 = MISSION_NAMES[info.mission_img1] * string(info.satellite_img1)
    s2 = MISSION_NAMES[info.mission_img2] * string(info.satellite_img2)
    return s1 == s2 ? s1 : "$s1 and $s2"
end

# Julia's `DateTime` is millisecond-precision, three orders coarser than the reference's microsecond
# timestamps; `date_center` loses precision this cannot recover. Unlike `time`/`date_created`, this
# attribute is not excluded from the golden comparison — a real, usually negligible, source of
# mismatch, not a choice made here.
function _date_dt_center(d0::DateTime, d1::DateTime)
    d1 >= d0 || throw(ArgumentError("acquisition_date_img1 must not be after acquisition_date_img2"))
    date_dt = Dates.value(d1 - d0) / 86_400_000.0
    date_center = d0 + (d1 - d0) ÷ 2
    return date_dt, date_center
end

_format_date(dt::DateTime) = Dates.format(dt, "yyyymmddTHH:MM:SS.sss")

# Landsat and Sentinel-2 both format `acquisition_date_img1`/`img2` the same way as `date_center` —
# `strftime(...).rstrip('0')` (`testautoRIFT.py:1527-1528` (S2), and the equivalent in the Landsat
# branch) — confirmed against a real captured `IMG_INFO_DICT`, not assumed. Sentinel-1's are sourced
# from a pre-formatted metadata string instead and may not follow this rule; unconfirmed against a
# real capture, since only an optical case has been checked so far.
_format_date_stripped(dt::DateTime) = rstrip(_format_date(dt), '0')

function _img_pair_info_attrs(info::ImagePairInfo, date_dt, date_center_str)
    d = Dict{String,Any}(
        "acquisition_date_img1" => _format_date_stripped(info.acquisition_date_img1),
        "acquisition_date_img2" => _format_date_stripped(info.acquisition_date_img2),
        "mission_img1" => info.mission_img1,
        "mission_img2" => info.mission_img2,
        "satellite_img1" => info.satellite_img1,
        "satellite_img2" => info.satellite_img2,
        "time_standard_img1" => "UTC",
        "time_standard_img2" => "UTC",
        "date_center" => date_center_str,
        "date_dt" => date_dt,
        "latitude" => info.latitude,
        "longitude" => info.longitude,
        "roi_valid_percentage" => info.roi_valid_percentage,
    )
    merge!(d, info.extra)
    return d
end

# ---------------------------------------------------------------------------
# `write_product`
# ---------------------------------------------------------------------------

function write_product(path::AbstractString, input::ItsLiveInput)
    input.pair_type in (:radar, :optical) ||
        throw(ArgumentError("pair_type must be :radar or :optical, got $(input.pair_type)"))
    is_radar = input.pair_type === :radar
    coeffs = input.coefficients
    if is_radar
        (coeffs.offset2vr === nothing || coeffs.offset2va === nothing) &&
            throw(ArgumentError("pair_type == :radar requires offset2vr and offset2va"))
        input.dt_seconds === nothing &&
            throw(ArgumentError("pair_type == :radar requires dt_seconds"))
    else
        (coeffs.offset2vr !== nothing || coeffs.offset2va !== nothing) &&
            throw(ArgumentError("pair_type == :optical requires offset2vr and offset2va to be nothing"))
        input.swath_bias === nothing ||
            throw(ArgumentError("pair_type == :optical requires swath_bias to be nothing — the " *
                                 "Sentinel-1 subswath-offset-bias correction has no meaning for an " *
                                 "optical pair"))
    end

    y, x = input.georef.y, input.georef.x
    size(input.dx) == size(input.dy) == (length(y), length(x)) ||
        throw(DimensionMismatch("dx, dy, and the x/y coordinate vectors must share one grid"))
    length(x) > 1 && !isapprox(x[2] - x[1], PIXEL_SIZE; atol = 1.0e-6) &&
        throw(ArgumentError("georef.x must be spaced at exactly $(PIXEL_SIZE) m, got $(x[2] - x[1])"))
    length(y) > 1 && !isapprox(y[1] - y[2], PIXEL_SIZE; atol = 1.0e-6) &&
        throw(ArgumentError("georef.y must be spaced at exactly -$(PIXEL_SIZE) m (descending), " *
                             "got step $(y[2] - y[1])"))

    dx, dy = input.dx, input.dy
    if input.swath_bias !== nothing
        dx, dy = _swath_offset_bias_correct(dx, dy, coeffs, input.swath_bias, input.scale_chip_size_y)
    end

    shift = _stable_shift_correct(dx, dy, coeffs, input.reference, input.stable_mask)
    vx, vy = shift.vx, shift.vy
    v = sqrt.(vx .^ 2 .+ vy .^ 2)
    nodata = isnan.(vx) .| isnan.(vy)

    date_dt, date_center_dt = _date_dt_center(input.img_pair_info.acquisition_date_img1,
                                               input.img_pair_info.acquisition_date_img2)

    # --- error/stable-shift attribute family: vx, vy, and (radar only) vr, va ---
    ssm, ssm1 = input.stable_mask, shift.ssm1
    stable_count, stable_count1 = shift.stable_count, shift.stable_count1
    applied = shift.stable_shift_applied

    vx_mean_shift = _velocity_mean_shift(coeffs.offset2vx_1, coeffs.offset2vx_2,
                                          shift.dx_mean_shift, shift.dy_mean_shift, ssm, stable_count != 0)
    vy_mean_shift = _velocity_mean_shift(coeffs.offset2vy_1, coeffs.offset2vy_2,
                                          shift.dx_mean_shift, shift.dy_mean_shift, ssm, stable_count != 0)
    vx_mean_shift1 = _velocity_mean_shift(coeffs.offset2vx_1, coeffs.offset2vx_2,
                                           shift.dx_mean_shift1, shift.dy_mean_shift1, ssm1, stable_count1 != 0)
    vy_mean_shift1 = _velocity_mean_shift(coeffs.offset2vy_1, coeffs.offset2vy_2,
                                           shift.dx_mean_shift1, shift.dy_mean_shift1, ssm1, stable_count1 != 0)

    vxref, vyref = input.reference.vx, input.reference.vy
    vx_error_model = is_radar ? ERROR_VECTOR_RADAR.vx : ERROR_VECTOR_OPTICAL.vx
    vy_error_model = is_radar ? ERROR_VECTOR_RADAR.vy : ERROR_VECTOR_OPTICAL.vy
    vx_stats = _error_family(vx, vxref, ssm, ssm1, stable_count, stable_count1, applied,
                              vx_mean_shift, vx_mean_shift1, vx_error_model, date_dt)
    vy_stats = _error_family(vy, vyref, ssm, ssm1, stable_count, stable_count1, applied,
                              vy_mean_shift, vy_mean_shift1, vy_error_model, date_dt)

    v_error_scalar = _v_error_monte_carlo(vx_stats.error, vy_stats.error)
    V_error = sqrt.((vx_stats.error .* vx ./ v) .^ 2 .+ (vy_stats.error .* vy ./ v) .^ 2)
    V_error[v .== 0] .= v_error_scalar

    vr = va = m11 = m12 = nothing
    dr_to_vr_factor = nothing
    vr_stats = va_stats = nothing
    pixel_size_y = input.georef.pixel_size_y
    if is_radar
        # `netcdf_output.py` explicitly casts `VRref`/`VAref` to float32 (`VRref = VRref.astype(np.float32)`)
        # even though `DXref`/`DYref`/`offset2vr` are float64 — confirmed against a real capture.
        vrref = Float32.(shift.dxref .* coeffs.offset2vr)
        varef = Float32.(shift.dyref .* coeffs.offset2va)
        vr = Float32.(shift.dxc .* coeffs.offset2vr)
        va = Float32.(shift.dyc .* coeffs.offset2va)

        vr_mean_shift = _radar_mean_shift(coeffs.offset2vr, shift.dx_mean_shift)
        va_mean_shift = _radar_mean_shift(coeffs.offset2va, shift.dy_mean_shift)
        vr_mean_shift1 = _radar_mean_shift(coeffs.offset2vr, shift.dx_mean_shift1)
        va_mean_shift1 = _radar_mean_shift(coeffs.offset2va, shift.dy_mean_shift1)

        vr_stats = _error_family(vr, vrref, ssm, ssm1, stable_count, stable_count1, applied,
                                  vr_mean_shift, vr_mean_shift1, ERROR_VECTOR_RADAR.vr, date_dt)
        va_stats = _error_family(va, varef, ssm, ssm1, stable_count, stable_count1, applied,
                                  va_mean_shift, va_mean_shift1, ERROR_VECTOR_RADAR.va, date_dt)

        dr_to_vr_factor = _finite_median(coeffs.offset2vr)
        denom = coeffs.offset2vx_1 .* coeffs.offset2vy_2 .- coeffs.offset2vx_2 .* coeffs.offset2vy_1
        m11 = Float32.(coeffs.offset2vy_2 ./ denom ./ coeffs.scale_factor_1)
        m12 = Float32.(.-coeffs.offset2vx_2 ./ denom ./ coeffs.scale_factor_1)

        seconds_per_year = 365.0 * 24.0 * 3600.0
        pixel_size_y = _finite_median(coeffs.offset2va) * input.dt_seconds / seconds_per_year
    end

    # --- chip sizes and interpolation mask, in metres/final form, nodata-zeroed ---
    chip_size_y_pixels = round.(input.chip_size_x .* input.scale_chip_size_y ./ 2) .* 2
    chip_size_width = Float64.(input.chip_size_x) .* input.georef.pixel_size_x
    chip_size_height = chip_size_y_pixels .* pixel_size_y
    chip_size_width[nodata] .= 0
    chip_size_height[nodata] .= 0
    interp_mask = copy(input.interp_mask)
    interp_mask[nodata] .= false

    # The reference's crop/no-crop choice is `process.py`'s literal filename check,
    # `not netcdf_file.name.endswith('_P000.nc')` — and that suffix is `floor(PPP)` where `PPP` is
    # `roi_valid_percentage`, not a fresh check of whether any pixel is valid. The two disagree
    # whenever `roi_valid_percentage` rounds down to zero without actually being zero — confirmed
    # against a real case with 19,994 valid pixels (0.39% of the grid) that the reference still
    # leaves uncropped, because its `roi_valid_percentage` is `0.9` and `floor(0.9) == 0`.
    align = floor(input.img_pair_info.roi_valid_percentage) == 0 ? nothing : _crop_and_align(x, y, v)

    img_pair_info_dict = _img_pair_info_attrs(input.img_pair_info, date_dt, "")   # date_center filled below
    mapping_dict = copy(input.georef.mapping_attrs)

    if align === nothing
        _write_file!(path, input, is_radar, x, y, vx, vy, v, V_error, vr, va, m11, m12, dr_to_vr_factor,
                     chip_size_width, chip_size_height, interp_mask, nodata,
                     vx_stats, vy_stats, vr_stats, va_stats, applied, stable_count, stable_count1,
                     pixel_size_y, img_pair_info_dict, mapping_dict, date_center_dt, nothing)
    else
        vxc = _crop_pad(vx, align, NaN32); vyc = _crop_pad(vy, align, NaN32)
        vcc = _crop_pad(v, align, NaN32); vec_ = _crop_pad(V_error, align, NaN32)
        vrc = vr === nothing ? nothing : _crop_pad(vr, align, NaN32)
        vac = va === nothing ? nothing : _crop_pad(va, align, NaN32)
        m11c = m11 === nothing ? nothing : _crop_pad(m11, align, NaN32)
        m12c = m12 === nothing ? nothing : _crop_pad(m12, align, NaN32)
        cszw = _crop_pad(chip_size_width, align, 0.0)
        cszh = _crop_pad(chip_size_height, align, 0.0)
        imask = _crop_pad(interp_mask, align, false)
        nodatac = _crop_pad(nodata, align, true)

        center_x = (align.xa[1] + align.xa[end]) / 2
        center_y = (align.ya[1] + align.ya[end]) / 2
        lon, lat = input.georef.lonlat(center_x, center_y)
        img_pair_info_dict["latitude"] = round(lat; digits = 2)
        img_pair_info_dict["longitude"] = round(lon; digits = 2)

        _write_file!(path, input, is_radar, align.xa, align.ya, vxc, vyc, vcc, vec_, vrc, vac, m11c, m12c,
                     dr_to_vr_factor, cszw, cszh, imask, nodatac,
                     vx_stats, vy_stats, vr_stats, va_stats, applied, stable_count, stable_count1,
                     pixel_size_y, img_pair_info_dict, mapping_dict, date_center_dt, align)
    end
    return path
end

function _v_error_monte_carlo(vx_error, vy_error; n = 1_000_000)
    vx_samples = vx_error .* randn(n)
    vy_samples = vy_error .* randn(n)
    return std(sqrt.(vx_samples .^ 2 .+ vy_samples .^ 2))
end

function _write_file!(path, input, is_radar, x, y, vx, vy, v, V_error, vr, va, m11, m12, dr_to_vr_factor,
                       chip_size_width, chip_size_height, interp_mask, nodata,
                       vx_stats, vy_stats, vr_stats, va_stats, stable_shift_applied,
                       stable_count, stable_count1, pixel_size_y,
                       img_pair_info_dict, mapping_dict, date_center_dt, align)
    cropped = align !== nothing
    ny, nx = length(y), length(x)

    if cropped
        dims2 = ("x", "y", "time")
        chunks2 = (CHUNK_SIZE, CHUNK_SIZE, 1)
        x_cell, y_cell = x[2] - x[1], y[2] - y[1]
    else
        dims2 = ("x", "y")
        chunk_lines = round(Int, min(ceil(8192 / ny) * 128, ny))
        chunks2 = (nx, chunk_lines)
        x_cell = length(x) > 1 ? x[2] - x[1] : PIXEL_SIZE
        y_cell = length(y) > 1 ? y[2] - y[1] : -PIXEL_SIZE
    end
    # The literal `0`/`0.0` rotation terms differ by path: `crop.py`'s f-string embeds a bare `0`,
    # but the uncropped path's is `str(0.0)` from a Python float list (`netcdf_output.py`'s own
    # `tran` array) — confirmed against a real uncropped product, not assumed identical to the
    # cropped path's convention.
    zero = cropped ? "0" : "0.0"
    mapping_dict["GeoTransform"] = "$(_geo_num(x[1])) $(_geo_num(x_cell)) $zero $(_geo_num(y[1])) $zero $(_geo_num(y_cell))"
    date_center_str = _format_date(date_center_dt)
    date_center_str = rstrip(date_center_str, '0')
    img_pair_info_dict["date_center"] = date_center_str

    _with_netcdf_create(path, cropped) do ds
        _setattr!(ds, "GDAL_AREA_OR_POINT", "Area")
        _setattr!(ds, "Conventions", "CF-1.8")
        _setattr!(ds, "date_created", Dates.format(Dates.now(), "dd-u-yyyy HH:MM:SS"))
        _setattr!(ds, "title", TITLE)
        _setattr!(ds, "autoRIFT_software_version", input.autorift_software_version)
        _setattr!(ds, "autoRIFT_parameter_file", input.parameter_file)
        _setattr!(ds, "scene_pair_type", String(input.pair_type))
        _setattr!(ds, "satellite", _satellite_attribute(input.img_pair_info))
        _setattr!(ds, "motion_detection_method", input.detection_method)
        _setattr!(ds, "motion_coordinates", input.coordinates)
        _setattr!(ds, "author", AUTHOR)
        _setattr!(ds, "institution", INSTITUTION)
        _setattr!(ds, "source", input.source)
        _setattr!(ds, "references", REFERENCES)

        # Declared in the reference's own dimension order (`time, y, x`) — cosmetic (a reader
        # accesses dimensions by name, not position), but free to match.
        cropped && NCDatasets.defDim(ds, "time", Inf)
        NCDatasets.defDim(ds, "y", ny)
        NCDatasets.defDim(ds, "x", nx)

        # `crop.py`'s explicit `encoding` dict only covers `ds.data_vars`, not the coordinate
        # variables — `x`/`y`/`time` fall through to xarray's default float64 encoding, which sets
        # `_FillValue = NaN` unless told otherwise. Confirmed against a real captured product: the
        # uncropped (`P000`) path is written directly by `netCDF_packaging`'s literal
        # `fill_value=None` instead, and genuinely has no fill value there.
        coord_fill = cropped ? NaN : nothing
        _defvar_coord!(ds, "x", x, "projection_x_coordinate", "x coordinate of projection"; fillvalue = coord_fill)
        _defvar_coord!(ds, "y", y, "projection_y_coordinate", "y coordinate of projection"; fillvalue = coord_fill)

        if cropped
            t = NCDatasets.defVar(ds, "time", Float64, ("time",); fillvalue = NaN)
            jitter_us = mod(hash(basename(path)), 1_000_000)
            seconds = (Dates.value(date_center_dt - GPS_EPOCH) * 1000 + jitter_us) / 1.0e6
            # `t[:] = [seconds]` on an unlimited dim at record-count 0 silently leaves netCDF's
            # default fill in place instead of growing the dimension; scalar record assignment does not.
            t[1] = seconds
            _setattr!(t, "standard_name", "time")
            _setattr!(t, "description", "mid-date between acquisition_date_img1 and acquisition_date_img2 " *
                                       "with microseconds added to ensure uniqueness.")
            _setattr!(t, "units", TIME_UNITS)
            _setattr!(t, "calendar", CALENDAR)
            _setattr!(t, "microseconds_added", jitter_us)
            _setattr!(t, "microseconds_added_description", "6-digit numeric hash of the filename.")
        end

        _defvar_scalar!(ds, "mapping", Int8(-127), mapping_dict)
        # `crop.py` resets `img_pair_info` to its fill value *before* `expand_dims(dim='time', ...)`,
        # and resets `mapping` *after` — so `img_pair_info` ends up dimensioned `(time)` and `mapping`
        # stays scalar, an accident of operation order rather than a schema decision. Confirmed
        # against a real product's `ncdump -h`; matched here rather than corrected, per this
        # package's rule of reproducing the reference exactly.
        img = _defvar_scalar!(ds, "img_pair_info", Int8(-127),
                               Dict{String,Any}("standard_name" => "image_pair_information"),
                               cropped ? ("time",) : ())
        for (k, val) in img_pair_info_dict
            _setattr!(img, k, val)
        end

        vx_desc = is_radar ? "velocity component in x direction from radar range and azimuth measurements" :
                             "velocity component in x direction"
        vy_desc = is_radar ? "velocity component in y direction from radar range and azimuth measurements" :
                             "velocity component in y direction"
        v_error_desc = is_radar ? "velocity magnitude error from radar range and azimuth measurements" :
                                  "velocity magnitude error"

        vxv = _defvar_2d!(ds, "vx", Int16, NODATA_I16, dims2, chunks2)
        _setattr!(vxv, "standard_name", "land_ice_surface_x_velocity")
        _setattr!(vxv, "description", vx_desc)
        _setattr!(vxv, "units", "meter/year")
        _setattr!(vxv, "grid_mapping", "mapping")
        _write_velocity_error_attribs!(vxv, vx_stats, "best estimate of x_velocity error: vx_error is " *
            "populated according to the approach used for the velocity bias correction as indicated in " *
            "\"stable_shift_flag\"", stable_shift_applied, stable_count, stable_count1)
        _put!(vxv, _quantize_i16(vx, nodata))

        vyv = _defvar_2d!(ds, "vy", Int16, NODATA_I16, dims2, chunks2)
        _setattr!(vyv, "standard_name", "land_ice_surface_y_velocity")
        _setattr!(vyv, "description", vy_desc)
        _setattr!(vyv, "units", "meter/year")
        _setattr!(vyv, "grid_mapping", "mapping")
        _write_velocity_error_attribs!(vyv, vy_stats, "best estimate of y_velocity error: vy_error is " *
            "populated according to the approach used for the velocity bias correction as indicated in " *
            "\"stable_shift_flag\"", stable_shift_applied, stable_count, stable_count1)
        _put!(vyv, _quantize_i16(vy, nodata))

        vv = _defvar_2d!(ds, "v", Int16, NODATA_I16, dims2, chunks2)
        _setattr!(vv, "standard_name", "land_ice_surface_velocity")
        _setattr!(vv, "description", "velocity magnitude")
        _setattr!(vv, "units", "meter/year")
        _setattr!(vv, "grid_mapping", "mapping")
        _put!(vv, _quantize_i16(v, nodata))

        vev = _defvar_2d!(ds, "v_error", Int16, NODATA_I16, dims2, chunks2)
        _setattr!(vev, "standard_name", "velocity_error")
        _setattr!(vev, "description", v_error_desc)
        _setattr!(vev, "units", "meter/year")
        _setattr!(vev, "grid_mapping", "mapping")
        _put!(vev, _quantize_i16(V_error, nodata))

        if is_radar
            vrv = _defvar_2d!(ds, "vr", Int16, NODATA_I16, dims2, chunks2)
            _setattr!(vrv, "standard_name", "range_velocity")
            _setattr!(vrv, "description", "velocity in radar range direction")
            _setattr!(vrv, "units", "meter/year")
            _setattr!(vrv, "grid_mapping", "mapping")
            _write_velocity_error_attribs!(vrv, vr_stats, "best estimate of range_velocity error: " *
                "vr_error is populated according to the approach used for the velocity bias correction " *
                "as indicated in \"stable_shift_flag\"", stable_shift_applied, stable_count, stable_count1)
            _put!(vrv, _quantize_i16(vr, nodata))

            vav = _defvar_2d!(ds, "va", Int16, NODATA_I16, dims2, chunks2)
            _setattr!(vav, "standard_name", "azimuth_velocity")
            _setattr!(vav, "description", "velocity in radar azimuth direction")
            _setattr!(vav, "units", "meter/year")
            _setattr!(vav, "grid_mapping", "mapping")
            _write_velocity_error_attribs!(vav, va_stats, "best estimate of azimuth_velocity error: " *
                "va_error is populated according to the approach used for the velocity bias correction " *
                "as indicated in \"stable_shift_flag\"", stable_shift_applied, stable_count, stable_count1)
            _put!(vav, _quantize_i16(va, nodata))

            m11v = _defvar_2d!(ds, "M11", Float32, NODATA_F32, dims2, chunks2)
            _setattr!(m11v, "standard_name", "conversion_matrix_element_11")
            _setattr!(m11v, "description", "conversion matrix element (1st row, 1st column) that can be " *
                "multiplied with vx to give range pixel displacement dr (see Eq. A18 in " *
                "https://www.mdpi.com/2072-4292/13/4/749)")
            _setattr!(m11v, "units", "pixel/(meter/year)")
            _setattr!(m11v, "grid_mapping", "mapping")
            _setattr!(m11v, "dr_to_vr_factor", dr_to_vr_factor)
            _setattr!(m11v, "dr_to_vr_factor_description", DR_TO_VR_FACTOR_DESCRIPTION)
            _put!(m11v, _fill_f32(m11, nodata))

            m12v = _defvar_2d!(ds, "M12", Float32, NODATA_F32, dims2, chunks2)
            _setattr!(m12v, "standard_name", "conversion_matrix_element_12")
            _setattr!(m12v, "description", "conversion matrix element (1st row, 2nd column) that can be " *
                "multiplied with vy to give range pixel displacement dr (see Eq. A18 in " *
                "https://www.mdpi.com/2072-4292/13/4/749)")
            _setattr!(m12v, "units", "pixel/(meter/year)")
            _setattr!(m12v, "grid_mapping", "mapping")
            _setattr!(m12v, "dr_to_vr_factor", dr_to_vr_factor)
            _setattr!(m12v, "dr_to_vr_factor_description", DR_TO_VR_FACTOR_DESCRIPTION)
            _put!(m12v, _fill_f32(m12, nodata))
        end

        cszwv = _defvar_2d!(ds, "chip_size_width", UInt16, UInt16(0), dims2, chunks2)
        _setattr!(cszwv, "standard_name", "chip_size_width")
        _setattr!(cszwv, "description", "width of search template (chip)")
        _setattr!(cszwv, "units", "m")
        _setattr!(cszwv, "grid_mapping", "mapping")
        if is_radar
            _setattr!(cszwv, "range_pixel_size", input.georef.pixel_size_x)
            _setattr!(cszwv, "chip_size_coordinates", CHIP_SIZE_COORDINATES.radar)
        else
            _setattr!(cszwv, "x_pixel_size", input.georef.pixel_size_x)
            _setattr!(cszwv, "chip_size_coordinates", CHIP_SIZE_COORDINATES.optical)
        end
        _put!(cszwv, _quantize_u16(chip_size_width))

        cszhv = _defvar_2d!(ds, "chip_size_height", UInt16, UInt16(0), dims2, chunks2)
        _setattr!(cszhv, "standard_name", "chip_size_height")
        _setattr!(cszhv, "description", "height of search template (chip)")
        _setattr!(cszhv, "units", "m")
        _setattr!(cszhv, "grid_mapping", "mapping")
        if is_radar
            _setattr!(cszhv, "azimuth_pixel_size", pixel_size_y)
            _setattr!(cszhv, "chip_size_coordinates", CHIP_SIZE_COORDINATES.radar)
        else
            _setattr!(cszhv, "y_pixel_size", pixel_size_y)
            _setattr!(cszhv, "chip_size_coordinates", CHIP_SIZE_COORDINATES.optical)
        end
        _put!(cszhv, _quantize_u16(chip_size_height))

        imv = _defvar_2d!(ds, "interp_mask", UInt8, UInt8(0), dims2, chunks2)
        _setattr!(imv, "standard_name", "interpolated_value_mask")
        _setattr!(imv, "description", "true where values have been interpolated")
        _setattr!(imv, "flag_values", UInt8[0, 1])
        _setattr!(imv, "flag_meanings", "measured interpolated")
        _setattr!(imv, "grid_mapping", "mapping")
        _put!(imv, UInt8.(interp_mask))
    end
    return path
end
