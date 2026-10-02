using Dates: DateTime

"""
    GeogridCoefficients

Per-pixel displacement-to-velocity conversion coefficients that Geogrid derives from the scene
geometry — not something `AutoRIFT.jl` computes. `scale_factor_1`/`scale_factor_2` are required for
every `pair_type`; only `offset2vr`/`offset2va` are radar-only (`nothing` for `pair_type ==
:optical`) — `scale_factor_1` additionally appears in the radar-only `M11`/`M12` formula.

Every field shares the shape of `ItsLiveInput.dx`.
"""
struct GeogridCoefficients
    offset2vx_1::Matrix{Float64}
    offset2vx_2::Matrix{Float64}
    offset2vy_1::Matrix{Float64}
    offset2vy_2::Matrix{Float64}
    scale_factor_1::Matrix{Float64}
    scale_factor_2::Matrix{Float64}
    offset2vr::Union{Matrix{Float64},Nothing}
    offset2va::Union{Matrix{Float64},Nothing}
end

"""
    coefficients(g) -> GeogridCoefficients

`g`'s displacement-to-velocity operator, scale factors, and (radar-only) direct axis-velocity band,
repackaged as a [`GeogridCoefficients`](@ref). Defined when `ImagePairGeometry` is loaded, for
`g::ImagePairGeometry.PairGeometry` — that package is what derives these from scene geometry; nothing
here could. See `geogrid.jl`.
"""
function coefficients end

"""
    image_location(g) -> (location_x, location_y)

`g`'s per-point image position — the range/azimuth (or x/y) pixel index each output grid point falls
at in the acquisition. Defined when `ImagePairGeometry` is loaded, for
`g::ImagePairGeometry.PairGeometry`. Feeds [`SwathOffsetBias`](@ref); nothing else in this package
needs it. See `geogrid.jl`.
"""
function image_location end

"""
    ReferenceVelocity

An external reference velocity field (e.g. a prior ITS_LIVE mosaic), used only to determine and
remove a stable-surface offset bias and to report the error attributes that bias correction implies.
Not an `AutoRIFT.jl` output.

`Float32`, matching the reference velocity raster's own on-disk type — confirmed against a real
capture, not assumed. It matters beyond precision: the `error`/`stable_shift` attribute family is
`std`/`median` of `vx .- reference.vx`, and the reference computes that difference in `Float32`
throughout; widening `reference` to `Float64` here would silently change those attributes' on-disk
type (and, since `Float32` and `Float64` rounding of the same quantity can round to different last
digits, occasionally their value) relative to what the reference actually writes.
"""
struct ReferenceVelocity
    vx::Matrix{Float32}
    vy::Matrix{Float32}
end

"""
    reference_velocity(inputs) -> (ReferenceVelocity, stable_mask::Matrix{Bool})

`inputs`' velocity and stable-surface bands, repackaged as a [`ReferenceVelocity`](@ref) and a
stable-surface mask. Defined when `ImagePairGeometry` is loaded, for
`inputs::ImagePairGeometry.GeometryInputs` — the same rasters `ImagePairGeometry.geometry_inputs`
fetches for the Geogrid computation itself, so building `ItsLiveInput.reference`/`stable_mask` this
way costs no second fetch. See `geogrid.jl`.
"""
function reference_velocity end

"""
    ItsLiveGeoref

The output grid and its CF grid-mapping metadata.

`x`/`y` are the full-resolution, uncropped, cell-center coordinates of every 2-D field in an
[`ItsLiveInput`](@ref); `x` increasing, `y` decreasing (north-up), matching the reference product's own
convention. `mapping_attrs` carries every CF grid-mapping attribute except `GeoTransform`, which
[`write_product`](@ref) always recomputes from whichever grid is actually written — the full grid when
[`ImagePairInfo`](@ref)'s `roi_valid_percentage` rounds down to `0`, the cropped-and-aligned grid
otherwise (this is the reference's own criterion, not a fresh check of whether any pixel is valid —
see `ItsLiveInput`). `lonlat` is called once, on the cropped grid's centroid, to fill in
`latitude`/`longitude` (skipped, along with cropping itself, in the uncropped case). Neither the
projection parameters nor the coordinate transform are this package's concern — both are supplied by
the caller, directly or via the `Rasters.Raster`-based constructor in `georef.jl` (see
[`cf_grid_mapping`](@ref)).

`pixel_size_x`/`pixel_size_y` are the scene's own pixel size in metres (range/azimuth for radar,
x/y for optical). For `pair_type == :radar`, [`write_product`](@ref) **ignores** `pixel_size_y` and recomputes
it from the median of `offset2va`, matching the reference exactly — a documented reference quirk, not
a choice made here.

`x`/`y` must already be spaced at exactly 120 m — the fixed ITS_LIVE product grid every schema uses,
not a resolution [`write_product`](@ref) resamples to. Cropping only crops and pads this existing grid to a
512-pixel/120 m-aligned tile; a caller whose native grid spacing differs must reproject/resample onto
a 120 m grid before calling [`write_product`](@ref). `write` checks the spacing and throws rather than
silently misaligning the output.
"""
struct ItsLiveGeoref
    x::Vector{Float64}
    y::Vector{Float64}
    mapping_attrs::Dict{String,Any}
    lonlat::Function
    pixel_size_x::Float64
    pixel_size_y::Float64
end

"""
    cf_grid_mapping(crs) -> Dict{String,Any}

The CF grid-mapping attributes GDAL's netCDF writer records for `crs` — everything
[`ItsLiveGeoref`](@ref)'s `mapping_attrs` needs except `GeoTransform`, which [`write_product`](@ref) always
recomputes itself. Defined when `Rasters`/`ArchGDAL` are loaded, for `crs` anything
`GeoFormatTypes.convert`ible to `WellKnownText` — a `Rasters.crs(raster)` result, a
`GeoFormatTypes.GeoFormat`, or an integer EPSG code. See `georef.jl`.

Supports only the two projections an ITS_LIVE product grid uses, polar stereographic and UTM
(transverse Mercator); any other throws.
"""
function cf_grid_mapping end

"""
    ImagePairInfo

The per-image-pair metadata the reference calls `IMG_INFO_DICT`.

`mission_img1`/`mission_img2` (`"L"`, `"S"`, or `"N"`) and `satellite_img1`/`satellite_img2` feed the
`satellite` global attribute. `satellite_img1`/`satellite_img2` are untyped because the reference's own
`IMG_INFO_DICT` carries them as whatever type the source metadata gives — a string for Landsat/Sentinel
(`"8"`, `"A"`) but a bare integer for NISAR — and both the type and the value are written through
verbatim to the per-variable `img_pair_info` attribute; only the derived `satellite` global attribute
stringifies it. `acquisition_date_img1`/`acquisition_date_img2` feed the computed
`date_dt`/`date_center` attributes. `latitude`/`longitude` are the pre-crop centroid — [`write_product`](@ref)
overwrites both after cropping, for every case except the uncropped (`P000`) one. `roi_valid_percentage`
does double duty beyond being a plain attribute: `write` skips cropping entirely whenever it rounds
down to `0`, matching the reference's own criterion — its filename ends in `_P000.nc` under the same
condition, and `process.py` decides whether to crop from that filename, not from a fresh check of
whether any pixel is valid. The two read differently only at the boundary, but they do: a real captured
case with 19,994 valid pixels (0.39% of the grid) has `roi_valid_percentage = 0.9`, floors to `0`, and
the reference leaves it uncropped despite the nonzero count. `extra` carries every other `img_pair_info`
attribute (`id_img1`, `sensor_img1`, orbit numbers, path/row, and so on) through verbatim; sourcing
those from raw sensor metadata (SAFE XML, MTL, HDF5, STAC) is outside this
package.

`DateTime` is millisecond-precision, one thousand times coarser than the reference's microsecond
timestamps; `date_center` therefore loses precision `write` cannot recover. `time` and `date_created`
are excluded from the golden-test comparison for an unrelated reason (an unseeded hash), but
`date_center` is not, so this is a real, if usually negligible, source of mismatch.
"""
struct ImagePairInfo
    acquisition_date_img1::DateTime
    acquisition_date_img2::DateTime
    mission_img1::String
    mission_img2::String
    satellite_img1
    satellite_img2
    roi_valid_percentage::Float64
    latitude::Float64
    longitude::Float64
    extra::Dict{String,Any}
end

"""
    roi_valid_percentage(chip_size_x, search_limit_x) -> Float64

The percentage of searched grid points the correlator resolved, which is [`ImagePairInfo`](@ref)'s
`roi_valid_percentage`.

The denominator counts every point with a nonzero search limit — the points a search was attempted at,
which is what "ROI" names here — and the numerator every point the correlator returned a chip size for.
Both arrays cover the whole grid before cropping, and both must already carry the no-data exclusions:
a point no search was attempted at belongs to neither count. `search_limit_x` is therefore the limit
the correlator was given with the dilated no-data mask applied over it, not the geogrid's raw band.

The ratio is quantized to three decimals before being scaled, so the result is a whole number of tenths
of a percent. It carries more than the attribute: a value that floors to `0` selects the uncropped
schema, which [`ImagePairInfo`](@ref) describes.

Throws when no point was searched, since a ratio over an empty ROI has no value to report.
"""
function roi_valid_percentage end

"""
    SwathOffsetBias

What the Sentinel-1 subswath-offset-bias correction (`cal_swath_offset_bias` in the reference) needs
beyond `dx`/`dy` themselves: `nothing` for an optical pair, for a same-platform S1 pair (the bias is a
difference between two spacecrafts' antenna patterns, so it has no meaning within one), or for a pair
whose product was not opened for all three subswaths.

`location_x`/`location_y` are the per-output-point image position `dx`/`dy` were measured at — see
[`image_location`](@ref) — needed because the correction masks by which subswath a point's *range*
position falls in, not by anything in `GeogridCoefficients`. `border12`/`border23`/`ncols` and
`same_platform`/`reference_platform` come from `SLCDatasets.subswath_borders`, unchanged; assembling
the two from a coregistered Sentinel-1 pair is that package's job, not this one's.

`grid_spacing_x` is the output grid's own spacing, in *original-image* pixels (`Params.grid_spacing.X`
— `AutoRIFT.params(g).grid_spacing.X`), used only by the correction's own coarse-grid smoothing step.
"""
struct SwathOffsetBias
    location_x::Matrix{Float64}
    location_y::Matrix{Float64}
    grid_spacing_x::Float64
    ncols::Int
    border12::Float64
    border23::Float64
    same_platform::Bool
    reference_platform::String
end

"""
    ItsLiveInput

Everything [`write_product`](@ref) needs to compute and package an ITS_LIVE product netCDF from raw pixel
offsets.

`dx`/`dy` are pixel offsets in the [`MultichipResult`](@ref) convention (secondary→reference, `dy`
row-positive); `chip_size_x` is `MultichipResult.chip_size`; `interp_mask` is
`MultichipResult.interpolated`. `AutoRIFT.jl` does not compute velocity, stable-shift correction, or
error estimates on its own — [`write_product`](@ref) does, from these fields plus the geogrid-derived
[`GeogridCoefficients`](@ref), the [`ReferenceVelocity`](@ref), and `stable_mask` (an external
stationary/slow-surface mask; the complementary "slowest 25%" mask is *not* a caller input — it is
derived from `reference` internally, matching the reference exactly).

`pair_type` selects the schema: `:radar` requires every radar-only field in `coefficients` to be given
(not `nothing`), and adds `vr`/`va`/`M11`/`M12` to the written product; `:optical` requires them to be
`nothing`. `dt_seconds` is the acquisition interval in seconds; it is used only for `pair_type ==
:radar`'s pixel-size recomputation (see [`ItsLiveGeoref`](@ref)) and has no other effect — the
reference itself passes `None` for it on an optical pair (confirmed against a real captured call, not
assumed), so it must be `nothing` there and non-`nothing` for `:radar`.

`swath_bias`, when not `nothing`, applies the Sentinel-1 subswath-offset-bias correction to `dx`/`dy`
before anything else — matching the reference's own order, where it runs on the raw displacement,
before stable-shift correction ever sees it (`testautoRIFT.py`'s call to `cal_swath_offset_bias`
precedes its stable-shift block).

Not included, and not computed here:

  - **Filename.** `write_product`'s first argument is the destination path.
  - **The `error_vector` dt-error model.** The reference hardcodes it (`[25.5, 25.5]` for optical; a
    fixed 2×6 table for radar, of which only the four columns feeding `vx`/`vy`/`vr`/`va` are ever
    read — see the comments in `write.jl`) rather than deriving it from
    anything in a run — `write_product` hardcodes the same constants.
"""
struct ItsLiveInput
    pair_type::Symbol
    detection_method::String
    coordinates::String
    autorift_software_version::String
    parameter_file::String
    source::String
    dt_seconds::Union{Float64,Nothing}

    georef::ItsLiveGeoref
    img_pair_info::ImagePairInfo
    coefficients::GeogridCoefficients
    reference::ReferenceVelocity
    stable_mask::Matrix{Bool}
    swath_bias::Union{SwathOffsetBias,Nothing}

    dx::Matrix{Float32}
    dy::Matrix{Float32}
    chip_size_x::Matrix{UInt16}
    scale_chip_size_y::Float64
    interp_mask::Matrix{Bool}
end

"""
    write_product(path, input::ItsLiveInput)

Compute velocity, stable-shift correction, error estimates, and the radar conversion matrix from
`input`, and write the result to `path` as an ITS_LIVE product netCDF, matching
`tools/golden/product.jl`'s reader. Three schemas: the 12-variable optical one, the 16-variable radar
one (`input.pair_type == :radar`), or the 11-variable uncropped one. `ImagePairInfo` says what selects
that last one.

The reference implementation (`hyp3-autorift`'s `testautoRIFT.py` + `netcdf_output.py` + `crop.py`)
computes this across three files and writes it in two passes — an uncropped file, then a full reopen
that crops to the valid extent, pads to a 512-pixel/120 m-aligned grid, and adds a `time` dimension.
This does all of it in one pass, from one struct.

"""
function write_product end
