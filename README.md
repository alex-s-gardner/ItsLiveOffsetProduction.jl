# ItsLiveOffsetProduction

ITS_LIVE offset products: the inputs a granule pair yields, and the netCDF they are written to.

[AutoRIFT.jl](https://github.com/alex-s-gardner/AutoRIFT.jl) correlates an image pair and returns
displacements. Turning those into an ITS_LIVE product is a separate job with its own conventions —
plugin and software version strings, scene naming, swath offset bias, the stable-shift correction, and
three netCDF schemas — and this package owns all of it.

`write_product` computes velocity, stable-shift correction, error estimates, and (for a radar pair) the
conversion matrix in a single pass, from one `ItsLiveInput`:

```julia
using ItsLiveOffsetProduction

input = ItsLiveInput(:optical, "feature", "map", autorift_version, parameter_file, source,
                     coefficients(geometry), reference_velocity(inputs)..., georef, info, dx, dy,
                     chip_size_x, interp_mask, scale_chip_size_y, dt_seconds, swath_bias)
write_product("S1A_..._X_..._G0120V02_P019.nc", input)
```

The reference implementation (`hyp3-autorift`'s `testautoRIFT.py`, `netcdf_output.py` and `crop.py`)
does this across three files and two writes — an uncropped file, then a full reopen that crops, pads to
a 512-pixel/120 m-aligned grid, and adds a `time` dimension. This does all of it in one pass.

## Building an `ItsLiveInput`

Each field of `ItsLiveInput` names an external input `write_product` does not compute; its docstring
says which package or caller responsibility supplies each one. The structs are `GeogridCoefficients`,
`ReferenceVelocity`, `ItsLiveGeoref`, `ImagePairInfo` and `SwathOffsetBias`.

Three helpers convert `ImagePairGeometry`'s own geometry into those structs. They are not exported, since
the names are generic enough to collide:

```julia
coeffs      = ItsLiveOffsetProduction.coefficients(g::PairGeometry)
refv, mask  = ItsLiveOffsetProduction.reference_velocity(inputs::GeometryInputs)
locx, locy  = ItsLiveOffsetProduction.image_location(g::PairGeometry)
```

`cf_grid_mapping(crs)` returns the CF grid-mapping attributes for a CRS, and `ItsLiveGeoref` has a
constructor taking a `Rasters.Raster` directly.

Note that `AutoRIFT.chip_size_scale` stays in AutoRIFT.jl: it configures the correlator's chip-size
pyramid, and `AutoRIFT.params` depends on it.

## Dependencies

Every dependency is a hard dependency. The netCDF writer, the geogrid repackaging and the CF grid
mapping are all needed on every run, so there is nothing for a package extension to defer.

`FastGeoProjections` is declared although no code here uses it: `ImagePairGeometry` requires a version
newer than the registered one, and Pkg resolves a `[sources]` entry only for a package also in `[deps]`.

## Testing

`Pkg.test()` covers the netCDF schemas, the metadata arithmetic and the radar corrections — everything
that runs offline. The granule-to-product pipeline and the comparison against the reference
implementation live in `tools/golden/` and are run by hand, since they need AWS credentials and
requester-pays egress.
