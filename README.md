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

## The `autorift_itslive` executable

`app/` builds this package and the rest of the chain (scene readers, geometry, AutoRIFT.jl) into one
self-contained executable that runs an image pair from its sources to the ITS_LIVE product netCDF.
[`app/RELEASE.md`](app/RELEASE.md) is the user's guide shipped in every bundle: inputs per sensor, credentials,
examples, and how to get a release. [`app/README.md`](app/README.md) covers building it.

```bash
# The product is written to the working directory under its ITS_LIVE name
autorift_itslive LC08_L1TP_021248_20210507_20210517_02_T1 LC08_L1TP_021248_20220611_20220617_02_T1

# A burst pair, its bursts listed as HyP3 lists them, uploaded and published
autorift_itslive \
  --reference S1_342239_IW3_20211202T161351_VV_3876-BURST S1_342240_IW3_20211202T161354_VV_3876-BURST \
  --secondary S1_342239_IW3_20211214T161351_VV_91C4-BURST S1_342240_IW3_20211214T161353_VV_91C4-BURST \
  --bucket <hyp3-bucket> --bucket-prefix <job-id> --publish-bucket its-live-data
```

Scenes may be given by scene ID, path or URL; the sensor is read from the reference. The bundles for
`linux-x86_64`, `linux-aarch64` and `macos-arm64` are attached to each tagged
[GitHub Release](https://github.com/alex-s-gardner/ItsLiveOffsetProduction.jl/releases).

### Docker

The Linux bundles also ship as the container image `ghcr.io/alex-s-gardner/itsliveoffsetproduction.jl`
(linux/amd64 and linux/arm64), its entry point the executable:

| tag | built from |
|---|---|
| `<tag>`, `latest` | a `v*` tag of this repository |
| `<version>`, `test` | `main`, when the Release workflow is run on it (Actions → Release → Run workflow); `<version>` is `git describe` |

```bash
docker run --rm -v "$PWD:/home/autorift/out" -w /home/autorift/out \
  -e AWS_ACCESS_KEY_ID -e AWS_SECRET_ACCESS_KEY -e AWS_SESSION_TOKEN \
  -e EARTHDATA_USERNAME -e EARTHDATA_PASSWORD -e PUBLISH_ACCESS_KEY_ID -e PUBLISH_SECRET_ACCESS_KEY \
  ghcr.io/alex-s-gardner/itsliveoffsetproduction.jl:latest \
  --reference <scene>... --secondary <scene>... --scratch /home/autorift/scratch \
  --bucket <bucket> --bucket-prefix <prefix> --publish-bucket <open-data-bucket>
```

The image runs as user 1000 in `/home/autorift`, and the product and its browse images are written to the
working directory. In AWS Batch, as HyP3 runs it, the job's arguments are the executable's: `--reference` and
`--secondary` take the space-separated lists HyP3 passes, and `--publish-bucket None` publishes nothing. The
secrets are the environment variables above, and the AWS credentials come from the job role. Give
`--scratch` a directory on disk with room for the pair (`app/RELEASE.md` lists the space each sensor needs).

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
