# autorift_itslive_production

One self-contained executable that runs the ITS_LIVE offset-production chain on a pair of scenes and writes
the product netCDF:

```
autorift_itslive_production [<sensor>] <reference> <secondary> [<output.nc>] [options]
```

The sensor, when omitted, is the reference's (`infer_sensor` in `src/production.jl`); the output, when
omitted, is the product's ITS_LIVE name (`<ref>_X_<sec>_G0120V02_P<valid>.nc`) in the working directory.
`--bucket` and `--bucket-prefix` upload the product to S3 when it is written, with its browse image and
thumbnail (`src/browse.jl`, pixel for pixel what `hyp3-autorift` draws); `--publish-bucket` publishes the three
to the ITS_LIVE Open Data bucket under `velocity_image_pair/<platform>/v02/<region>`.

The module is `AutoriftItsLiveProduction` (`src/`). It is built with `juliac --trim=safe`, so every value
whose type would depend on file contents is resolved by an explicit branch; `src/cli.jl`'s `USAGE` is the
full option list.

| sensor | inputs |
|---|---|
| `L4` `L5` `L7` `L8` `L9` | the scene's `_B2.TIF` (L4, L5) or `_B8.TIF` band with `_MTL.json` and `_ANG.txt` beside it, or its landsatlook STAC item |
| `S2` | the `B08.jp2` inside the product's `.SAFE`, or its Earth Search STAC item |
| `S1` | an IW SLC's `.SAFE` directory or `.zip` |
| `S1BURST` | the IW SLC `.zip`s the bursts are in, comma-separated, with `--subset`; the bursts' extractor URLs, comma-separated; or the SAFE `burst2safe` assembles from them |
| `NL1` `NL2` | a NISAR RSLC or GSLC `.h5` |

An input may also be the scene ID alone (`src/sceneid.jl`): `LC08_L1TP_021248_20210507_20210517_02_T1`,
`S2A_MSIL1C_20200624T181931_N0209_R127_T19XDJ_20200624T201221`,
`S1A_IW_SLC__1SDV_20211118T161558_20211118T161625_040628_04D1E5_D6DE`, a NISAR granule ID, or for `S1BURST`
burst IDs such as `S1_342239_IW3_20211202T161351_VV_3876-BURST`, comma-separated. Landsat IDs become their
landsatlook STAC item, Sentinel-2 IDs their `B08.jp2` on the Google Cloud mirror and Sentinel-1 SLC IDs their
zip at ASF. NISAR and burst IDs are looked up in ASF's catalog (CMR): a NISAR granule's directory carries its
collection's version, and a burst's record names the SLC that holds it and its index there, so the burst is
read out of that SLC's zip.

Each input is a local path, an `s3://` URI or an `https://` URL; scenes are read in place, a window at a
time. A Sentinel-2 band given as a URL is the exception: it is copied to the run's scratch directory first
(about 125 MB each, both at once), because a JPEG 2000 file read in place costs a request per tile-part header
and per tile. Landsat comes from requester-pays `s3://usgs-landsat`, so the AWS credential chain must name an account
that can pay. An Earthdata host (NISAR from ASF) is read with `EARTHDATA_TOKEN` or a token obtained with
`EARTHDATA_USERNAME` and `EARTHDATA_PASSWORD`, or else the `urs.earthdata.nasa.gov` entry of `~/.netrc`. Bursts come from ASF's burst extractor
(`https://sentinel1-burst.asf.alaska.edu/<SLC>/IW<n>/<POL>/<index>.tiff`, the URLs ASF's search lists),
authenticated with the same credentials, or out of their SLC zips (`src/s1subset.jl`: the bursts
both passes hold that a `--subset` polygon covers at least half of, matched by subswath and time after the
ascending node). Each burst's raster is downloaded, or inflated out of its zip up to the last burst needed,
into the run's scratch directory, and the product they make is the one `burst2safe` assembles from them. Sentinel-1 orbits are
fetched from the public `s1-orbits` bucket unless `--orbit-dir` holds them. A radar resample solves terrain
against the ITS_LIVE parameter region's DEM plus the EGM2008 geoid unless `--dem` gives a geographic DEM of
ellipsoidal heights.

## Building

```
app/build.sh [outdir]          # this machine's platform; default app/build
app/build_linux.sh <outdir>    # linux/x86_64, in an amd64 Julia container
```

`build.sh` uses JuliaC from the `@jcdrv` shared environment
(`julia +release --project=@jcdrv -e 'import Pkg; Pkg.add("JuliaC")'` once). The bundle under `outdir` is
relocatable as a whole.

`.github/workflows/Release.yml` builds the three bundles on GitHub's runners (`linux-x86_64`,
`linux-aarch64`, `macos-arm64`), each with `build.sh` on its own platform. A trimmed build with any
verifier error fails the job. Each bundle is packed with `tools/release/pack.sh`.

- A `v*` tag attaches the bundles, their `.sha256` and a `manifest.json` to that tag's GitHub Release, and
  pushes the container image (`Dockerfile`: the Linux bundles on `debian:bookworm-slim`, linux/amd64 and
  linux/arm64) to `ghcr.io/alex-s-gardner/itsliveoffsetproduction.jl` as `<tag>` and `latest`.
- Run workflow on GitHub builds any branch, keeping the bundles as workflow artifacts for 30 days. From
  `main` it also pushes the image as `<version>` (`git describe`) and `test`.
- The copy under `s3://its-live-data/production_code/image_pair/` is made by hand.

## Threads and FFTW wisdom

A `juliac` executable reads its thread count and AutoRIFT its FFTW wisdom directory only at start-up, so
the executable runs itself again when either is unset:

- with `JULIA_NUM_THREADS=auto`;
- with `AUTORIFT_FFTW_WISDOM_DIR` pointing at a temporary copy of the wisdom the bundle ships
  (`share/autorift_itslive_production/fftw-wisdom`, from `app/wisdom/`), deleted when the run ends.

Every run on a CPU model the bundle has wisdom for therefore plans its transforms alike: products are
reproducible run to run, and only bit-identical between runs that start from the same wisdom. Wisdom is
keyed by CPU model and FFTW version. `app/wisdom/` holds a file for Apple M2 Max, Graviton4 (`neoverse-v2`),
Graviton5 (`neoverse-v3`) and AMD EPYC 9R45; `RELEASE.md` says which runs grew each. On a CPU without a file
the run measures its plans and is slower. To give it one:

1. run the golden cases with `AUTORIFT_FFTW_WISDOM_DIR` set to one directory;
2. copy the file written there into `app/wisdom/`.

Set the variable yourself to keep and grow your own wisdom.

## Results

`tools/trim/README.md` records each phase's gate: products against the JIT chain and the Python reference,
wall time and peak memory per golden case.
