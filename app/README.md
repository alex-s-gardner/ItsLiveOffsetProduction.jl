# autorift_itslive_production

One self-contained executable that runs the ITS_LIVE offset-production chain on a pair of scenes and writes
the product netCDF:

```
autorift_itslive_production <sensor> <reference> <secondary> <output.nc> [options]
```

The module is `AutoriftItsLiveProduction` (`src/`). It is built with `juliac --trim=safe`, so every value
whose type would depend on file contents is resolved by an explicit branch; `src/cli.jl`'s `USAGE` is the
full option list.

| sensor | inputs |
|---|---|
| `L4` `L5` `L7` `L8` `L9` | the scene's `_B2.TIF` (L4, L5) or `_B8.TIF` band with `_MTL.json` and `_ANG.txt` beside it, or its landsatlook STAC item |
| `S2` | the `B08.jp2` inside the product's `.SAFE`, or its Earth Search STAC item |
| `S1` | an IW SLC's `.SAFE` directory or `.zip` |
| `S1BURST` | the bursts' extractor URLs, comma-separated, or the SAFE `burst2safe` assembles from them |
| `NL1` `NL2` | a NISAR RSLC or GSLC `.h5` |

Each input is a local path, an `s3://` URI or an `https://` URL; scenes are read in place, a window at a
time. Landsat comes from requester-pays `s3://usgs-landsat`, so the AWS credential chain must name an account
that can pay. An Earthdata host (NISAR from ASF) is read with `EARTHDATA_TOKEN` or a token obtained with the
`urs.earthdata.nasa.gov` entry of `~/.netrc`. Bursts come from ASF's burst extractor
(`https://sentinel1-burst.asf.alaska.edu/<SLC>/IW<n>/<POL>/<index>.tiff`, the URLs ASF's search lists),
authenticated with the same `~/.netrc` entry; each burst's raster is downloaded into the run's scratch
directory, and the product they make is the one `burst2safe` assembles from them. Sentinel-1 orbits are
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

## Threads and FFTW wisdom

A `juliac` executable reads its thread count and AutoRIFT its FFTW wisdom directory only at start-up, so
the executable runs itself again when either is unset:

- with `JULIA_NUM_THREADS=auto`;
- with `AUTORIFT_FFTW_WISDOM_DIR` pointing at a temporary copy of the wisdom the bundle ships
  (`share/autorift_itslive_production/fftw-wisdom`, from `app/wisdom/`), deleted when the run ends.

Every run on a CPU model the bundle has wisdom for therefore plans its transforms alike: products are
reproducible run to run, and only bit-identical between runs that start from the same wisdom. Wisdom is
keyed by CPU model and FFTW version. On a CPU without a file the run measures its plans and is slower.
To give it one:

1. run the golden cases with `AUTORIFT_FFTW_WISDOM_DIR` set to one directory;
2. copy the file written there into `app/wisdom/`.

Set the variable yourself to keep and grow your own wisdom.

## Results

`tools/trim/README.md` records each phase's gate: products against the JIT chain and the Python reference,
wall time and peak memory per golden case.
