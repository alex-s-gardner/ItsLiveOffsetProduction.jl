# autorift_itslive_production

A self-contained executable that runs the ITS_LIVE offset-production chain on a pair of scenes and writes the
product netCDF. It needs no Julia installation. This file describes a release and how to run it; it ships at
the root of every bundle and beside the releases in S3.

## The release files

Releases live under `s3://its-live-data/production_code/image_pair/`:

```
latest.json                                  the current version: {"version": ..., "manifest": ...}
README.md                                    this file
versions/<version>/
    autorift_itslive_production-<version>-<platform>.tar.zst         one bundle per platform
    autorift_itslive_production-<version>-<platform>.tar.zst.sha256
    manifest.json                            commit, build date, sizes, checksums, notes
```

- `<version>` is `git describe` of the commit the bundles were built from (`v0.1.0-11-ge53221e`: 11 commits
  after the `v0.1.0` tag, at commit `e53221e`). A version is written once and never overwritten;
  `latest.json` is the only object that changes.
- `<platform>` is `linux-aarch64` (AWS Graviton), `linux-x86_64` (Intel and AMD) or `macos-arm64` (Apple
  silicon). `manifest.json` records where each was built and tested.
- The `linux-x86_64` bundle is built for any x86_64 CPU (`generic`, `sandybridge`, `haswell` and `x86-64-v4`
  code paths, chosen at run time). It has run on an AMD EPYC 9R45 only.

Fetch, verify and unpack one:

```bash
aws s3 cp s3://its-live-data/production_code/image_pair/latest.json -          # the current version
V=v0.1.0-11-ge53221e; P=linux-aarch64; N=autorift_itslive_production-$V-$P
B=s3://its-live-data/production_code/image_pair/versions/$V
aws s3 cp $B/$N.tar.zst . && aws s3 cp $B/$N.tar.zst.sha256 .
sha256sum -c $N.tar.zst.sha256            # macOS: shasum -a 256 -c
zstd -dc $N.tar.zst | tar -xf -
$N/bin/autorift_itslive_production --help
```

## What is in a bundle

```
autorift_itslive_production-<version>-<platform>/
    bin/autorift_itslive_production     the executable
    lib/                                the runtime and its libraries (GDAL, PROJ, HDF5, NetCDF, FFTW, ...)
    share/                              data the libraries read (GDAL and PROJ data, the CA bundle, ...)
    share/autorift_itslive_production/fftw-wisdom/    FFTW wisdom, one file per CPU model (below)
    README.md                           this file
```

The directory is relocatable as a whole: move or copy it anywhere and run `bin/autorift_itslive_production`.

## Running it

```
autorift_itslive_production <sensor> <reference> <secondary> <output.nc> [options]
```

| sensor | inputs (each a local path, an `s3://` URI or an `https://` URL) |
|---|---|
| `L4` `L5` `L7` `L8` `L9` | the scene's band (`_B2.TIF` for L4 and L5, `_B8.TIF` for L7, L8 and L9) with `_MTL.json` and `_ANG.txt` beside it, or its landsatlook STAC item |
| `S2` | the `B08.jp2` inside the product's `.SAFE`, or its Earth Search STAC item |
| `S1` | an IW SLC's `.SAFE` directory or `.zip` |
| `S1BURST` | the bursts' extractor URLs, comma-separated, or the SAFE `burst2safe` assembles from them |
| `NL1` `NL2` | a NISAR RSLC (`NL1`) or GSLC (`NL2`) `.h5` |

```bash
# Sentinel-2 from the public Google Cloud mirror: needs no credentials
B=https://storage.googleapis.com/gcp-public-data-sentinel-2/tiles/07/V/EG
$N/bin/autorift_itslive_production S2 \
  $B/S2A_MSIL1C_20200626T204021_N0209_R014_T07VEG_20200626T235742.SAFE/GRANULE/L1C_T07VEG_A026180_20200626T204741/IMG_DATA/T07VEG_20200626T204021_B08.jp2 \
  $B/S2B_MSIL1C_20200701T204019_N0209_R014_T07VEG_20200701T223901.SAFE/GRANULE/L1C_T07VEG_A017343_20200701T204020/IMG_DATA/T07VEG_20200701T204019_B08.jp2 \
  s2_pair.nc --scratch /mnt/nvme/scratch
```

Options: `--scratch DIR`, `--cache-io pread|mmap`, `--filter-tile N`, `--heap-hint-gib N` (default 12),
`--orbit-dir DIR` (S1), `--dem PATH` (S1, NL1), `--quiet`. `--help` lists them in full.

**Output and exit codes.** Each stage's time goes to stderr (`download`, `geometry`, `imagery`, `correlate`,
`package`), then `wrote <path>`. Exit `0` is success, `1` a failed run (the reason follows `error:` on
stderr), `2` a usage error.

**Credentials.** The run reads them from the environment:

- Landsat is read from requester-pays `s3://usgs-landsat`: the AWS credential chain (environment,
  `AWS_PROFILE`, or the instance role) must name an account that can pay. Run in `us-west-2` to avoid transfer
  charges.
- NISAR granules, Sentinel-1 zips and Sentinel-1 bursts are read from ASF with Earthdata Login: set
  `EARTHDATA_TOKEN`, or put the `urs.earthdata.nasa.gov` entry in `~/.netrc`.
- Sentinel-2 (Google mirror) and Sentinel-1 orbits (the public `s1-orbits` bucket) need none.

**Threads.** If `JULIA_NUM_THREADS` is unset the executable runs itself again with `auto`, which is every core
the process may use. The measurements below use 4.

**Scratch.** `--scratch` must be on disk. The run refuses a RAM-backed `tmpfs` or `ramfs`, which includes `/tmp`
on Amazon Linux 2023. The run creates a directory of its own under it and removes it when it ends. A run killed
by a signal cannot remove its directory, so the next run on the same scratch does. A Sentinel-2 band given as
a URL is copied there first (about 125 MB each); the tile caches for a NISAR pair reach tens of GB, so use
local NVMe for the radar sensors.

**Resources.** Measured on a Graviton5 instance (m9gd.xlarge, 4 vCPU, 16 GiB, local NVMe) with 4 threads,
pairs drawn at random from the benchmark lists. Peak memory is the maximum resident anonymous memory.

| sensor | peak memory | peak scratch | time |
|---|---:|---:|---|
| `L4` `L5` `L7` `L8` `L9` | 5.2 GiB | 5 GB | 25 to 100 s |
| `S2` | 2.1 GiB | 2.5 GB | about 25 s |
| `S1` | 9.4 GiB | 30 GB | 2 to 3 min |
| `S1BURST` | 5.8 GiB | 18 GB | 1.5 to 4 min |
| `NL1` | 7.2 GiB (10 GiB for a large pair) | 67 GB | 1 to 17 min |
| `NL2` | 8.3 GiB (one pair) | about 50 GB* | about 25 min |

\* From an earlier run at 10 threads; the other scratch figures are the 4-thread maxima.

A 16 GiB instance with 100 GB of local scratch runs every sensor. Landsat and Sentinel-2 fit in 8 GiB. The
NISAR figures are not measured above 4 threads.

**NISAR frequency and polarization.** A NISAR product is read at frequency A, at HH if the reference lists it
and otherwise at VV. A granule with no frequency A (mode `NASV`, which acquires frequency B alone) is refused
with `<file> opens at frequency B, not A`.

**Failures that are not bugs.** A pair over no ice fails with `every search radius is zero`: the ITS_LIVE
parameter raster has no search range there. A Sentinel-1 pair whose subswaths differ in burst count is
refused. The ASF burst extractor sometimes answers a request with an HTTP or TLS error; run the pair again.

## FFTW wisdom

FFTW plans its transforms by measuring, and keeps what it learned as "wisdom". A run that starts from the
same wisdom plans every transform the same way, so its products are bit-identical to the last run's; a run that
has to measure may plan some transforms differently, and about 0.01 to 0.1% of the pixels then differ. The
executable therefore starts every run from the wisdom the bundle carries, in
`share/autorift_itslive_production/fftw-wisdom/`.

A file is named `<cpu>-fftw<version>.wisdom`, where `<cpu>` is the CPU model the OS reports (on Graviton, the
core name) and the version is FFTW's. A run reads only the file with its own CPU's name. It works on a
temporary copy, which is deleted when the run ends, so nothing a run measures is kept and runs do not change
the bundle.

| file | CPU | grown by |
|---|---|---|
| `neoverse-v3-fftw1.10.0.wisdom` | Graviton5 (tested on m9g, m9gd) | Landsat 4/5, 7 and 8/9, S2, S1, S1BURST, two NL1 pairs, NL2 |
| `neoverse-v2-fftw1.10.0.wisdom` | Graviton4 (tested on m8g) | Landsat 4/5, 7 and 8/9, S2, S1 |
| `AMD_EPYC_9R45-fftw1.10.0.wisdom` | AMD EPYC 9R45 (tested on m8a) | Landsat 4/5, 7 and 8/9, S2, S1, S1BURST |
| `Apple_M2_Max-fftw1.10.0.wisdom` | Apple M2 Max | the union of the wisdom the golden-case gates measured |

"Grown by" lists the runs that added plans to the file. It does not mean every pair of that sensor is
covered: plan sizes depend on a pair's chip and search sizes, and a pair the file has not seen measures its
missing plans, which can add tens of seconds. NISAR is covered on Graviton5 only.

**On any other CPU** (Graviton2 and 3, Intel, other AMD models) the bundle has no file. The run measures every
plan it needs: it is slower, and its products are only bit-reproducible against other runs that started from
the same wisdom, not against runs on another CPU. To give such a CPU a file:

1. On an instance of that type, copy the bundle's wisdom directory to a persistent one
   (`cp -r $N/share/autorift_itslive_production/fftw-wisdom ~/wisdom_grow`; make an empty one if you like) and
   run representative pairs, one after another, with `AUTORIFT_FFTW_WISDOM_DIR=~/wisdom_grow` set. A run with
   that variable set reads the directory and adds what it measures to it.
2. Copy the file written there (named for the CPU) into the bundle's `fftw-wisdom` directory on every instance
   of that type, or into a directory every run's `AUTORIFT_FFTW_WISDOM_DIR` points at. To make it part of the
   next release, add it to `app/wisdom/` in the repository.

Another CPU's file can be renamed to the new CPU's name to avoid measuring, but a plan tuned on one core's cache
hierarchy is a poor plan on another; treat that as a stopgap.

## Checking a bundle

`tools/release/smoke.sh <bundle>` in the repository runs one fast pair of each sensor through a bundle and
checks that each reaches its geometry stage (about a minute and a half in all). Run it on an instance of the
type you will deploy to, with the credentials above in place.
