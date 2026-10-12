# Trimmed-binary tooling and survey results

Tools for building `autorift_itslive_production` with `juliac --trim=safe` (Julia 1.13, JuliaC 0.3.9
from the `@jcdrv` shared environment).

- `probe.sh [--project ENV] [--threads N] <program.jl|package-dir> [args...]` — builds a trimmed,
  bundled executable, prints the verifier's errors grouped by the package holding each unresolved call
  (`verifier_report.jl`), and runs the result with `JULIA_DEPOT_PATH` pointing at an empty directory.
- `audit_closure.jl <project-dir>` — fails if the environment's dependency closure holds a package
  whose presence alone breaks a trimmed build (every package in the closure has its `__init__`
  compiled and run).

## What `--trim=safe` accepts on Julia 1.13.0

| construct | result |
|---|---|
| `Threads.@spawn` / `fetch` (with a type assertion), AutoRIFT `threaded = True()` | works; 4 threads bit-identical to 1 |
| Symbol keywords, `kw...` forwarding | works |
| `JSON3.read(s, T)` with StructTypes | works |
| `print(stdout, …)`, `println(stderr, …)` | fails: the globals are typed `IO`; use `Core.stdout` / `Core.stderr` |
| `sprint(showerror, err)` on an untyped exception | fails |
| FFTW.jl plan objects | fails (finalizer over `Vector{FFTWPlan}`); AutoRIFT's raw `Ptr` plans work |
| `repr(::Symbol)` in a message | fails (reaches the expression printer) |
| lazy `JSON3.read(s)[:k]` | fails |
| Downloads | fails: Base Timer/task error printers, `easy_hook` invokelatest |
| ArchGDAL, NCDatasets, HDF5.jl, AWS/AWSS3 | fail; replaced by GDAL.jl's C API, NetCDF_jll, HDF5_jll, own SigV4 |

### Loading versus calling

A package in the binary's closure has its `__init__` compiled and run. Built as a program that only
`import`s the package (`--trim=safe`, Julia 1.13), these load cleanly: **HTTP, Downloads, JSON, AWS,
ArchGDAL, NCDatasets**. These break the build by loading alone: **HDF5** (4 errors, its `@require`
blocks), **TiffImages** and **Rasters** (9 each, in Base's task and IO error printers). So only the
latter must leave a trimmed binary's closure (`audit_closure.jl`); for the rest, it is the *calls* the
binary makes that the verifier judges.

### Traps the verifier does not report

Each of these builds with 0 verifier errors and then fails when the binary runs:

- **A JLL library path captured in a `const`** (`const lib = HDF5_jll.libhdf5`): the binding is taken
  while the image is built, before the JLL's `__init__` sets it, so the binary calls into `""`. Name
  the JLL's global in every `ccall` — `ccall((:H5Fopen, HDF5_jll.libhdf5), …)` — so it is read at run
  time, where it points into the bundle. ChunkManifests' OpenJpeg extension and AutoRIFT's FFTW plans
  (`const LIBFFTW3F = libfftw3f_path`) both capture one.
- **`cglobal((:sym, lib))` with a library chosen at run time** raises a `MethodError`. Look the symbol
  up through the JLL's handle instead: `Libdl.dlsym(HDF5_jll.libhdf5_handle, :H5T_NATIVE_FLOAT_g)`.
- **Hiding the depot too thoroughly.** `JULIA_DEPOT_PATH=<empty dir>` also removes the bundle's own
  `share/julia`, where `--bundle` puts the artifacts; `JULIA_DEPOT_PATH=<empty dir>:` (trailing colon)
  hides only the user depot, which is what a deployed binary sees.

And two that the verifier does report, worth knowing because they look harmless: a variadic
`say(xs...)` called with mixed argument types, and `ErrorException.msg`, which is typed `AbstractString`
and must be narrowed (`m isa String`) before use.

## Library probes (`probes/`, macOS arm64)

Each probe is checked under the JIT first, then built with `probe.sh` and run with the user depot hidden.

| probe | what it exercises | verifier | run |
|---|---|---|---|
| `gdal.jl` | RasterIO window of a `/vsicurl/` prior COG; OGR point-in-polygon on the remote parameter shapefile; OSR WKT and semi-major axis; MEM→MEM `GDALWarp` | 0 | output identical to JIT; `GDAL_DATA`, `PROJ_LIB` resolve inside the bundle |
| `netcdf.jl` | netCDF-4: unlimited dim, deflate+shuffle, chunking, fill value, NC_CHAR/NC_STRING/Int64/double attributes, Int16 data | 0 | `ncdump` identical to the JIT-written file |
| `http.jl` | libcurl easy handles, four concurrent 64 KiB ranges from four threads; Earthdata redirect chain via `~/.netrc` + in-memory cookies; bundled `cert.pem` | 0 | bytes identical to `curl` |
| `s3.jl` | SigV4-signed ranged GET, credentials from env or `~/.aws/credentials`; requester-pays `usgs-landsat` and `its-live-data` | 0 | bytes identical to `aws s3api get-object` |
| `hdf5.jl` | HDF5_jll ccalls: fixed-length string, `H5Literate2`, filter pipeline, `H5Dchunk_iter`, compound complex hyperslab on a NISAR GSLC | 0 | identical to HDF5.jl |
| `stac.jl` | typed JSON3 + StructTypes parse of a landsatlook and an Earth Search item | 0 | identical to JIT |
| `xml.jl` | EzXML: S1 annotation XPath queries, burst list, orbit-file state vectors | 0 | identical to JIT |
| `zipindex.jl` | typed ZIP64 central directory; seekable-deflate index (Zlib_jll `Z_BLOCK`, `inflatePrime`, `inflateSetDictionary`) on a 1.47 GB S1 measurement member | 0 | windows identical to `unzip -p`; index 347 points / 11 MB built in 5.6 s; a 4 KiB window in ≤ 16 ms |
| `jp2.jl` | libopenjp2 decode of a 1024² area of an S2 B08 JP2 | 0 | identical to GDAL |
| `zarr.jl` | Zarr fork `readblock!` (`@async` producer, `Channel`, `@spawn` decoders) on a typed `ZArray`, 4 threads | 0 after the two fixes below (8 before) | identical to JIT |

**linux/amd64.** The same probes, built and run in the `julia:1.13` (1.13.1) image under Docker
Desktop's Rosetta emulation, give the same verifier counts and the same output, byte for byte; the
HTTP probe there takes its CA bundle from `/etc/ssl/certs`. Emulated Julia segfaults at a GC safepoint
when it runs multithreaded (Pkg's parallel artifact downloads crashed it), so the container runs with
`JULIA_NUM_THREADS=1 JULIA_PKG_CONCURRENT_DOWNLOADS=1`; that is enough for verifier counts and
correctness, not for timing or the threaded paths, which need a native x86_64 host. Data directories
must be mounted at a short container path (`-v …/golden_tests:/g`): mounted at their macOS path, they
appeared empty.

The Zarr probe's 8 errors before the fixes were two fork defects, not the task machinery: the scratch-buffer pool is keyed by
`Tuple{DataType,Tuple}` (`src/scratch.jl`), and `uncompress_raw!` passes a `Union{Nothing,T}`
`fill_value` as a keyword (`src/ZArray.jl:403`).

## Surveys

**Parameter regions** (`autorift_landice_0120m.shp`): 122 regions in three CRS families — UTM north
`326zz` and south `327zz` (all 120 zones), `NPS` EPSG:3413 and `SPS` EPSG:3031. Typed projections
therefore need UTM(zone, hemisphere), the two polar stereographics, and geographic.

**Encodings.**
- Landsat Collection 2 L1 bands: COG, DEFLATE, no predictor, 256² tiles, UInt16, nodata 0.
- ITS_LIVE prior COGs (`*_h.tif`, `*_vx.tif`, …): COG, LZW, 256² tiles; read through GDAL, not
  ChunkManifests.
- Sentinel-2 L1C: JPEG 2000.

**Landsat metadata from `_MTL.json`.** The product needs the time of day the product id lacks; STAC
`properties.datetime` (`2020-07-03T15:00:24.531814Z`) and MTL `DATE_ACQUIRED` + `SCENE_CENTER_TIME`
(`15:00:24.5318140Z`) are the same instant, and both truncate to the same millisecond `DateTime`.

**STAC items.**
- landsatlook (`landsat-c2l1`): `assets[<band>].alternate.s3.href` gives the requester-pays
  `s3://usgs-landsat` object; `MTL.json` and `ANG.txt` are assets of the same item.
- Earth Search (`sentinel-2-l1c`, ids like `S2A_7VEG_20200626_0_L1C`): `properties["s2:product_uri"]`
  names the `.SAFE`; asset hrefs point at `s3://sentinel-s2-l1c`, which the ITS_LIVE credentials cannot
  read (AccessDenied on list and get). An S2 item is therefore resolved through `s2:product_uri` to the
  anonymous Google mirror, as a path input is.

**Sentinel-1 SLC zips.** Every measurement TIFF in the golden cases' zips is *deflated*, not stored, so
no byte range of a remote or local zip maps to a byte range of the raster.

**ASF burst extractor.** A burst URL answers `202` while the burst is being generated, then redirects
(through Earthdata Login) to a presigned S3 URL valid for one hour that serves `Range` requests
(`206`, `Accept-Ranges: bytes`). Bursts can be range-read without a local copy; the presigned URL must
be resolved again when it expires. The golden S1-BURST case names are burst2safe-synthesized SAFE
names that ASF's search does not know, so those cases need their burst ids recorded separately.

## Parity baseline

The JIT app (`app/`) is the baseline for the binary: on LT05_L1GS_001013 two app runs agree exactly,
while two runs of the golden harness differ in coverage by ~1040 points in the same FFTW-wisdom state.
Of the 12 optical cases, the app matches the harness on 10, on S2A_07VEG in everything but `v_error`
(by 1, a field `compare.jl` treats as nondeterministic), and differs on LT05_001013 as above.

## Phase 1: the read path (done)

`probes/typedread.jl` is the gate: built with `--trim=safe` it has 0 verifier errors, and with the
user depot hidden it reads windows from a local Landsat COG, the same band from requester-pays S3, a
Sentinel-2 JPEG 2000 band locally and over HTTPS, and the Landsat band again through
`OpticalDatasets.open_optical(path, UInt16)`. Its sums match the JIT, and under the JIT the same
windows equal what `ChunkManifests.scan` reads, byte for byte.

What made it so, by package (all on local `trim` branches, uncommitted):

- **Zarr fork**: `filterdict` and `storageregexlist` are `const`; the scratch-buffer pool is keyed by a
  concrete type; `uncompress_raw!` passes a concrete fill-value keyword. Test counts unchanged
  (748 pass; the 15 errors are the Python-interop tests, failing identically before).
- **ChunkManifests**:
  - HTTP and S3 transports run on libcurl's multi interface (`src/transport/curl.jl`) instead of
    Downloads and AWS.jl. Requests in a batch run concurrently on one thread, which yields between
    10 ms GC-safe waits.
  - `S3Transport` signs its own SigV4 requests (checked against AWS's published example) with the SDK
    credential chain: environment, profile, ECS container endpoint, IMDSv2.
  - The HDF5 driver lives in the `ChunkManifestsHDF5Ext` extension, so HDF5.jl is no longer a hard
    dependency.
  - Typed opens `open_tiff(uri, T; transport)` / `open_jp2(uri, T; transport)` (`src/typed/typed.jl`)
    return `TypedArray{T}`, a `ZArray` over a concrete `TypedChunkStore` decoded by one `ChunkCodec`
    type. The general `scan` and its TiffImages extension are unchanged for interactive use.
  - Full suite 8287 pass, 0 fail.
- **OpticalDatasets**: `open_optical(path, T; transport)` reads through the typed opens, and the
  untyped `open_optical(path)` reads `T` from the file's header first. TiffImages is no longer a
  dependency. All tests pass, including the fixtures checked against GDAL.

## Phase 2: geometry and correlation (done)

`probes/correlate.jl` is the gate. It reads the inputs `probes/dump_correlate.jl` captures from the JIT
app — the two scene footprints, the parameter grid and its twelve windows, the two images and the
no-data mask — and runs coregistration, the blocked geogrid over a typed FastGeoProjections transform,
`AutoRIFT.params`/`pointset` from the geometry, and `autorift` with the pair's filter. Built with
`--trim=safe` it has 0 verifier errors. On five cases, each run from a copy of one frozen FFTW wisdom
file, the trimmed and JIT runs write identical geogrid bands and identical `dx`, `dy`, `chip_size` and
`interpolated`:

| case | grid → scene CRS | filter | levels | JIT s | trimmed s |
|---|---|---|---|---:|---:|
| LC09_L1GT_215109 | 3031 → 3031 (identity) | Highpass | 4 | 25.9 | 14.2 |
| LE07_L1TP_061018_2012 | 3413 → 32607 | WallisGapfill | 3 | 29.3 | 18.2 |
| LT04_L1TP_063018 | 3413 → 32607 | none (native `:fft`) | 4 | 22.4 | 12.2 |
| LC08_L1TP_009011 | 3413 → 32622 | Highpass | 4 | 109.7 | 99.6 |
| LE07_L1TP_063018_2004 | 3413 → 32607 | WallisGapfill | 4 | 91.3 | 78.6 |

Under the JIT the probe reproduces the app's own geogrid bands exactly; its correlation equals the app's
in the same wisdom state (LC09, where the dump and probe ran back to back). The `:fft` native filter's
destripe on raw plans is bitwise equal to FFTW.jl's on the LT04 scene, through the band-reject branch.

What made it so (local `trim` branches, uncommitted):

- **FastGeoProjections**: `NativeCRS` descriptors (`GeographicCRS`, `UTMCRS(zone, north)`,
  `PolarStereographicCRS(3413|3031)`, `native_crs(code)`) and `Transformation(source, target, T)` over
  them, whose type follows from the CRS kinds; bitwise equal to the EPSG constructor (`test/crs.jl`).
- **ImagePairGeometry**: `pairgeometry_blocked` specializes on its `transform`; windows are held as
  the concrete `GridWindow`; `_block_view` names every band; two error messages no longer print a
  struct or interpolate more than ten values. Tests: 56,211 pass.
- **AutoRIFT**: the geometry extension's `params` passes method objects and builds a four-entry
  subpixel ladder whose tail repeats the last level's factor, which `_check_subpixels` accepts, so the
  `Params` type does not depend on the level count; `_outliers` prints `nameof`; spawned tasks use
  `Threads.@spawn` (StableTasks' nested return-type inference left `_parallel_slices` uncompiled under a
  StableTasks task, and StableTasks is no longer a dependency); destripe runs on guru64 `ESTIMATE`
  plans planned as FFTW.jl plans them; every `ccall` names `FFTW_jll.libfftw3f` at run time. Tests:
  705,007 pass.

## Phase 3: the product writer (done)

`probes/writeproduct.jl` is the gate. `probes/capture_product.jl` runs a case through the JIT app with the
NCDatasets writer and captures the `ItsLiveInput` it was handed; the probe rebuilds that input and writes
it with the new writer, once as captured and once with a `roi_valid_percentage` that selects the
uncropped schema. Built with `--trim=safe` it has 0 verifier errors, and on `LC08_L1TP_009011` both
products equal the NCDatasets writer's: every line of `ncdump -h -s` (types, `NC_STRING`/`NC_CHAR`,
`LL`/`f` suffixes, chunking, deflate and shuffle, fill values) in the same order, apart from
`date_created`, and every variable's data (`probes/cmpnc.jl`). `cf_grid_mapping` on GDAL's OSR API
returns the same attributes and types as the ArchGDAL version for 3413, 3031 and six UTM zones.

ItsLiveOffsetProduction (`trim` branch, uncommitted):

- `write_product` writes through the netCDF C library (`src/netcdf.jl`); the cropped/uncropped string
  type is the file's `nc_string` flag rather than a `ScopedValue`.
- Attribute values are `AttrValue`, a closed union written by an explicit test of each member;
  `ImagePairInfo.satellite_img*` is `Union{String,Int64}`, `extra` and `ItsLiveGeoref.mapping_attrs`
  are `Dict{String,AttrValue}` (any dictionary of those values converts); `ItsLiveGeoref{F}` and
  `ItsLiveInput{F}` carry the `lonlat` callable's type.
- Inside `write_product`: fused broadcasts deep enough to have their type widened are loops; the error
  statistics and the written arrays are passed as two concrete structs with small-union fields; dates
  are formatted field by field (`Dates.format` pads through `textwidth`); `median` is the package's own
  (`Statistics.median!` prints the vector in its empty-input error).
- NCDatasets, Rasters, DimensionalData and ArchGDAL are no longer dependencies; `ItsLiveGeoref(::Raster)`
  is the `ItsLiveOffsetProductionRastersExt` extension. Tests: 106 pass.

## Phase 4: the optical binary (done)

`app/build.sh <dir>` builds `autorift_itslive_production` with `--trim=safe`: 0 verifier errors. The
gate ran each of the 12 optical golden cases through the JIT app (as it was before the trim rewrite) and
through the binary, as separate processes, each starting from a copy of one frozen FFTW wisdom file
(populated by one earlier run of every case), with `JULIA_NUM_THREADS=6,1` and the build depot hidden.

- **Products:** every data variable is identical in all 12 cases. The only header differences are the
  date attributes: the JIT app's NCDatasets writer stores a `SubString` attribute with the bytes of its
  parent string (`".080"`, `".000"`), the binary writes the substring (`".08"`, `"."`, as the Python
  S2 reference does).
- **Inputs:** an `s3://` path pair and a landsatlook STAC-item pair give the same LC08 product; an Earth
  Search STAC-item pair and a Google-mirror URL pair give the same S2 product.
- **Memory:** peak footprint at most 8.06 GiB (LT05_001013).
- **Time:** wall time from the first read to the written product, against `pylogs.py`:

| case | wall s | peak GiB | × Python |
|---|---:|---:|---:|
| LC08 009011 | 102.6 | 4.28 | 3.67 |
| LC08 060018 | 110.7 | 2.50 | 2.96 |
| LC08 062018 | 68.8 | 3.70 | 4.20 |
| LC09 215109 | 92.7 | 3.14 | 3.53 |
| LE07 061018 2012 | 65.5 | 2.24 | 4.14 |
| LE07 061018 2013 | 109.9 | 2.76 | 3.26 |
| LE07 063018 2004 | 106.0 | 3.46 | 3.16 |
| LT04 063018 | 41.6 | 6.26 | 4.68 |
| LT05 001013 | 46.8 | 8.06 | 4.06 |
| LT05 060018 | 39.3 | 4.78 | 3.00 |
| S2A 07VEG | 99.3 | 1.63 | 1.87 |
| S2B 22WEB | 99.5 | 1.87 | 1.88 |

The imagery stage, which reads the scenes over the network (about 20 MB/s here), is most of the wall
time: 73 s of LC08 060018's 111 s, 80–84 s of each S2 case's 99 s. From local JP2s the two S2 cases run
in 22.9 s and 26.4 s (about 8× and 7× Python, whose S2 runs read local files) with identical products.
The Landsat cases short of 4× are bound by the same link.

## Phase 5: NISAR (done, local and HTTP reads)

The binary takes `NL1` (RSLC) and `NL2` (GSLC) granules as paths or URLs: 0 verifier errors.

- **ChunkManifests** (`src/typed/h5.jl`): HDF5 metadata through `HDF5_jll`'s C API (`h5open_typed`,
  `h5_exists`, `h5_isgroup`, `h5_dims`, `h5_read`, `h5_read_elements`, `h5_read_strings`,
  `h5_attr_string`), with a range-request file driver for a remote file, and `h5_open_array`: a chunked
  dataset (shuffle and deflate) or a contiguous one (one chunk per row) as a `TypedArray`, its chunks
  located with `H5Dchunk_iter` and read through the transport. `fetchall` reads a whole small object
  without `Range`; `http_request` makes one bodyless request. Tests: `test/typedh5.jl` (local and HTTP),
  `test/http.jl`.
- **SLCDatasets**: the NISAR backend reads through `H5File`, opening the product for each read and closing
  it after, so a backend holds no file; HDF5.jl, HTTP, TiffImages and Base64 are no longer dependencies
  (the Sentinel-1 measurement raster opens with `open_tiff`; the Earthdata token request uses
  `http_request`). Tests: 3,171 pass.
- **ImagePairGeometry**: `rdr2geo` specializes on its `height` source (passed through unspecialized, it
  dispatched on every terrain sample), and `footprint_bounds` builds its sample tuples with a `Val` length.
- **App**: `nisar.jl` (GSLC crop, amplitude, identification) and `radar.jl` (RSLC setup, offset lattice,
  Doppler-demodulated eight-tap resample, `:radar` packaging).

Gate, from local granules, each run from a copy of one frozen wisdom file:

| case | wall s | peak GiB | × Python | product |
|---|---:|---:|---:|---|
| NL2 P093 | 552 | 10.6 | 5.1 | every data variable identical to the JIT app's |
| NL2 P093 over HTTP | 504 | 10.6 | — | identical to the local read |
| NL1 P094, `--dem` the reference `dem.tif` | 942 | 13.8 | 4.8 | against the Python product: `vx`, `vy`, `v`, `va` median \|Δ\| 1 m/yr, p95 8–10 m/yr; `vr` 1 and 3 |
| NL1 P094, ITS_LIVE DEM + EGM2008 | 909 | 13.7 | 4.9 | 99.7% of `vx`/`vy`/`v` identical to the `--dem` run |

Python: 2,817 s (NL2) and 4,476 s (NL1), both from local granules. The JIT app on NL2 took 522 s at a
12.0 GiB peak. The HTTP run reads from a local range server (`HTTP.jl`).

**From ASF.** NL2 read in place from `nisar.asf.earthdatacloud.nasa.gov` with an Earthdata Login token:
2,690 s, of which 2,242 s is the imagery stage at about 10 MB/s; peak 10.7 GiB; every data variable
identical to the local run. Each GSLC tile cache starts on a boundary of the granule's 512² chunks
(`TileFile.origin`), so its 1024-pixel tiles contain whole chunks and each chunk is fetched once; with the
cache aligned to the crop instead, a chunk was fetched by up to four tiles, and the run had read 34 GB when
it was stopped. An NL1 secondary is still read through the resample's support windows, which overlap.

- **The NL1 DEM.** The ITS_LIVE parameter DEM holds heights above the EGM2008 geoid; the reference's
  `dem.tif` holds heights above the ellipsoid. Over 1,450 ground points of P094 their difference is
  −30.9 m (median), and adding EGM2008's undulation (`us_nga_egm08_25.tif` from PROJ's CDN) leaves
  0.0 m median and 9.6 m at the 95th percentile of |difference|. `parameter_dem` therefore adds the
  undulation; `--dem` takes an ellipsoidal geographic DEM in its place.
- **The radar `dy` sign.** `testautoRIFT.py:785-786` negates a radar pair's `Dy` before packaging
  (`if optical_flag == 0: Dy = -Dy`), on top of the sign every pair's `dy` takes (`DY_SIGN`); a NISAR L2
  pair runs with `optical_flag = 1`. Without the second negation the NL1 product's `va` was the
  reference's negated (correlation −0.998, signed median difference −153 m/yr against a median |va| of
  77 m/yr), and so were `vx`/`vy`; with it the binary agrees with the reference as in the table. Before
  the change the binary's NL1 product equalled the harness's (`900/julia_native.nc`) in every data
  variable; the harness (`tools/golden/packaging.jl` `dy_sign`) now applies the same sign.

## Phase 6: Sentinel-1 (done)

The binary takes `S1` (an IW SLC's `.SAFE` or `.zip`, local or a URL) and `S1BURST` (the bursts' extractor
URLs, comma-separated, or the SAFE `burst2safe` assembles from them): 0 verifier errors. Orbits are found by sensing time in the
public `s1-orbits` bucket, or read from `--orbit-dir`; terrain heights come from the parameter region's DEM
plus EGM2008, or from `--dem`.

- **ChunkManifests** (`src/typed/zip.jl`): a zip archive's central directory (`zip_members`, ZIP64
  included) and its members by byte range through any transport: `read_member` inflates a small member
  whole, and `ZipEntry` indexes a deflated one — one streaming pass recording an access point (bit offset,
  input offset, 32 KiB window) every 4 MiB of output — so `read_entry` inflates a window from the nearest
  point. On the IW1 measurement member of `S1A_IW_SLC__1SSH_20151120T080202` (820 MB deflated, 1.14 GB
  raw) the index is 271 points, built in 5.3 s from a local file, and 40 random windows (1 B to 6 MiB,
  including the first and last bytes) equal `unzip -p`'s bytes. Tests: `test/zip.jl` (local and HTTP).
- **SLCDatasets**: a zipped product opens in place — annotations from the central directory, each
  subswath's measurement raster through a `ZippedMeasurement` over its member's index (built once per
  process) — and every subswath's burst window read from a zip equals the unpacked `.SAFE`'s. ZipArchives
  and Mmap are no longer dependencies; the S1 path no longer formats dates with `Dates.format`/`lpad` or
  prints vectors in its messages. Tests: 3,171 pass.
- **App** (`s1.jl`, `s1orbit.jl`): the `merge_bursts_in_swath`/`merge_swaths` mosaic, the TOPS carrier
  and its ISCE3 polynomial fit, the deramped eight-tap resample of each secondary burst through the offset
  lattice, the subswath bias inputs for a pair of two spacecraft, and the `:radar` packaging.
- **Harness**:
  - `ResampledSwath` reads each secondary burst with the secondary's own burst shape, which can differ from
    the reference's. P086's secondary bursts are two lines longer per burst, so burst `i` was read
    `2(i − 1)` lines early; `va` was off by 0, ±2, ±4 and ±6 pixels in burst-aligned bands.
  - `julia_e2e.jl --lattice` runs a full-SLC or burst pair through the lattice path rather than replaying
    COMPASS's offsets.
  - `s1_pair` reads the run's granules rather than the burst extractor.

Gate, each case from its run 400 (granules, orbits, DEM), one frozen wisdom file. On every case `vx`, `vy`,
`v`, `vr` and `va` equal the harness's lattice product exactly.

| case | sensor | wall s | peak GiB | × Python | against the reference: median \|Δ\| vr / va / v (m/yr) |
|---|---|---:|---:|---:|---|
| S1A 20150828 (P086) | S1 | 213 | 5.4 | 6.8 | 1 / 11 / 8 |
| S1A 20151120 (P002) | S1 | 275 | 8.9 | 6.7 | 27 / 53 / 46 |
| S1A×S1B 20170221 (P089) | S1 | 280 | 5.7 | 9.3 | 11 / 79 / 42 |
| S1B 20180809 | S1 | 275 | 5.5 | 9.7 | 4 / 27 / 17 |
| S1C 20250416 | S1 | 241 | 6.5 | 9.3 | 4 / 26 / 15 |
| S1A 20240618 025528 | S1BURST | 101 | 6.1 | 15.3 | 9 / 48 / 33 |
| S1A 20240618 025533 | S1BURST | 72 | 3.2 | 8.8 | 14 / 57 / 42 |
| S1C 20250416 010159 | S1BURST | 24 | 2.8 | 15.7 | 2 / 18 / 10 |

Python's times include downloading the granules; the binary's read local files. The differences from the
reference are the lattice coregistration's against COMPASS's, which resamples some reference bursts onto a
static layer's grid; P002 is 3% coherent.

**Burst-extractor URLs.** An `S1BURST` input may be the bursts' `sentinel1-burst.asf.alaska.edu` URLs, the
`url` ASF's search lists for each burst.

- **ChunkManifests**: `http_resolve(urls)` follows each URL's redirects, concurrently, and returns the final
  status and URL. An extractor URL leads through Earthdata Login, authenticated from `~/.netrc` (the
  extractor refuses a bearer token), to a pre-signed S3 URL.
- **SLCDatasets** (`src/extractor.jl`, `src/asf.jl`):
  - `AsfBurst(url)` reads a burst from its URL.
  - A burst's metadata is its parent product's whole annotation. It is fetched once per product and
    polarization and held in memory.
  - A raster is downloaded into the run's scratch directory in 16 MiB ranges, 8 in flight. All of a
    product's rasters are requested together, so the extractor generates them while earlier ones download;
    a `202` is asked again every 15 s, 20 times.
  - `Sentinel1Product(bursts; orbit)` is the product `burst2safe` assembles from the bursts, and the
    product's name is the one `burst2safe` gives it, with `0000` for the manifest checksum. Per subswath:
    - the bursts in azimuth order, refused if they are not consecutive;
    - scalar fields from the earliest burst's parent;
    - `samplesPerBurst` the widest burst's;
    - the FM-rate and Doppler-centroid lists of the parent products merged and cut to 500 s and 3 s of the
      bursts' span;
    - a `BurstStack` raster, each burst from its own file, zero past its own width.

  On all six products of the three golden burst cases, the assembled annotations equal those read from
  run 400's `burst2safe` SAFE field for field, and the stacked rasters equal the SAFE's measurement rasters
  on full-width bands at the first line, across the first seam and at the end. This includes P072, whose
  bursts come from two parent products. The name differs only in the checksum. `burst2safe` takes the
  checksum over a manifest stamped with its creation time, so it differs between runs: the S1C reference
  product names its inputs `6CB9`/`D395`, while run 400's SAFEs for the same bursts are `D7CD`/`F294`.
- **App**: `scene_source` names a burst-URL scene by its assembled product, and `open_s1` opens it with the
  bursts downloaded into the run's scratch directory.

Gate, each case with run 400's orbits and DEM, from a copy of the gate wisdom file. `vx`, `vy`, `v`, `vr`,
`va` and `v_error` equal the SAFE-input product's in every pixel. The only differences in the file are:
- `id_img*` and `product_unique_ID_img*` (`0000` for the SAFE's checksum);
- `time:microseconds_added`, a hash of the output file's name;
- `date_created`.

| case | bursts | wall s | peak GiB | × Python |
|---|---:|---:|---:|---:|
| S1A 20240618 025528 (P072) | 24 + 24 | 796 | 5.1 | 1.9 |
| S1A 20240618 025533 (P062) | 10 + 10 | 369 | 3.0 | 1.7 |
| S1C 20250416 010159 (P022) | 7 + 7 | 241 | 2.5 | 1.6 |

Both sets of times include downloading the bursts, so here the binary is under the 4× target. The run is
bound by the network. A burst raster is 120–150 MB, and the S1A 025528 pair downloads 6.7 GB. At the time of
the gate, one `curl` stream through the extractor reached 10 MB/s and 8 concurrent ranges 11 MiB/s. Earlier
the same day, the 8 ranges reached 26–28 MB/s. The SAFE-input times above are the compute alone.

## Binary against the JIT at 10 threads

Each case ran through the binary and through the same source under the JIT
(`julia --project=app -t 10,1`, packages precompiled). Each run had `JULIA_NUM_THREADS=10,1`, started from
a fresh copy of the shipped wisdom, and ran back to back with its counterpart in alternating order. Inputs
are as in each sensor's gate: optical scenes from S3, landsatlook and the Google mirror; S1 from local
zips; S1BURST from extractor URLs, downloads included; NISAR from local granules. `vx`, `vy` and `v` are
identical between binary and JIT in all eight cases. Times are seconds; the stage sum excludes start-up
and compilation.

| case | JIT wall | binary wall | JIT stages | binary stages | JIT / binary peak GiB |
|---|---:|---:|---:|---:|---|
| LC08 009011 | 110 | 69 | 73.5 | 67.5 | 6.2 / 4.3 |
| LE07 063018 2004 | 127 | 83 | 90.2 | 82.0 | 5.0 / 3.5 |
| S2B 22WEB | 114 | 76 | 77.5 | 75.7 | 3.3 / 1.9 |
| LT05 001013 | 86 | 46 | 48.9 | 44.8 | 8.0 / 7.6 |
| S1 P086 | 207 | 170 | 170.5 | 169.6 | 6.2 / 5.1 |
| S1BURST P022 (URLs) | 247 | 206 | 196.5 | 188.8 | 4.3 / 2.8 |
| NL2 P093 | 466 | 395 | 426.6 | 392.4 | 13.9 / 14.4 |
| NL1 P094 | 722 | 664 | 683.8 | 661.9 | 18.8 / 17.7 |

At 10 threads NL1's peak exceeds 16 GiB under both the binary and the JIT; at 6 threads it was 13.7–14.8 GiB.

## Remote reads from ASF

Every case below read its granules from ASF at 10 threads, with the bytes received on `en0` sampled every 2 s.
Times include the download. The link here carries 15–25 MiB/s, depending on the host and the hour.

| case | before | after | received | stored in what was read | link while reading |
|---|---:|---:|---:|---:|---|
| S1 P086 (ASF zip URLs) | not run: the window reader re-fetched each line once per tile column | 296–303 s | 4.4 GiB | 3.8 GiB zips | busy from the start; 20 MiB/s average |
| NL2 P093 | 2,200 s (imagery 1,821 s) | 1,877 s (imagery 1,447 s) | 21.8 GiB | 19.2 GiB of crop chunks | 13.6 MiB/s average against 15–16 |
| NL1 P094 | — | 1,910 s (imagery 1,287 s) | 23.0 GiB | 20.4 GiB of both rasters | 17–20 MiB/s |

The products equal the local-input products. A binary run and a JIT run of S1 P086 from the zip URLs take
303 s and 296 s, a difference within the link's variation. The S1C burst pair from extractor URLs takes
190 s.

What changed:

- **Zip members** (ChunkManifests `extract_member`, SLCDatasets `extract_measurements`). A zipped product's
  measurement members are inflated whole into the run's scratch directory, every subswath's at once. Each
  member keeps four 16 MiB ranges in flight while the previous batch inflates. A remote zip is first
  followed through Earthdata Login with `~/.netrc` (ASF's datapool answers a bearer token with `401`), and
  its central directory is kept for the session. The window reader this replaces re-inflated each line
  once per tile and, from a remote zip, re-fetched it. Scratch grows by the members' inflated size, about
  5.5 GiB for a pair.
- **Overlap.** The two S1 products open together. Their rasters download while the geometry is solved and
  the DEM read, the reference's first, so its mosaic is filled while the secondary's arrive
  (`prefetch_rasters`). An NL1 reference is read while the DEM is.
- **Tiles in flight** (`fill_tiles!(; inflight)`). A NISAR read keeps 4 tiles per thread in flight. With one
  per thread, the tiles at a granule's edges, whose chunks are mostly unwritten, left the link at 8–10 MiB/s.
- **NL1 secondary chunks** (ChunkManifests `cached`). The resample's support windows overlap their
  neighbours. Each chunk is fetched once into a file under the scratch directory, and a chunk another task
  is fetching is waited for rather than fetched twice.
- **HDF5 metadata.** The metadata bytes of a remote granule are kept for the session across opens, keyed
  by the credentials they were read with.
- **GDAL calls that wait.** Network opens, window reads and warps are GC-safe, so a collection that another
  thread starts meanwhile does not wait for them.
- **Errors.** A tile or fetch task that fails reports its own error rather than `TaskFailedException`.

Still latency-bound here: opening a NISAR granule. Each request goes through ASF's redirect, about 1.7–2.8 s
round trip on this link, and the first open makes about 35 of them, so geometry takes 40–65 s. In-region the
same requests take tens of milliseconds.

Scratch at 10 threads: NL2 49 GiB, NL1 47 GiB (the two tile caches and the chunk cache), S1 SLC about
15 GiB. A run killed with a signal leaves its scratch directory behind.

On macOS the first launch of a freshly built bundle is scanned before it starts, which took up to 40 s.
Every timing above was taken from a binary that had already been launched.

## Phase 7: release builds (macOS done; Linux built, not validated)

**Shipped FFTW wisdom.** `app/wisdom/` holds one file per CPU model; the bundle carries it under
`share/autorift_itslive_production/fftw-wisdom`. An executable started without `AUTORIFT_FFTW_WISDOM_DIR`
reruns itself with the variable pointing at a temporary copy, deleted when the run ends (`ensure_environment`
in `src/cli.jl`). The macOS file is the union of every wisdom file the gates measured (46 files). Two runs of
the final macOS bundle on LT05 060018 started that way gave identical data, equal to the Phase 4 product.

**macOS arm64**: `app/build/macos-arm64`, 0 verifier errors.

**Linux x86_64**: `app/build.sh` cannot cross-compile, so `app/build_linux.sh` builds in an amd64 Julia
1.13.1 container, with Julia's release multiversioning (`generic`, `sandybridge`, `haswell`, `x86-64-v4`):
`app/build/linux-x86_64`, 0 verifier errors. It has not run on x86_64 hardware. The `itslive` credentials
cannot launch EC2 instances: `ec2:RunInstances`, `ec2:DescribeImages`, `ec2:DescribeInstances` and
`ssm:GetParameter` are all denied.

Under Docker's x86_64 emulation on this Mac (Rosetta), the bundle runs LT05 060018 to completion (99 s, one
thread) but every `vx`/`vy` saturates.

- The fault is in FastGeoProjections' scalar polar-stereographic-to-UTM composition. At four test points
  it returns a different, non-invertible point when Julia compiles for `haswell` (AVX2/FMA). For
  `westmere`, `k8-sse3`, `--inline=no` or `-O0` it equals macOS bit for bit. The individual stages and
  the vector path equal macOS in every case.
- Under QEMU's emulator the same code gives the macOS result. A build for `generic` alone, which has no
  AVX2 clone, gives sensible products under Rosetta (median |Δv| 7 m/yr from macOS).
- So the evidence points at Rosetta's AVX2 emulation rather than the code, but only a run on real x86_64
  hardware can settle it, and no Linux product has been checked against the JIT chain or the reference.

**Not done**:
- the native Linux parity suite and the 22-case sweep, which need an x86_64 host;
- Linux wisdom files, which have to be measured on the deployment CPU models.
