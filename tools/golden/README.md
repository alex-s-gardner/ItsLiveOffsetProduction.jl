# The ITS_LIVE product golden tests

`s3://its-live-data/test-space/golden/` holds the acceptance products for the Python autoRIFT — 22
ITS_LIVE granules built by `hyp3_autorift` 0.28.4. These scripts build the same products with
`ItsLiveOffsetProduction` and compare them.

Run by hand, not in CI: they need `AWS_PROFILE` credentials, requester-pays egress, and Docker for the
reference container. `Pkg.test()` at the repository root covers everything that runs offline.

```bash
julia --project=tools/golden -e 'import Pkg; Pkg.instantiate()'

# The writer against a real captured `netCDF_packaging` call from the reference container.
julia --project=tools/golden tools/golden/validate_itslive_write.jl <case-fragment>
julia --project=tools/golden tools/golden/run_all_golden.jl              # every case

# Two raw granules to a product netCDF, with nothing taken from the reference. Optical granules are
# streamed from their buckets and NISAR granules from ASF (Earthdata Login from ~/.netrc).
AWS_PROFILE=itslive julia --project=tools/golden -t 8 \
    tools/golden/julia_e2e.jl LC08_L1TP_009011_20200703 --compare-intermediate

# Product-level comparison, and the reproducibility floor of the reference itself.
julia --project=tools/golden tools/golden/run.jl --status
julia --project=tools/golden tools/golden/run.jl --reproducibility S2B_MSIL1C_20200612
```

Data lives outside the repository, under `~/data/autorift/tests/golden_tests` — override with
`AUTORIFT_GOLDEN_CACHE`.

## What is here, and what stayed in AutoRIFT.jl

These scripts are the product half. The correlator half — the stage ladder (`stages.jl`,
`correlator.jl`), the block-size and profiling tools, the figures, `selftest.jl` and the gate ledger
in `dev/GATES.md` — stays in
[AutoRIFT.jl](https://github.com/alex-s-gardner/AutoRIFT.jl)`/tools/golden/`, because AutoRIFT's own
source comments cite it as the evidence for correlator behaviour.

The library files here are **copies** of that harness, not the originals: `manifest.jl` (the case
registry), `reference.jl` (the pinned container), `scenes.jl`, `parameters.jl`, `e2e.jl`, `e2e_run.jl`,
`mtl.jl`, `radar.jl`, `nisar.jl`, `tilecache.jl`, `stages.jl`, `correlator.jl`, `intermediate.jl`, and
`../ab/{xchg,memtrace}.jl`. They are duplicated because `julia_e2e.jl` reaches the granule through
`e2e_run.jl` → `e2e.jl`, and that chain is equally the on-ramp for AutoRIFT's correlator gate. A
change to the shared geometry or scene-access logic has to be applied in both places.

`product.jl` and `compare.jl` are the ITS_LIVE product reader and comparator and are also duplicated,
because AutoRIFT's `selftest.jl` gates them.

## Why names are qualified

Scripts here load both `AutoRIFT` and `ItsLiveOffsetProduction`, so the moved API is written
`ItsLiveOffsetProduction.ImagePairInfo` rather than bare. Writing it bare resolves only as long as
`AutoRIFT` does not also export the name, which it did until the extraction landed.

## Future change: fetching only the chunks a pair needs

Over sparse glacier coverage much of a scene is never correlated, so downloading all of it is waste.
`julia_e2e.jl` reads each scene a tile at a time through a `TileCache` (next section), but
`native_imagery` fills every tile of the overlap before the no-data mask is sampled, so the whole overlap
is still fetched. What is missing is a rule for which tiles to skip. This was **deferred, not rejected**,
blocked on a single thing that is now gone.

**The former blocker, cleared.** `AutoRIFT.bytescale` was a two-pass whole-array reduction whose global
mean and sample standard deviation set every quantized pixel, so a chunk never downloaded would have
changed the quantization of every pixel that was. `native_imagery`/`radar_imagery` no longer call it —
`dev/CORRECTNESS.md` item 4b on `AutoRIFT.jl` records the decision — so `highpass` (bounded local reach)
is now the only pass in this path, and nothing here still needs an image-wide statistic. This change is
unblocked; it has not been made.

**What is still required even then.** `_accumulate_zeros_banded!` samples the scene at grid points to build the
no-data mask, so a chunk holding grid points is still read. That is point-wise rather than global, which
is what makes skipping possible at all.

**How to decide what to fetch.** A chunk is needed if any *pre-imagery-searchable* point's
chip-plus-search footprint touches it, where searchability is the geogrid sentinel together with the
parameter rasters' own search limits. That is the only information available before downloading: the
pipeline's effective searchability is much narrower but imagery-derived, since `julia_e2e.jl` zeroes
search radii from the no-data mask after reading.

**Measured upside, on the two cases to hand:**

| case | pre-imagery searchable | chunks skippable | per band |
|---|---:|---:|---:|
| `LC08_L1TP_009011_20200703` | 67.6% | 630 of 4489 (14.0%) | 78.8 of 561 MiB |
| `LT05_L1GS_001013` | 82.3% | 34 of 1088 (3.1%) | 4.2 of 136 MiB |

**Both of those are high-coverage cases, so the case this is actually for is unmeasured.** A scene whose
ice occupies a fifth of the frame is where the strategy pays, and neither of these is that. Re-run the
method above on a sparse-coverage pair before investing in it.

**The hazard to measure against.** Per-chunk requests give up the batched, chunk-aligned slab reads that
took a 338 MiB band from 54.9 s to 32.3 s over `/vsis3` (`OpticalDatasets.jl`'s `src/raster.jl`,
`open_optical` docstring). A chunk is 0.12 MiB, so per-request latency can easily cost more than 14%
fewer bytes saves. The condition on this change is that it not degrade runtime.

## How `julia_e2e.jl` holds the imagery

**Optical pairs are lazy and correlated blocked.** Each scene is a chain of windowed arrays
(`streamed.jl`): the granule's object, the native filter `process.py` applies to it, the warp onto the
reference's projection for a pair that straddles a UTM zone, and the crop to the overlap. A `TileCache`
keeps each 1024² tile of the result on local scratch once derived, and the correlator filters a block at
a time (`preprocess` is the correlator filter, not `:none`), keeping each filtered 512² tile in a second
scratch file. Both caches are read and written with positioned I/O (`--cache-io pread`) rather than mapped,
so their pages stay in the kernel's file cache instead of the process's resident size; `--cache-io mmap`
maps them. Two parts are not windowed:

  * `:fft` (Landsat 4/5) — `Destripe`'s band-reject is chosen from a whole-scene spectrum, so the scene is
    streamed whole, filtered, and held in a memory-mapped scratch file.
  * The warp is GDAL's Lanczos warp of one window at a time. It reproduces `_warp_to`'s grid exactly and
    its values to within 1 DN at 99.8% of pixels on `LE07_L1TP_061018_20120428`; the rest differ where
    GDAL's approximate transformer, which interpolates across each warp chunk, lands a no-data edge on the
    other side of a pixel.

On axis-aligned grids blocking costs nothing in the answer: AutoRIFT's `dev/GATES.md` records the Landsat
sweep keeping `dx`/`dy`/`correlation` "bit-identical to an untiled run at every block size", and the
streamed, blocked `LC08_L1TP_009011_20200703` product is bit-identical to the staged untiled one.

**`WallisGapfill` noise is per window.** The native and correlator fills draw from a seeded generator per
window and per block, so filled values differ from a whole-scene call's; which pixels are filled does
not. The reference's own fill is unseeded, so it does not reproduce itself either (rung 5.4).

**Radar pairs are lazy and correlated blocked too.** The Sentinel-1 reference mosaic is assembled a window
at a time (`lazy_radar_mosaic`), the secondary is resampled a window at a time (`ResampledMosaic`,
`ResampledRSLC`), and each is behind a `TileCache`. On a rotated radar geogrid blocked and untiled runs do
not agree exactly: `dev/CORRECTNESS.md` on `AutoRIFT.jl` records 5.3% of an S1B run's points differing
(26,781 lost, 771 gained, 25,206 measured differently), traced to a coarse-pass radius-widening rule that
lets the untiled run search a scene corner no block's read window can reach.
