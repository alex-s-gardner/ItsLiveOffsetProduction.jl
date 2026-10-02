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

# Two raw granules to a product netCDF, with nothing taken from the reference.
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
`OpticalRaster` is already a `DiskArrays` array with a real chunk grid, so per-chunk reads are available;
what is missing is a cache and a consumer that asks for chunks rather than the whole overlap. This is
**deferred, not rejected**, and the blocker is a single thing that is on its way out.

**The blocker.** `AutoRIFT.bytescale` is a two-pass whole-array reduction whose global mean and sample
standard deviation set every quantized pixel, so a chunk never downloaded changes the quantization of
every pixel that was. It is the only image-wide statistic left in this path — `julia_e2e.jl` applies
`highpass` (bounded local reach) and passes `preprocess = :none` to the correlator — and it exists only
because the production driver sets `DataType = 0`. AutoRIFT's own register records that quantization as
matched-not-endorsed and states that production is moving to `Float32` (its
`tools/golden/README.md`, "Matched for agreement, not endorsed", the `UInt8` row), with the float path
bit-identical while the byte path fails on 1.7% of points by up to 36 px. **When `UInt8` goes, no
image-wide statistic remains and this becomes viable.**

**What is still required even then.** `_accumulate_zeros!` samples the scene at grid points to build the
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

## Future change: blocking the correlation to cut peak memory

Separate from the above, and the larger lever on memory. The L8 case peaks at 12.46 GiB resident, of
which the imagery is about 3.0 GiB — two full-scene `Float32` buffers at 1.09 GiB each, a `Bool`, and two
`UInt8`. The remainder is the untiled correlator, which `process_block_size` bounds and which needs no
lazy reading at all. Blocking has been measured elsewhere in this codebase at 36–56% lower peak.

The change is small: `julia_e2e.jl` correlates with `AutoRIFT.autorift(b2, b1, grid, p)`, and the blocked
form is `AutoRIFT.autorift(b2, b1, grid, p, (b.X, b.Y))` with
`b = AutoRIFT.block_size_for(grid, p, size(b1))` — the five-positional method already exists, and
`e2e_run.jl` does exactly this for the radar path. Neither `cache_budget` nor `filter_cache_tile` applies
here: the byte images are resident and the correlator is given `preprocess = :none`, so no per-block
filtering happens and the halo is the correlation reach alone.

**What blocking costs, by grid type.** On axis-aligned grids it costs nothing in the answer: AutoRIFT's
`dev/GATES.md` records the Landsat sweep keeping `dx`/`dy`/`correlation` "bit-identical to an untiled run
at every block size". On a rotated radar geogrid it is not yet safe — `dev/CORRECTNESS.md` records 5.3%
of an S1B run's points differing from untiled (26,781 lost, 771 gained, 25,206 measured differently),
traced to a coarse-pass radius-widening rule that lets the reference search a scene corner no block's
read window can reach, and names it the blocker for tiled processing. So this is adoptable for optical
before it is adoptable for radar.
