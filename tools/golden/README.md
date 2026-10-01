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
