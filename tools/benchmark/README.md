# Benchmark image pairs

`pairs/<group>.tsv` holds 1000 image pairs per benchmark group, one per line in the column order the EC2
runner reads: sensor code, product name, reference scene, secondary scene, reference input, secondary
input. The inputs are what `autorift_itslive_production` takes:

| Group | Sensor code | Inputs | Source of the pairs |
|---|---|---|---|
| `L4L5`, `L7`, `L8L9` | the reference's: L4, L5, L7, L8, L9 | landsatlook STAC items | ITS_LIVE v02 Landsat catalog |
| `S2` | S2 | the B08.jp2 on the Google Cloud mirror | ITS_LIVE v02 Sentinel-2 catalog |
| `S1` | S1 | the SLC zip on ASF's datapool | ITS_LIVE v02 Sentinel-1 catalog |
| `S1BURST` | S1BURST | the bursts' extractor URLs, comma-separated | ASF's burst search over NISAR frames |
| `NL1`, `NL2` | NL1, NL2 | the RSLC / GSLC `.h5` on ASF | NISAR beta products in CMR |

`pairs/<group>_meta.tsv` gives each pair's time separation, location and, for the catalog groups, its EPSG
code, `percent_valid_pix` and catalog directory; `pairs/counts.tsv` the number of pairs per group in the
catalogs.

Each group is a uniform random sample without replacement, reproducible: a pair's key is the first bytes of
`sha256(seed + product)`, and the 1000 smallest keys are taken.

- **Catalog groups.** `sample_pairs.py` scans every chunk of `s3://its-live-data/catalog_geojson/` (run it in
  us-west-2; the catalogs are about 70 GB). A Landsat pair is grouped by its reference scene; `LO08`/`LO09`
  references (OLI without TIRS) are counted but not sampled, because the executable's L8/L9 references are
  `LC08`/`LC09`. Scenes USGS has reprocessed since the catalogs were written have new IDs on landsatlook;
  `pairs/landsat_remap.tsv` lists the 160 replaced (old ID, current ID), matched by sensor, path/row and
  date. Column 2 keeps the ITS_LIVE product name.
- **NISAR groups.** `nisar_pairs.py` lists CMR's NISAR beta products and pairs each scene with the same
  frame (track, direction, frame, mode, polarization) 1 to 4 cycles (12-48 days) later; `nisar_ice.jl`
  keeps the pairs over land ice — any of a 5x5 grid of points spanning ±1° around the footprint centre on
  the ITS_LIVE parameter region's ROI raster; `nisar_draw.py` drops the pairs without frequency A, because
  the pipeline reads frequency A only, and takes the sample. Of the 13,546 NL1 and 13,466 NL2 land-ice
  candidates, 9,845 and 9,833 have frequency A; the 27% without it are acquired at frequency B alone (mode
  `NASV`, `NADV`). A pair's polarization is HH if its reference lists HH at frequency A, otherwise VV.
- **S1BURST group.** `s1_burst_pairs.py` selects bursts by NISAR frame, as production does: attempt i draws
  one of `nisar_ice.jl`'s land-ice frames and a start date from `sha256(seed + i)` and asks ASF's burst
  search for the IW bursts intersecting the frame's footprint in the next 30 days. A pass and the next pass
  on the same track 6 or 12 days later form a pair over the burst IDs both share, in HH if both passes have
  it, else VV. 1000 pairs took 1416 attempts; they hold 1 to 57 bursts (median 15), 577 HH and 423 VV.
- **Checked.** `check_inputs.sh` requested the first byte of every input (ASF through Earthdata Login): all
  answered, after the Landsat remap. The NISAR lists were redrawn for frequency A only: 285 NL1 and 255 NL2
  pairs replaced the ones without it, and their 537 and 478 inputs all answered. The S1BURST inputs are not all requested, because each request makes
  the extractor generate the burst; a uniform sample of 20 of the 35,253 burst URLs all answered 202
  (generating).
