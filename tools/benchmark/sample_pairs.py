#!/usr/bin/env python3
"""
A uniform random sample of ITS_LIVE v02 image pairs per benchmark group, from the catalogs in
s3://its-live-data/catalog_geojson, with the inputs autorift_itslive takes for each.

Every pair gets a key, the first 8 bytes of sha256(seed + filename); the K pairs with the smallest keys in
a group are a uniform sample of it without replacement, the same on every run, and a chunk's local K
smallest are all a merge needs.

    sample_pairs.py OUTDIR [K] [WORKERS]

Writes OUTDIR/<group>.tsv, one pair per line in run_cases order (sensor, product, reference, secondary,
reference input, secondary input), OUTDIR/<group>_meta.tsv (product, deldays, EPSG, lon, lat,
percent_valid_pix, catalog directory) and OUTDIR/counts.tsv (pairs per group in the catalogs).
"""
import concurrent.futures as cf
import hashlib
import heapq
import json
import os
import re
import sys
import urllib.request

BASE = "https://its-live-data.s3.amazonaws.com/catalog_geojson"
CATALOGS = {"landsatOLI": "v2_landsatOLI_catalogs.txt", "sentinel1": "v2_sentinel1_catalogs.txt",
            "sentinel2": "v2_sentinel2_catalogs.txt"}
SEED = b"itslive-benchmark-2026-10-12"
LANDSAT_SENSOR = {"LT04": "L4", "LT05": "L5", "LE07": "L7", "LC08": "L8", "LC09": "L9"}
LANDSAT_GROUP = {"LT04": "L4L5", "LT05": "L4L5", "LE07": "L7", "LC08": "L8L9", "LC09": "L8L9"}
STAC = "https://landsatlook.usgs.gov/stac-server/collections/landsat-c2l1/items/"
S2_MIRROR = "https://storage.googleapis.com/gcp-public-data-sentinel-2/tiles"


def get(url, tries=4):
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return r.read()
        except Exception:
            if k == tries - 1:
                raise


def scenes(filename):
    ref, sec = filename.split("_X_")
    return ref, re.sub(r"_G\d{4}V\d{2}_P\d{3}\.nc$", "", sec)


def group_of(mission, ref):
    if mission == "sentinel1":
        return "S1"
    if mission == "sentinel2":
        return "S2"
    # LO08/LO09 (OLI without TIRS) and anything unrecognized are counted but not sampled.
    return LANDSAT_GROUP.get(ref[:4], "Landsat_other_" + ref[:4])


def chunk(args):
    """The K smallest-key pairs per group in one catalog chunk, and the chunk's count per group."""
    mission, url, k = args
    feats = json.loads(get(url))["features"]
    best, counts = {}, {}
    for f in feats:
        p = f["properties"]
        fn = p["filename"]
        ref, sec = scenes(fn)
        g = group_of(mission, ref)
        counts[g] = counts.get(g, 0) + 1
        key = int.from_bytes(hashlib.sha256(SEED + fn.encode()).digest()[:8], "big")
        lon, lat = p.get("center_lonlat", [None, None])
        rec = (key, fn, ref, sec, p.get("deldays"), p.get("data_epsg"), lon, lat,
               p.get("percent_valid_pix"), p.get("directory"))
        h = best.setdefault(g, [])
        # A max-heap of the K smallest keys, as negated keys.
        if len(h) < k:
            heapq.heappush(h, (-key, rec))
        elif -h[0][0] > key:
            heapq.heapreplace(h, (-key, rec))
    return {g: [r for _, r in h] for g, h in best.items()}, counts


def s2_band_url(name):
    """The B08.jp2 of Sentinel-2 product `name` on the Google mirror, from its manifest.safe."""
    tile = name.split("_")[5]
    safe = f"{S2_MIRROR}/{tile[1:3]}/{tile[3]}/{tile[4:6]}/{name}.SAFE"
    manifest = get(safe + "/manifest.safe").decode()
    # Only the band itself: from processing baseline 04.00 the manifest also lists QI_DATA masks named
    # after B08.
    hits = sorted(set(re.findall(r"GRANULE/[^\"<>]*/IMG_DATA/[^\"<>]*_B08\.jp2", manifest)))
    # A product from before December 2016 holds every tile of its datatake: its own is the one named.
    if len(hits) > 1:
        hits = [h for h in hits if tile + "_B08" in h]
    if len(hits) != 1:
        raise ValueError(f"{name}: expected one B08 entry in manifest.safe, found {len(hits)}")
    return f"{safe}/{hits[0]}"


def inputs(group, ref, sec):
    if group in ("L4L5", "L7", "L8L9"):
        return LANDSAT_SENSOR[ref[:4]], STAC + ref, STAC + sec
    if group == "S1":
        return "S1", f"https://datapool.asf.alaska.edu/SLC/S{ref[2]}/{ref}.zip", \
               f"https://datapool.asf.alaska.edu/SLC/S{sec[2]}/{sec}.zip"
    return "S2", s2_band_url(ref), s2_band_url(sec)


def main():
    out = sys.argv[1]
    k = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
    workers = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    os.makedirs(out, exist_ok=True)
    jobs = []
    for mission, listing in CATALOGS.items():
        urls = get(f"{BASE}/{mission}/{listing}").decode().split()
        jobs += [(mission, u, k) for u in urls]
    print(f"{len(jobs)} catalog chunks", flush=True)
    best, counts = {}, {}
    with cf.ProcessPoolExecutor(workers) as ex:
        for n, (b, c) in enumerate(ex.map(chunk, jobs, chunksize=1), 1):
            for g, recs in b.items():
                best[g] = heapq.nsmallest(k, best.get(g, []) + recs, key=lambda r: r[0])
            for g, v in c.items():
                counts[g] = counts.get(g, 0) + v
            if n % 100 == 0:
                print(f"  {n}/{len(jobs)} chunks", flush=True)
    with open(os.path.join(out, "counts.tsv"), "w") as f:
        for g in sorted(counts):
            f.write(f"{g}\t{counts[g]}\n")
    for g in ("L4L5", "L7", "L8L9", "S1", "S2"):
        recs = sorted(best.get(g, []), key=lambda r: r[0])
        with cf.ThreadPoolExecutor(16) as ex:
            ins = list(ex.map(lambda r: _safe_inputs(g, r), recs))
        with open(os.path.join(out, f"{g}.tsv"), "w") as run, open(os.path.join(out, f"{g}_meta.tsv"), "w") as meta:
            meta.write("product\tdeldays\tepsg\tlon\tlat\tpercent_valid_pix\tdirectory\n")
            for r, i in zip(recs, ins):
                if i is None:
                    continue
                product = r[1][:-3]
                run.write("\t".join([i[0], product, r[2], r[3], i[1], i[2]]) + "\n")
                meta.write("\t".join(str(x) for x in (product,) + r[4:]) + "\n")
        print(f"{g}: {sum(1 for i in ins if i is not None)} of {len(recs)} written", flush=True)


def _safe_inputs(g, r):
    try:
        return inputs(g, r[2], r[3])
    except Exception as e:
        print(f"  skipped {r[1]}: {e}", flush=True)
        return None


if __name__ == "__main__":
    main()
