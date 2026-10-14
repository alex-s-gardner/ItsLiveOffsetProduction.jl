#!/usr/bin/env python3
"""
The S1BURST benchmark pairs: Sentinel-1 bursts selected by NISAR frames, as production selects them.

The frames are the NISAR frames over land ice (nisar_ice.jl's candidates), each with its footprint polygon
from CMR. Attempt i draws a frame and a start date, both from sha256(seed + i), and asks ASF's burst search
for the IW bursts intersecting the frame in the 30 days from that date. The bursts are grouped into passes
(S1 track and acquisition day); a pass and the next pass on the same track (6 or 12 days later) form a
pair, over the burst IDs both share, in HH if both passes have it, else VV. Attempts continue until K pairs.

    s1_burst_pairs.py DIR [K] [WORKERS]

Reads DIR/nisar_L1_ice.tsv; writes DIR/S1BURST.tsv (run_cases order, inputs the bursts' comma-separated
extractor URLs) and DIR/S1BURST_meta.tsv.
"""
import collections
import concurrent.futures as cf
import datetime as dt
import hashlib
import json
import time
import sys
import urllib.parse
import urllib.request

SEED = b"itslive-benchmark-2026-10-12"
CMR = "https://cmr.earthdata.nasa.gov/search/granules.json"
ASF = "https://api.daac.asf.alaska.edu/services/search/param"
# The Sentinel-1 IW archive the dates are drawn from.
FIRST, LAST = dt.date(2014, 10, 3), dt.date(2026, 9, 1)


def h(*parts):
    return int.from_bytes(hashlib.sha256(SEED + "|".join(map(str, parts)).encode()).digest()[:8], "big")


def frame_key(name):
    f = name.split("_")
    return f[5], f[6], f[7]                      # track, direction, frame


def frame_polygons(keys):
    """Each frame's footprint, as WKT, from one RSLC granule of it in CMR."""
    after, out = None, {}
    while True:
        req = urllib.request.Request(f"{CMR}?short_name=NISAR_L1_RSLC_BETA_V1&page_size=2000")
        if after:
            req.add_header("CMR-Search-After", after)
        for k in range(6):
            try:
                with urllib.request.urlopen(req, timeout=120) as r:
                    after = r.headers.get("CMR-Search-After")
                    entries = json.load(r)["feed"]["entry"]
                break
            except Exception:
                if k == 5:
                    raise
                time.sleep(2 ** k)
        for e in entries:
            k = frame_key(e["title"])
            if k in keys and k not in out and e.get("polygons"):
                v = [float(x) for x in e["polygons"][0][0].split()]
                ring = ", ".join(f"{lon:.4f} {lat:.4f}" for lat, lon in zip(v[0::2], v[1::2]))
                out[k] = f"POLYGON(({ring}))"
        if not entries or not after:
            return out


def bursts(wkt, start, end):
    q = {"dataset": "SLC-BURST", "intersectsWith": wkt, "start": start.isoformat() + "T00:00:00Z",
         "end": end.isoformat() + "T23:59:59Z", "beamMode": "IW", "maxResults": "5000", "output": "geojson"}
    data = urllib.parse.urlencode(q).encode()
    for k in range(6):
        try:
            with urllib.request.urlopen(urllib.request.Request(ASF, data=data), timeout=180) as r:
                return [f["properties"] for f in json.load(r)["features"]]
        except Exception:
            if k == 5:
                raise
            time.sleep(2 ** k)


def attempt(args):
    i, frames = args
    keys = sorted(frames)
    key = keys[h("frame", i) % len(keys)]
    start = FIRST + dt.timedelta(days=h("date", i) % ((LAST - FIRST).days))
    try:
        found = bursts(frames[key], start, start + dt.timedelta(days=30))
    except Exception as e:
        print(f"  attempt {i} skipped: {e}", flush=True)
        return None
    passes = collections.defaultdict(list)                 # (track, day) -> bursts
    for b in found:
        passes[(b["pathNumber"], b["startTime"][:10])].append(b)
    bytrack = collections.defaultdict(list)
    for (track, day) in passes:
        bytrack[track].append(day)
    candidates = []
    for track, days in bytrack.items():
        days = sorted(days)
        for d1, d2 in zip(days, days[1:]):
            gap = (dt.date.fromisoformat(d2) - dt.date.fromisoformat(d1)).days
            if gap not in (6, 12):
                continue
            p1, p2 = passes[(track, d1)], passes[(track, d2)]
            for pol in ("HH", "VV"):
                ids1 = {b["burst"]["fullBurstID"]: b for b in p1 if b["polarization"] == pol}
                ids2 = {b["burst"]["fullBurstID"]: b for b in p2 if b["polarization"] == pol}
                common = sorted(set(ids1) & set(ids2), key=lambda s: (s.split("_")[2], int(s.split("_")[1])))
                if common:
                    candidates.append((track, d1, d2, gap, pol, [ids1[c] for c in common], [ids2[c] for c in common]))
                    break
            break                                          # the first pair on a track
    if not candidates:
        return None
    track, d1, d2, gap, pol, b1, b2 = candidates[h("track", i) % len(candidates)]
    scene = lambda b: b["url"].split("/")[3]
    return {"frame": "_".join(key), "track": track, "d1": d1, "d2": d2, "gap": gap, "pol": pol,
            "n": len(b1), "ref": scene(b1[0]), "sec": scene(b2[0]),
            "rin": ",".join(b["url"] for b in b1), "sin": ",".join(b["url"] for b in b2)}


def main():
    d = sys.argv[1]
    k = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
    workers = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    keys = {frame_key(line.split("\t")[1]) for line in open(f"{d}/nisar_L1_ice.tsv")}
    frames = frame_polygons(keys)
    print(f"{len(frames)} NISAR frames over land ice with footprints", flush=True)
    got, seen, i = [], set(), 0
    with cf.ThreadPoolExecutor(workers) as ex:
        while len(got) < k:
            batch = list(ex.map(attempt, [(j, frames) for j in range(i, i + 4 * workers)]))
            i += 4 * workers
            for r in batch:
                if r is None or (r["ref"], r["sec"], r["rin"]) in seen or len(got) >= k:
                    continue
                seen.add((r["ref"], r["sec"], r["rin"]))
                got.append(r)
            print(f"  {len(got)} pairs after {i} attempts", flush=True)
    with open(f"{d}/S1BURST.tsv", "w") as run, open(f"{d}/S1BURST_meta.tsv", "w") as meta:
        meta.write("product\tnisar_frame\ts1_track\tdate1\tdate2\tdeldays\tpolarization\tbursts\n")
        for r in got:
            product = f"{r['ref']}_X_{r['sec']}_{r['frame']}"
            run.write("\t".join(["S1BURST", product, r["ref"], r["sec"], r["rin"], r["sin"]]) + "\n")
            meta.write("\t".join(map(str, [product, r["frame"], r["track"], r["d1"], r["d2"], r["gap"],
                                           r["pol"], r["n"]])) + "\n")
    print(f"S1BURST: {len(got)} pairs from {i} attempts", flush=True)


if __name__ == "__main__":
    main()
