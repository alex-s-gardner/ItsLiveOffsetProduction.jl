#!/usr/bin/env python3
"""
Candidate NISAR pairs from CMR, for the NL1 (RSLC) and NL2 (GSLC) benchmark groups.

A pair is two acquisitions of one frame — the same track, direction, frame, mode and polarization fields of
the granule name — 1 to MAXGAP cycles apart (12 days each), the later one the secondary. A scene with
several versions (composite release ID, counter) keeps the last by name. Writes OUTDIR/nisar_<level>_
candidates.tsv: level, reference, secondary, cycle gap, footprint centre lon, lat.

    nisar_pairs.py OUTDIR [MAXGAP]
"""
import collections
import json
import os
import sys
import urllib.request

CMR = "https://cmr.earthdata.nasa.gov/search/granules.json"
COLLECTIONS = {"L1": "NISAR_L1_RSLC_BETA_V1", "L2": "NISAR_L2_GSLC_BETA_V1"}


def granules(short_name):
    after, out = None, []
    while True:
        req = urllib.request.Request(f"{CMR}?short_name={short_name}&page_size=2000")
        if after:
            req.add_header("CMR-Search-After", after)
        with urllib.request.urlopen(req, timeout=120) as r:
            after = r.headers.get("CMR-Search-After")
            entries = json.load(r)["feed"]["entry"]
        out += entries
        if not entries or not after:
            return out


def centre(e):
    pts = []
    for ring in e.get("polygons", []):
        v = [float(x) for x in ring[0].split()]
        pts += list(zip(v[1::2], v[0::2]))           # CMR polygons are "lat lon lat lon ..."
    for b in e.get("boxes", []):
        s, w, n, ea = (float(x) for x in b.split())
        pts += [(w, s), (ea, n)]
    if not pts:
        return None
    # Longitude averaged on the unit circle, so a footprint across the antimeridian is centred on it.
    import math
    x = sum(math.cos(math.radians(lo)) for lo, _ in pts) / len(pts)
    y = sum(math.sin(math.radians(lo)) for lo, _ in pts) / len(pts)
    return math.degrees(math.atan2(y, x)), sum(la for _, la in pts) / len(pts)


def main():
    out = sys.argv[1]
    maxgap = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    os.makedirs(out, exist_ok=True)
    for level, short in COLLECTIONS.items():
        entries = granules(short)
        scenes = {}
        for e in entries:
            name = e["title"].removesuffix(".h5") if hasattr(str, "removesuffix") else e["title"]
            f = name.split("_")
            # NISAR_L1_PR_RSLC_<cycle>_<track>_<dir>_<frame>_<mode>_<pol>_..._<CRID>_..._<counter>
            key = (int(f[4]), f[5], f[6], f[7], f[8], f[9])
            if key not in scenes or name > scenes[key][0]:
                scenes[key] = (name, centre(e))
        frames = collections.defaultdict(dict)
        for (cycle, *frame), v in scenes.items():
            frames[tuple(frame)][cycle] = v
        n = 0
        with open(os.path.join(out, f"nisar_{level}_candidates.tsv"), "w") as fo:
            for frame, bycycle in frames.items():
                for c1, (r, cen) in bycycle.items():
                    for gap in range(1, maxgap + 1):
                        s = bycycle.get(c1 + gap)
                        if s is None or cen is None:
                            continue
                        fo.write(f"{level}\t{r}\t{s[0]}\t{gap}\t{cen[0]:.4f}\t{cen[1]:.4f}\n")
                        n += 1
        print(f"{level}: {len(entries)} granules, {len(scenes)} scenes, {len(frames)} frames, {n} candidate pairs",
              flush=True)


if __name__ == "__main__":
    main()
