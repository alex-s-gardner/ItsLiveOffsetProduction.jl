#!/usr/bin/env python3
"""
The NL1 and NL2 benchmark pairs: K of nisar_ice.jl's land-ice candidates that have frequency A, those with
the smallest sha256(seed + product) — a uniform sample, the same on every run — in run_cases order with ASF's
download URLs, and a _meta.tsv beside each.

    nisar_draw.py DIR [K]
"""
import hashlib
import sys

SEED = b"itslive-benchmark-2026-10-12"


def has_frequency_a(name):
    """The pipeline reads frequency A only. A granule's polarization code (field 9 of its name) is two letters
    for frequency A then two for B, and NA marks a frequency the mode does not acquire: NASV is B alone."""
    return name.split("_")[9][:2] != "NA"


def url(level, name):
    product = "RSLC" if level == "L1" else "GSLC"
    return f"https://nisar.asf.earthdatacloud.nasa.gov/NISAR/NISAR_{level}_{product}_BETA_V1/{name}/{name}.h5"


def main():
    d = sys.argv[1]
    k = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
    for level in ("L1", "L2"):
        rows = [line.rstrip("\n").split("\t") for line in open(f"{d}/nisar_{level}_ice.tsv")]
        total = len(rows)
        rows = [r for r in rows if has_frequency_a(r[1]) and has_frequency_a(r[2])]
        rows.sort(key=lambda r: hashlib.sha256(SEED + (r[1] + "_X_" + r[2]).encode()).digest())
        group = "NL" + level[1]
        with open(f"{d}/{group}.tsv", "w") as run, open(f"{d}/{group}_meta.tsv", "w") as meta:
            meta.write("product\tcycle_gap\tdeldays\tlon\tlat\tregion\troi_fraction\n")
            for r in rows[:k]:
                product = r[1] + "_X_" + r[2]
                run.write("\t".join([group, product, r[1], r[2], url(level, r[1]), url(level, r[2])]) + "\n")
                meta.write("\t".join([product, r[3], str(12 * int(r[3])), r[4], r[5], r[6], r[7]]) + "\n")
        print(f"{group}: {min(k, len(rows))} of {len(rows)} candidates with frequency A ({total} over land ice) written")


if __name__ == "__main__":
    main()
