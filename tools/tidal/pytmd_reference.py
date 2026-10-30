"""Reference tide heights from pyTMD, for `app/test/tides.jl`.

    uv venv tidal-env && VIRTUAL_ENV=tidal-env uv pip install pyTMD
    tidal-env/bin/python tools/tidal/pytmd_reference.py synthetic app/test/data/tide_pytmd_synthetic.csv
    tidal-env/bin/python tools/tidal/pytmd_reference.py cats <CATS2008_v2023.nc> tools/tidal/tide_pytmd_cats.csv
    tidal-env/bin/python tools/tidal/pytmd_reference.py gr1km <h_Gr1kmTM_v1> <grid_Gr1kmTM_v1> tools/tidal/tide_pytmd_gr1km.csv

`synthetic` predicts from random complex constants of the ten CATS2008 constituents, so it checks the
harmonic sum, the nodal corrections and the inferred minor constituents, independent of any model file.
`cats` predicts from CATS2008_v2023 itself at points on Antarctic ice shelves, with the model's flexure
applied, so it also checks the reading of the converted model and the spatial interpolation.
`gr1km` predicts from Gr1kmTM v1's OTIS files at points on six Greenland ice tongues: pyTMD's OTIS path,
corrections "OTIS", which for these constituents are the TMD3 ones.
All use delta-T 0 and infer the minor constituents.
"""

import sys
import numpy as np
import pandas as pd
import xarray as xr
import timescale
import pyTMD.io
import pyTMD.io.OTIS
import pyTMD.io.dataset  # noqa: F401  registers the `tmd` accessor

CONSTITUENTS = ["m2", "s2", "n2", "k2", "k1", "o1", "p1", "q1", "mf", "mm"]

TIMES = np.array(
    [
        "2014-10-03T05:12:41.250",
        "2015-01-01T00:00:00",
        "2016-02-29T23:59:59",
        "2017-06-15T12:30:00",
        "2018-03-11T03:00:09.5",
        "2019-01-14T08:48:00",
        "2020-07-01T18:00:00",
        "2021-11-30T21:21:21",
        "2022-04-18T02:02:02",
        "2023-09-09T09:09:09",
        "2024-12-31T23:00:00",
        "2025-05-05T05:05:05",
        "2026-10-19T14:44:44",
    ],
    dtype="datetime64[ms]",
)

# (name, lon, lat): free-floating shelf interiors and points a few km seaward of grounding lines.
POINTS = [
    ("ronne_south", -60.0, -80.5),
    ("ronne_mid", -55.0, -78.0),
    ("filchner", -40.0, -79.0),
    ("ross_east", -160.0, -81.0),
    ("ross_west", 175.0, -80.0),
    ("ross_front", -175.0, -78.3),
    ("larsen_c", -62.0, -67.5),
    ("amery", 71.0, -70.5),
    ("pine_island", -100.5, -75.1),
    ("totten", 116.3, -67.0),
    ("brunt", -26.0, -75.5),
    ("getz", -125.0, -74.2),
]


# (name, lon, lat) on Greenland's floating tongues, from BedMachine Greenland v7's floating ice.
GREENLAND_POINTS = [
    ("79n", -20.9, 79.53),
    ("petermann", -60.9, 80.82),
    ("zachariae", -19.85, 78.95),
    ("ryder", -50.49, 81.70),
    ("storstrommen", -22.55, 76.71),
    ("steensby", -54.62, 81.53),
]

GR1KM_CRS = "+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +datum=WGS84 +units=km"


def tide_days(times):
    return np.asarray(timescale.from_datetime(times).tide, dtype=np.float64)


def predict(ds, t, corrections="TMD3"):
    major = ds.tmd.predict(t, deltat=0.0, corrections=corrections)
    minor = ds.tmd.infer(t, deltat=0.0, corrections=corrections)
    return major + minor


def synthetic(out):
    rng = np.random.default_rng(20261029)
    npoint = 6
    re = rng.normal(scale=0.4, size=(len(CONSTITUENTS), npoint))
    im = rng.normal(scale=0.4, size=(len(CONSTITUENTS), npoint))
    ds = xr.Dataset(
        {c: (("y", "x"), (re[i] + 1j * im[i])[None, :]) for i, c in enumerate(CONSTITUENTS)},
        coords=dict(x=np.arange(npoint, dtype=float), y=[0.0]),
    )
    t = tide_days(TIMES)
    h = predict(ds, t).values[0]  # (x, time)
    rows = []
    for p in range(npoint):
        for k, time in enumerate(TIMES):
            row = {"point": p, "time": str(time), "tide_days": t[k], "height": h[p, k]}
            for i, c in enumerate(CONSTITUENTS):
                row[c + "_re"] = re[i, p]
                row[c + "_im"] = im[i, p]
            rows.append(row)
    pd.DataFrame(rows).to_csv(out, index=False, float_format="%.17g")


def open_cats(model_file):
    # The reader `pyTMD.io.model.open_dataset` uses for a TMD3 file; it takes the CRS from the file.
    ds = pyTMD.io.OTIS.open_tmd3_dataset(model_file, group="z")
    for c in ds.tmd.constituents:
        ds[c] *= ds["flexure"]
    return ds


def open_gr1km(h_file, grid_file):
    ds = pyTMD.io.OTIS.open_otis_dataset(h_file, grid_file, group="z", crs=GR1KM_CRS)
    ds["flexure"] = xr.ones_like(ds["bathymetry"])
    return ds


def reference(ds, points, corrections, out):
    t = tide_days(TIMES)
    rows = []
    for name, lo, la in points:
        x, y = ds.tmd.transform_as(np.array([lo]), np.array([la]), crs=4326)
        local = ds.tmd.interp(xr.DataArray(x, dims="p"), xr.DataArray(y, dims="p"), method="linear")
        h = np.asarray(predict(local, t, corrections=corrections)).reshape(-1)
        flex = float(np.asarray(local["flexure"]).reshape(-1)[0])
        rows.extend(
            {"point": name, "lon": lo, "lat": la, "x_km": float(x[0]), "y_km": float(y[0]), "flexure": flex,
             "time": str(time), "tide_days": t[k], "height": h[k]}
            for k, time in enumerate(TIMES)
        )
    pd.DataFrame(rows).to_csv(out, index=False, float_format="%.17g")


if __name__ == "__main__":
    if sys.argv[1] == "synthetic":
        synthetic(sys.argv[2])
    elif sys.argv[1] == "cats":
        reference(open_cats(sys.argv[2]), POINTS, "TMD3", sys.argv[3])
    elif sys.argv[1] == "gr1km":
        reference(open_gr1km(sys.argv[2], sys.argv[3]), GREENLAND_POINTS, "OTIS", sys.argv[4])
    else:
        raise SystemExit("usage: pytmd_reference.py synthetic <out.csv> | cats <model.nc> <out.csv> | "
                         "gr1km <h_file> <grid_file> <out.csv>")
