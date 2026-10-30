"""Check a tide-corrected Sentinel-1 product against an uncorrected run of the same pair, independently.

    tidal-env/bin/python tools/tidal/check_product.py <corrected.nc> <uncorrected.nc> <reference.SAFE> \
        <secondary.SAFE> <CATS2008_v2023.nc>

Nothing here is shared with the Julia correction. The tide comes from pyTMD and the original model netCDF;
the incidence angle comes from ESA's own geolocation grid in the SAFE annotations, interpolated to each
pixel by position; the inverse barometer is left out, so it shows as the residual. Only the range
velocity `vr` and the range-pixel factor `dr_to_vr_factor` are read from the products.

For every corrected pixel the expected change in `vr` is

    Δvr = dr_to_vr_factor · w · Δh_tide · cos θ / dr

with `w` the model's flexure. The script reports how the observed change compares: its slope against the
expected one, their correlation, and the residual. Pixels inside the 2 km grounding-line taper are
excluded, since this script does not rebuild the taper.
"""

import sys
import glob
import numpy as np
import netCDF4
import pyproj
import timescale
import xml.etree.ElementTree as ET
from scipy.interpolate import LinearNDInterpolator
import pyTMD.io.OTIS
import pyTMD.io.dataset  # noqa: F401  registers the `tmd` accessor
import xarray as xr


def geolocation_grid(safe):
    pts = []
    for f in sorted(glob.glob(safe + "/annotation/s1*-iw*-slc-*.xml")):
        root = ET.parse(f).getroot()
        for p in root.iter("geolocationGridPoint"):
            pts.append([float(p.find(k).text) for k in ("longitude", "latitude", "incidenceAngle")])
        dr = float(root.find(".//rangePixelSpacing").text)
    return np.array(pts), dr


def mid_time(safe):
    starts, stops = [], []
    for f in glob.glob(safe + "/annotation/s1*-iw*-slc-*.xml"):
        root = ET.parse(f).getroot()
        starts.append(np.datetime64(root.find(".//productFirstLineUtcTime").text))
        stops.append(np.datetime64(root.find(".//productLastLineUtcTime").text))
    return min(starts) + (max(stops) - min(starts)) / 2


def main(corrected, uncorrected, ref_safe, sec_safe, model):
    a, b = netCDF4.Dataset(corrected), netCDF4.Dataset(uncorrected)
    vr_a = np.ma.filled(a["vr"][0].astype(float), np.nan)
    vr_b = np.ma.filled(b["vr"][0].astype(float), np.nan)
    k = float(a["M11"].dr_to_vr_factor)
    x, y = a["x"][:], a["y"][:]
    epsg = int(a["mapping"].spatial_epsg)
    changed = np.isfinite(vr_a) & np.isfinite(vr_b) & (vr_a != vr_b)
    iy, ix = np.nonzero(changed)
    print(f"pixels changed: {iy.size}")

    lon, lat = pyproj.Transformer.from_crs(epsg, 4326, always_xy=True).transform(x[ix], y[iy])
    grid, dr = geolocation_grid(ref_safe)
    inc = LinearNDInterpolator(grid[:, :2], grid[:, 2])(lon, lat)
    cos_t = np.cos(np.radians(inc))

    ds = pyTMD.io.OTIS.open_tmd3_dataset(model, group="z")
    flex = ds["flexure"]
    for c in ds.tmd.constituents:
        ds[c] *= ds["flexure"]
    mx, my = ds.tmd.transform_as(lon, lat, crs=4326)
    local = ds.tmd.interp(xr.DataArray(mx, dims="p"), xr.DataArray(my, dims="p"), method="linear")
    t = np.array([mid_time(ref_safe), mid_time(sec_safe)], dtype="datetime64[ms]")
    tide = timescale.from_datetime(t).tide
    h = (local.tmd.predict(tide, deltat=0.0, corrections="TMD3") +
         local.tmd.infer(tide, deltat=0.0, corrections="TMD3")).values
    dh = h[:, 1] - h[:, 0]
    w = np.asarray(local["flexure"])
    expected = k * dh * cos_t / dr
    observed = vr_a[iy, ix] - vr_b[iy, ix]

    # vr is stored as whole m/yr, so each observed change carries up to ±1 m/yr of rounding.
    full = np.isfinite(expected) & (w > 0.999)
    e, o = expected[full], observed[full]
    slope = np.sum(e * o) / np.sum(e * e)
    print(f"freely floating pixels compared: {full.sum()}")
    print(f"Δh_tide: median {np.median(dh[full]):+.3f} m, range {dh[full].min():+.3f} … {dh[full].max():+.3f} m")
    print(f"θ from ESA's geolocation grid: {inc[full].min():.2f}–{inc[full].max():.2f}°, dr {dr:.6f} m")
    print(f"expected Δvr (tide only): median {np.median(e):+.2f} m/yr")
    print(f"observed Δvr:             median {np.median(o):+.2f} m/yr")
    print(f"slope observed/expected {slope:.4f}, correlation {np.corrcoef(e, o)[0, 1]:.4f}")
    r = o - e
    print(f"residual (inverse barometer + rounding): median {np.median(r):+.2f}, "
          f"5–95% {np.percentile(r, 5):+.2f} … {np.percentile(r, 95):+.2f} m/yr")


if __name__ == "__main__":
    main(*sys.argv[1:6])
