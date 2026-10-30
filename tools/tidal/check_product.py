"""Check a tide-corrected Sentinel-1 product against an uncorrected run of the same pair, independently.

    tidal-env/bin/python tools/tidal/check_product.py <corrected.nc> <uncorrected.nc> <reference.SAFE> \
        <secondary.SAFE> <CATS2008_v2023.nc>
    tidal-env/bin/python tools/tidal/check_product.py <corrected.nc> <uncorrected.nc> <reference.SAFE> \
        <secondary.SAFE> <h_Gr1kmTM_v1> <grid_Gr1kmTM_v1> [--ib]

With one model file the model is CATS2008_v2023's netCDF and its flexure is applied; with two it is
Gr1kmTM's OTIS heights and grid, which have no flexure. Pixels within the 2 km taper are excluded: for
CATS, where the flexure is below 1; for Gr1kmTM, within 2 km of a BedMachine Greenland v7 pixel that is
neither ocean nor floating ice (`BEDMACHINE_MASK`).

Nothing here is shared with the Julia correction. The tide comes from pyTMD and the original model netCDF;
the incidence angle comes from ESA's own geolocation grid in the SAFE annotations, interpolated to each
pixel by position. The inverse barometer is reported separately: ERA5 mean sea level pressure read with
xarray from Google's ARCO-ERA5 Zarr store when `--ib` is given, left in the residual otherwise. Only the
range velocity `vr` and the range-pixel factor `dr_to_vr_factor` are read from the products.

For every corrected pixel the expected change in `vr` is

    Δvr = dr_to_vr_factor · w · (Δh_tide + Δh_IB) · cos θ / dr

with `w` the model's flexure and Δh_IB = −Δp / (ρ g), ρ = 1030 kg/m³. The script reports how the observed change compares: its slope against the
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


BEDMACHINE_MASK = "https://its-live-data.s3.amazonaws.com/tidal/BedMachineGreenland_v7_mask_0120m.tif"


def away_from_land(x, y, metres=2000.0):
    """Whether each point (EPSG:3413 metres) lies farther than `metres` from grounded ice or land."""
    import rasterio
    from rasterio.windows import from_bounds
    from scipy import ndimage
    with rasterio.open(BEDMACHINE_MASK) as r:
        pad = metres + 240
        w = from_bounds(x.min() - pad, y.min() - pad, x.max() + pad, y.max() + pad, r.transform)
        m = r.read(1, window=w)
        t = r.window_transform(w)
    land = (m == 1) | (m == 2)
    d = ndimage.distance_transform_edt(~land) * 120.0
    col = ((x - t.c) / t.a).astype(int)
    row = ((y - t.f) / t.e).astype(int)
    return d[row, col] > metres


def open_model(files):
    from pytmd_reference import open_cats, open_gr1km
    return (open_cats(files[0]), "TMD3") if len(files) == 1 else (open_gr1km(*files), "OTIS")


def ib_dh(lon, lat, t):
    import gcsfs  # noqa: F401  registers gs:// for xarray
    ds = xr.open_zarr("gs://gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3",
                      storage_options={"token": "anon"}, chunks=None)
    msl = ds["mean_sea_level_pressure"]
    pts = dict(longitude=xr.DataArray(np.mod(lon, 360.0), dims="p"), latitude=xr.DataArray(lat, dims="p"))
    p = [msl.sel(time=slice(np.datetime64(ti, "h") - np.timedelta64(1, "h"),
                            np.datetime64(ti, "h") + np.timedelta64(2, "h"))).interp(time=ti, **pts).values
         for ti in t]
    return -(p[1] - p[0]) / (1030.0 * 9.81)


def main(corrected, uncorrected, ref_safe, sec_safe, *model):
    use_ib = "--ib" in model
    model = tuple(m for m in model if m != "--ib")
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

    ds, corrections = open_model(model)
    mx, my = ds.tmd.transform_as(lon, lat, crs=4326)
    local = ds.tmd.interp(xr.DataArray(mx, dims="p"), xr.DataArray(my, dims="p"), method="linear")
    t = np.array([mid_time(ref_safe), mid_time(sec_safe)], dtype="datetime64[ms]")
    tide = timescale.from_datetime(t).tide
    h = (local.tmd.predict(tide, deltat=0.0, corrections=corrections) +
         local.tmd.infer(tide, deltat=0.0, corrections=corrections)).values
    dh = h[:, 1] - h[:, 0]
    w = np.asarray(local["flexure"])
    expected = k * dh * cos_t / dr
    observed = vr_a[iy, ix] - vr_b[iy, ix]

    # vr is stored as whole m/yr, so each observed change carries up to ±1 m/yr of rounding.
    full = np.isfinite(expected) & (w > 0.999)
    if corrections == "OTIS":
        full &= away_from_land(x[ix], y[iy])
    if use_ib:
        dib = ib_dh(lon[full], lat[full], t)
        print(f"Δh_IB (xarray ERA5): median {np.median(dib):+.3f} m")
        expected[full] += k * w[full] * dib * cos_t[full] / dr
    e, o = expected[full], observed[full]
    slope = np.sum(e * o) / np.sum(e * e)
    print(f"freely floating pixels compared: {full.sum()}")
    print(f"Δh_tide: median {np.median(dh[full]):+.3f} m, range {dh[full].min():+.3f} … {dh[full].max():+.3f} m")
    print(f"θ from ESA's geolocation grid: {inc[full].min():.2f}–{inc[full].max():.2f}°, dr {dr:.6f} m")
    print(f"expected Δvr ({'tide + IB' if use_ib else 'tide only'}): median {np.median(e):+.2f} m/yr")
    print(f"observed Δvr:             median {np.median(o):+.2f} m/yr")
    print(f"slope observed/expected {slope:.4f}, correlation {np.corrcoef(e, o)[0, 1]:.4f}")
    r = o - e
    print(f"residual ({'rounding' if use_ib else 'inverse barometer + rounding'}): median {np.median(r):+.2f}, "
          f"5–95% {np.percentile(r, 5):+.2f} … {np.percentile(r, 95):+.2f} m/yr")


if __name__ == "__main__":
    main(*sys.argv[1:])
