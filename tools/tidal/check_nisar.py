"""Check tide-corrected NISAR products against uncorrected runs of the same pair, independently.

    tidal-env/bin/python tools/tidal/check_nisar.py <nl2_tide.nc> <nl2_none.nc> <early GSLC.h5 (local or
        https)> <late GSLC.h5> <CATS2008_v2023.nc> [<nl1_tide.nc> <nl1_none.nc>]

NISAR L2 (GSLC). The expected change in the map velocities comes from JPL's own geometry: the GSLC's
`radarGrid` cube of line-of-sight unit vectors and incidence angles (ISCE3 `metadataCubes`), at the
ellipsoid-height slice nearest the surface. Its line of sight is the native-Doppler one, which squints
by about 1° along track; the zero-Doppler direction the image is focused in is that vector with its
along-track part removed, which is what is compared. A surface that rose by Δh is geocoded Δh·cot θ along
the horizontal line of sight toward the sensor, so

    Δ(vx, vy) = −Δh · cos θ / sin² θ · (los_x', los_y') / Δt_years

with (los_x', los_y') the zero-Doppler horizontal line of sight, rotated from east/north into the
product's polar stereographic x/y by the meridian convergence, and scaled by the projection's scale factor,
since the GSLC grid and the product's velocities are in grid metres. The tide is pyTMD's on the original
model with its flexure; the inverse barometer is ERA5 mean sea level pressure read with xarray from
Google's ARCO-ERA5 Zarr store, −Δp/(ρ g) with ρ = 1030 kg/m³, also scaled by the flexure.

NISAR L1 (RSLC), when given: `Δvr` against `dr_to_vr_factor · Δh · cos θ / dr` with θ from the same cube,
as `check_product.py` does for Sentinel-1.

Pixels inside the 2 km grounding-line taper are excluded; this script does not rebuild it.
"""

import sys
import numpy as np
import netCDF4
import h5py
import fsspec
import pyproj
import timescale
import xarray as xr
import pyTMD.io.OTIS
import pyTMD.io.dataset  # noqa: F401  registers the `tmd` accessor


def band(path, v):
    ds = netCDF4.Dataset(path)
    ds.set_auto_mask(False)
    a = ds[v][0].astype(float)
    a[a == -32767] = np.nan
    return a, ds


def open_h5(uri):
    if uri.startswith("http"):
        import netrc
        login, _, password = netrc.netrc().authenticators("urs.earthdata.nasa.gov")
        fs = fsspec.filesystem("https", client_kwargs={"auth": __import__("aiohttp").BasicAuth(login, password)})
        return h5py.File(fs.open(uri, block_size=8 << 20), "r")
    return h5py.File(uri, "r")


def mid_time(h):
    i = h["science/LSAR/identification"]
    t0 = np.datetime64(i["zeroDopplerStartTime"][()].decode()[:26])
    t1 = np.datetime64(i["zeroDopplerEndTime"][()].decode()[:26])
    return t0 + (t1 - t0) / 2


def cube_at(h, x, y):
    """JPL's zero-Doppler horizontal LOS (east, north) and cos θ, bilinear at product x, y (EPSG of the cube)."""
    g = h["science/LSAR/GSLC/metadata/radarGrid"]
    hs = g["heightAboveEllipsoid"][:]
    k = int(np.argmin(np.abs(hs)))
    cx, cy = g["xCoordinates"][:], g["yCoordinates"][:]
    lx, ly = g["losUnitVectorX"][k], g["losUnitVectorY"][k]
    ax, ay = g["alongTrackUnitVectorX"][k], g["alongTrackUnitVectorY"][k]
    inc = g["incidenceAngle"][k]
    fi = np.interp(y, cy[::-1], np.arange(len(cy))[::-1])
    fj = np.interp(x, cx, np.arange(len(cx)))

    def bil(a):
        i0 = np.clip(np.floor(fi).astype(int), 0, len(cy) - 2)
        j0 = np.clip(np.floor(fj).astype(int), 0, len(cx) - 2)
        di, dj = fi - i0, fj - j0
        return ((1 - di) * (1 - dj) * a[i0, j0] + (1 - di) * dj * a[i0, j0 + 1] +
                di * (1 - dj) * a[i0 + 1, j0] + di * dj * a[i0 + 1, j0 + 1])

    le, ln, c = bil(lx), bil(ly), np.cos(np.radians(bil(inc)))
    te, tn = bil(ax), bil(ay)
    tn_norm = np.hypot(te, tn)
    te, tn = te / tn_norm, tn / tn_norm
    # Remove the along-track component of the LOS (horizontal part; the along-track vector is horizontal).
    s = le * te + ln * tn
    le, ln = le - s * te, ln - s * tn
    # Renormalize the 3-vector (le, ln, c) to unit length.
    norm = np.sqrt(le**2 + ln**2 + c**2)
    return le / norm, ln / norm, c / norm


def tide_dh(model, lon, lat, t):
    ds = pyTMD.io.OTIS.open_tmd3_dataset(model, group="z")
    for c in ds.tmd.constituents:
        ds[c] *= ds["flexure"]
    mx, my = ds.tmd.transform_as(lon, lat, crs=4326)
    local = ds.tmd.interp(xr.DataArray(mx, dims="p"), xr.DataArray(my, dims="p"), method="linear")
    tide = timescale.from_datetime(np.array(t, dtype="datetime64[ms]")).tide
    h = (local.tmd.predict(tide, deltat=0.0, corrections="TMD3") +
         local.tmd.infer(tide, deltat=0.0, corrections="TMD3")).values
    return h[:, 1] - h[:, 0], np.asarray(local["flexure"])


def ib_dh(lon, lat, t, flex):
    import gcsfs  # noqa: F401  registers gs:// for xarray
    ds = xr.open_zarr("gs://gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3",
                      storage_options={"token": "anon"}, chunks=None)
    msl = ds["mean_sea_level_pressure"]
    pts = dict(longitude=xr.DataArray(np.mod(lon, 360.0), dims="p"), latitude=xr.DataArray(lat, dims="p"))
    p = [msl.sel(time=slice(np.datetime64(ti, "h") - np.timedelta64(1, "h"), np.datetime64(ti, "h") + np.timedelta64(2, "h")))
            .interp(time=ti, **pts).values for ti in t]
    return -(p[1] - p[0]) / (1030.0 * 9.81) * flex


def report(name, e, o):
    slope = np.sum(e * o) / np.sum(e * e)
    r = o - e
    # Means, not medians or slopes, are what integer rounding leaves unbiased when the expected values span
    # well over 1 m/yr; their difference is the systematic error.
    print(f"  {name}: expected mean {np.mean(e):+.3f} (5–95% {np.percentile(e, 5):+.2f} … {np.percentile(e, 95):+.2f}), "
          f"observed mean {np.mean(o):+.3f} m/yr, difference {np.mean(o) - np.mean(e):+.3f}")
    print(f"  {name}: expected median {np.median(e):+.2f}, observed median {np.median(o):+.2f} m/yr; "
          f"slope {slope:.4f}; residual median {np.median(r):+.2f}, 5–95% {np.percentile(r, 5):+.2f} … "
          f"{np.percentile(r, 95):+.2f} m/yr (rounding to 1 m/yr)")


def main(nl2_tide, nl2_none, early, late, model, nl1_tide=None, nl1_none=None):
    he, hl = open_h5(early), open_h5(late)
    t = [mid_time(he), mid_time(hl)]
    dt_years = (t[1] - t[0]) / np.timedelta64(1, "s") / (365 * 86400)
    print(f"mid times {t[0]} → {t[1]}, Δt {dt_years * 365:.4f} d")

    vxa, a = band(nl2_tide, "vx")
    vxb, b = band(nl2_none, "vx")
    vya, _ = band(nl2_tide, "vy")
    vyb, _ = band(nl2_none, "vy")
    changed = np.isfinite(vxa) & np.isfinite(vxb) & ((vxa != vxb) | (vya != vyb))
    iy, ix = np.nonzero(changed)
    x, y = a["x"][:][ix], a["y"][:][iy]
    epsg = int(a["mapping"].spatial_epsg)
    lon, lat = pyproj.Transformer.from_crs(epsg, 4326, always_xy=True).transform(x, y)
    print(f"NL2: {iy.size} pixels changed; grid EPSG:{epsg}")
    dh, w = tide_dh(model, lon, lat, t)
    dh = dh + ib_dh(lon, lat, t, w)
    full = w > 0.999
    le, ln, c = cube_at(he, x, y)
    s2 = le**2 + ln**2
    # East/north → polar stereographic x/y: rotate by the meridian convergence (EPSG:3031, lon_0 = 0).
    gam = np.radians(lon) if epsg == 3031 else np.zeros_like(lon)
    gx = le * np.cos(gam) + ln * np.sin(gam)
    gy = -le * np.sin(gam) + ln * np.cos(gam)
    # The GSLC is posted on the map grid, so the correlator measures grid metres: a true ground shift
    # appears scaled by the projection's scale factor (0.978 at 81°S in EPSG:3031), and so do the
    # product's velocities.
    kmap = pyproj.Proj(f"EPSG:{epsg}").get_factors(lon, lat).parallel_scale
    k = -dh * c / s2 / dt_years * kmap
    ex, ey = k * gx, k * gy
    ox, oy = (vxa - vxb)[iy, ix], (vya - vyb)[iy, ix]
    print(f"  freely floating pixels compared: {full.sum()}; Δh (tide + IB) median {np.median(dh[full]):+.3f} m; "
          f"θ {np.degrees(np.arccos(c[full])).min():.2f}–{np.degrees(np.arccos(c[full])).max():.2f}°")
    report("Δvx", ex[full], ox[full])
    report("Δvy", ey[full], oy[full])

    if nl1_tide:
        vra, a1 = band(nl1_tide, "vr")
        vrb, _ = band(nl1_none, "vr")
        kf = float(a1["M11"].dr_to_vr_factor)
        dr = float(he["science/LSAR/GSLC/metadata/sourceData/swaths/frequencyA/slantRangeSpacing"][()])
        changed = np.isfinite(vra) & np.isfinite(vrb) & (vra != vrb)
        iy, ix = np.nonzero(changed)
        x, y = a1["x"][:][ix], a1["y"][:][iy]
        lon, lat = pyproj.Transformer.from_crs(int(a1["mapping"].spatial_epsg), 4326, always_xy=True).transform(x, y)
        print(f"NL1: {iy.size} pixels changed; dr {dr:.4f} m, dr_to_vr_factor {kf:.4f}")
        dh, w = tide_dh(model, lon, lat, t)
        dh = dh + ib_dh(lon, lat, t, w)
        _, _, c = cube_at(he, x, y)
        full = w > 0.999
        report("Δvr", (kf * dh * c / dr)[full], (vra - vrb)[iy, ix][full])


if __name__ == "__main__":
    main(*sys.argv[1:])
