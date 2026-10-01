"""Capture the arguments the reference hands to `netCDF_packaging`, from inside the container.

`netCDF_packaging` is `hyp3_autorift.vend.netcdf_output`'s entry point for the uncropped product
write — the boundary `AutoRIFT.jl`'s `ItsLiveAutoRIFT` extension is meant to replace. Its own inputs
(`VX`, `VY`, the geogrid conversion coefficients, the reference velocity, the stable-surface mask, the
already-determined stable-shift scalars, `IMG_INFO_DICT`) are computed inside `testautoRIFT.py` and
never written to disk on their own — the uncropped file `netCDF_packaging` produces is deleted by
`process.py` once `crop.py` has cropped it, so there is no artifact to read after a normal run. This
captures them at the call itself, the same discipline `capture.py` uses one level down at
`runAutorift`.

Arrays are copied before the original function runs, since several arguments (`VX`, `VY`, `CHIPSIZEX`,
...) are mutated in place inside it — capturing a reference rather than a copy would record the
post-mutation state under the name of the input.
"""

import json
import os
import sys
from pathlib import Path

import numpy as np

import xchg

OUT = Path(os.environ.get('CAPTURE_DIR', '/home/ubuntu/work/capture_packaging'))

# Every 2-D array argument `netCDF_packaging` takes. `offset2vr`/`offset2va` are `None` for an
# optical pair and skipped; everything else is always a real array.
ARRAY_ARGS = (
    'VX', 'VY', 'DX', 'DY', 'INTERPMASK', 'CHIPSIZEX', 'CHIPSIZEY', 'SSM', 'SSM1', 'SX', 'SY',
    'offset2vx_1', 'offset2vx_2', 'offset2vy_1', 'offset2vy_2', 'offset2vr', 'offset2va',
    'scale_factor_1', 'scale_factor_2', 'MM', 'VXref', 'VYref', 'DXref', 'DYref',
)


def _scalar(v):
    if v is None or isinstance(v, (bool, int, float, str)):
        return v
    if isinstance(v, np.generic):
        return v.item()
    if isinstance(v, np.ndarray):
        return v.tolist()
    return str(v)


def _write(name, a):
    a = np.asarray(a)
    if a.dtype == np.bool_:
        a = a.view(np.uint8)
    if a.dtype.str not in xchg.TAGS:
        return {'skipped_dtype': a.dtype.str}
    xchg.write(str(OUT / name), a)
    return {'file': f'{name}.abx', 'dtype': a.dtype.str, 'shape': list(a.shape)}


def install():
    OUT.mkdir(parents=True, exist_ok=True)
    import hyp3_autorift.vend.netcdf_output as no

    original = no.netCDF_packaging
    state = {'n': 0}

    def patched(VX, VY, DX, DY, INTERPMASK, CHIPSIZEX, CHIPSIZEY, SSM, SSM1, SX, SY,
                offset2vx_1, offset2vx_2, offset2vy_1, offset2vy_2, offset2vr, offset2va,
                scale_factor_1, scale_factor_2, MM, VXref, VYref, DXref, DYref,
                rangePixelSize, azimuthPixelSize, dt, epsg, srs, tran, out_nc_filename,
                pair_type, detection_method, coordinates, IMG_INFO_DICT,
                stable_count, stable_count1, stable_shift_applied,
                dx_mean_shift, dy_mean_shift, dx_mean_shift1, dy_mean_shift1,
                error_vector, parameter_file):
        state['n'] += 1
        n = state['n']
        manifest = {'call': n, 'arrays': {}}

        args = dict(locals())
        for name in ARRAY_ARGS:
            v = args.get(name)
            if v is None:
                continue
            v = np.array(v, copy=True)
            if v.ndim != 2:
                continue
            manifest['arrays'][name] = _write(name, v)

        manifest['scalars'] = {
            'rangePixelSize': _scalar(rangePixelSize),
            'azimuthPixelSize': _scalar(azimuthPixelSize),
            'dt': _scalar(dt),
            'epsg': _scalar(epsg),
            'pair_type': pair_type,
            'detection_method': detection_method,
            'coordinates': coordinates,
            'stable_count': _scalar(stable_count),
            'stable_count1': _scalar(stable_count1),
            'stable_shift_applied': _scalar(stable_shift_applied),
            'dx_mean_shift': _scalar(dx_mean_shift),
            'dy_mean_shift': _scalar(dy_mean_shift),
            'dx_mean_shift1': _scalar(dx_mean_shift1),
            'dy_mean_shift1': _scalar(dy_mean_shift1),
            'error_vector': np.asarray(error_vector).tolist(),
            'parameter_file': parameter_file,
            'tran': [float(v) for v in tran],
            'srs_wkt': srs.ExportToWkt(),
            'srs_proj4': srs.ExportToProj4(),
            'srs_projection': srs.GetAttrValue('PROJECTION'),
            'out_nc_filename': out_nc_filename,
        }
        manifest['img_pair_info'] = {k: _scalar(v) for k, v in IMG_INFO_DICT.items()}

        with open(OUT / f'call{n}.json', 'w') as f:
            json.dump(manifest, f, indent=2)
        print(f'[capture_packaging] call {n}: {len(manifest["arrays"])} arrays -> {OUT}', flush=True)

        return original(VX, VY, DX, DY, INTERPMASK, CHIPSIZEX, CHIPSIZEY, SSM, SSM1, SX, SY,
                         offset2vx_1, offset2vx_2, offset2vy_1, offset2vy_2, offset2vr, offset2va,
                         scale_factor_1, scale_factor_2, MM, VXref, VYref, DXref, DYref,
                         rangePixelSize, azimuthPixelSize, dt, epsg, srs, tran, out_nc_filename,
                         pair_type, detection_method, coordinates, IMG_INFO_DICT,
                         stable_count, stable_count1, stable_shift_applied,
                         dx_mean_shift, dy_mean_shift, dx_mean_shift1, dy_mean_shift1,
                         error_vector, parameter_file)

    no.netCDF_packaging = patched
    print(f'[capture_packaging] netCDF_packaging patched; output -> {OUT}', flush=True)


def main():
    install()
    from hyp3_autorift.process import main as process_main

    sys.argv = ['hyp3_autorift'] + sys.argv[1:]
    process_main()


if __name__ == '__main__':
    main()
