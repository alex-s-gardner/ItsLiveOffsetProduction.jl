"""Capture the arrays the reference correlator is actually handed, from inside the container.

Driven by `tools/golden/intermediate.jl`, which mounts this file and `tools/ab/xchg.py` into the
container and runs it in place of the ordinary entry point.

**Why intercept rather than reconstruct.** The correlator's inputs are not the scene files. By the
time `runAutorift` is called the pair has been read, high-pass filtered, nodata-filled and cast, and
the grid has been snapped to the even-chip convention. Rebuilding that chain in Julia to feed
AutoRIFT.jl the "same" arrays would mean reimplementing the code under test, so a preprocessing
difference would surface as a correlator difference. Taking the arrays at the call boundary makes the
comparison a comparison of the correlator alone — the same discipline `tools/ab` uses in handing the
reference arrays AutoRIFT.jl already filtered.

Arrays are written through `tools/ab/xchg.py`, so the element type and both dimensions travel in the
file and the Julia reader takes the layout from the header rather than from a caller-supplied shape.
These grids are square, and a square array read with the wrong convention is silently transposed.

**The argument-order trap.** `arImgDisp_s(a, b)` cuts its chip from `b` and its search window from
`a`: `a` binds to `I1`, but the body calls the C++ as `(I2.ravel(), I1.ravel())` and the C++ binds its
first array to `sec_img`, where the chip is taken. Two swaps that compose rather than cancel. Arrays
here are named for the attribute they came from (`I1`, `I2`), so a reader must apply that mapping
itself; `tools/ab/README.md` states it and `tools/ab/stage2_python.py` asserts it.
"""

import json
import os
import sys
import types
from pathlib import Path

import numpy as np

import xchg

OUT = Path(os.environ.get('CAPTURE_DIR', '/home/ubuntu/work/capture'))

# Scalars that change the answer. Captured so the Julia side configures `Params` from what the
# reference actually used rather than from what the driver is believed to set — `OverSampleRatio`
# alone is a per-chip-size dictionary assembled at run time (`testautoRIFT.py:488-510`).
SCALARS = (
    'ChipSize0X', 'ChipSizeMinX', 'ChipSizeMaxX', 'ScaleChipSizeY',
    'GridSpacingX', 'SkipSampleX', 'SkipSampleY',
    'OverSampleRatio', 'WallisFilterWidth', 'MultiThread',
    'BuffDistanceC', 'CoarseCorCutoff', 'sparseSearchSampleRate',
    # `DataType` decides which branch of `uniform_data_type` ran and therefore which correlator the
    # pyramid reached: 0 quantizes to 256 levels and dispatches to `arImgDisp_u`, 1 keeps the filtered
    # `Float32` field and dispatches to `arImgDisp_s`. Recorded so a reader takes the path from the
    # capture rather than from the directory it happened to find it in — the two captures are otherwise
    # indistinguishable except by the element type of `in_I1`.
    'DataType',
    'DataTypeInput', 'ChipSizeMaxXInput', 'preproc_filt_width',
    # `minSearch` is the floor the level loop applies to every nonzero search limit
    # (`autoRIFT.py:598-602`), so the radius the correlator receives is not the radius this capture
    # records in `in_SearchLimitX` — a point asking for 1 searches at 6. Without this scalar the Julia
    # side cannot reconstruct what the correlator was actually handed.
    'minSearch',
    # The outlier filter and the hole fill. `autorift()` derives `DispFiltC`/`DispFiltF` from
    # `FracValid`, `FracSearch`, `FiltWidth`, `Iter` and `MadScalar` (`autoRIFT.py:484-505`), and the
    # three-pass fill uses `fillFiltWidth` for its median (`:775-790`). These are the stages that own
    # the coverage difference, so comparing them against a *believed* parameter set would attribute a
    # parameter mismatch to a reducer or to the merge.
    'FracValid', 'FracSearch', 'FiltWidth', 'Iter', 'MadScalar',
    'fillFiltWidth', 'colfiltChunkSize', 'StandardDeviationCutoff',
)

# Arrays to take before the call: the filtered pair, the grid, the priors, the per-point limits.
INPUTS = ('I1', 'I2', 'xGrid', 'yGrid', 'Dx0', 'Dy0',
          'SearchLimitX', 'SearchLimitY', 'ChipSizeMinX', 'ChipSizeMaxX', 'zeroMask')

# And after: the answer.
OUTPUTS = ('Dx', 'Dy', 'InterpMask', 'ChipSizeX')


def _write(name, a):
    xchg.write(str(OUT / name), a)
    return {'file': f'{name}.abx', 'dtype': a.dtype.str, 'shape': list(a.shape)}


def _scalar(v):
    """A JSON-representable form of an autoRIFT attribute, or None for an array."""
    if v is None or isinstance(v, (bool, int, float, str)):
        return v
    if isinstance(v, dict):
        return {str(k): _scalar(x) for k, x in v.items()}
    if isinstance(v, np.generic):
        return v.item()
    if isinstance(v, np.ndarray):
        return None
    return str(v)


def _dump(obj, names, prefix, manifest):
    for name in names:
        v = getattr(obj, name, None)
        if not isinstance(v, np.ndarray) or v.ndim != 2:
            continue
        # `InterpMask` and `zeroMask` are boolean. NumPy stores those one byte per element, so a view
        # as uint8 is the same bytes with a name xchg's type table already has — no value changes.
        # Every other dtype outside the table is recorded as skipped rather than cast, since a cast
        # that loses precision would silently change what is being compared.
        if v.dtype == np.bool_:
            v = v.view(np.uint8)
        if v.dtype.str not in xchg.TAGS:
            manifest['skipped'][f'{prefix}{name}'] = v.dtype.str
            continue
        manifest['arrays'][f'{prefix}{name}'] = _write(f'{prefix}{name}', v)


def record_filtered_scenes(manifest):
    """Record the filtered scenes `process.py` wrote, so a reader takes them from the manifest.

    **Landsat 4/5 and 7 are not filtered inside `autorift()`.** `process.py:312-336` dispatches on the
    platform and applies the filter to the *native* scenes before geogrid runs — `apply_fft_filter` for
    L4/5 and `apply_wallis_nodata_fill_filter` for L7 and any L8 paired with an L7 — writing `Float32`
    GeoTIFFs to `Path.cwd()/'filtered'` (`create_filtered_filepath`, `:252`). The `FIXME` there says why:
    the FFT filter locates the scene corners, and geogrid rounds and chops corners when subsetting to the
    common overlap, so the filter has to see the native footprint.

    That means the byte arrays this capture takes at the `runAutorift` boundary are **downstream of a
    filter whose output is a file on disk**, and comparing against them alone cannot separate a filter
    difference from a correlator one. Recording the paths makes those scenes the input a Julia-side
    comparison can be fed, which is what turns the L4/5 and L7 cases from a product diff into a staged one.

    Only the properties go in the manifest, not the arrays: they are full scenes at native resolution, and
    the GeoTIFFs are already on disk beside the capture. Shape, dtype and finite count are enough for a
    reader to assert it opened the file the capture meant.
    """
    filtered = manifest.setdefault('filtered', {})
    directory = Path.cwd() / 'filtered'
    if not directory.is_dir():
        # Not an error: an `hps` pair is filtered inside `autorift()` and writes nothing here. The
        # absence is the fact, and a reader distinguishes "no filtered scenes" from "not looked for".
        manifest['filtered_dir'] = None
        print('[capture] no filtered/ directory: this pair is filtered inside autorift()', flush=True)
        return

    manifest['filtered_dir'] = str(directory)
    try:
        from osgeo import gdal

        gdal.UseExceptions()
    except ImportError:
        gdal = None

    for path in sorted(directory.iterdir()):
        if path.suffix.lower() not in ('.tif', '.tiff'):
            continue
        rec = {'file': path.name, 'bytes': path.stat().st_size}
        if gdal is not None:
            ds = gdal.Open(str(path))
            band = ds.GetRasterBand(1)
            a = band.ReadAsArray()
            rec.update(shape=[int(ds.RasterYSize), int(ds.RasterXSize)],
                       dtype=str(a.dtype),
                       finite=int(np.count_nonzero(np.isfinite(a))),
                       nonzero=int(np.count_nonzero(a)),
                       geotransform=[float(v) for v in ds.GetGeoTransform()])
            ds = None
        filtered[path.stem] = rec
        print(f'[capture] filtered scene {path.name}: {rec.get("shape")} {rec.get("dtype")}', flush=True)


def reference_module():
    """The `autoRIFT.autoRIFT` module, which is where `arImgDisp_*`, `DISP_FILT` and `colfilt` live.

    `from autoRIFT import autoRIFT` does **not** give this. The package's `__init__` binds that name to
    the `autoRIFT` *class*, so the import returns a class that has none of these symbols — and
    `setattr` on it then patches an attribute nothing calls, while reporting success. Every capture
    taken that way carries `levels: []`, which reads as "the reference resolved no levels" rather than
    as a harness fault.

    `sys.modules` is the reliable route: the module is already imported by the time anything here runs,
    and its dotted name is unambiguous where the attribute lookup is not.
    """
    mod = sys.modules.get('autoRIFT.autoRIFT')
    if mod is None:
        import autoRIFT.autoRIFT  # noqa: F401  -- imported for its side effect on sys.modules

        mod = sys.modules['autoRIFT.autoRIFT']
    if not isinstance(mod, types.ModuleType):
        raise TypeError(f'autoRIFT.autoRIFT resolved to {type(mod)}, not a module')
    return mod


def install_levels(module, manifest):
    """Patch the per-level correlator and filter calls to dump what each pyramid level decided.

    `runAutorift`'s own inputs and outputs describe the *merged* answer, so which level a
    disagreement came from has to be inferred from the reported `ChipSizeX`. That inference is
    weakest exactly where it matters — a point the two implementations assign to different levels
    is the case under investigation, and it is the case where the merged output cannot say which
    level's decision diverged.

    The pyramid's per-level state is local to `runAutorift`, so it is unreachable from outside. What
    *is* reachable is the calls that loop makes: `arImgDisp_u`/`arImgDisp_s` once for the coarse pass
    and once for the fine pass of each level, and `DISP_FILT.filtDisp` once per pass. Wrapping those
    records each level's raw measurement and each level's rejection mask without reimplementing the
    loop, so the recorded values are the reference's own rather than a reconstruction.

    Calls are numbered in the order they happen. `ChipSizeX` is an argument to the correlator, so the
    level is recorded rather than deduced, and the coarse and fine passes are distinguished by
    `SubPixFlag` — the coarse pass runs with it `False`.

    `module` must be the `autoRIFT.autoRIFT` *module*, which is what `autorift()` resolves these names
    against. See `reference_module`.
    """
    seq = {'n': 0}
    levels = manifest.setdefault('levels', [])

    def wrap_corr(name):
        # `raise`, not a silent return. A missing symbol means this patch cannot observe anything, and
        # the run that follows still writes a product, a log and an exit status of 0 — so a quiet
        # skip is indistinguishable from a level that genuinely never ran, which is exactly the
        # reading a wrong `module` produces.
        if not hasattr(module, name):
            raise AttributeError(
                f'{module.__name__} has no {name}; the per-level patch cannot be installed. '
                'This is the symptom of patching the autoRIFT *class* rather than the module '
                'that defines these functions — see `reference_module`.')
        original = getattr(module, name)

        def patched(I1, I2, xGrid, yGrid, ChipSizeX, ChipSizeY, SearchLimitX, SearchLimitY,
                    Dx0, Dy0, SubPixFlag, overSampleRatio, *rest):
            seq['n'] += 1
            n = seq['n']
            dx, dy = original(I1, I2, xGrid, yGrid, ChipSizeX, ChipSizeY, SearchLimitX,
                              SearchLimitY, Dx0, Dy0, SubPixFlag, overSampleRatio, *rest)
            rec = {
                'seq': n,
                'kind': 'fine' if SubPixFlag else 'coarse',
                'chip_size_x': float(ChipSizeX),
                'chip_size_y': float(ChipSizeY),
                'oversample': float(overSampleRatio),
                'grid_shape': list(np.shape(xGrid)),
                'measured': int(np.count_nonzero(~np.isnan(dx))),
                'arrays': {},
            }
            # The grid too: a level's coarse pass runs on a decimated grid, and the decimation is
            # what places a coarse estimate over the fine points it stands for.
            #
            # **The prior as well, and it is not redundant with `in_Dx0`.** A level is handed `Dx00`, the
            # cell-mean of `Dx0` over `1/Scale` cells resized to the level's lattice
            # (`autoRIFT.py:161-179`) — a different array from the full-resolution `in_Dx0` the outer
            # capture records. The prior decides which window each point searches, so replaying a level's
            # pass without it searches somewhere else: reconstructing it instead of recording it is
            # possible, since `self.Dx0` is never mutated in the loop, but it puts the reference's
            # `colfilt` and two `cv2.resize` calls in the replay's path, where a reconstruction error is
            # indistinguishable from a correlator difference. Recording what the correlator was actually
            # handed removes that whole class of doubt.
            #
            # Note the sign: this is `Dx00`/`Dy00` as passed, so `Dy00` is still cartesian-Y — the
            # correlator's own `Dy0 = -Dy0` happens inside the call, after this wrapper sees it. A replay
            # must negate it exactly as `pointset_from_capture` negates `in_Dy0`.
            for label, arr in (('dx', dx), ('dy', dy), ('xgrid', xGrid), ('ygrid', yGrid),
                               ('searchx', SearchLimitX), ('searchy', SearchLimitY),
                               ('dx0', Dx0), ('dy0', Dy0)):
                a = np.asarray(arr)
                if a.ndim != 2 or a.dtype.str not in xchg.TAGS:
                    continue
                rec['arrays'][label] = _write(f'lvl{n}_{label}', a)
            levels.append(rec)
            print(f'[capture] level call {n}: {rec["kind"]} chip {ChipSizeX} '
                  f'grid {rec["grid_shape"]} measured {rec["measured"]}', flush=True)
            return dx, dy

        setattr(module, name, patched)

    wrap_corr('arImgDisp_u')
    wrap_corr('arImgDisp_s')

    # The rejection mask each pass keeps. Paired with the raw `dx` above, this separates "this level
    # never measured the point" from "this level measured it and threw it away", which the merged
    # `ChipSizeX` cannot distinguish.
    if not hasattr(module, 'DISP_FILT'):
        raise AttributeError(f'{module.__name__} has no DISP_FILT; see `reference_module`.')
    disp_filt = module.DISP_FILT
    original_filt = disp_filt.filtDisp

    def patched_filt(self, Dx, Dy, SearchLimitX, SearchLimitY, M, OverSampleRatio):
        seq['n'] += 1
        n = seq['n']
        kept = original_filt(self, Dx, Dy, SearchLimitX, SearchLimitY, M, OverSampleRatio)
        levels.append({
            'seq': n,
            'kind': 'filtDisp',
            'filt_width': int(self.FiltWidth),
            'frac_valid': float(self.FracValid),
            'iterations': int(self.Iter),
            'oversample': float(OverSampleRatio),
            'grid_shape': list(np.shape(Dx)),
            'in_mask': int(np.count_nonzero(M)),
            'kept': int(np.count_nonzero(kept)),
            'arrays': {'kept': _write(f'lvl{n}_kept', np.asarray(kept).view(np.uint8))},
        })
        print(f'[capture] level call {n}: filtDisp width {self.FiltWidth} '
              f'frac {self.FracValid:.4f} kept {np.count_nonzero(kept)} '
              f'of {np.count_nonzero(M)}', flush=True)
        return kept

    disp_filt.filtDisp = patched_filt


def install():
    """Patch `autoRIFT.runAutorift` to dump its inputs and outputs around the real call."""
    OUT.mkdir(parents=True, exist_ok=True)
    from autoRIFT.autoRIFT import autoRIFT

    original = autoRIFT.runAutorift
    state = {'n': 0}

    def patched(self):
        state['n'] += 1
        call = state['n']
        manifest = {'call': call, 'arrays': {}, 'scalars': {}, 'skipped': {}}

        # Per-level patches go on before the call, since they record what happens inside it. They
        # write into this call's manifest, so a driver that runs the correlator more than once keeps
        # the levels attributed to the right run.
        install_levels(reference_module(), manifest)
        # Off unless asked for: a line-level trace over a 2344x2336 grid run costs real time, and the
        # ordinary capture does not need it.
        if os.environ.get('CAPTURE_STAGES'):
            install_stage_trace(manifest, level=int(os.environ.get('CAPTURE_STAGE_LEVEL', '1')))

        # `CAPTURE_FLOAT32` correlates the float field instead of the 256-level quantization of it.
        #
        # **`DataType` is not the lever, because `uniform_data_type` never runs.** The class has it
        # (`autoRIFT.py:356-404`, with `DataType == 1` keeping a `Float32` field), but the container's
        # vendored driver *inlines* the `DataType == 0` arithmetic instead of calling the method — the same
        # rescale-and-round, written out at `vend/testautoRIFT.py:449-481` — so setting the attribute
        # changes nothing and `grep uniform_data_type` finds only the definition. Setting `DataType` and
        # reporting success is exactly the silent no-op this harness has been bitten by before, so the
        # arrays are widened here instead, where their state can be checked rather than assumed.
        #
        # Widening the *quantized* bytes rather than reaching back for the pre-quantization field: the
        # driver has already overwritten `obj.I1`, so the float field is gone by the time any patch on this
        # method can see it. What that buys is still the measurement wanted — the two runs then differ only
        # in which C++ template correlates the same values, `arImgDisp_u` on bytes against `arImgDisp_s` on
        # floats, which is the reference disagreeing with itself and the floor every other rung is read
        # against. It does *not* measure what the quantization cost, since both runs see quantized values;
        # `tools/ab` stage 1 measures that on a windowed float field.
        if os.environ.get('CAPTURE_FLOAT32'):
            if self.I1.dtype != np.uint8:
                raise RuntimeError(
                    f'expected the driver to have quantized to uint8; found {self.I1.dtype}. '
                    'The float path is selected by widening those bytes, so a different input dtype '
                    'means this is no longer the comparison it claims to be.')
            self.I1 = self.I1.astype(np.float32)
            self.I2 = self.I2.astype(np.float32)
            self.DataType = 1
            print('[capture] widened I1/I2 to Float32: the pyramid now reaches arImgDisp_s '
                  'rather than arImgDisp_u on the same values', flush=True)

        original(self)
        # The patch has to have *fired*, and the element type is what says so. A run that reported the
        # widening and then correlated bytes would look identical in every other respect.
        if os.environ.get('CAPTURE_FLOAT32') and self.I1.dtype != np.float32:
            raise RuntimeError(f'CAPTURE_FLOAT32 was set but I1 is {self.I1.dtype} after the call')
        sys.settrace(None)
        # That the patch was *installed* does not establish that it *fired* — the same distinction the
        # stale-`autoRIFT_intermediate.nc` trap taught, applied to the wrapper rather than to the run.
        # `autorift()` calls the correlator at least twice per resolved level, so an empty list means
        # the wrapper never ran, and the only way that happens after the checks in `install_levels` is
        # if the loop `continue`d out of every level. Either way it is a fault to surface, not a
        # measurement to record: a downstream reader seeing `levels: []` cannot tell the two apart.
        if not manifest['levels']:
            raise RuntimeError(
                'runAutorift returned with no per-level records. The correlator wrapper was '
                'installed but never called, so nothing inside the pyramid was observed.')
        # Inputs are taken *after* the call, not before. `runAutorift` rewrites them as its first
        # action — `self.xGrid = np.round(self.xGrid[0:rlim, 0:clim]) + 0.5` and the same for `yGrid`,
        # then truncates `Dx0`, `Dy0`, `SearchLimit*` and the chip bounds to that same window
        # (`autoRIFT.py:883-905`) — and it is the rewritten arrays the correlator sees. Dumping them
        # first captures an integer grid missing the half pixel, which puts every search centre half a
        # pixel from where the reference put it: a residual that is zero under uniform motion and
        # grows with the velocity gradient, so it hides in the median and shows up only as a
        # gradient-correlated difference map.
        #
        # `self.I1`/`self.I2` are also rewritten, by the uniform-data-type conversion, so taking them
        # after is right for the same reason.
        _dump(self, INPUTS, 'in_', manifest)
        _dump(self, OUTPUTS, 'out_', manifest)

        record_filtered_scenes(manifest)

        for name in SCALARS:
            if hasattr(self, name):
                s = _scalar(getattr(self, name))
                if s is not None:
                    manifest['scalars'][name] = s

        with open(OUT / f'call{call}.json', 'w') as f:
            json.dump(manifest, f, indent=2)
        print(f'[capture] call {call}: {len(manifest["arrays"])} arrays -> {OUT}', flush=True)
        if manifest['skipped']:
            print(f'[capture] skipped (dtype not in xchg): {manifest["skipped"]}', flush=True)

    autoRIFT.runAutorift = patched
    print(f'[capture] runAutorift patched; output -> {OUT}', flush=True)


def main():
    install()
    # The rest of the command line goes to the ordinary entry point, so this runs the real pipeline
    # rather than a reconstruction of it.
    from hyp3_autorift.process import main as process_main

    sys.argv = ['hyp3_autorift'] + sys.argv[1:]
    process_main()


#
# `DxF0`/`DyF0` are the previous level's merged field brought onto this level's grid — `colfilt` over a
# `(Scale+1)²` window, then `INTER_AREA` (`autoRIFT.py:823-834`) — and they are what fills this level's
# holes at `:847`. On a coarse level that is most of the array: the raw fine field is 90.6% NaN at chip
# 384 and 92.6% at chip 768, so the array reaching the `INTER_CUBIC` upsample is mostly `DxF0` rather
# than anything this level measured. Without them the fill path is the one part of the merge that cannot
# be attributed to a step.
STAGE_LOCALS = (
    'xGrid0', 'yGrid0', 'M0', 'SearchLimitX0', 'SearchLimitY0', 'Dx00', 'Dy00',
    'xGrid0C', 'yGrid0C', 'SearchLimitX0C', 'SearchLimitY0C', 'Dx0C', 'Dy0C',
    'DxC', 'DyC', 'M0C', 'MC', 'MC2',
    'DxF', 'DyF', 'DxF0', 'DyF0', 'DxFM', 'DyFM', 'MF', 'MM',
    'Dx', 'Dy', 'ChipSizeX', 'InterpMask',
)


def install_stage_trace(manifest, level=1):
    """Dump each named local of `runAutorift` as this chip-size level produces it.

    `level` is the chip-size loop's own index `i`, so it is **0-based**: `level=0` is the base chip
    size, `level=1` the next one up. One level per run, because a line trace over a whole-scene grid
    costs real time and dumping four multiplies the output for no gain.

    **A name bound in this level's frame is not necessarily this level's value.** `autorift`'s loop
    body rebinds `xGrid0`, `M0`, `DxC` and the rest each iteration, so at the first line of an
    iteration every one of them still holds the *previous* level's array — and dumping the first value
    seen while the level index matches records that leftover instead. On the golden Landsat case that
    put a 2344x2336 base-level `xGrid0` in the files for level 1, where the resized grid is
    1172x1168: the shape is the only thing that gives it away, and a reader comparing against a
    level-1 array of the right shape would find no counterpart at all while one comparing loosely
    would diff two different levels.

    So the change detector is *seeded* at the level boundary with whatever is bound there, rather than
    cleared. A name is then dumped when its bytes first differ from what this level inherited, which is
    the first value this level actually computed.

    Names in `REDUMP` are dumped at every subsequent change too, numbered in order. `SearchLimitX0` is
    built at `:599`, read by the coarse pass at `:629`, then zeroed against the coarse mask at `:724`
    before the fine pass reads it — three distinct values under one name, and dumping only the first
    hides the handoff the trace exists to check. `Dx`, `Dy`, `ChipSizeX` and `InterpMask` are the
    cross-level accumulators, where the same reasoning applies for a different reason: they are
    *supposed* to carry in from the previous level, and seeding is what makes `rev0` this level's
    contribution rather than the state it started from.
    """
    state = {'level': None}
    stages = manifest.setdefault('stages', {})

    # **Every state of every name, not the first.** Almost every local in this loop is rebound at least
    # once, and which state a downstream stage consumes is exactly what a stage comparison has to
    # establish — so dumping one state per name is what makes a comparison silently compare the wrong
    # quantity. Five names cost a rung each before this became the default:
    #
    #   `xGrid0`/`yGrid0`  `cv2.resize` builds them, then the even-chip snap replaces them with
    #                      `round(x + 0.5) - 0.5` (`:509-530`). The correlator reads the snapped grid;
    #                      the resize's values sit at quarter-fractions where the snapped ones are all
    #                      half-integers, so the first state reads as an interpolation mismatch when the
    #                      interpolation agrees to the bit.
    #   `Dx00`/`Dy00`      built by `colfilt`, then resized and **rounded** at `:585-586`. The prior
    #                      displaces the chip and `chip_bounds` floors the result, so an unrounded prior
    #                      moves the chip by a pixel.
    #   `SearchLimitX0`    built at `:599`, read by the coarse pass at `:629`, zeroed against the coarse
    #                      mask at `:724` — and raised to `minSearch` in between. Four states at a
    #                      coarse level, two of them before the rewrite.
    #   `MF`               created all-zero at `:790`; only the fill loop sets anything, so one dump
    #                      says the reference filled nothing.
    #   `DxF`/`DyF`        bound by the correlator, then nulled against the rejection mask at
    #                      `:772-773`. Medianing the un-nulled field sees neighbours the reference has
    #                      already discarded.
    #
    # The cost of dumping all of them is disk, and the cost of dumping one was five wrong measurements.
    # A consumer names the state it wants — by a property that identifies it, not by an index, since the
    # count varies with the level.
    revs = {}

    def _sig(v):
        # Over the raw bytes rather than the values: a displacement array is full of `NaN`, so a
        # numeric sum both throws on conversion and compares unequal to itself. Hashing detects any
        # change, including one a sum would cancel, and does not care what the values mean.
        return hash(v.tobytes())

    def tracer(frame, event, arg):
        if event == 'call':
            return tracer if frame.f_code.co_name == 'autorift' else None
        if event != 'line':
            return tracer
        loc = frame.f_locals
        # `i` is the chip-size loop variable; it advances once per level.
        if 'i' in loc and isinstance(loc['i'], (int, np.integer)):
            new = int(loc['i'])
            if new != state['level']:
                state['level'] = new
                # Seed, not clear: every array bound right now belongs to the level that just ended,
                # so recording its signature is what makes the next dump of that name a value this
                # level computed.
                revs.clear()
                for nm in STAGE_LOCALS:
                    vv = loc.get(nm)
                    if isinstance(vv, np.ndarray) and vv.ndim == 2:
                        revs[nm] = _sig(vv)
        if state['level'] != level:
            return tracer
        for name in STAGE_LOCALS:
            v = loc.get(name)
            if not isinstance(v, np.ndarray) or v.ndim != 2:
                continue
            sig = _sig(v)
            if revs.get(name) == sig:
                continue
            revs[name] = sig
            n = sum(1 for kk in stages
                    if kk.startswith(name + '_rev') and kk.endswith('_L%d' % level))
            key = '%s_rev%d_L%d' % (name, n, level)
            a = v.view(np.uint8) if v.dtype == np.bool_ else v
            if a.dtype.str not in xchg.TAGS:
                continue
            # `.copy()`, not `ascontiguousarray`: that returns the *same* buffer for an array which is
            # already contiguous, so a later in-place write reaches the dumped bytes. `SearchLimitX0`
            # is mutated at `autoRIFT.py:724` — zeroed wherever the coarse mask rejected — long after
            # this dump, and without a copy the file ends up holding the post-mutation state while
            # claiming to be the value the coarse pass consumed. That misattributed 420 coarse radii to
            # a reducer difference when the arrays were simply from different moments.
            stages[key] = _write('stage_' + key, np.ascontiguousarray(a).copy())
            stages[key]['level'] = level
            print('[capture] stage %s %s %s' % (key, a.dtype.str, a.shape), flush=True)
        return tracer

    sys.settrace(tracer)
    return tracer


if __name__ == '__main__':
    main()


# Ordered stage trace: the local arrays `runAutorift` builds between its calls.
#
# The wrappers above record what crosses a function boundary, which is not enough to say *where*
# agreement is lost. Between `arImgDisp_u` and `filtDisp` the reference builds a dozen intermediates —
# the resized grid, the decimated search limits, the coarse mask, its distance-transform dilation, the
# expanded mask that gates the fine search — and a disagreement can enter at any of them. Comparing
# only the endpoints leaves the interior unmeasured, which is how five successive hypotheses about
# this pipeline each turned out to be wrong: each was consistent with the endpoints and none was
# tested against the intermediate that would have refuted it.
#
# `sys.settrace` on the frame is what makes them reachable: they are local to a loop body, so nothing
# short of a line-level trace can see them without rewriting the function.
