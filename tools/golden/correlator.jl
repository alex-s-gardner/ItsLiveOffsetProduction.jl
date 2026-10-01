# Phase 1: AutoRIFT.jl against the reference correlator, on production imagery.
#
#   julia --project=tools/golden -t 8 tools/golden/correlator.jl S2B_MSIL1C_20200612
#
# A diagnostic, not the gate. The gate is the product comparison, and when that disagrees this says
# whether the correlator or the packaging is responsible — which is worth having before the packaging
# exists, because it is the difference between one question and two.
#
# Both sides are handed the *same* arrays: the filtered pair, grid, priors, per-point search limits
# and chip bounds that `capture.py` took at the reference's own `runAutorift` boundary. So a
# preprocessing difference cannot appear here as a correlator difference, and the comparison is of
# the correlator alone. This is the `tools/ab` stage-2 discipline on production imagery instead of a
# hand-cut window.
#
# Four conventions have to be right, each documented and asserted in `tools/ab/README.md` and each
# capable of producing a plausible-looking wrong answer rather than an obvious failure:
#
#   **Index base.** The reference's grid is 0-based pixel indices; `PointSet` is 1-based. So every
#   coordinate gains 1.
#
#   **The half pixel.** `arImgDisp_s` adds `+0.5` to the grid internally, and `runAutorift` snaps an
#   even chip's grid to `round(x + 0.5) - 0.5` before that. AutoRIFT.jl adds the same `0.5` in
#   `_shift_points`, so it must *not* be applied here — doing it twice moves every search centre.
#
#   **The `Dy` sign.** `arImgDisp_s` converts its answer from matrix-Y to cartesian-Y before
#   returning — "Y from down being positive to up being positive"
#   (`autoRIFT.py:1142-1143`, and again at `:1308-1309` for the unsigned entry point). So the
#   reference's `Dy` is up-positive while AutoRIFT.jl's `dy` is row-positive, and comparing them
#   needs one negation. This is the correlator's own convention and applies whatever `optical_flag`
#   says; the separate `if optical_flag == 0: Dy = -Dy` in `testautoRIFT.py` is a *second* flip that
#   the driver applies to radar afterwards, and it acts on the already-cartesian value.
#
#   The sign is *measured* rather than asserted, the way `tools/ab/bench_figures.jl` does it: all four
#   combinations are scored and the comparison reports which won. A hardcoded flip is only correct
#   while the writer's convention holds and fails silently when it changes, whereas a near-perfect
#   *anti*-correlation is a flipped sign and nothing else.
#
#   **Array layout.** Handled by `xchg`, which puts the element type and both dimensions in the file.
#   These grids are square, where a wrong-convention read is a silent transpose.
#
# Both are compared where both measured. The reference reports a search-window corner for a
# degenerate chip where AutoRIFT.jl reports nothing (`REFERENCE.md`), so a coverage difference is
# expected and is counted rather than hidden — it is the mechanism that can later shift
# `stable_shift`.

include("manifest.jl")
include("reference.jl")
include("intermediate.jl")

using AutoRIFT
using AutoRIFT: PointSet, params, rebuild
using Printf, Statistics

"""
    pointset_from_capture(k::Capture) -> PointSet{2}

The reference's search grid as a `PointSet`.

Coordinates become 1-based, and a point the reference skipped — search limit zero — keeps a zero
radius, which is how `PointSet` marks a point to skip.

**The grid needs `+1`, not `+0.5`, and this was measured rather than reasoned.** `runAutorift` sets
`xGrid = round(xGrid) + 0.5` before correlating (`autoRIFT.py:890-891`) and `capture.py` dumps the
arrays after that, so the values arriving here are half-integers in 0-based pixel coordinates. The
tempting move is to subtract that half pixel while adding the index base, on the grounds that
`_shift_points` adds AutoRIFT.jl's own `0.5` for the even-chip centroid — but scanning the offset says
otherwise: at `+1.0` exact agreement is **85.0%** with a median residual of 0, and at `+0.5` it is
49.2% with a median of 1/32. The two half pixels do not cancel; they are the same convention counted
once on each side.

Getting this wrong is invisible in a median over the whole scene, because a half-pixel offset produces
no residual under uniform motion and one proportional to the local velocity gradient. Before the fix,
exact agreement fell from 66.5% in the flattest gradient decile to 5.6% in the steepest, and the
difference map showed structure only along the fast-flow margins. Scan the offset before believing any
argument about which convention applies — including this one.
"""
function pointset_from_capture(k::Capture)
    xg = k.arrays["in_xGrid"]
    yg = k.arrays["in_yGrid"]

    # The grid must be the one the correlator saw, which is `round(xGrid) + 0.5` — half-integer and
    # `Float32`. An integer grid means the capture was taken *before* `runAutorift` rewrote it, and
    # comparing against it puts every search centre half a pixel from where the reference put it.
    #
    # That failure is worth an error rather than a warning because it is invisible in the result: the
    # residual is zero under uniform motion and grows with the velocity gradient, so the median stays
    # small and only a heatmap shows the red/blue dipoles along fast flow. One such capture scored
    # 22.7% exact where a correct one on the same sensor and glacier scores above 80%.
    let nz = filter(!iszero, vec(Float64.(xg)))
        isempty(nz) && error("captured grid is entirely zero for this run")
        all(≈(0.5), nz .- floor.(nz)) || error(
            "captured grid is not on the half-integer convention the correlator uses: " *
            "fractional parts $(unique(nz .- floor.(nz))), eltype $(eltype(xg)). " *
            "This capture predates the fix that dumps inputs after `runAutorift` rewrites them " *
            "(`capture.py`); redo it with `intermediate.jl <case> --force`.")
    end
    srx = k.arrays["in_SearchLimitX"]
    sry = k.arrays["in_SearchLimitY"]
    csmin = k.arrays["in_ChipSizeMinX"]
    csmax = k.arrays["in_ChipSizeMaxX"]
    dx0 = k.arrays["in_Dx0"]
    dy0 = k.arrays["in_Dy0"]

    chip0 = Int(k.scalars["ChipSize0X"])
    scale_y = Float64(k.scalars["ScaleChipSizeY"])

    rx, ry = _level_search_limits(srx, sry, k)

    # **`Dy0` is negated, because the captured value is not the one the correlator uses.**
    # `arImgDisp_u`/`arImgDisp_s` do `Dy0 = -Dy0` as their first act (`autoRIFT.py:1058`,
    # and `:1231` for the signed entry point), converting the prior from cartesian-Y to
    # matrix-Y before any chip is cut; the capture records `self.Dy0`, which is pre-flip. So
    # the reference offsets its chip by `+Dy0` where AutoRIFT.jl's `chip_bounds` offsets by
    # `-dy_prior`, and handing the captured value through unchanged places the chip `2 * Dy0`
    # rows from where the reference put it.
    #
    # This is the same convention as the `Dy` on the way out, counted once on each side: the
    # correlator flips the prior going in and flips the answer coming out, so a comparison
    # must undo both. Undoing only the output — which is what a reader checking `Dy` against
    # `dy` naturally does — leaves the input flip in place.
    #
    # Measured on the worst 64x64 block of the golden Landsat case, where `Dy0` reaches 14:
    # `dx` exact goes 6.10% -> **99.61%** and `dy` 8.15% -> **99.55%**, with the mean residual
    # falling from +0.90 px to -0.0004. The failure is invisible where the prior is small and
    # grows with it, so it concentrates on the fast-flow tongues and reads as a velocity-
    # dependent bias rather than as a sign error. It is also invisible in `dy` alone: the
    # misplaced chip biases `dx` while `dy` stays near zero, because the wrong rows still
    # correlate best at a similar vertical offset.
    #
    # `positioned` marks the points geogrid placed outside the image. The capture holds the reference's
    # own post-rewrite grid, where such a point is at `0.5`: `runAutorift` overwrites the `-32767`
    # sentinel with zero (`testautoRIFT.py:390-391`) and then stores `round(xGrid) + 0.5`. A point with a
    # real position is at `1.5` or beyond, so the threshold sits below every one of them rather than being
    # tuned. Without this a blocked run has no way to tell a placeholder coordinate from a position and
    # loses 5% of a Sentinel-1 granule — see `dev/plan-16gib.md` and `block_gate.jl`.
    return PointSet(
        Float64.(xg) .+ 1, Float64.(yg) .+ 1,
        rx, ry,
        Float64.(dx0), .-Float64.(dy0),
        fill(chip0, size(xg)), fill(round(Int, chip0 * scale_y), size(xg)),
        Int.(csmin), Int.(csmax),
        [Float64(v) > 1 for v in xg],
    )
end

"""
    _thin(grid::PointSet{2}, stride::Integer; block = 128) -> PointSet{2}

`grid` with all but a scattered `1/stride^2` of its points marked skipped — `block`-sized square tiles
on a `stride * block` lattice — for a comparison that costs a fraction of the whole scene.

**Tiles, not a point lattice, because the correlator is not pointwise.** `filtDisp` and the
level merge consult each point's neighbors, so a point whose neighbors are all skipped loses its
base-level measurement to an interpolated one. Thinning to every 16th *point* on the L1 RSLC grid
drives `exact` to **0.00%** against 73.90% on the whole grid and moves `dx` `bias_core` from +0.051 to
+0.109. A tile keeps each kept point's neighborhood.

**Tiles rather than one window** because a NISAR grid is a rotated radar footprint on a map grid: the
chip sizes, search radii and levels a point resolves at vary across it, and one window carries
whichever mix its corner holds. Scattered tiles keep the mix roughly in proportion.

**A thinned run's answers are its own, whatever the tiling.** A level's coarse grid is this grid
decimated by 1, 2, 4, 8, and a thinned one can fall below the width its filter needs, at which point
the level silently produces nothing — so the two sides resolve different level sets and the effect is
not confined to a tile's border. On L1 at `stride = 4`, `exact` runs 16.81% at a 128-px tile, 9.41% at
256 and 0.00% at 512, non-monotonically, on unchanged code. Compare a thinned run only against another
at the same `stride` *and* `block`; `GATES.md` records the tiling every threshold was drawn at.

The shape is preserved and points are dropped by zeroing their radius — the marker
`pointset_from_capture` already uses for a point the reference skipped. Cropping the arrays instead
would renumber the grid, and every coordinate here is one convention away from a plausible wrong
answer.

**Coverage is not comparable across a thinned run.** A zeroed point produces `NaN`, which
`_axis_stats` cannot distinguish from a genuine no-measurement, so `only_reference` counts the
thinned-out points too and rises by construction. Read `exact`, `bias_core` and `correlation` from a
thinned run; read coverage from the whole grid only.
"""
function _thin(grid::PointSet{2}, stride::Integer; block::Integer = 128)
    stride >= 1 || throw(ArgumentError("stride must be >= 1, got $stride"))
    block >= 1 || throw(ArgumentError("block must be >= 1, got $block"))
    # One `block`-wide tile per `stride * block` period along each axis. Tiles start from the array's
    # own first index, so an offset axis tiles the way a 1-based one does.
    tiles(d) = [s:min(s + block - 1, lastindex(grid.radius_x, d))
                for s in firstindex(grid.radius_x, d):(stride * block):lastindex(grid.radius_x, d)]
    keep = falses(axes(grid.radius_x))
    for rows in tiles(1), cols in tiles(2)
        keep[rows, cols] .= true
    end
    z = zero(eltype(grid.radius_x))
    return rebuild(grid; radius_x = ifelse.(keep, grid.radius_x, z),
                   radius_y = ifelse.(keep, grid.radius_y, z))
end

"""
    _level_search_limits(srx, sry, k) -> (rx, ry)

The search radii the correlator is handed, which are **not** the ones the capture records.

`runAutorift` rewrites them at the top of every level (`autoRIFT.py:598-602`) and it is the rewritten
array the correlator sees:

```python
idxZero = (SearchLimitX0 <= 0) | (SearchLimitY0 <= 0)
SearchLimitX0[idxZero] = 0
SearchLimitY0[idxZero] = 0
SearchLimitX0[~idxZero & (SearchLimitX0 < minSearch)] = minSearch
```

Two rules, and each changes which window is searched:

  * **Either axis zero zeroes both.** A point wanting 5 across and 0 down is skipped entirely, not
    searched as a horizontal line.
  * **A nonzero radius below `minSearch` is raised to it.** A point asking for 1 searches at 6.

Skipping this is not a small effect: 944,036 points on the golden Landsat case carry a radius that
differs from the captured one, by up to 23. And it biases a comparison in a way that looks like the
opposite of its cause — the rewritten points are the *small*-radius ones, so agreement appears best
where the radius is smallest, which reads as the large radii being at fault.

`minSearch` comes from the capture when present. A capture taken before it was recorded falls back to
the reference's own default with a warning, since silently using a wrong floor reproduces exactly the
bug this function exists to fix.
"""
function _level_search_limits(srx, sry, k::Capture)
    minsearch = if haskey(k.scalars, "minSearch")
        Int(k.scalars["minSearch"])
    else
        @warn "capture predates `minSearch` being recorded; using the reference's default of 6 " *
              "(`autoRIFT.py:946`). Redo the capture with `--force` to take it from the run."
        6
    end

    rx = Int.(srx)
    ry = Int.(sry)
    out_x = similar(rx)
    out_y = similar(ry)
    @inbounds for i in eachindex(rx, ry)
        x, y = rx[i], ry[i]
        if x <= 0 || y <= 0
            out_x[i] = 0
            out_y[i] = 0
        else
            out_x[i] = max(x, minsearch)
            out_y[i] = max(y, minsearch)
        end
    end
    return out_x, out_y
end

"""
    subpixel_from_capture(k::Capture) -> Tuple{Vararg{PyramidRefine}}

The subpixel methods the reference used, one per chip-size level, finest first.

`OverSampleRatio` is a per-chip-size dictionary the driver assembles at run time
(`testautoRIFT.py:488-510`) rather than a scalar, and `autoRIFT.py:652` looks it up per level. So the
quantization step is a property of the level: 1/16 px at the base chip size and 1/32 or 1/64 above
it for optical input. A scalar is also handled, because the driver sets one when the chip bounds are
absent.

Entries are ordered by chip size and truncated to the levels this run will use, since a tuple longer
than the level list is a configuration error rather than something to ignore.
"""
function subpixel_from_capture(k::Capture)
    osr = k.scalars["OverSampleRatio"]
    chip0 = Int(k.scalars["ChipSize0X"])
    maxchip = Int(maximum(k.arrays["in_ChipSizeMaxX"]))
    # The levels this run will actually correlate: chip0 * 2^j up to the largest bound present.
    chips = [chip0 << j for j in 0:floor(Int, log2(maxchip / chip0))]

    osr isa Number && return (PyramidRefine(Int(osr)),)
    # JSON object keys arrive as symbols, and the driver keys them by chip size.
    return Tuple(PyramidRefine(Int(osr[Symbol(c)])) for c in chips)
end

"""
    kwargs_from_capture(k::Capture) -> NamedTuple

Correlator settings taken from what the reference used, rather than from what the driver is believed
to set.

`autorift(a, b, ::PointSet)` takes keywords rather than a `Params`, so these are passed through as
such. The per-point fields — coordinates, priors, search radii, chip bounds — travel in the
`PointSet` instead; only the scene-wide settings are here.

`preprocess = :none` because the captured pair is already filtered. Filtering again would compare two
different images.
"""
function kwargs_from_capture(k::Capture)
    chip0 = Int(k.scalars["ChipSize0X"])
    scale_y = Float64(k.scalars["ScaleChipSizeY"])
    spacing = Int(k.scalars["GridSpacingX"])
    maxchip = Int(maximum(k.arrays["in_ChipSizeMaxX"]))

    # **Both chip bounds carry `ScaleChipSizeY`, not just the minimum.** `in_ChipSizeMaxX` is the X
    # axis alone — the array's name says so — and the pyramid doubles both axes together, so Y's
    # maximum is Y's minimum times the same number of doublings. Scaling only `chip_size` leaves the
    # two axes reaching their maxima after different doublings, which `_check_levels` rejects: on a
    # Sentinel-1 pair at `ScaleChipSizeY = 0.25` that is `8` in X against `32` in Y.
    #
    # Invisible wherever the pixel is square, which is where the scale is 1.0 and the two forms
    # coincide. That is a property of the pixel, not of the platform: `ScaleChipSizeY` is the ratio of
    # the ITS_LIVE parameter chip sizes *in metres* (`vend/testautoRIFT.py:377-378`), so it departs
    # from 1.0 exactly when the y and x pixel sizes differ, and it varies per acquisition — 0.2353,
    # 0.25 and 0.2857 across the eight golden radar cases. Testing only on cases with a square pixel
    # exercises neither form.
    # `threaded` follows the thread count the process was started with, since `Params` cannot know it:
    # `threaded = true` on one thread pays the parallel machinery for no workers, and `false` on eight
    # leaves seven idle. It parallelises over grid points within a pass (`src/params.jl:392`), so it
    # changes how long a pass takes and not what it answers — measured bit-identical on a 201x201
    # window of the L1 RSLC grid, 38,983 points either way, 148.9 s against 22.3 s.
    #
    # A golden run is one pair at a time, which is the case this is for. The documentation's advice to
    # leave it `false` and run one pair per worker is about batch throughput across many pairs, where
    # the outer parallelism is the better axis; a single whole-scene comparison has no outer axis.
    return (; chip_size = (X = chip0, Y = round(Int, chip0 * scale_y)),
            chip_size_max = (X = maxchip, Y = round(Int, maxchip * scale_y)),
            grid_spacing = (X = spacing, Y = spacing),
            subpixel = subpixel_from_capture(k),
            threaded = Threads.nthreads() > 1,
            preprocess = :none)
end

"""
    compare_correlator(c::GoldenCase; n = 100) -> NamedTuple

Run AutoRIFT.jl on the reference's own captured inputs and diff `dx`/`dy` against its `Dx`/`Dy`.

Returns the per-axis statistics plus the coverage split, since a point one side answered and the
other did not is a different finding from a point they both answered differently.
"""
function compare_correlator(c::GoldenCase; n::Integer = 100, stride::Integer = 1,
                            block::Integer = 128)
    k = read_capture(c; n)

    grid = pointset_from_capture(k)
    stride > 1 && (grid = _thin(grid, stride; block))
    kw = kwargs_from_capture(k)

    a = k.arrays["in_I1"]
    b = k.arrays["in_I2"]

    # Argument order. `arImgDisp_s(a, b)` cuts its chip from `b` and its window from `a`; the
    # reference calls it as `arImgDisp_s(self.I2, self.I1)`, so `I1` supplies the chip. AutoRIFT.jl's
    # `autorift(reference, secondary)` cuts its chip from `secondary`, so `I1` binds to `secondary`.
    @info "correlating" scene=size(a) npoints=length(grid.x) subpixel=kw.subpixel
    t = @elapsed out = autorift(b, a, grid; kw...)

    rdx = k.arrays["out_Dx"]
    rdy = k.arrays["out_Dy"]

    # The reference's output grid is one point smaller in each direction than its input grid, so the
    # comparison is over the overlap, stated rather than assumed.
    ny = min(size(out.dx, 1), size(rdx, 1))
    nx = min(size(out.dx, 2), size(rdx, 2))
    @info "grid overlap" julia=size(out.dx) reference=size(rdx) compared=(ny, nx)

    stats = map((:dx, :dy)) do axis
        jul = getproperty(out, axis)[1:ny, 1:nx]
        ref = (axis === :dx ? rdx : rdy)[1:ny, 1:nx]
        # Score both signs and keep the better. `dy` needs the flip and `dx` does not, but measuring
        # says so rather than a comment asserting it — and a flip that stops being needed shows up
        # here as a changed `sign` instead of as a silent bias.
        pos = _axis_stats(String(axis), jul, ref, +1)
        neg = _axis_stats(String(axis), jul, ref, -1)
        # **Chosen on correlation, not on `exact`.** `exact` is zero for *both* signs whenever no point
        # is quantized, which is every pair whose base chip size was skipped — on the golden L7 pair the
        # reference resolves only chips 32 and 64, where both implementations replace the measurement
        # with a bicubic resize and 0.03% of values land on any 1/N grid. The comparison then ties at
        # zero, `>=` keeps `+1`, and the report shows `dy` at sign `+` with a −0.85 correlation and a
        # 0.76 px bias: a flipped axis presented as a measured choice.
        #
        # Correlation is the right discriminant because it is what a sign error actually destroys — a
        # near-perfect *anti*-correlation is a flipped sign and nothing else, which is why `_axis_stats`
        # computes it. It is also defined wherever two points vary, so it does not collapse on an
        # unquantized level the way `exact` does.
        pick(a, b) = (isnan(a.correlation) ? -Inf : a.correlation) >=
                     (isnan(b.correlation) ? -Inf : b.correlation) ? a : b
        pick(pos, neg)
    end

    return (; time = t, overlap = (ny, nx), dx = stats[1], dy = stats[2], stride,
            julia_size = size(out.dx), reference_size = size(rdx), result = out, capture = k)
end

"""
    _axis_stats(name, jul, ref, sign) -> NamedTuple

Compare `jul` against `sign * ref` on one axis.

`NaN` is the correlator's no-measurement marker on both sides, so a point one side answered alone is
counted as coverage rather than folded into the value statistics: the reference reports a
search-window corner for a degenerate chip where AutoRIFT.jl reports nothing (`REFERENCE.md`), and
that difference has a different cause from a disagreement about a measured value.

`correlation` is over both-measured points and is the diagnostic for a sign error specifically — a
near-perfect negative correlation is a flipped sign and nothing else.
"""
function _axis_stats(name, jul, ref, sign::Int)
    both = 0; only_j = 0; only_r = 0; exact = 0
    d = Float64[]
    va = Float64[]
    vb = Float64[]
    for i in eachindex(jul, ref)
        mj = isnan(jul[i]); mr = isnan(ref[i])
        if mj && mr
            continue
        elseif mr
            only_j += 1
        elseif mj
            only_r += 1
        else
            both += 1
            a = Float64(jul[i]); b = sign * Float64(ref[i])
            push!(va, a); push!(vb, b)
            a == b ? (exact += 1) : push!(d, a - b)
        end
    end
    ad = abs.(d)
    # **Two biases, because one mean cannot answer the question.** A mean over every both-measured
    # point is not robust to a heavy tail, and these comparisons have one: on the Sentinel-1 SLC pair
    # 187 of 60,611 points disagree by up to 45.9 px and drag `bias` to +0.027 px while the 78% that
    # agree within a pixel sit at −0.0008. Those 187 are two-sided (99 positive, 88 negative), so they
    # are the peak picked on the other side of a nearly flat surface — this pair correlates at a median
    # of 0.148, SAR speckle over a 24-day repeat — and not a systematic offset that a mean should be
    # reporting as one.
    #
    # `bias` stays the mean over everything, because a tail that grew would otherwise go unnoticed.
    # `bias_core` is the same mean over points agreeing within one pixel: that is the population a
    # *systematic* error lives in, and it is the number to read when asking whether the correlator is
    # offset. Reporting only one of the two invites the wrong conclusion in whichever direction that
    # comparison happens to be tailed.
    core = d[ad .<= 1.0]
    return (; name, sign, both, only_julia = only_j, only_reference = only_r, exact,
            exact_fraction = both == 0 ? 1.0 : exact / both,
            max_abs = isempty(ad) ? 0.0 : maximum(ad),
            p99 = isempty(ad) ? 0.0 : quantile(ad, 0.99),
            median = isempty(ad) ? 0.0 : median(ad),
            bias = isempty(d) ? 0.0 : mean(d),
            bias_core = isempty(core) ? 0.0 : mean(core),
            n_core = length(core), n_tail = count(>(10.0), ad),
            correlation = length(va) < 2 ? NaN : cor(va, vb))
end

function report(r)
    @printf("\ncorrelated in %.1f s; grid julia %s, reference %s, compared %s\n",
            r.time, r.julia_size, r.reference_size, r.overlap)
    # A thinned run's coverage columns are not the scene's, and they read as a coverage collapse to
    # anyone who does not know the stride. Stated here rather than left to the reader: `_thin` marks the
    # points it drops the same way the reference marks one it skipped, so they land in `only ref`.
    if r.stride > 1
        @printf("STRIDE %d: 1 point in %d searched. `only ref` counts the thinned-out points, so read\n\
                 exact, bias core and corr here; read coverage from a whole-grid run.\n",
                r.stride, r.stride^2)
    end
    # The exact **count** beside the fraction, because the fraction is not comparable between runs whose
    # coverage differs — and coverage is one of the things being fixed. Measured on the S2A case: two
    # changes moved `exact` from 92.66% to 96.74% while the count moved by 37 points of 543,071, because
    # the 24,680 points that left `both` were disproportionately ones AutoRIFT.jl had got wrong. Dropping
    # them from the denominator raises the fraction and improves nothing.
    #
    # The two exclusive sets are printed for the same reason and read together: a change that shrinks
    # both is unambiguously better, and one that grows either is trading coverage for a percentage.
    @printf("\n%-5s %5s %9s %8s %8s %8s %10s %9s %9s %9s %8s\n",
            "axis", "sign", "both", "only jl", "only ref", "exact", "exact n", "median", "p99", "max", "corr")
    for s in (r.dx, r.dy)
        @printf("%-5s %5s %9d %8d %8d %7.2f%% %10d %9.4g %9.4g %9.4g %+8.5f\n",
                s.name, s.sign > 0 ? "+" : "-", s.both, s.only_julia, s.only_reference,
                100 * s.exact_fraction, s.exact, s.median, s.p99, s.max_abs, s.correlation)
    end
    @printf("\nmeasured by julia %d, by reference %d\n",
            r.dx.both + r.dx.only_julia, r.dx.both + r.dx.only_reference)
    # Both biases, and the tail count that explains any gap between them. See `_axis_stats`: the mean
    # over everything is not robust to a heavy tail, so a reader given only that number reads a
    # cancelling tail as a systematic offset.
    @printf("\nbias:      dx %+.6g, dy %+.6g   (mean over all both-measured points)\n",
            r.dx.bias, r.dy.bias)
    @printf("bias core: dx %+.6g, dy %+.6g   (within 1 px: %d and %d points)\n",
            r.dx.bias_core, r.dy.bias_core, r.dx.n_core, r.dy.n_core)
    @printf("tail >10px: dx %d, dy %d of %d\n", r.dx.n_tail, r.dy.n_tail, r.dx.both)
    return nothing
end

function main(args)
    isempty(args) &&
        error("usage: correlator.jl <product-name-fragment> [--run N] [--stride S] [--block B]")
    c = only(cases(args[1]))
    n = 100
    i = findfirst(==("--run"), args); i === nothing || (n = parse(Int, args[i + 1]))
    stride = 1
    j = findfirst(==("--stride"), args); j === nothing || (stride = parse(Int, args[j + 1]))
    block = 128
    b = findfirst(==("--block"), args); b === nothing || (block = parse(Int, args[b + 1]))
    report(compare_correlator(c; n, stride, block))
    return nothing
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS)
