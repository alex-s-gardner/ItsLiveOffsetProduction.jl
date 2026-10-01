# The stage ladder: AutoRIFT.jl against the reference, one intermediate at a time.
#
#   CAPTURE_STAGES=1 CAPTURE_STAGE_LEVEL=0 julia --project=tools/golden \
#       tools/golden/intermediate.jl LC08_L1TP_009011 --force
#   julia --project=tools/golden -t 8 tools/golden/stages.jl LC08_L1TP_009011
#
# **Why one stage at a time.** `autorift()` builds about two dozen intermediates inside its
# chip-size loop (`autoRIFT.py:407-874`) and every one of them is a local, so the function boundaries
# are the only place the two implementations have ever been compared. A comparison there is a
# comparison of the whole composition: when it disagrees, any of the twenty-four steps could be
# responsible, and `tools/golden/README.md` records five successive hypotheses that each survived an
# endpoint check and then turned out to be wrong. What refutes a hypothesis about step 12 is the
# array step 12 produced.
#
# **Every stage is fed the reference's own input.** A stage reads the reference's dumped input, runs
# AutoRIFT.jl's counterpart on exactly those bytes, and is diffed against the reference's dumped
# output. Julia output is never chained into the next Julia stage — that would rebuild the composed
# comparison this exists to take apart, and a difference entering at step 3 would then show up at
# every step after it.
#
# **The gate depends on what the stage produces.** A mask or an integer grid is a *decision*: it is
# exact or it is a semantic difference, and there is no rounding to hide behind. A displacement at the
# base chip size is quantized to the level's upsampling step, so it is also exact. A displacement
# above the base level is not comparable that way at all, because both sides replace the measurement
# with a bicubic resize (`autoRIFT.py:856-866`) and neither field is quantized — there, bias and the
# within-one-step fraction are the statistics and `exact` is meaningless. `tools/golden/README.md`
# records the measurement that closed that question.
#
# **Every rung that runs a correlator runs it on both element types, and the reference against itself is
# the floor.** The reference has two correlators — `arImgDisp_u` on bytes, `arImgDisp_s` on floats,
# separate C++ templates — and production reaches the byte one, because `uniform_data_type` quantizes to
# 256 levels before `runAutorift` is called (`autoRIFT.py:359-384`). Running the pair separates a
# difference in what the correlator *computes* from a tie the quantization created in the surface it
# computes *on*: collapsing a filtered float field onto 256 levels makes plateaus the float field does
# not have, and a plateau broken differently puts the peak whole pixels away rather than one step away.
#
# Handing the two templates the same information — the captured bytes, and those bytes widened to
# `Float32` — gives a disagreement that belongs to neither implementation and bounds what any rung can
# be asked to achieve. Measured at chip 96 on the S2A case: **98.54% exact, maximum 7 px**, the
# reference against itself. A rung reporting that figure has found nothing; one reporting worse has.
#
# Stages run in order and the ladder stops at the first red one, because a stage fed a correct input
# tells you about itself, while the stages after a red one are being asked a question whose premise
# has already failed.

include("manifest.jl")
include("reference.jl")
include("intermediate.jl")
include("correlator.jl")

using AutoRIFT
using AutoRIFT: PointSet, windowmax, windowmean, windowrange, windowmedian, sanitize!,
                rebuild, params, extent, rescale, relax, window, DisplacementField, Params,
                resample, Nearest, Area, measure_at, subpixel_at
using Printf, Statistics

# ---------------------------------------------------------------------------
# Comparing one stage
# ---------------------------------------------------------------------------

"""
    StageResult

One rung of the ladder: what was compared, how it was judged, and whether it passed.

`detail` carries the numbers behind the verdict so a red rung says *how* it differs rather than only
that it does — a stage off by one quantization step everywhere and a stage transposed both report
"not equal".
"""
struct StageResult
    name::String
    reference::String
    gate::String
    passed::Bool
    n::Int
    detail::String
end

"""
    exact_stage(name, ref_name, jl, ref; nanok = true) -> StageResult

Compare a decision — a mask, a grid, a search radius, a chip size — for exact equality.

`NaN` counts as equal to `NaN`, since it is the no-value marker on both sides rather than a value.
Nothing else is tolerated: these arrays are selected or copied from existing values, so a difference
of any size is a difference in which value was selected.

The failure detail names the *first* disagreeing position and the two values there. A count alone
cannot distinguish a one-column shift from a scattered difference, and the position is what a heatmap
would show.
"""
function exact_stage(name, ref_name, jl, ref)
    if size(jl) != size(ref)
        return StageResult(name, ref_name, "exact", false, 0,
                           "shape $(size(jl)) against reference $(size(ref))")
    end
    bad = 0
    first_bad = nothing
    for i in eachindex(IndexCartesian(), jl)
        a, b = jl[i], ref[i]
        same = (a isa AbstractFloat && isnan(a) && b isa AbstractFloat && isnan(b)) || a == b
        same && continue
        bad += 1
        first_bad === nothing && (first_bad = (Tuple(i), a, b))
    end
    n = length(ref)
    detail = if bad == 0
        "all $n equal"
    else
        pos, a, b = first_bad
        @sprintf("%d of %d differ (%.4f%%), first at %s: julia %s, reference %s",
                 bad, n, 100 * bad / n, pos, a, b)
    end
    return StageResult(name, ref_name, "exact", bad == 0, n, detail)
end

"""
    quantized_stage(name, ref_name, jl, ref, step) -> StageResult

Compare two displacement fields that are both quantized to `step`.

This is the base chip size, where a residual below one upsampling step is not a disagreement about
position but about which of two adjacent representable values a peak rounded to. `exact` is the
headline and coverage is reported beside it, because a point one side measured and the other did not
is a different finding from a point they both measured differently.

Passing is judged against the pre-existing L8/L9 benchmark rather than against 100%: below a real
correlation peak both implementations are choosing from noise, and no two implementations agree about
that. `tools/ab/README.md` records the benchmark at 77.4% exact and 97.3% within one step.
"""
function quantized_stage(name, ref_name, jl, ref, step; min_exact = 0.774, min_within = 0.973)
    if size(jl) != size(ref)
        return StageResult(name, ref_name, "quantized", false, 0,
                           "shape $(size(jl)) against reference $(size(ref))")
    end
    both = 0; only_j = 0; only_r = 0; exact = 0; within = 0
    d = Float64[]
    for i in eachindex(jl, ref)
        mj = isnan(jl[i]); mr = isnan(ref[i])
        mj && mr && continue
        if mr
            only_j += 1
        elseif mj
            only_r += 1
        else
            both += 1
            δ = Float64(jl[i]) - Float64(ref[i])
            push!(d, δ)
            δ == 0 && (exact += 1)
            abs(δ) <= step && (within += 1)
        end
    end
    fe = both == 0 ? 0.0 : exact / both
    fw = both == 0 ? 0.0 : within / both
    ad = abs.(d)
    detail = @sprintf("both %d, only jl %d, only ref %d; exact %.2f%%, within one step %.2f%%, median %.4g, p99 %.4g, bias %+.4g",
                      both, only_j, only_r, 100fe, 100fw,
                      isempty(ad) ? 0.0 : median(ad), isempty(ad) ? 0.0 : quantile(ad, 0.99),
                      isempty(d) ? 0.0 : mean(d))
    gate = string("exact>=", round(100min_exact; digits = 1), "% within>=",
                  round(100min_within; digits = 1), "%")
    return StageResult(name, ref_name, gate, fe >= min_exact && fw >= min_within, both, detail)
end

"""
    unquantized_stage(name, ref_name, jl, ref) -> StageResult

Compare two displacement fields that are **not** quantized, which is every level above the base.

Both sides overwrite their own measurements with a bicubic resize of a decimated field
(`autoRIFT.py:811,856-866`; `_undecimate_level`), so neither field's values sit on any quantization
grid and exact agreement is unreachable by construction rather than by defect. Measured: multiples of
1/16, 1/32, 1/64 and 1/128 account for 0.01–0.03% of *either* side's chip-48 values, against 99.4% at
the base level.

So the gate is bias and dispersion. A bias means the two disagree about position systematically, which
is a bug; a symmetric spread is the accumulated Float32 difference between two implementations of a
chain that decimates, median-filters, area-resizes, hole-fills and bicubic-resizes.
"""
function unquantized_stage(name, ref_name, jl, ref; max_bias = 0.01, min_within_tenth = 0.45)
    if size(jl) != size(ref)
        return StageResult(name, ref_name, "unquantized", false, 0,
                           "shape $(size(jl)) against reference $(size(ref))")
    end
    d = Float64[]
    for i in eachindex(jl, ref)
        (isnan(jl[i]) || isnan(ref[i])) && continue
        push!(d, Float64(jl[i]) - Float64(ref[i]))
    end
    isempty(d) && return StageResult(name, ref_name, "unquantized", false, 0,
                                    "no point measured on both sides")
    ad = abs.(d)
    bias = mean(d)
    tenth = count(<=(0.1), ad) / length(ad)
    # The core reported beside the mean, and the *mean* gated. `endpoint_stage` gates the core because a
    # speckle pair's two-sided tail drags the mean, and the same could be true here — so it is reported,
    # to say whether a bias is a tail's or the field's. On the one case where the two are far enough apart
    # to matter they are not: −0.0394 core against −0.0446 mean, so the tail is not what fails it and
    # gating the core would move the threshold without answering anything.
    core = [x for x in d if abs(x) <= 1.0]
    detail = @sprintf("both %d; bias %+.5f (core %+.5f of %d), median %.4g, p95 %.4g, \
                       within 0.1 px %.1f%%",
                      length(d), bias, isempty(core) ? bias : mean(core), length(core),
                      median(ad), quantile(ad, 0.95), 100tenth)
    gate = string("|bias|<=", max_bias, " within0.1>=", round(Int, 100min_within_tenth), "%")
    return StageResult(name, ref_name, gate,
                       abs(bias) <= max_bias && tenth >= min_within_tenth, length(d), detail)
end

# ---------------------------------------------------------------------------
# The ladder
# ---------------------------------------------------------------------------

# The stage trace records one chip-size level per run, and which one it was is in the keys.
function traced_level(k::Capture)
    levels = Set{Int}()
    for name in keys(k.stages)
        m = match(r"_L(\d+)$", name)
        m === nothing || push!(levels, parse(Int, m.captures[1]))
    end
    isempty(levels) && error(
        "this capture has no stage trace. Re-take it with " *
        "`CAPTURE_STAGES=1 CAPTURE_STAGE_LEVEL=<i> intermediate.jl <case> --force`, " *
        "where `<i>` is the chip-size loop's own 0-based index.")
    length(levels) == 1 || error("stage trace spans levels $(sort!(collect(levels))); one per run")
    return only(levels)
end

"""
    stage(k, name, level) -> Matrix

One traced array, by its bare name. Errors rather than returning `nothing` when it is absent, because
a stage silently skipped is a rung the ladder reports as green.
"""
function stage(k::Capture, name::AbstractString, level::Integer)
    revs = _revisions(k, name, level)
    isempty(revs) && error("no stage `$(name)_L$(level)` in this capture; traced: " *
                           join(sort!(collect(keys(k.stages))), ", "))
    # One state, so there is nothing to choose. More than one and the choice *is* the question — which
    # state a downstream stage consumes — so this refuses rather than guessing, and names them.
    length(revs) == 1 && return k.stages[only(revs)]
    error("`$name` has $(length(revs)) states at level $level: " * join(revs, ", ") *
          ". Ask for the one the consuming stage reads, by the property that identifies it.")
end

stage_rev(k::Capture, name::AbstractString, rev::Integer, level::Integer) =
    k.stages["$(name)_rev$(rev)_L$(level)"]

"""
    ladder(c::GoldenCase; n = 100) -> Vector{StageResult}

Walk the stages of one chip-size level in order, stopping at the first that disagrees.

The level is whichever one the capture traced. Level 0 is the base chip size, where the grid passes
through unresized and the interesting stages are the search-limit rewrite onward; a level above it
additionally exercises the grid resize, the mask dilation and the radius widening.
"""
function ladder(c::GoldenCase; n::Integer = 100, stop_on_red::Bool = true)
    k = read_capture(c; n)
    L = traced_level(k)
    chip0 = Int(k.scalars["ChipSize0X"])
    chip = chip0 << L
    spacing = Int(k.scalars["GridSpacingX"])
    minsearch = Int(k.scalars["minSearch"])
    @info "stage ladder" case=c.product level=L chip spacing minsearch

    # One `Params` for the whole ladder. Every rung configures AutoRIFT.jl from the same capture, and
    # `kwargs_from_capture` reduces `in_ChipSizeMaxX` twice per call — a multi-megapixel array on a
    # production grid — so building it per rung reduces the same array once per rung to reach the same
    # answer.
    p = params(; kwargs_from_capture(k)...)

    out = StageResult[]
    push!(out, rung_grid(k, p, L, chip, chip0))
    push!(out, rung_search_rewrite(k, L, minsearch))
    push!(out, rung_priors(k, p, L, chip, chip0, spacing))
    append!(out, rungs_coarse_sampling(k, p, L, chip))
    append!(out, rungs_coarse_correlation(k, p, L, chip))
    append!(out, rungs_filter_params(k, p, L))
    append!(out, rungs_coarse_mask(k, p, L, chip))
    append!(out, rungs_fill(k, p, L))
    append!(out, rungs_readback(k, p, L, chip, chip0))
    append!(out, rungs_merge(k, L, chip, chip0))

    if stop_on_red
        i = findfirst(r -> !r.passed, out)
        i === nothing || (out = out[1:i])
    end
    return out
end

# The grid state the correlator actually consumed, selected by the property that identifies it rather
# than by a revision number.
#
# Above the base chip size the reference rebinds `xGrid0` twice in three lines — `cv2.resize` builds it
# and the even-chip snap replaces it with `round(x + 0.5) - 0.5` (`autoRIFT.py:509-530`) — and it is the
# snapped grid the coarse and fine passes read. The two are distinguishable without knowing the order
# they were written in: **every non-zero coordinate of the snapped grid is a half-integer**, where the
# resize leaves quarter-fractions. `correlator.jl` already errors on a grid that fails this test, for the
# same reason and against the same convention.
#
# Selecting by the property rather than by `rev1` is what keeps this correct when the trace's dump order
# changes. Comparing against the pre-snap state instead reports 47.8% of the grid differing while the
# interpolation agrees to the bit — the difference is entirely `20.75` against `20.5`.
function _consumed_grid(k::Capture, name::AbstractString, L::Int)
    cands = String[]
    key = "$(name)_L$L"
    haskey(k.stages, key) && push!(cands, key)
    append!(cands, sort!([kk for kk in keys(k.stages)
                         if occursin(Regex("^$(name)_rev\\d+_L$L\$"), kk)],
                        by = kk -> parse(Int, match(r"rev(\d+)_", kk).captures[1])))
    isempty(cands) && error("no `$name` at level $L in this capture")
    for kk in cands
        nz = filter(!iszero, vec(Float64.(k.stages[kk])))
        isempty(nz) && continue
        all(≈(0.5), nz .- floor.(nz)) && return k.stages[kk]
    end
    # No candidate is on the half-integer convention. At the base chip size that is expected only if the
    # capture predates the grid rewrite, which `pointset_from_capture` already rejects; above it, it
    # means the snap did not run and the reader should know rather than silently compare.
    error("no state of `$name` at level $L is on the half-integer convention the correlator uses; " *
          "candidates " * join(cands, ", ") * ". Fractional parts of the last: " *
          string(sort(unique(let v = filter(!iszero, vec(Float64.(k.stages[last(cands)])))
                                 v .- floor.(v)
                             end))))
end

# The stride level `L` runs on, from the pyramid's own two functions rather than from a rule restated
# here. A rung that recomputes the reference's `chip ÷ chip0` compares the reference against a
# reimplementation of the reference, which agrees however production behaves; the oracle is the
# reference's traced array, so the stride should come from the code under test.
_level_stride(p::Params, L::Int) = AutoRIFT._level_decimation(p, AutoRIFT.chip_sizes(p)[L + 1])

# ---------------------------------------------------------------------------
# 3.1 — the level's grid
# ---------------------------------------------------------------------------
#
# At the base chip size `xGrid0` is `self.xGrid` copied (`autoRIFT.py:587-593`), so this rung asserts
# the capture is self-consistent: the traced grid must be the one `runAutorift` was handed. That is
# worth a rung of its own rather than an assumption, because a capture taken before `runAutorift`
# rewrote the grid holds an integer grid rather than `round(x) + 0.5`, and comparing against it moves
# every search centre half a pixel — a residual that is zero under uniform motion and grows with the
# velocity gradient, so it hides in every summary statistic.
#
# Above the base size the reference resizes with `INTER_AREA` and snaps an even chip's grid to
# `round(x + 0.5) - 0.5` (`autoRIFT.py:509-530`), where AutoRIFT.jl decimates by taking every
# `stride`-th point and taking the cell mean (`_decimate_level`, `_cell_means`). Those two
# reach the same place by different routes, so the rung compares the positions rather than the method.
function rung_grid(k::Capture, p::Params, L::Int, chip::Int, chip0::Int)
    xg0 = _consumed_grid(k, "xGrid0", L)
    if L == 0
        return exact_stage("3.1 grid, base level", "xGrid0", k.arrays["in_xGrid"], xg0)
    end
    # **AutoRIFT.jl decimates where the reference resizes, and the two reach different positions by
    # design.** The reference builds its grid with `INTER_AREA` and snaps to `round(x + 0.5) - 0.5`,
    # which averages `stride` coordinates and lands on the cell's centre; `_decimate_level` takes every
    # `stride`-th point and takes the cell mean (`_cell_means`). Those agree in the interior
    # and cannot agree at the grid's zero margin, where the reference averages real coordinates with the
    # zeros the driver wrote and AutoRIFT.jl does not.
    #
    # So the gate is the interior, and the margin is counted rather than tolerated: a position
    # difference there changes which pixels a chip covers, and `src/multichip.jl` records that
    # self-consistency between correlation position and read-back is what the accuracy depends on —
    # not agreement with either half of the reference separately.
    full = pointset_from_capture(k)
    stride = _level_stride(p, L)
    sub = AutoRIFT._decimate_level(full, trues(size(full)), stride)
    sub === nothing && return StageResult("3.1 grid, level $L", "xGrid0", "exact", false, 0,
                                          "AutoRIFT.jl decimated to nothing at stride $stride")
    jl = Float32.(sub.grid.x .- 1)
    # **Gated, not reported.** A level whose grid is a different size from the reference's is not a
    # position difference to characterise — it means the two are correlating a different number of
    # points, and every downstream rung on this level is then comparing arrays that do not correspond.
    size(jl) == size(xg0) || return StageResult("3.1 grid, level $L", "xGrid0", "exact", false, 0,
        "shape $(size(jl)) against reference $(size(xg0)) at stride $stride — AutoRIFT.jl's " *
        "decimation and the reference's resize disagree about the level's grid size")
    # **Above the base chip size the two constructions differ by design, and this rung reports the
    # difference rather than gating on it.** The reference resizes with `INTER_AREA` and snaps to
    # `round(x + 0.5) - 0.5`; `_decimate_level` takes every `stride`-th point and `_cell_means` averages
    # it to the cell centre. Three things separate them and only the first is arithmetic noise:
    #
    #   * On a **rotated** grid — `x` varies by 1 px per row on this scene — the reference's block mean is
    #     not the x-centre of the column pair, and the snap moves it a further half pixel. Measured at
    #     level 1: the cell centre is 3992.5 on both sides, the block mean exactly 3993.0, the snapped
    #     node 3993.5. 38% of nodes carry that one-pixel offset and 62% carry none.
    #   * At the **nodata margin** the reference averages the zeros the driver wrote
    #     (`testautoRIFT.py:394-403`) with real coordinates, so its node is a *fraction* of the
    #     coordinate — 112.5 where the cell holds 224.5. Not a position error on either side; the two
    #     are simply not describing the same thing there.
    #   * The offset **grows with stride**, because both effects scale with the cell.
    #
    # `src/multichip.jl` records why AutoRIFT.jl does not follow: `_undecimate_level` reads a coarse node
    # back from the cell centre, so moving the correlation position without moving the read-back measures
    # the field in one place and attributes it to another, and matching one of the reference's two halves
    # alone measured worse than matching neither. Self-consistency is the property that carries the
    # accuracy, and the coarse-level statistics that matter are gated at rung 3.17 instead.
    #
    # So this reports the median offset in cells, which is the number a reader should watch: it is stable
    # at half a cell by construction, and a *change* in it means the decimation or the read-back moved.
    offs = Float64[]
    for i in eachindex(jl, xg0)
        iszero(xg0[i]) && continue
        push!(offs, abs(Float64(xg0[i]) - Float64(jl[i])))
    end
    cell = stride * _grid_step_of(full)
    med = isempty(offs) ? 0.0 : median(offs)
    exact = count(iszero, offs)
    # Reported, not gated: the difference is a documented decision, so a pass/fail here would either
    # always fail or encode a tolerance nobody derived. Rung 3.17 gates the coarse levels on what is
    # actually comparable.
    return StageResult("3.1 grid, level $L", "xGrid0", "reported (see the comment)", true, length(offs),
                       @sprintf("median offset %.4g px = %.3f of a %.0f px cell; exact at %d of %d nodes",
                                med, med / cell, cell, exact, length(offs)))
end

# ---------------------------------------------------------------------------
# 3.5 — the search-limit rewrite
# ---------------------------------------------------------------------------
#
# `autoRIFT.py:596-602` rewrites the limits at the top of every level, and it is the rewritten array
# the correlator receives:
#
#     idxZero = (SearchLimitX0 <= 0) | (SearchLimitY0 <= 0)
#     SearchLimitX0[idxZero] = 0; SearchLimitY0[idxZero] = 0
#     SearchLimitX0[~idxZero & (SearchLimitX0 < minSearch)] = minSearch
#
# Two rules, each of which changes the window searched: either axis zero zeroes both, and a nonzero
# radius below `minSearch` is raised to it. `sanitize!` is AutoRIFT.jl's form of the same rewrite.
#
# The rung compares `rev1` against `rev0` — the value after the rewrite against the value before —
# because that pair isolates the rewrite from whatever produced its input. Reading `rev0` as the
# correlator's radius is the error this rung exists to prevent: it costs 944,036 points on this case,
# by up to 23 pixels, and biases the comparison opposite to its cause, since the rewritten points are
# the small-radius ones.
function rung_search_rewrite(k::Capture, L::Int, minsearch::Int)
    # The two states this rung is about, chosen by *shape* rather than by revision number. Above the
    # base chip size the resize rebinds the name, so `rev0` still holds the previous level's full-grid
    # array and the first value this level computed is a later revision — the level's own grid size is
    # what identifies them.
    # **The two states are identified by content, not by position.** The rewrite raises every nonzero
    # radius to at least `minSearch`, so the state *before* it has nonzero values below `minSearch` and
    # the state *after* does not — and the pair to compare is the last pre-rewrite state with the first
    # post-rewrite one. Selecting by revision index instead picks whichever the trace happened to write
    # in those slots: at level 1 the sequence is rev1 and rev2 pre-rewrite, rev3 and rev4 after, so
    # `(rev1, rev2)` compares two pre-rewrite states and reports 22.6% differing.
    want = size(_consumed_grid(k, "xGrid0", L))
    revs = [kk for kk in _revisions(k, "SearchLimitX0", L) if size(k.stages[kk]) == want]
    below(kk) = any(v -> 0 < v < minsearch, k.stages[kk])
    pre = findlast(below, revs)
    post = pre === nothing ? nothing : findfirst(i -> !below(revs[i]), (pre + 1):length(revs))
    (pre === nothing || post === nothing) && return StageResult(
        "3.5 search-limit rewrite", "SearchLimitX0", "exact", false, 0,
        "no pre/post-`minSearch` pair among " *
        join(["$kk min-nonzero $(minimum(filter(!iszero, vec(k.stages[kk])); init = 0.0f0))"
              for kk in revs], ", "))
    after_key = revs[pre + post]
    before_x = k.stages[revs[pre]]
    after_x = k.stages[after_key]
    # `SearchLimitY0` is rewritten in the same statement, so its matching state is the one at the same
    # index — but its revision count can differ, since the trace numbers each name independently.
    yrevs = [kk for kk in _revisions(k, "SearchLimitY0", L) if size(k.stages[kk]) == want]
    ybelow(kk) = any(v -> 0 < v < minsearch, k.stages[kk])
    ypre = findlast(ybelow, yrevs)
    before_y = k.stages[ypre === nothing ? first(yrevs) : yrevs[ypre]]

    pts = PointSet(zeros(want), zeros(want), Int.(before_x), Int.(before_y),
                   zeros(want), zeros(want), fill(1, want), fill(1, want),
                   zeros(Int, want), zeros(Int, want))
    sanitize!(pts, minsearch)
    return exact_stage("3.5 search-limit rewrite", after_key, Float32.(pts.radius_x), after_x)
end

# ---------------------------------------------------------------------------
# 3.3 / 3.4 — a decimated cell's radius and prior
# ---------------------------------------------------------------------------
#
# A decimated point stands for its whole cell, so its window must cover every fine point in it: both
# the widest radius any of them asked for and the spread of their priors, since two fine points with
# different priors search around different centres. The reference adds exactly those two terms over
# `1 / Scale` cells (`autoRIFT.py:540-586`):
#
#     SearchLimitX0 = colfilt(SearchLimitX, (1/Scale, 1/Scale), 0)   # max
#                   + colfilt(Dx0,          (1/Scale, 1/Scale), 4)   # range
#     Dx00          = colfilt(Dx0,          (1/Scale, 1/Scale), 2)   # mean
#
# At the base chip size there is no decimation and the three are copies, so this rung asserts that —
# a level that widened its own base radius would be searching a window the reference does not.
function rung_priors(k::Capture, p::Params, L::Int, chip::Int, chip0::Int, spacing::Int)
    dx00 = _integral_state(k, "Dx00", L)
    if L == 0
        return exact_stage("3.4 prior, base level", "Dx00", k.arrays["in_Dx0"], dx00)
    end
    # The reference reduces with `colfilt(Dx0, (1/Scale, 1/Scale), 2)` — a NaN-aware mean over the cell
    # — then resizes with `INTER_NEAREST` and rounds (`autoRIFT.py:568-586`). So the comparison is the
    # reduced field sampled at the level's nodes, and `windowmean` is the reducer Gate 2 pinned against
    # `colfilt` option 2.
    # The cell's first point, which is where `INTER_NEAREST` lands for this size ratio at all but one
    # destination on this level. `resample(..., Nearest())` is *not* interchangeable here — it disagrees
    # with OpenCV's mapping at 18,422 of 1,368,896 destinations against this slice's 1, so the slice is
    # the better model of `cv2.resize`'s nearest rule and the residual point is reported rather than
    # absorbed.
    # The pyramid's own stride, so this reduces over the cell it actually decimates to.
    stride = _level_stride(p, L)
    mx = windowmean(k.arrays["in_Dx0"], stride)
    rows = 1:stride:size(mx, 1)
    cols = 1:stride:size(mx, 2)
    jl = Float32.(round.([mx[i, j] for i in rows, j in cols]))
    size(jl) == size(dx00) || return StageResult("3.4 prior, level $L", "Dx00", "exact", false, 0,
        "shape $(size(jl)) against reference $(size(dx00)) — the trace recorded `Dx00` on a " *
        "different grid than this level's, so the two are not the same quantity")
    # Gated on a small share rather than on zero, and on the *size* of the difference as well as its
    # count: the residual is `cv2.resize`'s nearest-pixel mapping differing from a first-of-cell slice
    # at a handful of destinations, which moves the prior by exactly one pixel there. One pixel of prior
    # displaces the chip by one pixel, which the search radius covers; a larger difference, or a
    # systematic one, would not be this.
    r = exact_stage("3.4 prior, level $L", "Dx00", jl, dx00)
    bad = 0
    worst = 0.0
    for i in eachindex(jl, dx00)
        jl[i] == dx00[i] && continue
        bad += 1
        worst = max(worst, abs(Float64(jl[i]) - Float64(dx00[i])))
    end
    return StageResult(r.name, r.reference, "<=1e-4 of points, each <=1 px",
                       bad / length(jl) <= 1e-4 && worst <= 1.0, length(jl),
                       @sprintf("%d of %d differ (%.5f%%), largest by %.4g px",
                                bad, length(jl), 100bad / length(jl), worst))
end

# The grid's own spacing, as a mode over adjacent steps, so a
# production grid's zeroed nodata margin cannot make it read as zero.
_grid_step_of(pts::PointSet{2}) = AutoRIFT._grid_step(pts.x, 2)

# The `PointSet` this level's coarse pass runs on, built from the level's own traced arrays.
#
# At the base chip size that is the capture's own grid. Above it the reference has resized everything by
# `ChipSize0X / chip` (`autoRIFT.py:509-586`), and taking the full grid instead compares a 2344-wide
# array against an 1172-wide one — a `BoundsError` if you are lucky and a wrong answer if the shapes
# happen to admit the index.
#
# The radii come from the first revision whose shape is this level's, not from `rev0`: the resize
# rebinds the name, so at a level above the base `rev0` still holds the previous level's full-grid
# array. Selecting by shape states that requirement rather than encoding a revision count that shifts
# with the level.
function _level_pointset(k::Capture, L::Int, chip::Int, p::Params)
    full = pointset_from_capture(k)
    L == 0 && return AutoRIFT._level_points(full, p, extent(chip), trues(size(full)))

    xg = _consumed_grid(k, "xGrid0", L)
    yg = _consumed_grid(k, "yGrid0", L)
    # The post-`minSearch` state, which is what the coarse pass's `colfilt` reduces and what the
    # correlator receives. Identified by content — no nonzero value below `minSearch` — rather than by
    # revision index, for the reason rung 3.5 records.
    minsearch = Int(k.scalars["minSearch"])
    sx = _post_rewrite(k, "SearchLimitX0", L, size(xg), minsearch)
    sy = _post_rewrite(k, "SearchLimitY0", L, size(xg), minsearch)
    # The **rounded** prior, which is what the coarse slice and the fine pass both read
    # (`autoRIFT.py:585-586`). Identified by its values being integral: the pre-round state carries the
    # cell mean's quarter-fractions. The prior displaces the chip and `chip_bounds` floors the displaced
    # centre, so a half-pixel prior moves the chip a whole pixel.
    dx0 = _integral_state(k, "Dx00", L)
    dy0 = _integral_state(k, "Dy00", L)
    # `Dx00` is traced on the full grid here, so it is reduced onto the level's lattice the way the
    # reference builds it: a NaN-aware cell mean, `INTER_NEAREST` back to the level's shape, and
    # **rounded** (`autoRIFT.py:568-586`). The rounding is not cosmetic — the prior displaces the chip,
    # and `chip_bounds` floors the displaced centre, so a half-pixel prior moves the chip by a pixel.
    # Omitting it puts 360 of 21,316 coarse points a half pixel from where the reference put them.
    if size(dx0) != size(xg)
        stride = size(dx0, 1) ÷ size(xg, 1)
        rows = 1:stride:size(dx0, 1)
        cols = 1:stride:size(dx0, 2)
        dx0 = round.(windowmean(dx0, stride)[rows, cols][1:size(xg, 1), 1:size(xg, 2)])
        dy0 = round.(windowmean(dy0, stride)[rows, cols][1:size(xg, 1), 1:size(xg, 2)])
    end
    n = size(xg)
    # `dy0` is negated for the reason `pointset_from_capture` negates it: the trace records the prior
    # in cartesian-Y, and `arImgDisp_u`/`arImgDisp_s` flip it to matrix-Y before cutting a chip
    # (`autoRIFT.py:1058`, `:1231`). A `PointSet` carries the matrix-Y form, so passing the traced
    # value through unchanged places every chip `2 * Dy0` rows from where the reference placed it.
    return PointSet(
        Float64.(xg) .+ 1, Float64.(yg) .+ 1,
        Int.(sx), Int.(sy),
        Float64.(dx0), .-Float64.(dy0),
        fill(chip, n), fill(chip, n),
        zeros(Int, n), zeros(Int, n),
    )
end

# The state of `name` at this level whose shape is `want`, or `nothing`. Shape identifies a state
# wherever the reference decimates or resizes under one name — `SearchLimitX0C` is the `colfilt` output
# and then the coarse slice of it.
function _state_shaped(k::Capture, name::AbstractString, L::Int, want::Tuple{Int,Int})
    for kk in _revisions(k, name, L)
        size(k.stages[kk]) == want && return k.stages[kk]
    end
    return nothing
end

# The state of `name` at this level whose values are all integers, which is the reference's `round`ed
# form. `Dx00` is built by a cell mean and then rounded (`autoRIFT.py:585-586`), so the two states are
# told apart by their fractional parts rather than by which was written first.
function _integral_state(k::Capture, name::AbstractString, L::Int)
    revs = _revisions(k, name, L)
    isempty(revs) && error("no `$name` at level $L")
    for kk in revs
        all(v -> isnan(v) || v == round(v), k.stages[kk]) && return k.stages[kk]
    end
    # At the base chip size the prior is copied unchanged and may legitimately be fractional, so the
    # last state is the answer there rather than an error.
    return k.stages[last(revs)]
end

# The first state of `name` at this level, on grid `want`, that the `minSearch` rewrite has already
# been applied to: no nonzero value below `minSearch` (`autoRIFT.py:596-602`). That is the array the
# coarse reduction and the correlator both read.
function _post_rewrite(k::Capture, name::AbstractString, L::Int, want::Tuple{Int,Int}, minsearch::Int)
    revs = [kk for kk in _revisions(k, name, L) if size(k.stages[kk]) == want]
    isempty(revs) && error("no state of `$name` at level $L on grid $want")
    for kk in revs
        any(v -> 0 < v < minsearch, k.stages[kk]) || return k.stages[kk]
    end
    error("every state of `$name` at level $L still has a nonzero radius below minSearch=$minsearch; " *
          "the rewrite was expected to have run before the coarse pass")
end

# The first revision of `name` at this level whose shape is `want`, or the unversioned array.
function _level_revision(k::Capture, name::AbstractString, L::Int, want::Tuple{Int,Int})
    key = "$(name)_L$L"
    haskey(k.stages, key) && size(k.stages[key]) == want && return k.stages[key]
    revs = sort!([kk for kk in keys(k.stages) if occursin(Regex("^$(name)_rev\\d+_L$L\$"), kk)],
                 by = kk -> parse(Int, match(r"rev(\d+)_", kk).captures[1]))
    for kk in revs
        size(k.stages[kk]) == want && return k.stages[kk]
    end
    error("no revision of `$name` at level $L has shape $want; have " *
          join(["$kk $(size(k.stages[kk]))" for kk in revs], ", "))
end

# ---------------------------------------------------------------------------
# 3.6 — the coarse pass's sample lattice, its radii, and its priors
# ---------------------------------------------------------------------------
#
# The coarse pass samples every `sparseSearchSampleRate * ratio`-th point of the level's grid
# (`autoRIFT.py:603-614`), and each sampled point stands in for the cell around it — so its radius is a
# maximum over a window rather than its own value, and its prior is the cell's mean.
#
# **The window is one wider than the step when the step is even.** `filtWidth = stride + 1` for an even
# stride and `stride` for an odd one (`autoRIFT.py:618-626`), so the reduction is symmetric about the
# node it is sampled at. Reducing over the stride instead under-covers a coarse point's window, and
# always downward: at stride 8 it cost 1,809 of 85,556 coarse points on this case, by up to 152 pixels.
# A radius too small searches a narrower window than the reference did, so it rails out or misses the
# peak at exactly the points where the prior was doing work — which reads as a correlator disagreement.
#
# This rung is *not* crop-safe: which fine points the lattice lands on depends on the grid's extent, so
# a sub-window samples a different set and the comparison would be against the wrong nodes.
function rungs_coarse_sampling(k::Capture, p::Params, L::Int, chip::Int)
    # **This level's grid, not the full one.** Above the base chip size the reference resizes by
    # `ChipSize0X / chip` and everything after that runs on the resized grid — at chip 32 the grid is
    # 1172x1168 against the full 2344x2336. So the points handed to the coarse setup have to be this
    # level's, and the level's own traced arrays are what supply them: the resized grid, and the radii
    # and priors as the level rewrote them.
    #
    # `SearchLimitX0_rev0` is the *leftover* full-grid value at a level above the base, because the
    # resize rebinds the name — `rev1` is the first value this level computed. Selecting by shape rather
    # than by revision number states the requirement instead of encoding a count that changes with the
    # level.
    lp = _level_pointset(k, L, chip, p)
    setup = AutoRIFT._coarse_points(lp, p, extent(chip))
    setup === nothing && return [StageResult("3.6 coarse sampling", "xGrid0C", "exact", false, 0,
                                            "AutoRIFT.jl found no coarse grid at this level")]

    out = StageResult[]
    # The lattice, by the positions it lands on. `+1`/`-1` is the index base and nothing else.
    push!(out, exact_stage("3.6a coarse lattice", "xGrid0C",
                           Float32.(setup.coarse.x .- 1), stage(k, "xGrid0C", L)))
    # The radii, which are the reduction the width rule above governs. The reference dumps the
    # *undecimated* `colfilt` output under this name, so the comparison decimates it on the lattice
    # just verified rather than assuming the two lattices coincide.
    # `SearchLimitX0C` has two states: the undecimated `colfilt` output at `:629`, then the decimated
    # slice at `:640` that the coarse pass actually receives. The consumer wants the second, and it is
    # identified by shape — the coarse grid's — rather than by which was written first.
    for (axis, refname, jl) in (("x", "SearchLimitX0C", setup.coarse.radius_x),
                                ("y", "SearchLimitY0C", setup.coarse.radius_y))
        ref = _state_shaped(k, refname, L, size(jl))
        if ref === nothing
            # Only the undecimated state was traced, so decimate it on the lattice rung 3.6a verified.
            full = last(_revisions(k, refname, L))
            ref = k.stages[full][setup.rows, setup.cols]
        end
        push!(out, exact_stage("3.6b coarse radius $axis", refname, Float32.(jl), ref))
    end
    # The priors, sampled at the node rather than reduced: the reference slices `Dx00` with the same
    # `rIdxC` it slices the grid with (`autoRIFT.py:643-644`), so no window is involved here.
    #
    # **The y prior is compared in cartesian-Y, which is the convention the trace is in.** A `PointSet`
    # holds the matrix-Y form, so the sign is put back rather than the two being compared across the
    # flip — an axis whose only difference is the convention would otherwise report every nonzero point
    # as a disagreement while the chips it places are identical.
    for (axis, refname, jl) in (("x", "Dx0C", setup.coarse.dx_prior),
                                ("y", "Dy0C", .-setup.coarse.dy_prior))
        ref = _state_shaped(k, refname, L, size(jl))
        ref === nothing && continue
        push!(out, exact_stage("3.6c coarse prior $axis", refname, Float32.(jl), ref))
    end
    return out
end

# ---------------------------------------------------------------------------
# 3.8 / 3.13 — the outlier filter's parameters
# ---------------------------------------------------------------------------
#
# `autorift()` derives two filters from one parameter set (`autoRIFT.py:484-505`): `DispFiltC` for the
# coarse pass, loosened for the decimation, and `DispFiltF` for the fine pass. AutoRIFT.jl's `relax`
# and `rescale` reproduce those derivations, and this rung checks the *results* against the reference's
# own `filtDisp` records rather than re-deriving the formula on this side.
#
# Worth a rung because a parameter mismatch and a reducer mismatch are indistinguishable in the mask
# they produce, and a comparison of masks alone would attribute one to the other. The reference reports
# `FiltWidth`, `FracValid` and `Iter` per call, so this is a direct read.
function rungs_filter_params(k::Capture, p::Params, L::Int)
    out = StageResult[]
    # **A level's `filtDisp` calls are found by grid shape, not by counting two per level.** A level
    # whose coarse pass falls below `CoarseCorCutoff` `continue`s out (`autoRIFT.py:704-706`) and never
    # runs a fine pass, so it contributes *one* record — and every later level's index shifts. On the
    # golden L7 pair the base level does exactly that: 102 of 20,525 coarse points survive, 0.50%
    # against a 1% cutoff, and both implementations agree to skip it. Indexing at `2L + 1` then compares
    # the *coarse* filter against the fine filter's parameters and reports a difference that is entirely
    # the harness's.
    #
    # The coarse record is on the level's decimated coarse grid and the fine record on its level grid, so
    # the shapes name them. `nothing` for a pass that did not run, and the rung says so rather than
    # inventing a comparison.
    # **Shape alone is ambiguous across levels, so the search is bounded to this level's own records.**
    # A coarse grid at one level can have the same shape as a *fine* grid several levels up: on the golden
    # Landsat case level 0's coarse grid is 293x292 and level 3's fine grid is 293x292 too, since both are
    # the full grid divided by 8. Searching the whole list picks whichever matched last and compares the
    # coarse filter against fine parameters.
    #
    # `k.levels` is in call order, so this level's records are the ones between the correlator calls that
    # carry its chip size — and the coarse record precedes the fine one within that span.
    coarse_shape = size(_state_coarse(k, "xGrid0C", L))
    level_shape = size(_consumed_grid(k, "xGrid0", L))
    chip_here = Float64(Int(k.scalars["ChipSize0X"]) << L)
    span = findall(r -> r.kind != "filtDisp" && r.chip_size[1] == chip_here, k.levels)
    lo = isempty(span) ? 1 : first(span)
    # Up to the next level's first correlator call, so a shape shared with a later level cannot be reached.
    nxt = findfirst(i -> i > lo && k.levels[i].kind != "filtDisp" &&
                         k.levels[i].chip_size[1] != chip_here, eachindex(k.levels))
    hi = nxt === nothing ? length(k.levels) : nxt - 1
    mine = k.levels[lo:hi]
    calls = [r for r in mine if r.kind == "filtDisp"]
    coarse_call = findfirst(r -> r.grid_shape == coarse_shape, calls)
    fine_call = findfirst(r -> r.grid_shape == level_shape, calls)
    coarse_ratio = AutoRIFT._oversample(p)
    if fine_call === nothing
        push!(out, StageResult("3.8 filter parameters, fine", "filtDisp", "level ran no fine pass",
                               true, 0,
                               "the coarse pass fell below CoarseCorCutoff, so this level was skipped " *
                               "on both sides — rung 3.9 is where that agreement is gated"))
    end
    for (label, ci, filt) in
        (("coarse", coarse_call, rescale(relax(p.outliers), coarse_ratio, p.coarse_stride)),
         ("fine", fine_call, rescale(p.outliers, coarse_ratio)))
        ci === nothing && continue
        call = calls[ci]
        got = (window(filt), filt.min_agree_fraction, filt.iterations)
        want = (call.counts["filt_width"], call.oversample, call.counts["iterations"])
        ref_frac = call.frac_valid
        ok = got[1] == want[1] && got[3] == want[3] && isapprox(got[2], ref_frac; atol = 1e-6)
        push!(out, StageResult("3.8 filter parameters, $label", "filtDisp", "exact", ok, 1,
                               @sprintf("width %d/%d, frac %.4f/%.4f, iterations %d/%d (julia/reference)",
                                        got[1], want[1], got[2], ref_frac, got[3], want[3])))
    end
    return out
end

# ---------------------------------------------------------------------------
# 3.14 / 3.15 — the median field and the three-pass hole fill
# ---------------------------------------------------------------------------
#
# The reference nulls its rejected points, takes one `fillFiltWidth`-wide median of what remains, and
# then fills for three passes from *that* field (`autoRIFT.py:764-808`):
#
#     DxF[~M0] = nan
#     DxFM = colfilt(DxF, (fillFiltWidth, fillFiltWidth), 3)      # once, before the loop
#     MM = ~isnan(DxFM)
#     for j in range(3):
#         foo  = MF | M0
#         foo1 = (filter2D(foo, ones(3,3)) >= 6) | foo             # area closing
#         fillIdx = ~bwareaopen(~foo1, 5) & ~foo & MM              # or a small component
#         MF[fillIdx] = True; DxF[fillIdx] = DxFM[fillIdx]
#
# Both rungs feed AutoRIFT.jl the reference's own inputs: `DxF` from the level record and `M0` from the
# `filtDisp` record, so the rejection is taken from the reference and only the fill is under test.
#
# **`DxFM` is computed once, outside the loop.** Every pass therefore fills from the median of the
# *original* field, and a point filled in pass 1 does not change what pass 2 fills with. `_fill_holes!`
# recomputes its median per pass over the progressively filled field, so a second-pass fill there sees
# first-pass values. That is a real difference in the values filled — not in which points are filled —
# and the rung measures it rather than asserting either behaviour.
function rungs_fill(k::Capture, p::Params, L::Int)
    out = StageResult[]

    # `DxF`, `DyF` and `MF` are mutated in place by the fill, so the trace numbers their revisions and a
    # consumer has to name the one it wants. `rev0` is the state before the fill: the raw correlator
    # output for `DxF`, and all-zero for `MF`.
    # The raw correlator output on this level's grid, which is the first state: the fill overwrites
    # entries in place afterwards, and above the base level the merge resizes it to the full grid.
    dxf_states = [kk for kk in _revisions(k, "DxF", L)]
    isempty(dxf_states) && return out
    dxf_key = first(dxf_states)

    # **`rev0` of `DxF` is the value *before* `DxF[~M0] = nan`.** `DxF` is bound by the correlator at
    # `:735`, so the first state on disk is the raw fine measurement, while the reference medians the
    # field *after* `filtDisp` has nulled the rejected points (`autoRIFT.py:764-773`). Medianing the raw
    # field compares two different inputs: it reports 23.6% of `DxFM` differing, entirely because
    # AutoRIFT.jl sees neighbours the reference has already thrown away.
    #
    # The `filtDisp` record supplies the mask that closes the gap. Two calls per resolved level, coarse
    # first, so the fine pass is the second — and its `kept` count is checked against the nulling to
    # confirm the right record was taken rather than assumed.
    dxf = copy(k.stages[dxf_key])
    dyf = copy(k.stages[replace(dxf_key, "Dx" => "Dy")])

    # The `filtDisp` record for *this level's fine pass*, chosen by grid shape rather than by a
    # positional index. A level that `continue`d out before its fine pass contributes fewer records than
    # the two a resolved level does, so `2L + 2` is not reliably its index — and the fine pass is the one
    # whose grid matches this level's `DxF`.
    fine = [r for r in k.levels if r.kind == "filtDisp" && r.grid_shape == size(dxf)]
    isempty(fine) && return [StageResult("3.13 fine rejection", "filtDisp", "exact", false, 0,
        "no filtDisp record on this level's $(size(dxf)) grid; shapes present: " *
        join(unique(r.grid_shape for r in k.levels if r.kind == "filtDisp"), ", "))]
    fine_filt = last(fine)
    kept = fine_filt.arrays["kept"] .!= 0
    dxf[.!kept] .= NaN32
    dyf[.!kept] .= NaN32
    push!(out, StageResult("3.13 fine rejection", "filtDisp kept", "count matches",
                           count(!isnan, dxf) == fine_filt.counts["kept"],
                           count(kept),
                           @sprintf("kept %d of %d; nulled field has %d valid",
                                    fine_filt.counts["kept"], fine_filt.counts["in_mask"],
                                    count(!isnan, dxf))))

    # The reference's `MM`: where its one-shot median exists. This gates every fill, so a difference
    # here bounds everything downstream.
    fillw = Int(get(k.scalars, "fillFiltWidth", 3))
    jl_dxfm = windowmedian(dxf, fillw)
    # **The *first* `DxFM` on this level's grid.** Above the base chip size the name is rebound twice
    # more during the merge: `:843` recomputes it at width 5 to fill the level's remaining holes, and
    # `:851` overwrites entries from the previous levels' answer. The fill loop reads the `fillFiltWidth`
    # median built at `:776`, which is the first.
    dxfm_key = first([kk for kk in _revisions(k, "DxFM", L) if size(k.stages[kk]) == size(dxf)])
    mm_key = first([kk for kk in _revisions(k, "MM", L) if size(k.stages[kk]) == size(dxf)])
    push!(out, exact_stage("3.14 fill median (width $fillw)", dxfm_key, jl_dxfm, k.stages[dxfm_key]))
    push!(out, exact_stage("3.14 fill median gate", mm_key,
                           UInt8.(map(!isnan, jl_dxfm)), k.stages[mm_key]))

    # The fill itself, on the reference's own rejected field. `_fill_holes!` mutates, so it gets a copy.
    d = DisplacementField(copy(dxf), copy(dyf),
                          fill(NaN32, size(dxf)), fill(NaN32, size(dxf)),
                          map(!isnan, dxf))
    filled = AutoRIFT._fill_holes!(d, p)
    jl_mf = falses(size(dxf))
    for i in filled
        jl_mf[i] = true
    end
    # The *last* `MF` revision, which is the state after the three fill passes. `rev0` is the all-zero
    # initialization at `:790`, and comparing against it reports that the reference filled nothing.
    # The last `MF` state *on this level's grid*. `rev0` is the all-zero initialization at `:790`, and
    # above the base chip size the final revision is the merge's `INTER_NEAREST` resize back to the full
    # grid (`:861`) — a different quantity from the mask the fill produced.
    mine = [kk for kk in _revisions(k, "MF", L) if size(k.stages[kk]) == size(dxf)]
    isempty(mine) && return out
    push!(out, exact_stage("3.15 fill mask", last(mine), UInt8.(jl_mf), k.stages[last(mine)]))
    return out
end

# ---------------------------------------------------------------------------
# 3.7 — the coarse correlation, on both element types
# ---------------------------------------------------------------------------
#
# **Every rung that runs a correlator runs it twice, on `UInt8` and on `Float32`, and the pair is the
# measurement rather than either half.** The reference has two correlators — `arImgDisp_u` on bytes and
# `arImgDisp_s` on floats, separate C++ templates — and production reaches the byte one, because
# `uniform_data_type` rescales each scene by its own mean and standard deviation and quantizes to 256
# levels before `runAutorift` is called (`autoRIFT.py:359-384`).
#
# What the pair separates is a difference in what the correlator *computes* from a tie the quantization
# created in the surface it computes *on*. Collapsing a filtered float field onto 256 levels makes
# plateaus the float field does not have, and a plateau broken differently puts the peak far away rather
# than one step away — which is why the byte path's disagreements are measured in whole pixels while the
# float path's are bit-exact. A byte-only comparison reports the sum of the two and cannot apportion it.
#
# **The reference against itself is the floor.** Handing the two templates the same information — the
# captured bytes, and those same bytes widened to `Float32` — gives a disagreement that belongs to
# neither implementation. At chip 96 on the S2A case that floor is 98.54% exact with a maximum of 7 px,
# and AutoRIFT.jl sits *on* it against both paths. A rung cannot be asked to do better than the
# reference reproducing itself, so the floor is what it is gated against.
#
# The grid is the reference's own captured coarse grid, not one rebuilt by `_coarse_points`: rebuilding
# it measures the setup as well, which rung 3.6 already does separately, and on this case that
# conflation read as 85% where the correlator alone reads 98.5%.
function rungs_coarse_correlation(k::Capture, p::Params, L::Int, chip::Int)
    out = StageResult[]
    xg = _state_coarse(k, "xGrid0C", L)
    xg === nothing && return out
    yg = _state_coarse(k, "yGrid0C", L)
    srx = _state_coarse(k, "SearchLimitX0C", L, size(xg))
    sry = _state_coarse(k, "SearchLimitY0C", L, size(xg))
    dx0 = _state_coarse(k, "Dx0C", L, size(xg))
    refdx = _state_coarse(k, "DxC", L, size(xg))
    any(isnothing, (yg, srx, sry, dx0, refdx)) && return out
    dy0 = _state_coarse(k, "Dy0C", L, size(xg))

    n = size(xg)
    # `Dy0C` is cartesian-Y, as traced; a `PointSet` carries matrix-Y. See `_level_pointset`.
    pts = PointSet(Float64.(xg) .+ 1, Float64.(yg) .+ 1, Int.(srx), Int.(sry),
                   Float64.(dx0), .-Float64.(dy0), fill(chip, n), fill(chip, n),
                   zeros(Int, n), zeros(Int, n))
    a = k.arrays["in_I1"]
    b = k.arrays["in_I2"]
    level = findfirst(r -> r.kind == "coarse" && r.chip_size[1] == Float64(chip), k.levels)

    for (label, aa, bb) in (("UInt8", a, b), ("Float32", Float32.(a), Float32.(b)))
        pair = AutoRIFT._prepare(ImagePair(bb, aa), p)
        cd = AutoRIFT.run_pass(AutoRIFT.WholeScene(pair), deepcopy(pts), p,
                               measure_at(p, L + 1), AutoRIFT.NoRefine())
        both = 0; oj = 0; orf = 0; ex = 0
        d = Float64[]
        for i in eachindex(cd.dx, refdx)
            mj = isnan(cd.dx[i]); mr = isnan(refdx[i])
            mj && mr && continue
            if mr
                oj += 1
            elseif mj
                orf += 1
            else
                both += 1
                δ = Float64(cd.dx[i]) - Float64(refdx[i])
                push!(d, δ)
                δ == 0 && (ex += 1)
            end
        end
        ad = abs.(d)
        frac = both == 0 ? 0.0 : ex / both
        # Coverage is the gate: a coarse node the reference measured and this did not is a decision, and
        # `exact` is reported beside it because the byte path's ties are not a position disagreement.
        # The reference's own record says how many it measured, so a silent shape or grid error cannot
        # pass as agreement.
        refn = level === nothing ? count(!isnan, refdx) : k.levels[level].counts["measured"]
        # **A reference-only node whose value is exactly `-radius_x` is the degenerate-chip corner, not a
        # measurement AutoRIFT.jl failed to make.** A constant chip carries no information about
        # displacement, and the reference returns the search window's corner there — `dx = -radius_x`,
        # `dy = +radius_y` — where AutoRIFT.jl reports nothing. `tools/golden/README.md` registers that as
        # matched-for-agreement-not-endorsed, in the direction where AutoRIFT.jl deliberately differs.
        #
        # Counted separately rather than tolerated: it is a fabricated value over masked and low-texture
        # ground, so folding it into the coverage gate would make the gate reward reproducing it. At chip
        # 32 on the golden Landsat case 10 of 13 reference-only nodes are exactly `-radius_x` and the other
        # 3 are one off it, at a neighbouring radius.
        fabricated = 0
        for i in eachindex(cd.dx, refdx)
            isnan(cd.dx[i]) && !isnan(refdx[i]) || continue
            r = Float64(pts.radius_x[i])
            abs(Float64(refdx[i]) + r) <= 1 && (fabricated += 1)
        end
        push!(out, StageResult("3.7 coarse correlation, $label",
                               "DxC (reference measured $refn)",
                               "coverage identical but for degenerate chips",
                               oj == 0 && orf == fabricated && count(!isnan, refdx) == refn, both,
                               @sprintf("both %d, only jl %d, only ref %d (%d the degenerate-chip corner); exact %.2f%%, median %.4g, p99 %.4g, max %.4g, bias %+.4g",
                                        both, oj, orf, fabricated, 100frac,
                                        isempty(ad) ? 0.0 : median(ad),
                                        isempty(ad) ? 0.0 : quantile(ad, 0.99),
                                        isempty(ad) ? 0.0 : maximum(ad),
                                        isempty(d) ? 0.0 : mean(d))))
    end
    return out
end

# A coarse-grid state of `name` at this level: the one on the coarse grid, which is the smallest shape
# the level traced. `want` pins it when the name has states on more than one grid — `SearchLimitX0C` is
# dumped both undecimated and decimated.
function _state_coarse(k::Capture, name::AbstractString, L::Int, want = nothing)
    revs = _revisions(k, name, L)
    isempty(revs) && return nothing
    want === nothing || return _state_shaped(k, name, L, want)
    return k.stages[argmin(kk -> prod(size(k.stages[kk])), revs)]
end

# ---------------------------------------------------------------------------
# 3.9 / 3.10 / 3.11 — the coarse mask, and the fine search it restricts
# ---------------------------------------------------------------------------
#
# This is the stage that decides *where the fine pass is allowed to look*, and therefore the one that
# owns a coverage difference the correlator cannot explain. Three steps (`autoRIFT.py:702-724`):
#
#     ROIC = SearchLimitX0C > 0
#     CoarseCorValidFac = sum(MC[ROIC]) / sum(M0C[ROIC])     # the level's go/no-go gate
#     if CoarseCorValidFac < CoarseCorCutoff: continue
#     MC2 = distance_transform_edt(!MC) < BuffDistanceC      # dilate on the *coarse* grid
#     MC2 = resize(MC2, coarse_shape * stride, INTER_NEAREST)
#     ... pad by edge replication if that falls short of the fine grid ...
#     SearchLimitX0[!MC2] = 0                                # and this is what the fine pass sees
#
# The last line is the one that matters: a point outside `MC2` has its radius zeroed and is never
# searched, so **every measurement lost here is lost before the correlator runs**. `M0C` and `MC` are
# already compared at rungs 3.7 and 3.8; what this adds is the dilation, the expansion, and the
# resulting radius — the three places a matching mask can still restrict a different region.
function rungs_coarse_mask(k::Capture, p::Params, L::Int, chip::Int)
    out = StageResult[]
    haskey(k.stages, "MC_rev0_L$L") || return out
    mc = _last_revision(k, "MC", L)
    mc === nothing && return out
    keep = k.stages[mc] .!= 0

    # The go/no-go gate, on the reference's own `MC` and `M0C` so only the ratio is under test. A level
    # the reference dropped and AutoRIFT.jl kept — or the reverse — is a whole chip size of coverage.
    m0c = _last_revision(k, "M0C", L)
    if m0c !== nothing
        measured = k.stages[m0c] .!= 0
        srx = _state_shaped(k, "SearchLimitX0C", L, size(keep))
        if srx !== nothing
            roic = srx .> 0
            denom = count(roic .& measured)
            numer = count(roic .& keep)
            frac = denom == 0 ? 0.0 : numer / denom
            cutoff = p.min_coarse_valid_fraction
            push!(out, StageResult("3.9 coarse valid fraction", mc,
                                   string("same verdict at cutoff ", cutoff),
                                   (frac >= cutoff) == (frac >= Float64(k.scalars["CoarseCorCutoff"])),
                                   denom,
                                   @sprintf("%d of %d coarse points kept = %.4f; cutoff julia %.4g, reference %.4g; both %s",
                                            numer, denom, frac, cutoff,
                                            Float64(k.scalars["CoarseCorCutoff"]),
                                            frac >= cutoff ? "continue" : "skip the level")))
        end
    end

    # The dilation, on the coarse grid. `dilate_within` is AutoRIFT.jl's form of
    # `distance_transform_edt(!MC) < BuffDistanceC`, and `test/fixtures/disttransform` pins the transform
    # itself — so a difference here is the threshold or the buffer, not the metric.
    grown = AutoRIFT.dilate_within(keep, p.coarse_buffer)
    ref_mc2 = _state_shaped(k, "MC2", L, size(keep))
    if ref_mc2 !== nothing
        push!(out, exact_stage("3.10 coarse dilation", "MC2 (coarse)",
                               UInt8.(grown), ref_mc2))
    end

    # The expansion back to the fine grid, and then the radius the fine pass actually receives. The
    # radius is the consequential comparison: `_expand_coarse_mask` and the reference's
    # `INTER_NEAREST` + edge-replication padding are two independent derivations of one correspondence,
    # and `src/multichip.jl` records that an offset between them puts the mask a row off the evidence
    # that produced it.
    # **The reference's expansion is offset from its own evidence, and AutoRIFT.jl's is not.** Its
    # `INTER_NEAREST` resize is left-aligned — coarse cell `k` covers fine `(k-1)*stride+1 .. k*stride`,
    # so the node at fine index `k*stride` is the *last* point of its own cell — while the radius it
    # reduces for that node is a centred `filtWidth`-wide window, fine `k*stride-4 .. k*stride+4` at the
    # default stride. So the reference gathers coherence evidence from fine 4..12 and applies it to fine
    # 1..8: offset by 3. `_expand_coarse_mask` inverts `_cell_max_radius!`'s own cell assignment instead,
    # so the mask lands on exactly the points the evidence was gathered over.
    #
    # This rung therefore *reports*. Reproducing the reference here would mean restricting the fine
    # search using evidence from three cells away, and `src/multichip.jl` records the measurement that
    # matching one of the reference's two halves alone is worse than matching neither. What is worth
    # watching is the count: it is the number of points the two disagree about searching, and rung 3.11c
    # gives the totals it resolves to.
    want = size(_consumed_grid(k, "xGrid0", L))
    stride = AutoRIFT._sparse_stride(p)
    expanded = AutoRIFT._expand_coarse_mask(grown, want, stride)
    full_mc2 = _state_shaped(k, "MC2", L, want)
    if full_mc2 !== nothing
        ref = full_mc2 .!= 0
        push!(out, StageResult("3.11a coarse mask on the fine grid", "MC2 (fine)", "reported", true,
                               length(ref),
                               @sprintf("julia searches %d, reference %d, differ %d of %d (%.4f%%); julia-only %d, reference-only %d",
                                        count(expanded), count(ref), count(expanded .!= ref),
                                        length(ref), 100count(expanded .!= ref) / length(ref),
                                        count(expanded .& .!ref), count(ref .& .!expanded))))
    end

    # `SearchLimitX0` after the zeroing: the array the fine correlator is handed. Identified as the
    # post-rewrite state whose zero count exceeds the pre-mask one — the mask only ever zeroes more.
    minsearch = Int(k.scalars["minSearch"])
    revs = [kk for kk in _revisions(k, "SearchLimitX0", L) if size(k.stages[kk]) == want]
    post = [kk for kk in revs if !any(v -> 0 < v < minsearch, k.stages[kk])]
    if length(post) >= 2
        before = k.stages[post[1]]
        after = k.stages[last(post)]
        jl = copy(before)
        jl[.!expanded] .= 0.0f0
        # Inherits 3.11a's offset by construction, so this reports the same disagreement in the units
        # that matter: how many points each side ends up searching.
        push!(out, StageResult("3.11b fine search radius", last(post), "reported", true,
                               length(jl),
                               @sprintf("differ %d of %d (%.4f%%); julia searches %d, reference %d",
                                        count(i -> jl[i] != after[i], eachindex(jl)), length(jl),
                                        100count(i -> jl[i] != after[i], eachindex(jl)) / length(jl),
                                        count(>(0), jl), count(>(0), after))))
        push!(out, StageResult("3.11c points the mask removes", last(post), "reported", true,
                               count(>(0), before),
                               @sprintf("reference searches %d of %d points after the mask (%d zeroed); julia %d",
                                        count(>(0), after), count(>(0), before),
                                        count(>(0), before) - count(>(0), after), count(>(0), jl))))
    end
    return out
end

# ---------------------------------------------------------------------------
# 3.12 — the coarse read-back: from a level's own grid to the full one
# ---------------------------------------------------------------------------
#
# Above the base chip size a level measures on a decimated grid and its answer has to be carried back to
# the full one. The reference does it in four steps (`autoRIFT.py:820-856`), and every one is traced:
#
#     DxF0 = colfilt(Dx, (Scale+1, Scale+1), 2)          # the finer levels' answer, NaN-aware mean
#     DxF0 = resize(DxF0, level_shape, INTER_AREA)        # reduced onto this level's grid
#     DxFM = colfilt(DxF, (5,5), 3)                       # this level's own holes, median-filled
#     DxFM[isnan(DxF) & !isnan(DxF0)] = DxF0[...]          # prefer the finer answer where it exists
#     DxF[isnan(DxF) & !isnan(DxFM)] = DxFM[...]           # then the median
#     DxF = resize(DxF, full_shape, INTER_CUBIC)           # and up to the full grid
#
# `_undecimate_level` and `_fill_level_holes` are AutoRIFT.jl's form of the same chain, and the ordering
# is load-bearing on both sides: the prior first, then the local median, then the nearest finite value.
# Filling nearest-first lets one distant value stand in for a whole gap, which at the edge of a fast
# feature means a value from across the discontinuity.
#
# **This is the stage that decides how much of a level's measurement survives**, so it is where a level
# that measures plenty and posts little has to be looked at. The comparison is per step rather than
# end-to-end: a bicubic resize of a slightly different field differs everywhere, so only the step that
# first diverges says anything.
function rungs_readback(k::Capture, p::Params, L::Int, chip::Int, chip0::Int)
    out = StageResult[]
    L == 0 && return out                      # the base level is assigned, not resized
    # This level's stride, which is the factor `_undecimate_level`'s `reduce_prior` uses — `step(rows) + 1`
    # for the mean and `step(rows)` for the area resize — so rung 3.12b models that step only if the two
    # are the same number.
    scale = _level_stride(p, L)

    # `DxFM` at a coarse level has three states: the `fillFiltWidth` median the fill loop reads, the
    # width-5 median built here, and that same array with `DxF0` written into it. The second is the one
    # this step produces, and it is identified by being the first state whose window is 5 rather than
    # `fillFiltWidth` — so it is taken by position among the states on this level's grid instead, since
    # both medians have the same shape.
    lvlshape = size(_consumed_grid(k, "xGrid0", L))
    dxf_states = [kk for kk in _revisions(k, "DxF", L) if size(k.stages[kk]) == lvlshape]
    dxfm_states = [kk for kk in _revisions(k, "DxFM", L) if size(k.stages[kk]) == lvlshape]
    (isempty(dxf_states) || length(dxfm_states) < 2) && return out

    # **The `DxF` this step reads is the one the fill loop left, not the last state on the grid.** `DxF`
    # is rebound once more at `:854` by the strong interpolation that *follows* this median, so the last
    # state is downstream of the very step being checked — medianing it reports 25.7% differing where the
    # right input reports none.
    #
    # The two are told apart by role rather than by index: the median's input is the state whose own
    # width-5 median is the array under test, which is a property this can check rather than a count it
    # has to know. On this level that selects `rev4` of five states, and it is exact.
    ref_dxfm = k.stages[dxfm_states[2]]
    match = findfirst(dxf_states) do kk
        m = windowmedian(k.stages[kk], 5)
        all(i -> (isnan(m[i]) && isnan(ref_dxfm[i])) || m[i] == ref_dxfm[i], eachindex(m))
    end
    dxf_key = match === nothing ? last(dxf_states) : dxf_states[match]
    push!(out, exact_stage("3.12a level median (width 5), from $dxf_key", dxfm_states[2],
                           windowmedian(k.stages[dxf_key], 5), ref_dxfm))

    # The finer levels' answer reduced onto this level's grid: a `Scale + 1` wide NaN-aware mean, then an
    # `INTER_AREA` resize. The mean is not redundant with the area average — its window is one wider than
    # a cell, so it reaches across each cell boundary and closes a partly-missing neighbourhood before the
    # area step weights what remains. `_undecimate_level`'s `reduce_prior` is this pair.
    dxf = k.stages[dxf_key]
    prev = _last_revision(k, "Dx", L - 1)
    if prev !== nothing && size(k.stages[prev]) != lvlshape
        jl = resample(windowmean(k.stages[prev], scale + 1), lvlshape, Area();
                      scale = (Float64(scale), Float64(scale)))
        # Reported: `DxF0` is a local and the trace does not carry it, so this states the reconstruction
        # against the input it was built from rather than against the reference's own array.
        push!(out, StageResult("3.12b prior reduced onto the level", "DxF0 (not traced)", "reported",
                               true, length(jl),
                               @sprintf("from %s: %d of %d finite; median %.4g",
                                        prev, count(!isnan, jl), length(jl),
                                        (v = filter(!isnan, vec(jl)); isempty(v) ? 0.0 : median(v)))))
    end
    return out
end

# ---------------------------------------------------------------------------
# 3.16 / 3.17 — the merge, and what this level contributed to it
# ---------------------------------------------------------------------------
#
# At the base chip size the merge is an assignment (`autoRIFT.py:810-816`): `Dx = DxF` wholesale, and
# `ChipSizeX[M0 | MF] = chip`. Above it the level's field is resized back up and written only where no
# finer level has claimed the point (`:818-866`):
#
#     idxRaw  = M0 & (ChipSizeX == 0)
#     idxFill = MF & (ChipSizeX == 0)
#     ChipSizeX[idxRaw | idxFill] = chip
#     Dx[idxRaw | idxFill] = DxF[idxRaw | idxFill]
#
# `ChipSizeX` is the rung that matters, and it is a *decision* rather than a measurement — which level
# owns each point — so it is gated on exact equality. The displacement is gated differently by level:
# quantized at the base chip size, where both sides land on the upsampling grid, and on bias above it,
# where both replace their measurement with a bicubic resize and neither field is quantized.
#
# The reference's `ChipSizeX` accumulates across levels, so the state after this level is its last
# revision at this level. Comparing against `rev0` would compare against what this level inherited.
function rungs_merge(k::Capture, L::Int, chip::Int, chip0::Int)
    out = StageResult[]
    cs_key = _last_revision(k, "ChipSizeX", L)
    cs_key === nothing && return out
    ref_cs = k.stages[cs_key]

    # Which points this level owns in the reference's answer: those it set to this chip size.
    ref_owned = ref_cs .== Float32(chip)
    push!(out, StageResult("3.16 level ownership", cs_key, "reported", true, count(ref_owned),
                           @sprintf("reference assigns %d points to chip %d; grid %s",
                                    count(ref_owned), chip, size(ref_cs))))

    # The displacement this level contributed, against the reference's own accumulated `Dx` restricted
    # to the points it owns. Both sides are on the full grid here, so no resampling is involved in the
    # comparison itself.
    dx_key = _last_revision(k, "Dx", L)
    dx_key === nothing && return out
    ref_dx = k.stages[dx_key]
    step = 1 / Float64(_level_upsampling(k, L))
    owned_ref = [ref_dx[i] for i in eachindex(ref_dx) if ref_owned[i]]
    onstep = count(v -> !isnan(v) && abs(v * (1 / step) - round(v * (1 / step))) < 1e-4, owned_ref)
    nval = count(!isnan, owned_ref)
    # Whether the level's values are quantized decides which statistic is meaningful downstream, so it
    # is measured here rather than assumed from the level index.
    push!(out, StageResult("3.17 quantization of this level", dx_key,
                           L == 0 ? "expected quantized" : "expected unquantized",
                           true, nval,
                           @sprintf("%d of %d owned values on the 1/%d grid (%.3f%%)",
                                    onstep, nval, round(Int, 1 / step),
                                    100onstep / max(nval, 1))))
    return out
end

# The upsampling factor this level used, from the reference's own per-level record.
function _level_upsampling(k::Capture, L::Int)
    corr = [r for r in k.levels if r.kind in ("coarse", "fine")]
    for r in corr
        r.chip_size[1] == Float64(Int(k.scalars["ChipSize0X"]) << L) && return r.oversample
    end
    return 16.0
end

# Every state of `name` at this level, in the order the trace wrote them, with the unversioned array
# first if there is one. The states of a name the pyramid rewrites are what a rung has to choose
# between, and choosing needs the list.
function _revisions(k::Capture, name::AbstractString, L::Int)
    out = String[]
    key = "$(name)_L$L"
    haskey(k.stages, key) && push!(out, key)
    append!(out, sort!([kk for kk in keys(k.stages)
                        if occursin(Regex("^$(name)_rev\\d+_L$L\$"), kk)],
                       by = kk -> parse(Int, match(r"rev(\d+)_", kk).captures[1])))
    return out
end

# The highest-numbered revision of `name` at this level, or `nothing`.
function _last_revision(k::Capture, name::AbstractString, L::Int)
    revs = sort!([kk for kk in keys(k.stages) if occursin(Regex("^$(name)_rev\\d+_L$L\$"), kk)],
                 by = kk -> parse(Int, match(r"rev(\d+)_", kk).captures[1]))
    isempty(revs) || return last(revs)
    key = "$(name)_L$L"
    return haskey(k.stages, key) ? key : nothing
end

# ---------------------------------------------------------------------------

function report(rs::Vector{StageResult})
    @printf("\n%-28s %-24s %-28s %-6s %10s\n", "stage", "reference", "gate", "state", "n")
    for r in rs
        @printf("%-28s %-24s %-28s %-6s %10d\n", r.name, r.reference, r.gate,
                r.passed ? "GREEN" : "RED", r.n)
        @printf("%s%s\n", " "^4, r.detail)
    end
    red = count(r -> !r.passed, rs)
    @printf("\n%d rung%s, %d green, %d red\n", length(rs), length(rs) == 1 ? "" : "s",
            length(rs) - red, red)
    return red == 0
end

"""
    dtype_pair(c; byte, float) -> Vector{StageResult}

The same level of the same case captured on both of the reference's data paths, compared against each
other.

`byte` is an ordinary capture (`DataType = 0`, quantized to 256 levels, `arImgDisp_u`) and `float` one
taken with `CAPTURE_FLOAT32=1` (`DataType = 1`, the filtered field kept, `arImgDisp_s`). Everything else
about the two runs is identical, so the difference between them **is** the quantization's contribution —
measured on production imagery rather than inferred from a windowed test.

This is what makes a byte-path residual attributable. A residual that is also present between the
reference's two own paths is the quantizer's; one that is not is the pipeline's. The comparison is on the
reference's own `DxC` arrays, so AutoRIFT.jl is not involved and the result is a property of the reference
alone — which is exactly what a floor has to be.
"""
function dtype_pair(c::GoldenCase; byte::Integer, float::Integer)
    kb = read_capture(c; n = byte)
    kf = read_capture(c; n = float)
    out = StageResult[]

    # The captures must differ in the data path and nothing else. **The element type of `in_I1` is the
    # authority, not the `DataType` scalar**: the scalar was added later, so a capture taken before it
    # simply lacks the key — and the correlator dispatches on the array's dtype rather than on the
    # attribute anyway (`autoRIFT.py:655,745`), so the array is what decides which template ran.
    tb = eltype(kb.arrays["in_I1"])
    tf = eltype(kf.arrays["in_I1"])
    push!(out, StageResult("dtype: the two captures are the two paths", "in_I1 eltype",
                           "UInt8 and Float32", tb === UInt8 && tf === Float32, 2,
                           "byte run $byte: $tb (DataType $(get(kb.scalars, "DataType", "absent"))); " *
                           "float run $float: $tf (DataType $(get(kf.scalars, "DataType", "absent")))"))
    (tb === UInt8 && tf === Float32) || return out

    L = traced_level(kb)
    L == traced_level(kf) || return push!(out, StageResult(
        "dtype: same level", "stage trace", "equal", false, 0,
        "byte capture traced level $L, float capture traced level $(traced_level(kf))"))

    # The grid the two runs correlated. It is built before `uniform_data_type` matters, so a difference
    # here means the two runs are not the same problem and nothing below them is comparable.
    gb = _consumed_grid(kb, "xGrid0C", L)
    gf = _consumed_grid(kf, "xGrid0C", L)
    push!(out, exact_stage("dtype: same coarse grid", "xGrid0C", gb, gf))
    gb == gf || return out

    for (label, name) in (("dx", "DxC"), ("dy", "DyC"))
        rb = _state_coarse(kb, name, L, size(gb))
        rf = _state_coarse(kf, name, L, size(gb))
        (rb === nothing || rf === nothing) && continue
        both = 0; ob = 0; of = 0; ex = 0
        d = Float64[]
        for i in eachindex(rb, rf)
            mb = isnan(rb[i]); mf = isnan(rf[i])
            mb && mf && continue
            if mf
                ob += 1
            elseif mb
                of += 1
            else
                both += 1
                δ = Float64(rb[i]) - Float64(rf[i])
                push!(d, δ)
                δ == 0 && (ex += 1)
            end
        end
        ad = abs.(d)
        # Reported, never gated: this is the reference disagreeing with itself, so there is no verdict to
        # reach — only a number every other rung on this level has to be read against.
        push!(out, StageResult("dtype: reference byte vs float, $label", name, "reported", true, both,
                               @sprintf("both %d, byte only %d, float only %d; exact %.2f%%, median %.4g, p99 %.4g, max %.4g, bias %+.4g",
                                        both, ob, of, both == 0 ? 0.0 : 100ex / both,
                                        isempty(ad) ? 0.0 : median(ad),
                                        isempty(ad) ? 0.0 : quantile(ad, 0.99),
                                        isempty(ad) ? 0.0 : maximum(ad),
                                        isempty(d) ? 0.0 : mean(d))))
    end
    return out
end

function main(args)
    isempty(args) && error("usage: stages.jl <product> [--run N] [--all] [--dtype-pair BYTE,FLOAT]")
    c = only(cases(args[1]))
    n = 100
    i = findfirst(==("--run"), args); i === nothing || (n = parse(Int, args[i + 1]))
    j = findfirst(==("--dtype-pair"), args)
    if j !== nothing
        b, f = parse.(Int, split(args[j + 1], ","))
        ok = report(dtype_pair(c; byte = b, float = f))
        ok || exit(1)
        return nothing
    end
    ok = report(ladder(c; n, stop_on_red = !("--all" in args)))
    ok || exit(1)
    return nothing
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS)
