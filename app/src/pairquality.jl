# Whether a pair's two images can be correlated at all, read from a sample of each before correlation, and
# the diagnostics given when no point of the pair resolves.
#
# Only conditions under which correlation cannot succeed are warned about: a sampled image with no data, or
# with the same value everywhere, which leaves a chip nothing to match. Brightness is reported but never
# tested: good pairs differ by tens of dB (Landsat sun angle and gain, wet snow at C band), and the
# correlation of the two images' block means is below 0.1 on some good pairs, so neither separates good
# pairs from bad. Whether a pair correlated is decided by its result instead.

# Sample windows: a 3 × 3 lattice of squares this many pixels on a side.
const QUALITY_WINDOW = 512

"""
    ImageSample

One image over the sample windows: `median` of its samples with data (NaN when there are none),
`valid` the fraction of samples with data, and `flat` whether every sample with data has one value.
"""
struct ImageSample
    median::Float64
    valid::Float64
    flat::Bool
end

"""
    PairQuality

The two images over the same sample windows, and `db` the late image's brightness relative to the early
one: `20 log10` of the ratio of medians, the same in amplitude and in digital number. NaN unless both
medians are positive, as a Wallis-filtered scene's need not be.
"""
struct PairQuality
    early::ImageSample
    late::ImageSample
    db::Float64
end

# Zero is no data in every scene the chain reads.
_has_data(v) = isfinite(v) && v != 0

function _sample(img, corners::Vector{NTuple{2,Int}}, w::Int)
    vals = Float64[]
    n = 0
    for (r0, c0) in corners
        a = img[r0:r0 + w - 1, c0:c0 + w - 1]
        n += length(a)
        for v in a
            _has_data(v) && push!(vals, v)
        end
    end
    isempty(vals) && return ImageSample(NaN, 0.0, false)
    return ImageSample(_median(vals), length(vals) / n, minimum(vals) == maximum(vals))
end

"""
    pair_quality(img::Imagery) -> PairQuality

[`PairQuality`](@ref) of the two images over a 3 × 3 lattice of windows across their common extent. The
images are the run's local tile caches, so this reads a few megabytes from disk.
"""
function pair_quality(img::Imagery)
    nr, nc = size(img.early)
    size(img.late) == (nr, nc) || throw(DimensionMismatch("the two images are not on one grid"))
    w = min(QUALITY_WINDOW, nr, nc)
    corners = [(clamp(round(Int, fr * nr / 4 - w / 2), 1, nr - w + 1), clamp(round(Int, fc * nc / 4 - w / 2), 1, nc - w + 1))
               for fr in (1, 2, 3) for fc in (1, 2, 3)]
    e, l = _sample(img.early, corners, w), _sample(img.late, corners, w)
    return PairQuality(e, l, _db(e.median, l.median))
end

# The brightness difference, which only an amplitude has; see `PairQuality`.
_db(early::Float64, late::Float64) = early > 0 && late > 0 ? 20 * log10(late / early) : NaN

# The median of a non-empty vector holding no NaN.
_median(v::Vector{Float64}) = (s = sort(v); m = length(s); isodd(m) ? s[(m + 1) ÷ 2] : (s[m ÷ 2] + s[m ÷ 2 + 1]) / 2)

_g3(x::Float64) = isnan(x) ? "n/a" : string(round(x; sigdigits = 3))
_f1(x::Float64) = isnan(x) ? "n/a" : string(round(x; digits = 1))
_percent(x::Float64) = string(round(100x; digits = 1)) * "%"

_describe(s::ImageSample) = "median " * _g3(s.median) * ", data in " * _percent(s.valid) * " of samples"

"""
    quality_report(early, late, q) -> String

The pair's diagnostics on one line: each scene with its median and how much of its sample has data, and
the brightness difference in dB.
"""
quality_report(early::String, late::String, q::PairQuality) =
    "early " * early * " " * _describe(q.early) * "; late " * late * " " * _describe(q.late) *
    "; brightness difference " * _f1(q.db) * " dB"

"""
    quality_warnings(q) -> Vector{String}

What about `q` stops the pair from correlating: an image with no data in any sample window, or with one
value throughout.
"""
function quality_warnings(q::PairQuality)
    w = String[]
    for (name, s) in (("early", q.early), ("late", q.late))
        s.valid == 0 && push!(w, "the " * name * " image has no data in any sample window")
        s.flat && push!(w, "the " * name * " image has one value throughout the sample windows")
    end
    return w
end

"""
    warn_quality(early, late, q)

Report on stderr, in one line with the pair's diagnostics, whatever [`quality_warnings`](@ref) finds.
"""
function warn_quality(early::String, late::String, q::PairQuality)
    w = quality_warnings(q)
    isempty(w) || say("warning: " * join(w, "; ") * ". " * quality_report(early, late, q))
    return nothing
end

"""
    check_resolved(c::Correlation, early, late, q)

Throw, naming the pair and its diagnostics, when no point of the correlation resolved: a product with no
velocity at all is never a result, whatever the reason.
"""
function check_resolved(c::Correlation, early::String, late::String, q::PairQuality)
    any(!iszero, c.chip_size) && return nothing
    w = quality_warnings(q)
    throw(ErrorException("no point of the pair resolved, so no product is written" *
                         (isempty(w) ? "" : "; " * join(w, "; ")) * ". " * quality_report(early, late, q)))
end
