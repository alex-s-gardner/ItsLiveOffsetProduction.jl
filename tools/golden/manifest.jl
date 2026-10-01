# The golden test cases: the 22 ITS_LIVE products in
# `s3://its-live-data/test-space/golden/`, and the job that produced each one.
#
# `manifest.json` is generated from ASFHyP3/hyp3-testing's `autorift_golden.json.j2` and committed,
# rather than parsed from the template at run time: the template is Jinja, so reading it needs
# comment and expression stripping, and the job-to-product mapping is not a naming rule that can be
# applied blindly. Two things make it a lookup rather than a derivation:
#
#   - A product is named `<img1>_X_<img2>_G0120V02_P<nn>` in **acquisition order**, which is not the
#     job's reference/secondary order. Jobs 7 and 15 of the template are reversed relative to their
#     product, and `img_pair_info:id_img1` follows the product, not the job. So `reference` must not
#     be assumed to be `img1`.
#   - A scene may appear in more than one pair. `LC08_L1TP_060018_20130330_20200912_02_T1` is in
#     two, so matching on one scene name or on a date is ambiguous; both halves are needed.
#
# Burst jobs are the exception: a product never names its bursts, so the first reference burst's
# acquisition timestamp identifies `img1`.

using JSON3

const MANIFEST_PATH = joinpath(@__DIR__, "manifest.json")

"""
    GoldenCase

One golden test case: a product in the golden bucket and the HyP3 job that made it.

`reference` and `secondary` are the job's inputs, in job order — see this file's header for why
that is not necessarily the order the product names them in. `phase` is the phase of
`tools/golden/README.md` that this case belongs to, which tracks what AutoRIFT.jl can currently
process rather than anything about the case itself.
"""
struct GoldenCase
    product::String
    platform::String
    phase::Int
    note::String
    reference::Vector{String}
    secondary::Vector{String}
    frame_id::Union{String,Nothing}
end

"""
    cases() -> Vector{GoldenCase}

Every golden test case, ordered by phase then product name.
"""
function cases()
    raw = JSON3.read(read(MANIFEST_PATH, String))
    return [GoldenCase(c.product, c.platform, c.phase, c.note,
                       collect(String, c.reference), collect(String, c.secondary),
                       get(c, :frame_id, nothing)) for c in raw]
end

"""
    cases(phase::Integer) -> Vector{GoldenCase}
    cases(name::AbstractString) -> Vector{GoldenCase}

The cases in one phase, or the single case whose product name contains `name` — enough of the name
to be unambiguous, since a full product name is 100 characters.
"""
cases(phase::Integer) = filter(c -> c.phase == phase, cases())

function cases(name::AbstractString)
    hits = filter(c -> occursin(name, c.product), cases())
    isempty(hits) && throw(ArgumentError("no golden case matches \"$name\""))
    length(hits) == 1 || throw(ArgumentError(
        "\"$name\" matches $(length(hits)) cases; be more specific:\n" *
        join(("  " * c.product for c in hits), "\n")))
    return hits
end

# Where downloaded data lives. Outside the repository: the products alone are 145 MB and the scenes
# they were made from are gigabytes. Mirrors `AUTORIFT_TESTDATA` in `test/utils.jl`.
const CACHE = get(ENV, "AUTORIFT_GOLDEN_CACHE",
                  joinpath(homedir(), "data", "autorift", "tests", "golden_tests"))

const BUCKET = "s3://its-live-data/test-space/golden"

# The sidecars beside each product. `.head` is ncdump output, so a variable or attribute can be
# checked without opening the product; `.stac.json` carries the acquisition times and footprint,
# which is what `img_pair_info` is built from — so the metadata path is testable without querying
# STAC again.
#
# Two naming conventions, which is why these are separate lists rather than one: `.head`, `.premet`
# and `.spatial` hang off the *file* (`<product>.nc.head`), while `.stac.json` hangs off the
# *product* (`<product>.stac.json`).
const FILE_SIDECARS = (".head", ".premet", ".spatial")
const PRODUCT_SIDECARS = (".stac.json",)

golden_dir() = joinpath(CACHE, "products")
golden_path(c::GoldenCase) = joinpath(golden_dir(), c.product * ".nc")

"""
    sidecar_path(c::GoldenCase, ext) -> String

Where sidecar `ext` of `c` lives in the cache. `ext` must be one of [`FILE_SIDECARS`](@ref) or
[`PRODUCT_SIDECARS`](@ref); the two differ in whether the extension follows `.nc`.
"""
function sidecar_path(c::GoldenCase, ext::AbstractString)
    ext in FILE_SIDECARS && return joinpath(golden_dir(), c.product * ".nc" * ext)
    ext in PRODUCT_SIDECARS && return joinpath(golden_dir(), c.product * ext)
    throw(ArgumentError("unknown sidecar \"$ext\"; expected one of " *
                        join((FILE_SIDECARS..., PRODUCT_SIDECARS...), ", ")))
end

"""
    have_golden(c::GoldenCase) -> Bool

Whether `c`'s product is in the cache. The comparison is skipped rather than failed when it is not,
so a fresh checkout runs whatever it can.
"""
have_golden(c::GoldenCase) = isfile(golden_path(c))
