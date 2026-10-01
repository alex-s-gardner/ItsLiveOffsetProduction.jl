# Comparing two ITS_LIVE products, per variable.
#
# The comparison is deliberately strict by default. Everything in a product is deterministic except
# three fields, so a tolerance is a claim that needs evidence — and the evidence is what
# `reference.jl` measures by running the reference twice. Until a number has been measured, the
# honest default is exact equality, which is why `EXACT` is the fallback rather than some small
# epsilon that would quietly absorb a real difference.
#
# The three exceptions, all properties of the reference rather than of either implementation:
#
#   `v_error`   where `v == 0`, `netcdf_output.py::v_error_cal` draws 10^6 samples from an
#               *unseeded* `default_rng()`, so this moves between two runs of the reference itself.
#   `time`      `crop.py::numeric_hash` jitters the coordinate by `hash(filename) % 10^6`
#               microseconds, and `PYTHONHASHSEED` is unset, so Python salts it per process.
#   global      `date_created` is wall clock; the version strings record the build.
#
# Nodata is compared as a value, not skipped. Where one product measured a pixel and the other did
# not, that is a disagreement about coverage — often the most informative kind, since it points at
# the degenerate-chip and stable-shift divergences rather than at arithmetic.

using Printf, Statistics

"""
    VarDiff

How one variable of two products compares.

`n_both` counts pixels both products measured, and `only_a`/`only_b` the pixels one measured alone —
kept separate from the value statistics because a coverage difference and a value difference have
different causes and different fixes. `max_abs` and `p99` are over `n_both` only, since a difference
against nodata is not a number.
"""
struct VarDiff
    name::String
    n_both::Int
    only_a::Int
    only_b::Int
    n_exact::Int
    max_abs::Float64
    p99::Float64
    bias::Float64
end

agrees(d::VarDiff; tol = 0.0) =
    d.only_a == 0 && d.only_b == 0 && (d.n_exact == d.n_both || d.max_abs <= tol)

exact_fraction(d::VarDiff) = d.n_both == 0 ? 1.0 : d.n_exact / d.n_both

"""
    compare_plane(name, a, b) -> VarDiff

Compare one plane of two products. `a` and `b` must have identical axes; a shape difference is an
error rather than a finding, since nothing below it would be meaningful.
"""
function compare_plane(name::AbstractString, a::AbstractMatrix, b::AbstractMatrix)
    axes(a) == axes(b) || throw(DimensionMismatch(
        "$name: $(size(a)) vs $(size(b)) — products are on different grids, so no per-pixel " *
        "comparison is possible; compare their `x`/`y` coordinates first"))

    both = 0; onlya = 0; onlyb = 0; exact = 0
    diffs = Float64[]
    for i in eachindex(a, b)
        ma, mb = ismissing(a[i]), ismissing(b[i])
        if ma && mb
            continue
        elseif mb
            onlya += 1
        elseif ma
            onlyb += 1
        else
            both += 1
            va, vb = Float64(a[i]), Float64(b[i])
            va == vb ? (exact += 1) : push!(diffs, abs(va - vb))
        end
    end

    # `bias` is over both-measured pixels including the equal ones, so it is the mean signed error
    # rather than the mean over disagreeing pixels — a bias computed only where the two differ would
    # overstate it by the fraction that agree.
    signed = Float64[]
    for i in eachindex(a, b)
        (ismissing(a[i]) || ismissing(b[i])) && continue
        push!(signed, Float64(a[i]) - Float64(b[i]))
    end

    return VarDiff(String(name), both, onlya, onlyb, exact,
                   isempty(diffs) ? 0.0 : maximum(diffs),
                   isempty(diffs) ? 0.0 : quantile(diffs, 0.99),
                   isempty(signed) ? 0.0 : mean(signed))
end

"""
    ProductDiff

The whole comparison: per-variable differences, coordinate agreement, and the attributes that
differ.

`attrib_diffs` maps `"variable.attribute"` to an `(a, b)` pair. Attributes carry as much of the
answer as the pixels — `stable_shift`, the four `error` estimates, `stable_count` — so a product
whose planes matched and whose attributes did not has not matched.
"""
struct ProductDiff
    a::String
    b::String
    vars::Vector{VarDiff}
    x_max::Float64
    y_max::Float64
    time_delta::Union{Float64,Nothing}
    attrib_diffs::Dict{String,Tuple{Any,Any}}
    missing_vars::Vector{String}
end

# Attributes that cannot match, per this file's header. Excluded from `attrib_diffs` rather than
# reported and ignored, so that what remains is a list every entry of which is a real finding.
const VOLATILE_GLOBAL = ("date_created",)

"""
    compare_products(a::Product, b::Product; skip = ()) -> ProductDiff

Compare `a` against `b` variable by variable, plus coordinates and attributes.

Variables in `skip` are not compared. A variable present in one product and not the other is
recorded in `missing_vars` rather than silently dropped: the optical and radar schemas differ by four
variables, and comparing across schemas is a mistake worth surfacing.
"""
function compare_products(a::Product, b::Product; skip = ())
    shared = sort!(collect(String, intersect(keys(a.planes), keys(b.planes))))
    missed = sort!(collect(String, symdiff(keys(a.planes), keys(b.planes))))

    vars = VarDiff[]
    for name in shared
        name in skip && continue
        push!(vars, compare_plane(name, a.planes[name], b.planes[name]))
    end

    xmax = length(a.x) == length(b.x) ? maximum(abs, a.x .- b.x) : Inf
    ymax = length(a.y) == length(b.y) ? maximum(abs, a.y .- b.y) : Inf

    td = if a.time === nothing || b.time === nothing
        a.time === b.time ? nothing : Inf
    else
        abs((a.time - b.time).value) / 1000  # seconds
    end

    ad = Dict{String,Tuple{Any,Any}}()
    for var in sort!(collect(String, intersect(keys(a.attribs), keys(b.attribs))))
        aa, bb = a.attribs[var], b.attribs[var]
        for k in sort!(collect(String, intersect(keys(aa), keys(bb))))
            isequal(aa[k], bb[k]) || (ad["$var.$k"] = (aa[k], bb[k]))
        end
    end
    for k in sort!(collect(String, intersect(keys(a.global_attribs), keys(b.global_attribs))))
        k in VOLATILE_GLOBAL && continue
        isequal(a.global_attribs[k], b.global_attribs[k]) ||
            (ad["global.$k"] = (a.global_attribs[k], b.global_attribs[k]))
    end

    return ProductDiff(a.name, b.name, vars, xmax, ymax, td, ad, collect(String, missed))
end

"""
    Base.show(io, d::ProductDiff)

A per-variable table. `exact` is the fraction of both-measured pixels that agree bit for bit, and
`only a`/`only b` the coverage difference — the column to read first, since a coverage difference
makes the value columns a comparison of different pixel sets.
"""
function Base.show(io::IO, d::ProductDiff)
    println(io, "a: ", d.a)
    println(io, "b: ", d.b)
    isempty(d.missing_vars) || println(io, "\nvariables in one product only: ",
                                       join(d.missing_vars, ", "))
    @printf(io, "\n%-18s %10s %8s %8s %8s %10s %10s %10s\n",
            "variable", "both", "only a", "only b", "exact", "max", "p99", "bias")
    for v in d.vars
        @printf(io, "%-18s %10d %8d %8d %7.2f%% %10.4g %10.4g %10.4g\n",
                v.name, v.n_both, v.only_a, v.only_b, 100 * exact_fraction(v),
                v.max_abs, v.p99, v.bias)
    end
    @printf(io, "\ncoordinates: max |dx| = %.6g, max |dy| = %.6g\n", d.x_max, d.y_max)
    d.time_delta === nothing || @printf(io, "time: %.6g s apart\n", d.time_delta)
    if isempty(d.attrib_diffs)
        println(io, "attributes: all shared attributes equal")
    else
        println(io, "\nattributes differing (", length(d.attrib_diffs), "):")
        for k in sort(collect(keys(d.attrib_diffs)))
            va, vb = d.attrib_diffs[k]
            println(io, "  ", rpad(k, 34), " ", repr(va), "  vs  ", repr(vb))
        end
    end
    return nothing
end

"""
    identical(d::ProductDiff) -> Bool

Whether the two products agree completely: every shared variable bit-exact with no coverage
difference, identical coordinates, identical time, no differing attribute, and no variable present in
only one.

This is what a product compared against *itself* must satisfy. It deliberately does not accept a
tolerance — a tolerance belongs to a specific measured comparison, not to the definition of
agreement. Two separate runs of the reference do not satisfy it, because the time coordinate is
jittered per process; [`agrees_on_data`](@ref) is the question to ask of those.
"""
identical(d::ProductDiff) =
    isempty(d.missing_vars) && all(agrees, d.vars) &&
    d.x_max == 0 && d.y_max == 0 && isempty(d.attrib_diffs) &&
    (d.time_delta === nothing || d.time_delta == 0)

# The time coordinate and the attribute recording its own jitter. Both are `crop.py::numeric_hash`
# of the filename under a per-process hash salt, so they differ between any two runs — including two
# runs of the reference on one granule.
const VOLATILE_TIME = ("time.microseconds_added",)

"""
    agrees_on_data(d::ProductDiff) -> Bool

Whether the two products agree on everything a correlator controls: every plane bit-exact with no
coverage difference, identical coordinates, and every attribute equal except the time jitter.

This is the gate for comparing two *different* runs or two implementations. It differs from
[`identical`](@ref) only in ignoring the time coordinate, which no implementation can reproduce.
"""
agrees_on_data(d::ProductDiff) =
    isempty(d.missing_vars) && all(agrees, d.vars) &&
    d.x_max == 0 && d.y_max == 0 &&
    all(k -> k in VOLATILE_TIME, keys(d.attrib_diffs))
