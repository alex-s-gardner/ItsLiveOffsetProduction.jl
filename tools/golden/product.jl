# Reading an ITS_LIVE product into something two products can be compared through.
#
# A product is not one schema. Three variants appear in the golden set, and a reader that assumes
# any one of them silently mis-reads the others:
#
#   optical  12 variables, `(x, y, time)`      — Landsat, Sentinel-2, and NISAR L2 GSLC
#   radar    16 variables, `(x, y, time)`      — adds `vr`, `va`, `M11`, `M12`
#   P000      11 variables, `(x, y)`, no `time` — `process.py` skips cropping when no pixel is
#                                                 valid, and cropping is what adds the time axis
#
# So `time` and the radar quartet are optional, and the code asks rather than assumes. The `P000`
# case is not a degenerate input to skip: it is the one product in the set that exercises the
# uncropped path, and its absence of a time axis is the observable difference.
#
# **Dimension order.** NCDatasets reports a netCDF `(time, y, x)` variable as `(x, y, time)`, because
# netCDF is row-major and Julia column-major and the library preserves memory layout rather than
# index order. So `ds["vx"]` is indexed `[x, y, time]` — column first — while AutoRIFT indexes
# `[row, col]`. This reader drops the singleton time axis and *keeps* the file's `[x, y]` order,
# stating it rather than transposing, so that a comparison between two products is a comparison of
# identically-shaped arrays and no transpose can hide in it. Anything handing these to the
# correlator must transpose, as `tools/realdata/prepare.jl` does.

using NCDatasets, Dates, Statistics

# `_FillValue` is set on every 2-D variable, so NCDatasets returns `Union{Missing,T}`. `missing` is
# the product's nodata and has to survive into the comparison — it is a value both sides must agree
# about, not an absence to be filled.
const Plane{T} = AbstractMatrix{Union{Missing,T}}

"""
    Product

One ITS_LIVE product: its data planes, its coordinates, its per-variable attributes and its global
attributes.

`planes` is keyed by variable name and holds each 2-D array in the file's own `[x, y]` order, with
the singleton time axis dropped. `attribs` holds every variable's attributes, including the
`stable_shift`/`error` family that carries as much of the answer as the pixels do — a product whose
`vx` matched but whose `stable_shift` did not would be wrong in a way only the attributes reveal.
"""
struct Product
    name::String
    planes::Dict{String,Any}
    x::Vector{Float64}
    y::Vector{Float64}
    time::Union{DateTime,Nothing}
    attribs::Dict{String,Dict{String,Any}}
    global_attribs::Dict{String,Any}
end

# Variables that are coordinates or scalar metadata rather than data planes.
const NON_PLANE = ("x", "y", "time", "mapping", "img_pair_info")

"""
    read_product(path) -> Product

Read the product at `path`.

Every 2-D variable becomes a plane; `mapping` and `img_pair_info` are read for their attributes
only, since both are single fill bytes whose entire content is metadata.
"""
function read_product(path::AbstractString)
    NCDataset(path) do ds
        planes = Dict{String,Any}()
        attribs = Dict{String,Dict{String,Any}}()

        for (name, var) in ds
            attribs[name] = Dict{String,Any}(String(k) => v for (k, v) in var.attrib)
            name in NON_PLANE && continue
            a = Array(var)
            # Drop the singleton time axis that cropping adds. `P000` products never had one.
            planes[name] = ndims(a) == 3 ? dropdims(a; dims = 3) : a
        end

        t = haskey(ds, "time") ? first(skipmissing(Array(ds["time"]))) : nothing
        g = Dict{String,Any}(String(k) => v for (k, v) in ds.attrib)

        return Product(basename(path), planes,
                       collect(Float64, skipmissing(Array(ds["x"]))),
                       collect(Float64, skipmissing(Array(ds["y"]))),
                       t, attribs, g)
    end
end

read_product(c::GoldenCase) = read_product(golden_path(c))

"""
    valid_fraction(p::Product, var = "vx") -> Float64

Fraction of `var` that is not nodata. The product name's `P<nn>` records this as a percentage, so it
is checkable against the name — which is how a reader that dropped the fill value is caught.
"""
function valid_fraction(p::Product, var::AbstractString = "vx")
    a = p.planes[var]
    return count(!ismissing, a) / length(a)
end

"""
    schema(p::Product) -> Symbol

`:radar` if `p` carries the range/azimuth quartet, `:optical` otherwise.

Read from the variables present rather than from the product name: NISAR L2 GSLC is a radar sensor
whose product has the optical schema, because the geocoded product measures displacement on a map
grid rather than in range and azimuth.
"""
schema(p::Product) = all(haskey(p.planes, v) for v in ("vr", "va", "M11", "M12")) ? :radar : :optical

"""
    cropped(p::Product) -> Bool

Whether `p` went through `crop.py`. An uncropped product has no time coordinate, since the time axis
is added during cropping; `process.py` skips cropping when no pixel is valid.
"""
cropped(p::Product) = p.time !== nothing
