# Rebuild the `ItsLiveInput` `capture_product.jl` captured and write it with `write_product`, once as
# captured and once with a `roi_valid_percentage` that selects the uncropped schema.
#
#   writeproduct <capturedir> <out.nc> <out_p000.nc>

using Dates: DateTime
using ItsLiveOffsetProduction
const P = ItsLiveOffsetProduction
import FastGeoProjections as FGP

include(joinpath(@__DIR__, "rawio.jl"))

say(s::String) = (write(Core.stderr, s); write(Core.stderr, "\n"); nothing)

readtext(d::String, k::String) = read(joinpath(d, k * ".txt"), String)

# `2020-07-03T15:00:24.531`, as `string(::DateTime)` writes it, parsed field by field: a `DateFormat`
# parse builds its error messages at run time, which a trimmed program cannot.
function parse_datetime(s::String)
    date, time = split(s, 'T')
    y, m, d = (parse(Int, p) for p in split(date, '-'))
    hms = split(time, '.')
    H, M, S = (parse(Int, p) for p in split(hms[1], ':'))
    ms = 0
    if length(hms) == 2
        f = hms[2]
        ms = length(f) >= 3 ? parse(Int, f[1:3]) : parse(Int, f) * 10^(3 - length(f))
    end
    return DateTime(y, m, d, H, M, S, ms)
end

function parse_attr(t::AbstractString, v::AbstractString)::P.AttrValue
    t == "String" && return String(v)
    t == "Float64" && return parse(Float64, v)
    t == "Float32" && return parse(Float32, v)
    t == "Int64" && return parse(Int64, v)
    t == "Int32" && return parse(Int32, v)
    throw(ArgumentError("no attribute rule for type " * String(t)))
end

function readattrs(path::String)
    d = Dict{String,P.AttrValue}()
    for line in eachline(path)
        isempty(line) && continue
        k, t, v = split(line, '\t'; limit = 3)
        d[String(k)] = parse_attr(t, v)
    end
    return d
end

satellite(t::String, v::String)::Union{String,Int64} = t == "Int64" ? parse(Int64, v) : v

function build(d::String, lonlat::F, roi_override::Union{Nothing,Float64}) where {F}
    x = vec(readraw(d, "x", Float64))
    y = vec(readraw(d, "y", Float64))
    px = meta_floats(d, "pixel_size")
    georef = ItsLiveGeoref(x, y, readattrs(joinpath(d, "mapping_attrs.tsv")), lonlat, px[1], px[2])
    dates = meta_strings(d, "acquisition_dates")
    missions = meta_strings(d, "missions")
    sats = meta_strings(d, "satellites")
    rll = meta_floats(d, "roi_lat_lon")
    info = ImagePairInfo(parse_datetime(dates[1]), parse_datetime(dates[2]), missions[1], missions[2],
                         satellite(sats[1], sats[2]), satellite(sats[3], sats[4]),
                         roi_override === nothing ? rll[1] : roi_override, rll[2], rll[3],
                         readattrs(joinpath(d, "extra.tsv")))
    r(k) = readraw(d, k, Float64)
    coeffs = GeogridCoefficients(r("c_offset2vx_1"), r("c_offset2vx_2"), r("c_offset2vy_1"),
                                 r("c_offset2vy_2"), r("c_scale_factor_1"), r("c_scale_factor_2"),
                                 nothing, nothing)
    refv = ReferenceVelocity(readraw(d, "ref_vx", Float32), readraw(d, "ref_vy", Float32))
    meta_string(d, "pair_type") == "optical" || throw(ArgumentError("only optical captures are rebuilt"))
    meta_string(d, "dt_seconds") == "nothing" || throw(ArgumentError("an optical capture has no dt_seconds"))
    return ItsLiveInput(:optical, readtext(d, "detection_method"), readtext(d, "coordinates"),
                        readtext(d, "autorift_software_version"), readtext(d, "parameter_file"),
                        readtext(d, "source"), nothing, georef, info, coeffs, refv,
                        readraw(d, "stable_mask", Bool), nothing,
                        readraw(d, "dx", Float32), readraw(d, "dy", Float32),
                        readraw(d, "chip_size_x", UInt16), meta_float(d, "scale_chip_size_y"),
                        readraw(d, "interp_mask", Bool))
end

function write_both(d::String, out::String, out_p000::String, lonlat::F) where {F}
    write_product(out, build(d, lonlat, nothing))
    write_product(out_p000, build(d, lonlat, 0.5))
    return nothing
end

# The grid CRS's own transform to longitude and latitude, one concrete type per branch.
function run(d::String, out::String, out_p000::String)
    g = FGP.native_crs(meta_int(d, "grid_epsg"))
    if g isa FGP.UTMCRS
        write_both(d, out, out_p000, FGP.Transformation(g, FGP.GeographicCRS()))
    elseif g isa FGP.PolarStereographicCRS
        write_both(d, out, out_p000, FGP.Transformation(g, FGP.GeographicCRS()))
    else
        throw(ArgumentError("a product grid is UTM or polar stereographic"))
    end
    return nothing
end

function (@main)(args::Vector{String})::Cint
    length(args) == 3 || (say("usage: writeproduct <capturedir> <out.nc> <out_p000.nc>"); return Cint(2))
    try
        run(args[1], args[2], args[3])
    catch err
        m = err isa ArgumentError ? err.msg : err isa ErrorException ? err.msg : nothing
        say("error: " * (m isa String ? m : string(nameof(typeof(err)))))
        return Cint(1)
    end
    say("wrote " * args[2] * " and " * args[3])
    return Cint(0)
end
