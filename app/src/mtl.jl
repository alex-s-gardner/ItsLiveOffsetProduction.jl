# A Landsat scene's scan geometry from its `_ANG.txt` ephemeris, which orients the Landsat 4/5 band-reject
# filter along the scan.

"""
    ang_fields(text) -> Dict{String,Vector{Float64}}

Every numeric field of a Landsat `_ANG.txt`, scalars as one-element vectors.

The fields that matter are parenthesised lists spanning several lines, `EPHEMERIS_ECEF_X = (-2513196.897534,
-2519426.030063, ...`, which a line-at-a-time `KEY = VALUE` reader truncates to its first row.
"""
function ang_fields(text::String)
    out = Dict{String,Vector{Float64}}()
    key = ""
    buf = Float64[]
    open_list = false
    for raw in split(text, '\n')
        s = strip(raw)
        if open_list
            # A continuation row of the list opened above; the closing parenthesis ends it.
            done = occursin(')', s)
            append!(buf, _numbers(s))
            if done
                out[key] = copy(buf)
                open_list = false
            end
            continue
        end
        i = findfirst('=', s)
        i === nothing && continue
        k = strip(s[1:prevind(s, i)])
        v = strip(s[nextind(s, i):end])
        (startswith(k, "GROUP") || startswith(k, "END_GROUP")) && continue
        if startswith(v, "(")
            key = String(k)
            buf = _numbers(v)
            occursin(')', v) ? (out[key] = copy(buf)) : (open_list = true)
        else
            n = tryparse(Float64, v)
            n === nothing || (out[String(k)] = [n])
        end
    end
    isempty(out) && throw(ArgumentError("an _ANG.txt parsed to no numeric fields"))
    return out
end

_numbers(s::AbstractString) =
    Float64[parse(Float64, m.match) for m in eachmatch(r"-?\d+\.?\d*(?:[eE][-+]?\d+)?", s)]

"""
    ang_ephemeris(text) -> NamedTuple

The satellite's trajectory from `_ANG.txt`: `time` in seconds of day and `x`, `y`, `z` in ECEF metres.
`NUMBER_OF_POINTS` is checked against all four vectors, because a truncated list would otherwise give a
shorter trajectory and a bearing taken from the wrong pair of points.
"""
function ang_ephemeris(text::String)
    f = ang_fields(text)
    get1(k::String) = haskey(f, k) ? f[k] : throw(ArgumentError(
        "the _ANG.txt has no " * k * "; the ephemeris is what the scan geometry is derived from"))
    t, x, y, z = get1("EPHEMERIS_TIME"), get1("EPHEMERIS_ECEF_X"), get1("EPHEMERIS_ECEF_Y"),
                 get1("EPHEMERIS_ECEF_Z")
    n = haskey(f, "NUMBER_OF_POINTS") ? Int(first(f["NUMBER_OF_POINTS"])) : length(t)
    (length(t) == n && length(x) == n && length(y) == n && length(z) == n) || throw(ArgumentError(
        "the _ANG.txt declares " * string(n) * " ephemeris points but its lists differ in length"))
    return (; time = t, x, y, z)
end

"""
    orbit_scan_angles(ang_text, epsg; spacing) -> (along_track, cross_track)

The scan geometry from the orbit, in the reference's own angle convention: two consecutive ephemeris
positions at the middle of the ephemeris, transformed into the raster's CRS, give the ground-track
direction. The projection step matters — a slope on a UTM raster is a grid bearing, which departs from true
north by the meridian convergence. The two angles are perpendicular by construction. `spacing` undoes pixel
anisotropy: an angle on a grid sampled differently in x and y is not the angle in metres.
"""
function orbit_scan_angles(ang_text::String, epsg::Int; spacing::Tuple{Float64,Float64})
    e = ang_ephemeris(ang_text)
    k = length(e.time) ÷ 2
    # Geocentric ECEF is EPSG:4978; traditional axis order keeps the projected result (x, y).
    p1, p2 = with_transform(4978, epsg) do ct
        (transform_point(ct, e.x[k], e.y[k], e.z[k]), transform_point(ct, e.x[k + 1], e.y[k + 1], e.z[k + 1]))
    end
    # Magnitudes, so the frame stays y-up: dividing by a north-up raster's signed `-30` would flip the
    # bearing into a row-down frame and negate both results.
    sx, sy = abs(spacing[1]), abs(spacing[2])
    dx = (p2[1] - p1[1]) / sx
    dy = (p2[2] - p1[2]) / sy
    # `atan` of a ratio, in degrees, so the result is in (-90, 90) and carries no quadrant. The ground track
    # is the along-track direction and the cross-track its perpendicular, folded back into the same interval.
    along = rad2deg(atan(dy / dx))
    cross = along + (along > 0 ? -90.0 : 90.0)
    return (along, cross)
end
