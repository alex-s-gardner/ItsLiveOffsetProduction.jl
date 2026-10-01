# The scan geometry of a Landsat scene, from its `_MTL.txt`.
#
# **Temporary, and deliberately harness-only.** The durable home for this is `ImagePairGeometry`, which
# already owns acquisition geometry — `ImageFootprint`, `ProjectedCoordinate`, `Orbit`, `incidence_angle` —
# and has the Rasters extension to attach it to a raster at load time. Parallax will want the same fields,
# so a parser here would be duplicated there. This exists so the L4/5 destripe gate is not blocked on that
# package, and it is deleted when the real loader lands: nothing outside `tools/golden` calls it, and the
# gate it feeds compares angles against the reference's own logged values rather than against this parse.
#
# **The corner fields do not carry the scan geometry, and that is measured.** `_fft_filter` recovers the
# footprint from pixels — `connectedComponentsWithStats`, `findContours`, `moments`, `minAreaRect`, a
# `warpAffine`d quadrant map, four `distanceTransform` argmaxes (`autoRIFT.py:153-206`) — and reading the
# MTL corners instead does *not* reproduce it: an L1 product is `ORIENTATION = NORTH_UP`, so both its
# projected and its lat/lon corner sets are the axis-aligned bounding box and every slope from them is 0
# or ±90, against the 71.60°/−20.25° the reference logs.
#
# What does carry it is the **orbit**: the `_ANG.txt` beside the MTL holds `EPHEMERIS_ECEF_{X,Y,Z}` at 1 s
# spacing, and two consecutive positions transformed into the raster's own CRS give the ground-track
# direction. The projection step is not optional — a slope on a UTM raster is a grid bearing, and grid north
# departs from true north by the meridian convergence, 3.2° at `LT05_L1TP_060018` and more near the poles.
#
# **Measured over all six scenes of the three L4/5 pairs** (`orbit_angle_check`, gate `5.orbit`): the
# orbit's cross-track agrees with the reference's to **0.10° at worst**, while its along-track sits
# **1.45° to 1.89° below** it. The error is the reference's and it is all in one axis — the orbit's two
# directions are perpendicular to the last digit by construction, and the reference's are 91.42° to 91.86°
# apart, never 90°. Its own non-perpendicularity accounts for the along-track gap case by case, which is
# what identifies `nanmax` over two edge slopes as the cause: the worse-conditioned edge wins and only the
# along-track pair is affected.
#
# `tools/golden/README.md` registers the orbit route as the more correct alternative, to be adopted once the
# L4/5 pairs agree. `mtl_corners` and `scan_angles` below are kept for that comparison, not because the
# corners are the answer.

using Printf

"""
    mtl_fields(path) -> Dict{String,String}

Every `KEY = VALUE` pair in a Landsat `_MTL.txt`, as strings.

The format is flat `KEY = VALUE` inside `GROUP`/`END_GROUP` blocks, and keys are unique across groups in
every Landsat collection, so the groups are ignored rather than parsed. Values keep their quotes stripped
and are left as text: a caller that wants a number says so, which keeps a malformed field an error at the
point of use rather than a silent zero here.
"""
function mtl_fields(path::AbstractString)
    isfile(path) || error("no MTL at $path")
    out = Dict{String,String}()
    for line in eachline(path)
        s = strip(line)
        (isempty(s) || s == "END") && continue
        i = findfirst('=', s)
        i === nothing && continue
        key = strip(s[1:prevind(s, i)])
        val = strip(s[nextind(s, i):end])
        (startswith(key, "GROUP") || startswith(key, "END_GROUP")) && continue
        out[String(key)] = String(strip(val, ['"']))
    end
    isempty(out) && error("$path parsed to no fields; is it an MTL?")
    return out
end

"""
    mtl_corners(fields) -> NamedTuple

The four scene corners in projected metres, as `(ul, ur, ll, lr)` of `(x, y)`.

The `CORNER_*_PROJECTION_{X,Y}_PRODUCT` fields. On a north-up product — which every Landsat L1 is — these
are the **axis-aligned bounding box** of the raster and not the imaged swath, so slopes taken from them are
0 and ±90 rather than the scan directions. Kept for that comparison, since establishing what a field does
*not* carry is worth as much as establishing what it does.

A scene missing them is an error rather than a default: a `Float64` zero would put the footprint at the
origin and produce two plausible angles from nothing.
"""
function mtl_corners(fields::AbstractDict)
    get2(c) = let kx = "CORNER_$(c)_PROJECTION_X_PRODUCT", ky = "CORNER_$(c)_PROJECTION_Y_PRODUCT"
        haskey(fields, kx) && haskey(fields, ky) ||
            error("MTL has no $kx / $ky; these are the projected product corners the scan " *
                  "geometry is read from, and a scene without them cannot be destriped from metadata")
        (parse(Float64, fields[kx]), parse(Float64, fields[ky]))
    end
    return (; ul = get2("UL"), ur = get2("UR"), ll = get2("LL"), lr = get2("LR"))
end

"""
    scan_angles(corners; spacing) -> (along_track, cross_track)

The along- and cross-track angles in degrees, on the reference's own convention.

Reproduces `_get_slopes` (`autoRIFT.py:138-151`) exactly, and the exactness matters because it is a
*choice* rather than a derivation: the reference takes `nanmax` of two slopes per axis, so a footprint whose
opposite edges differ slightly resolves to the larger angle rather than to their mean. Reading corners from
metadata changes where the slopes come from, not which of the two is picked.

Its `_calculate_slope` is `atan((y1 - y2) / (x1 - x2))` in degrees — note `atan` of a ratio and not
`atan2`, so the result is in `(-90, 90)` and a vertical edge gives `±90` by the division overflowing rather
than by a branch. Along-track is from the bottom and top edges, cross-track from the two side edges.

**Corners arrive in projected metres and the reference's are in pixels.** A slope is a ratio, so a uniform
scale cancels — but `spacing` is required rather than defaulted because a non-square pixel does *not*
cancel, and silently assuming square is the kind of thing that produces a plausible wrong angle.
"""
function scan_angles(corners; spacing::Tuple{Real,Real})
    sx, sy = Float64(spacing[1]), Float64(spacing[2])
    # Into pixel units, where the reference's slopes are taken. The y sign is irrelevant to a `nanmax`
    # over both edges of an axis, and is left as the projection has it.
    px(p) = (p[1] / sx, p[2] / sy)
    ul, ur, ll, lr = px(corners.ul), px(corners.ur), px(corners.ll), px(corners.lr)
    slope(a, b) = rad2deg(atan((a[2] - b[2]) / (a[1] - b[1])))
    # `_get_slopes(tl, tr, bl, br)`: along-track from the bottom and top edges, cross-track from the sides.
    along = maximum(skip_nan((slope(ll, lr), slope(ul, ur))))
    cross = maximum(skip_nan((slope(lr, ur), slope(ll, ul))))
    return (along, cross)
end

# `skip_nan`: the reference uses `np.nanmax`, which ignores a `NaN` slope rather than propagating it. A `NaN`
# arises when two corners coincide, which a degenerate footprint can produce.
skip_nan(t) = (v = filter(!isnan, collect(t)); isempty(v) ? [NaN] : v)

"""
    mtl_scan_angles(path; spacing) -> (along_track, cross_track)

[`scan_angles`](@ref) of [`mtl_corners`](@ref) of [`mtl_fields`](@ref), for the common case.
"""
mtl_scan_angles(path::AbstractString; spacing) =
    scan_angles(mtl_corners(mtl_fields(path)); spacing)

"""
    reference_scan_angles(log) -> Vector{NamedTuple}

The along- and cross-track angles the reference *printed*, one entry per filtered scene.

`_fft_filter` logs both before it decides anything (`autoRIFT.py:201-202`), so the container log carries the
values its own pixel-derived route produced. That is what the metadata route is checked against — the gate
is agreement with the reference, and its log is the only place its intermediate geometry is visible.
"""
function reference_scan_angles(log::AbstractString)
    out = NamedTuple{(:along, :cross),Tuple{Float64,Float64}}[]
    along = nothing
    for line in eachline(log)
        m = match(r"Along track angle is\s+(-?[\d.]+)\s+degrees", line)
        m === nothing || (along = parse(Float64, m.captures[1]); continue)
        m = match(r"Cross track angle is\s+(-?[\d.]+)\s+degrees", line)
        if m !== nothing && along !== nothing
            push!(out, (; along, cross = parse(Float64, m.captures[1])))
            along = nothing
        end
    end
    return out
end

"""
    ang_fields(path) -> Dict{String,Vector{Float64}}

Every numeric field of a Landsat `_ANG.txt`, scalars as one-element vectors.

A separate parser from [`mtl_fields`](@ref) because the format is not the same one: the fields that
matter here are **parenthesised lists spanning several lines**, `EPHEMERIS_ECEF_X = (-2513196.897534,
-2519426.030063, ...`, which a line-at-a-time `KEY = VALUE` reader truncates to its first row without
complaining. Numeric rather than textual for the same reason `mtl_fields` is textual: every consumer
here wants numbers, and a malformed ephemeris should fail at the parse rather than reach a bearing.
"""
function ang_fields(path::AbstractString)
    isfile(path) || error("no ANG at $path")
    out = Dict{String,Vector{Float64}}()
    text = read(path, String)
    key = nothing
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
    isempty(out) && error("$path parsed to no numeric fields; is it an ANG?")
    return out
end

_numbers(s::AbstractString) =
    [parse(Float64, m.match) for m in eachmatch(r"-?\d+\.?\d*(?:[eE][-+]?\d+)?", s)]

"""
    ang_ephemeris(path) -> NamedTuple

The satellite's trajectory from `_ANG.txt`: `time` in seconds of day and `x`, `y`, `z` in ECEF metres.

`GROUP = EPHEMERIS` carries `NUMBER_OF_POINTS` samples at 1 s spacing, centred on the acquisition. The
count is checked against all four vectors rather than trusted, because a truncated list would otherwise
produce a shorter trajectory and a bearing taken from the wrong pair of points.
"""
function ang_ephemeris(path::AbstractString)
    f = ang_fields(path)
    get1(k) = haskey(f, k) ? f[k] :
              error("$path has no $k; the ephemeris is what the scan geometry is derived from")
    t, x, y, z = get1("EPHEMERIS_TIME"), get1("EPHEMERIS_ECEF_X"),
                 get1("EPHEMERIS_ECEF_Y"), get1("EPHEMERIS_ECEF_Z")
    n = haskey(f, "NUMBER_OF_POINTS") ? Int(first(f["NUMBER_OF_POINTS"])) : length(t)
    all(v -> length(v) == n, (t, x, y, z)) || error(
        "$path declares $n ephemeris points but parsed " *
        "$(length(t))/$(length(x))/$(length(y))/$(length(z)); the list parse is wrong")
    return (; time = t, x, y, z)
end

"""
    orbit_scan_angles(ang_path, epsg; spacing) -> (along_track, cross_track)

The scan geometry from the **orbit**, in the reference's own angle convention.

Two consecutive ephemeris positions transformed into the raster's CRS give the ground-track direction
directly, which is what `_fft_filter` is trying to recover from the valid-data region's shape. The pair
is taken at the middle of the ephemeris, which is where the acquisition is centred.

**The projection step is not optional.** A slope on a UTM raster is a *grid* bearing, and grid north
departs from true north by the meridian convergence — 3.2° at `LT05_L1TP_060018`, more near the poles.
Taking the bearing from ECEF or from lat/lon and calling it a raster angle would be wrong by exactly
that.

The two returned angles are perpendicular by construction, which is the property the reference's pair
does *not* have: it logs 71.60° and −20.25°, which are 91.85° apart. `spacing` is required for the
reason it is in [`scan_angles`](@ref) — a non-square pixel does not cancel out of a slope.
"""
function orbit_scan_angles(ang_path::AbstractString, epsg::Integer; spacing::Tuple{Real,Real})
    e = ang_ephemeris(ang_path)
    k = length(e.time) ÷ 2
    # Geocentric ECEF is EPSG:4978. `order = :trad` keeps the projected result (x, y).
    p = Vector{Tuple{Float64,Float64}}(undef, 2)
    ArchGDAL.crs2transform(EPSG(4978), EPSG(Int(epsg)); order = :trad) do tf
        for (j, i) in enumerate((k, k + 1))
            q = ArchGDAL.createpoint(e.x[i], e.y[i], e.z[i])
            ArchGDAL.transform!(q, tf)
            p[j] = (ArchGDAL.getx(q, 0), ArchGDAL.gety(q, 0))
        end
    end
    # **Magnitudes, so the frame stays y-up.** `spacing` is only here to undo pixel anisotropy — an
    # angle on a grid sampled differently in x and y is not the angle in metres — and dividing by a
    # north-up raster's *signed* `-30` would flip the bearing into a row-down frame and negate both
    # results. Measured on both scenes of `LT05_L1TP_060018`: y-up reproduces the reference's signs,
    # row-down inverts them.
    sx, sy = abs(Float64(spacing[1])), abs(Float64(spacing[2]))
    dx = (p[2][1] - p[1][1]) / sx
    dy = (p[2][2] - p[1][2]) / sy
    # The reference's convention: `atan` of a ratio, in degrees, so the result is in (-90, 90) and
    # carries no quadrant. The ground track is the **along**-track direction and the cross-track is its
    # perpendicular, folded back into the same open interval.
    along = rad2deg(atan(dy / dx))
    cross = along + (along > 0 ? -90.0 : 90.0)
    return (along, cross)
end

"""
    ang_path(scene_path) -> String

The `_ANG.txt` beside a Landsat band raster, on whatever route reached the band.

Derived from the band's own path rather than resolved separately, so a scene read from `/vsis3` finds
its ephemeris in the same bucket and one read from disk finds it in the same directory.
"""
ang_path(scene::AbstractString) =
    replace(scene, r"_B\d+\.TIF$"i => "_ANG.txt")

"""
    vsi_text(path) -> String

A whole text object read through GDAL's virtual filesystem.

The route exists because an `_ANG.txt` is only reachable over requester-pays S3: the STAC item's plain
`https` href redirects to an HTML landing page rather than the object, so `Downloads.download` returns a
login form that parses to no fields. `/vsis3` with `AWS_REQUEST_PAYER` is the same access this harness
already uses for the band rasters, and `ArchGDAL.GDAL` re-exports the VSI calls, so this needs no
dependency `ArchGDAL` does not already bring.
"""
function vsi_text(path::AbstractString)
    G = ArchGDAL.GDAL
    ArchGDAL.setconfigoption("AWS_REQUEST_PAYER", "requester")
    h = G.vsifopenl(path, "rb")
    h == C_NULL && error("cannot open $path through GDAL's virtual filesystem; a requester-pays " *
                         "object needs AWS_PROFILE to name credentials that can pay")
    try
        G.vsifseekl(h, 0, 2)                 # SEEK_END
        n = Int(G.vsiftelll(h))
        G.vsifseekl(h, 0, 0)
        buf = Vector{UInt8}(undef, n)
        got = G.vsifreadl(pointer(buf), 1, n, h)
        Int(got) == n || error("short read of $path: $got of $n bytes")
        return String(buf)
    finally
        G.vsifclosel(h)
    end
end

"""
    scene_ang(name, cache) -> String

A local copy of `name`'s `_ANG.txt`, fetched once through the STAC item and cached.

Cached because the gate reads it on every run and the object is requester-pays: 34 KiB is not the cost,
the round trip is. The name is the granule's, so a cache entry is unambiguous.
"""
function scene_ang(name::AbstractString, cache::AbstractString)
    local_path = joinpath(cache, name * "_ANG.txt")
    isfile(local_path) && return local_path
    mkpath(cache)
    url = "https://landsatlook.usgs.gov/stac-server/collections/landsat-c2l1/items/$name"
    item = JSON3.read(String(take!(Downloads.download(url, IOBuffer(); timeout = 60))))
    asset = get(item.assets, Symbol("ANG.txt"), nothing)
    asset === nothing && error("STAC item $name has no `ANG.txt` asset")
    href = asset.alternate.s3.href
    text = vsi_text("/vsis3/" * href[6:end])
    # Through a temporary, so an interrupted fetch does not leave a truncated file that parses.
    tmp = local_path * ".partial"
    write(tmp, text)
    mv(tmp, local_path; force = true)
    return local_path
end

"""
    orbit_angle_check(c::GoldenCase, run, cache) -> Vector{NamedTuple}

Each filtered scene's orbit-derived scan angles beside the reference's own logged pair.

The comparison the register turns on: whether the `_ANG.txt` ephemeris reproduces what `_fft_filter`
recovers from the valid-data region's shape, and where it does not, whether the orbit or the pixels are
the better answer. Scenes are taken in **job order**, since `apply_landsat_filtering(reference,
secondary)` filters and logs them in that order rather than in acquisition order.

**Each scene's EPSG and spacing come from the scene itself**, one per scene: the two scenes of a
cross-zone pair are in different projections, and a bearing is a grid bearing, so using one CRS for both
would put the second scene's angles several degrees out. The granule rather than the run's `filtered/`
copy, which carries the same native grid and which [`prune_run`](@ref) deletes as regenerable.
"""
function orbit_angle_check(c::GoldenCase, run::AbstractString, cache::AbstractString)
    logged = reference_scan_angles(joinpath(run, "capture.log"))
    names = [first(c.reference), first(c.secondary)]
    length(logged) >= length(names) || error(
        "$(run)/capture.log logs $(length(logged)) angle pairs for $(length(names)) filtered " *
        "scenes; the log is truncated or this is not an L4/L5 pair")
    # `names` is in job order and `scene_path` takes acquisition order, so the two are paired by name
    # rather than by position.
    early, late = acquisition_order(c)
    paths = Dict(early => scene_path(c, :reference), late => scene_path(c, :secondary))
    out = NamedTuple[]
    for (i, name) in enumerate(names)
        haskey(paths, name) || error("\"$name\" is neither acquisition of $(c.product)")
        ds = ArchGDAL.read(paths[name])
        epsg = parse(Int, ArchGDAL.toEPSG(ArchGDAL.importWKT(ArchGDAL.getproj(ds))) |> string)
        gt = ArchGDAL.getgeotransform(ds)
        along, cross = orbit_scan_angles(scene_ang(name, cache), epsg; spacing = (gt[2], gt[6]))
        r = logged[i]
        push!(out, (; name, epsg, along, cross, ref_along = r.along, ref_cross = r.cross,
                    d_along = along - r.along, d_cross = cross - r.cross,
                    ours_apart = abs(along - cross), ref_apart = abs(r.along - r.cross)))
    end
    return out
end

"""
    reference_banding(log) -> Vector{Bool}

Whether the reference's band-reject **fired**, one entry per filtered scene.

`_fft_filter` prints its two band powers unconditionally and adds a "No banding filter applied" line
only when it declines (`autoRIFT.py:211-226`), so the decision is readable from the log and does not
have to be inferred from the output — which would be circular, since this is what a rung comparing that
output needs to know.

The decision is a **binary branch on a ratio**, and it can be marginal: on
`LT05_L1GS_001013_19920425` the powers are 1588 and 3279, clearing the `>= 2` test by 3.2%. A scene
that close can be decided differently by two implementations whose input fields differ slightly, and
then its whole output differs — a declined reject returns the clamped input, a fired one returns the
band-rejected field.
"""
function reference_banding(log::AbstractString)
    out = Bool[]
    pending = false
    for line in eachline(log)
        if occursin(r"Cross track power is", line)
            pending && push!(out, true)
            pending = true
        elseif occursin("No banding filter applied", line)
            pending && (push!(out, false); pending = false)
        end
    end
    pending && push!(out, true)
    return out
end
