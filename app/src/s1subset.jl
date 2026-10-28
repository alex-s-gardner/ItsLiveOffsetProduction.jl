# A Sentinel-1 burst pair given as SLC zips and a polygon: the bursts of the zips that the polygon covers,
# read out of the zips rather than from ASF's burst extractor.
#
# Each scene is one or more SLC zips of a datatake, comma-separated. The bursts kept are those at least half
# of whose footprint — from the annotation's geolocation grid — lies inside the polygon, in either scene,
# and that both scenes hold. Half, because consecutive bursts overlap in azimuth: a polygon drawn around a
# set of bursts, such as an OPERA frame's, also touches the edges of the bursts on either side. A burst is
# the same ground on every pass of its track, starting the same time after the ascending node to well under
# a burst's length (`azimuthAnxTime`), which is what its burst ID counts, so the two scenes' bursts are
# matched by subswath and that time.
#
# Footprints are compared in longitude and latitude, in one frame of continuous longitudes centred on the
# polygon. The polygon is read as written, its edges straight in longitude and latitude, so one crossing
# the antimeridian is given split at ±180° into a `MULTIPOLYGON`, as OPERA stores such frames; each piece
# is moved by whole turns into the frame, which rejoins them. A burst footprint is a fine grid of points,
# so a step of more than half a turn between neighbours is a crossing: its longitudes are made continuous,
# then moved into the frame. Two slices of one datatake overlap and
# can list the same burst; it is taken once, from the zip named first.

"""
    is_zip_list(input) -> Bool

Whether a scene's input is SLC zips, comma-separated, rather than one product or extractor URLs.
"""
is_zip_list(input::String) = all(u -> endswith(strip(u), ".zip"), split(input, ','; keepempty = false))

# The fraction of a burst's footprint the polygon must cover for the burst to be kept.
const S1_SUBSET_COVER = 0.5

# Bursts of one track and subswath on two passes start within this much of the same time after the
# ascending node. A burst lasts about 2.76 s, and the ESA burst ID grid is that wide; the orbit is
# controlled to a fraction of a second.
const S1_ANX_MATCH = 1.0

"""
    SubsetBurst

One burst a scene's zips hold: where it is, its subswath and index, its time after the ascending node and
its footprint ring of `(longitude, latitude)` vertices.
"""
struct SubsetBurst
    zip::String
    swath::Int
    burst::Int
    anx_time::Float64
    ring::Ring
end

# `lon` moved by whole turns to within half a turn of `center`.
_near(lon::Float64, center::Float64) = lon - 360 * round((lon - center) / 360)

"""
    unwrap(ring, center) -> Ring

`ring` with its longitudes made continuous — each vertex within half a turn of the one before — and then
shifted by whole turns so that its first vertex is within half a turn of `center`.
"""
function unwrap(ring::Ring, center::Float64)
    out = similar(ring)
    prev = _near(first(ring)[1], center)
    for (i, (x, y)) in enumerate(ring)
        prev = _near(x, prev)
        out[i] = (prev, y)
    end
    return out
end

"""
    subset_area(wkt) -> (Vector{Vector{Ring}}, Float64)

The polygons `wkt` describes, each moved by whole turns of longitude into one frame, and that frame's
centre: the circular mean of the polygons' longitude midpoints. The pieces of a polygon split at the
antimeridian then share edges and join.

Refused: a polygon wider than half a turn of longitude, which reads differently according to which way
round it is taken, and a set of polygons that together span more than half a turn, which no frame of
continuous longitudes holds.
"""
function subset_area(wkt::String)
    polys = wkt_polygons(wkt)
    sx = sy = 0.0
    for p in polys
        lo, hi = extrema(first, first(p))
        hi - lo <= 180 || throw(ArgumentError(
            "\"" * wkt * "\" has a polygon " * string(round(hi - lo; digits = 1)) * "° wide in longitude; split " *
            "one crossing the antimeridian at ±180° into a MULTIPOLYGON"))
        sx += cosd((lo + hi) / 2)
        sy += sind((lo + hi) / 2)
    end
    center = atand(sy, sx)
    out = [Ring[[(x - 360 * round(((lo + hi) / 2 - center) / 360), y) for (x, y) in r] for r in p]
           for p in polys for (lo, hi) in (extrema(first, first(p)),)]
    lo, hi = extrema(x for p in out for (x, _) in first(p))
    hi - lo <= 180 || throw(ArgumentError(
        "\"" * wkt * "\" spans " * string(round(hi - lo; digits = 1)) * "° of longitude; a subset spans at most 180°"))
    return out, center
end

# Every burst of every zip of a scene in `polarization`, a burst listed by two overlapping zips once, from
# the one listed first.
function scene_bursts(input::String, polarization::String)
    zips = String[String(strip(u)) for u in split(input, ','; keepempty = false)]
    out = SubsetBurst[]
    for z in zips, sw in S1_SWATHS
        for (k, f) in enumerate(SLCDatasets.burst_footprints(z, sw, polarization))
            any(b -> b.swath == sw && abs(b.anx_time - f.anx_time) < S1_ANX_MATCH, out) && continue
            push!(out, SubsetBurst(z, sw, k, f.anx_time, f.ring))
        end
    end
    return out
end

"""
    subset_bursts(reference, secondary, wkt, polarization) -> (Vector{SubsetBurst}, Vector{SubsetBurst})

The bursts of each scene's zips a burst pair over polygon `wkt` (geographic, `(longitude latitude)`
vertices) correlates: those `S1_SUBSET_COVER` of whose footprint the polygon covers in either
scene, held by both, matched burst for burst and in the same order in the two vectors.
"""
function subset_bursts(reference::String, secondary::String, wkt::String, polarization::String)
    polys, center = subset_area(wkt)
    rs, ss = scene_bursts(reference, polarization), scene_bursts(secondary, polarization)
    pairs = Tuple{SubsetBurst,SubsetBurst}[]
    cover(area, b) = ring_overlap(area, unwrap(b.ring, center))
    with_polygons(polys) do area
        for r in rs
            k = findfirst(s -> s.swath == r.swath && abs(s.anx_time - r.anx_time) < S1_ANX_MATCH, ss)
            k === nothing && continue
            s = ss[k]
            max(cover(area, r), cover(area, s)) >= S1_SUBSET_COVER && push!(pairs, (r, s))
        end
    end
    isempty(pairs) && throw(ArgumentError(
        "no burst both scenes hold lies inside the subset polygon; check that the polygon is in longitude " *
        "and latitude and that the zips are of the same track"))
    return SubsetBurst[p[1] for p in pairs], SubsetBurst[p[2] for p in pairs]
end

"""
    subset_spec(bursts, polarization) -> String

The bursts as the scene input a burst pair is run from, each `<zip>#IW<n>/<POL>/<index>` with the index
counted from zero, comma-separated.
"""
subset_spec(bursts::Vector{SubsetBurst}, polarization::String) =
    join([SLCDatasets.zip_burst_spec(SLCDatasets.ZipBurst(b.zip, b.swath, polarization, b.burst)) for b in bursts], ",")

"""
    subset_inputs(reference, secondary, wkt) -> (String, String)

The two scenes' inputs as the bursts of their SLC zips ([`subset_spec`](@ref)) polygon `wkt` selects
([`subset_bursts`](@ref)), in the like-polarized channel both zips' names carry.
"""
function subset_inputs(reference::String, secondary::String, wkt::String)
    name(input) = String(chop(basename(String(strip(first(split(input, ','))))); tail = 4))
    pol = uppercase(s1_polarization(name(reference)))
    pol == uppercase(s1_polarization(name(secondary))) || throw(ArgumentError(
        "the two scenes' zips carry different like-polarized channels, so their bursts make no pair"))
    rb, sb = subset_bursts(reference, secondary, wkt, pol)
    return subset_spec(rb, pol), subset_spec(sb, pol)
end

"""
    is_zip_burst_list(input) -> Bool

Whether a burst pair's input is bursts of SLC zips, as [`subset_spec`](@ref) writes them.
"""
is_zip_burst_list(input::String) = SLCDatasets.is_zip_burst_spec(strip(first(split(input, ','))))

"""
    zip_burst_list(input, dir) -> Vector{SLCDatasets.ZipBurst}

The bursts a comma-separated list of zip bursts names, each raster inflated into `dir` when it is first
read.
"""
zip_burst_list(input::String, dir::Union{Nothing,String}) =
    [SLCDatasets.ZipBurst(strip(u); dir) for u in split(input, ','; keepempty = false)]
