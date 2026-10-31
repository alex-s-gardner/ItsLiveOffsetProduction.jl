# The pass an optical scene was imaged on, from the satellite's own trajectory at the scene: `"ascending"`
# when it moves north there, `"descending"` when it moves south. A radar scene's comes from its own
# metadata (`s1.jl`, `nisar.jl`).
#
# Landsat's `_ANG.txt` holds the scene's ephemeris. A Sentinel-2 product's datastrip metadata
# (`MTD_DS.xml`) holds the GPS fixes of the whole datastrip, which can run over the pole; the direction is
# the fix nearest the tile's own sensing time (`MTD_TL.xml`), not the datastrip's
# `SENSING_ORBIT_DIRECTION`, which names its start.

# The datastrip metadata is several megabytes, most of it attitude and quality records after the GPS
# fixes; this much of its start holds them in every product measured (at most about 650 kB). The tile
# metadata gives its sensing time within its first kilobyte.
const S2_DATASTRIP_HEAD_BYTES = 2 << 20
const S2_TILE_HEAD_BYTES = 16 << 10

_pass(vz::Float64) = vz > 0 ? "ascending" : "descending"

"""
    landsat_flight_direction(ang_text) -> String

The pass of a Landsat scene from the ephemeris in its `_ANG.txt`: the sign of the change in ECEF `z` over the
middle step.
"""
function landsat_flight_direction(ang_text::String)
    e = ang_ephemeris(ang_text)
    k = length(e.time) ÷ 2
    k >= 1 || throw(ArgumentError("the _ANG.txt ephemeris has fewer than two points"))
    return _pass(e.z[k + 1] - e.z[k])
end

"""
    s2_flight_direction(datastrip_head, sensing_time) -> String

The pass of a Sentinel-2 tile from its datastrip's GPS fixes (`datastrip_head`, the start of `MTD_DS.xml`):
the sign of the ECEF `z` velocity at the fix nearest `sensing_time`. Throws when no fix is within a minute.
"""
function s2_flight_direction(datastrip_head::String, sensing_time::DateTime)
    best, dt = 0.0, typemax(Int)
    for m in eachmatch(r"<GPS_Point>(.*?)</GPS_Point>"s, datastrip_head)
        p = m.captures[1]
        t = match(r"<GPS_TIME>([^<]+)</GPS_TIME>", p)
        v = match(r"<VELOCITY_VALUES[^>]*>([^<]+)</VELOCITY_VALUES>", p)
        (t === nothing || v === nothing) && continue
        d = abs(Dates.value(parse_utc(String(t.captures[1])) - sensing_time))
        d < dt || continue
        vel = _numbers(String(v.captures[1]))
        length(vel) == 3 || continue
        best, dt = vel[3], d
    end
    dt <= 60_000 || throw(ArgumentError(
        "the datastrip has no GPS fix within a minute of the tile's sensing time; its pass is unknown"))
    return _pass(best)
end

# `<SENSING_TIME>` of a Sentinel-2 tile's `MTD_TL.xml`.
function _s2_sensing_time(tile_xml::String)
    m = match(r"<SENSING_TIME[^>]*>([^<]+)</SENSING_TIME>", tile_xml)
    m === nothing && throw(ArgumentError("the tile metadata has no SENSING_TIME"))
    return parse_utc(String(m.captures[1]))
end

"""
    flight_direction(scene::SceneSource) -> String

The pass of an optical scene: [`landsat_flight_direction`](@ref) of a Landsat scene's `_ANG.txt`, or
[`s2_flight_direction`](@ref) of a Sentinel-2 product's metadata on the public mirror, found through the
product's `manifest.safe` as its band is.
"""
function flight_direction(scene::SceneSource)
    startswith(scene.name, "S2") || return landsat_flight_direction(read_text(scene.ang))
    safe = s2_safe_url(scene.name)
    manifest = read_text(safe * "/manifest.safe")
    tile = _s2_tile(scene.name)
    ds = _manifest_entry(manifest, r"DATASTRIP/[^\"<>]*/MTD_DS\.xml", scene.name)
    tl = _manifest_entry(manifest, Regex("GRANULE/[^\"<>]*" * tile * "[^\"<>]*/MTD_TL\\.xml"), scene.name)
    # The two heads are independent reads.
    dsh = Threads.@spawn _object_head(safe * "/" * ds, S2_DATASTRIP_HEAD_BYTES)
    tlh = _object_head(safe * "/" * tl, S2_TILE_HEAD_BYTES)
    return s2_flight_direction(_await(dsh)::String, _s2_sensing_time(tlh))
end

# The first `n` bytes of an object, or all of it when it is shorter.
function _object_head(url::String, n::Int)
    t = ChunkManifests.HTTPTransport()
    m = min(UInt64(n), ChunkManifests.objectsize(t, url))
    return String(ChunkManifests.fetchrange(t, url, ChunkManifests.ByteRange(UInt64(0), m)))
end

function _manifest_entry(manifest::String, pattern::Regex, name::String)
    hits = unique(String[m.match for m in eachmatch(pattern, manifest)])
    length(hits) == 1 || throw(ArgumentError("expected one " * pattern.pattern * " in " * name * "'s manifest.safe, found " *
                                             string(length(hits))))
    return hits[1]
end
