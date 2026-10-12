# Sentinel-1 orbit files: a product carries no state vectors of its own, so its precise orbit (POEORB) is
# found by acquisition time in the public `s1-orbits` bucket, a mirror of the Copernicus orbit products.

const S1_ORBITS = "https://s1-orbits.s3.us-west-2.amazonaws.com"

# A product name's mission and sensing interval: `S1A_IW_SLC__1SSH_20150828T162412_20150828T162431_…`.
function s1_name_times(name::String)
    m = match(r"^(S1[ABCD])_\w\w_\w{4}_\w{4}_(\d{8}T\d{6})_(\d{8}T\d{6})_", name)
    m === nothing && throw(ArgumentError("cannot read a mission and sensing times out of " * name))
    c = m.captures
    return String(something(c[1])), _ccsds(String(something(c[2]))), _ccsds(String(something(c[3])))
end

# `YYYYMMDDTHHMMSS`.
function _ccsds(s::String)
    length(s) == 15 && s[9] == 'T' || throw(ArgumentError("not a YYYYMMDDTHHMMSS time: " * s))
    return DateTime(parse(Int, s[1:4]), parse(Int, s[5:6]), parse(Int, s[7:8]), parse(Int, s[10:11]),
                    parse(Int, s[12:13]), parse(Int, s[14:15]))
end

# The characters a query value may hold unescaped (RFC 3986 unreserved).
_unreserved(c::UInt8) = (UInt8('A') <= c <= UInt8('Z')) || (UInt8('a') <= c <= UInt8('z')) ||
                        (UInt8('0') <= c <= UInt8('9')) || c in (UInt8('-'), UInt8('_'), UInt8('.'), UInt8('~'))

function _query_escape(s::String)
    out = UInt8[]
    hex = codeunits("0123456789ABCDEF")
    for c in codeunits(s)
        if _unreserved(c)
            push!(out, c)
        else
            push!(out, UInt8('%'), hex[(c >> 4) + 1], hex[(c & 0x0f) + 1])
        end
    end
    return String(out)
end

"""
    s1_orbit_url(name) -> String

The URL of the precise orbit file for the Sentinel-1 product `name`: of the `AUX_POEORB` files of its
mission whose validity spans the sensing interval with a minute to spare at each end, the most recently
created. The mission's files are listed a page at a time, a thousand keys per request.
"""
function s1_orbit_url(name::String)
    mission, start, stop = s1_name_times(name)
    prefix = "AUX_POEORB/" * mission * "_OPER_AUX_POEORB_OPOD_"
    want0, want1 = start - Dates.Minute(1), stop + Dates.Minute(1)
    best = ""
    best_created = DateTime(0)
    token = ""
    t = ChunkManifests.HTTPTransport()
    while true
        url = S1_ORBITS * "/?list-type=2&prefix=" * _query_escape(prefix) *
              (isempty(token) ? "" : "&continuation-token=" * _query_escape(token))
        page = String(ChunkManifests.fetchall(t, url))
        for k in eachmatch(r"<Key>([^<]+)</Key>", page)
            key = String(something(k.captures[1]))
            m = match(r"_OPOD_(\d{8}T\d{6})_V(\d{8}T\d{6})_(\d{8}T\d{6})\.EOF$", key)
            m === nothing && continue
            c = m.captures
            created, v0, v1 = _ccsds(String(something(c[1]))), _ccsds(String(something(c[2]))),
                              _ccsds(String(something(c[3])))
            (v0 <= want0 && want1 <= v1 && created > best_created) || continue
            best, best_created = key, created
        end
        n = match(r"<NextContinuationToken>([^<]+)</NextContinuationToken>", page)
        n === nothing && break
        token = String(something(n.captures[1]))
    end
    isempty(best) && throw(ArgumentError("s1-orbits holds no precise orbit spanning " * name))
    return S1_ORBITS * "/" * best
end

"""
    fetch_s1_orbit(name, dir) -> String

The precise orbit file for product `name`, downloaded into `dir`; its path.
"""
function fetch_s1_orbit(name::String, dir::String)
    url = s1_orbit_url(name)
    path = joinpath(dir, basename(url))
    isfile(path) && return path
    write(path, ChunkManifests.fetchall(ChunkManifests.HTTPTransport(), url))
    return path
end
