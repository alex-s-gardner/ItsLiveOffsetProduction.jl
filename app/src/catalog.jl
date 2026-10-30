# A product's catalog records, as the ITS_LIVE metadata step (`itslive-metadata`'s
# `cryoforge/generate.py`, which `hyp3_itslive_metadata` runs after the product is published) writes
# them: its STAC item (`<product>.stac.json`) and NSIDC's premet and spatial files (`<product>.nc.premet`,
# `<product>.nc.spatial`). Every value is the reference's arithmetic on what the product records of
# itself: its footprint through PROJ to longitude and latitude, rounded as `numpy.round` rounds; its
# dates as pandas parses and prints them; the rest copied from `img_pair_info`.

# The ITS_LIVE data version the reference writes into every record.
const ITS_LIVE_DATA_VERSION = "002"

# ---------------------------------------------------------------------------
# Numbers and dates as the reference prints them
# ---------------------------------------------------------------------------

# `numpy.round(x, decimals = d)`: `x` scaled, rounded half to even, and scaled back.
_np_round(x::Float64, d::Int) = (s = 10.0^d; round(x * s) / s)

# Python's `round(x, d)`, correctly rounded from `x`'s exact binary value.
_py_round(x::Float64, d::Int) = parse(Float64, Base.Ryu.writefixed(x, d))

"""
    PyDateTime

An instant to the microsecond, as pandas holds the dates a product records: Julia's `DateTime` keeps only
milliseconds.
"""
struct PyDateTime
    date::Date
    hour::Int
    minute::Int
    second::Int
    microsecond::Int
end

# `YYYYMMDDTHH:MM:SS[.f…]`, as `img_pair_info` writes its dates: `pd.to_datetime` of one.
function parse_pair_date(s::String)
    m = match(r"^(\d{4})(\d{2})(\d{2})T(\d{2}):(\d{2}):(\d{2})(?:\.(\d*))?$", s)
    m === nothing && throw(ArgumentError("\"" * s * "\" is not a YYYYMMDDTHH:MM:SS.f date"))
    c = String[something(x, "") for x in m.captures]
    frac = c[7]
    length(frac) <= 6 || throw(ArgumentError("\"" * s * "\" holds a date finer than a microsecond"))
    us = isempty(frac) ? 0 : parse(Int, frac) * 10^(6 - length(frac))
    return PyDateTime(Date(parse(Int, c[1]), parse(Int, c[2]), parse(Int, c[3])), parse(Int, c[4]),
                      parse(Int, c[5]), parse(Int, c[6]), us)
end

# `DD-Mon-YYYY HH:MM:SS`, as a product writes `date_created`: `pd.to_datetime` of one.
function parse_created_date(s::String)
    m = match(r"^(\d{2})-([A-Z][a-z]{2})-(\d{4}) (\d{2}):(\d{2}):(\d{2})$", s)
    m === nothing && throw(ArgumentError("\"" * s * "\" is not a DD-Mon-YYYY HH:MM:SS date"))
    c = String[something(x, "") for x in m.captures]
    month = findfirst(==(c[2]), ItsLiveOffsetProduction._MONTH_ABBREVIATIONS)
    month === nothing && throw(ArgumentError("\"" * s * "\" names no month"))
    return PyDateTime(Date(parse(Int, c[3]), month, parse(Int, c[1])), parse(Int, c[4]), parse(Int, c[5]),
                      parse(Int, c[6]), 0)
end

_two(n::Int) = string(n; pad = 2)

# `Timestamp.isoformat()` of a UTC instant with `+00:00` written `Z`: the microseconds only when not zero.
function iso_utc(t::PyDateTime)
    d = t.date
    s = string(Dates.year(d); pad = 4) * "-" * _two(Dates.month(d)) * "-" * _two(Dates.day(d)) * "T" *
        _two(t.hour) * ":" * _two(t.minute) * ":" * _two(t.second)
    t.microsecond == 0 || (s *= "." * string(t.microsecond; pad = 6))
    return s * "Z"
end

_ymd(t::PyDateTime) = string(Dates.year(t.date); pad = 4) * "-" * _two(Dates.month(t.date)) * "-" *
                      _two(Dates.day(t.date))
# `strftime('%H:%M:%S')` and the milliseconds, truncated.
_hms_ms(t::PyDateTime) = _two(t.hour) * ":" * _two(t.minute) * ":" * _two(t.second) * "." *
                         string(t.microsecond ÷ 1000; pad = 3)

# ---------------------------------------------------------------------------
# The footprint
# ---------------------------------------------------------------------------

"""
    Footprint4326

A product's footprint in longitude and latitude, as `cryoforge`'s `get_geom` computes it: the grid's
outer edge walked counterclockwise from its lower-left corner, four segments to an edge, each vertex
projected and rounded to four decimals; the bounding box of those vertices; the grid's centre, rounded
likewise; and the four corners from the upper left, clockwise.
"""
struct Footprint4326
    ring::Vector{Tuple{Float64,Float64}}
    bbox::NTuple{4,Float64}
    center::Tuple{Float64,Float64}
    corners::NTuple{4,Tuple{Float64,Float64}}
end

function footprint_4326(s::ItsLiveOffsetProduction.ProductSummary)
    gt = split(s.geotransform)
    length(gt) == 6 || throw(ArgumentError("GeoTransform \"" * s.geotransform * "\" is not six numbers"))
    px, py = parse(Float64, gt[2]), parse(Float64, gt[6])
    # Cell centres to the grid's outer edge; `py` is negative.
    minx, maxx = s.x[1] - px / 2.0, s.x[end] + px / 2.0
    miny, maxy = s.y[end] + py / 2.0, s.y[1] - py / 2.0
    with_transform(s.epsg, 4326) do ct
        lonlat(x, y) = (lon_lat = transform_point(ct, x, y);
                        (_np_round(lon_lat[1], 4), _np_round(lon_lat[2], 4)))
        ll, lr, ur, ul = lonlat(minx, miny), lonlat(maxx, miny), lonlat(maxx, maxy), lonlat(minx, maxy)
        center = lonlat((s.x[1] + s.x[end]) / 2.0, (s.y[1] + s.y[end]) / 2.0)
        ring = Tuple{Float64,Float64}[]
        # Each edge from its first corner, the three points a quarter, half and three quarters along it.
        for (corner, x0, y0, dx, dy) in ((ll, minx, miny, maxx - minx, 0.0), (lr, maxx, miny, 0.0, maxy - miny),
                                         (ur, maxx, maxy, minx - maxx, 0.0), (ul, minx, maxy, 0.0, miny - maxy))
            push!(ring, corner)
            for frac in (0.25, 0.5, 0.75)
                push!(ring, lonlat(x0 + frac * dx, y0 + frac * dy))
            end
        end
        push!(ring, ll)
        bbox = (minimum(first, ring), minimum(last, ring), maximum(first, ring), maximum(last, ring))
        return Footprint4326(ring, bbox, center, (ul, ur, lr, ll))
    end::Footprint4326
end

# ---------------------------------------------------------------------------
# The STAC item
# ---------------------------------------------------------------------------

const STAC_EXTENSIONS = ["https://stac-extensions.github.io/projection/v2.0.0/schema.json",
                         "https://stac-extensions.github.io/alternate-assets/v1.2.0/schema.json",
                         "https://stac-extensions.github.io/version/v1.2.0/schema.json"]

# The sat extension, declared only by an item carrying `sat:orbit_state`: its schema requires at least one
# `sat:` field, and `sat:orbit_state` must be "ascending", "descending" or "geostationary".
const STAC_SAT_EXTENSION = "https://stac-extensions.github.io/sat/v1.1.0/schema.json"

# `s3://bucket/key` as the bucket's public HTTPS address.
function s3_https(uri::String)
    startswith(uri, "s3://") || throw(ArgumentError("\"" * uri * "\" is not an s3:// URI"))
    bucket, key = split(uri[6:end], '/'; limit = 2)
    return "https://" * bucket * ".s3.amazonaws.com/" * key
end

# The frames of a pair's scenes, as the reference reads them: a Landsat scene's path and row, a Sentinel-2
# product's tile, a NISAR granule's relative orbit, direction and frame, and a Sentinel-1 product's
# `frame_img1`/`frame_img2`. For Landsat, Sentinel-2 and NISAR the second frame is the first scene's: the
# reference splits the first scene's ID for both.
function scene_frames(s::ItsLiveOffsetProduction.ProductSummary, mission::String)
    startswith(mission, "S1") && return _na(s.frame_img1), _na(s.frame_img2)
    parts = split(s.id_img1, '_')
    f = startswith(mission, "L") ? String(parts[3]) : startswith(mission, "S2") ? String(parts[6]) :
        startswith(mission, "N") ? join(parts[6:8], "_") : "N/A"
    return f, f
end

# An attribute the product does not carry, as the reference writes it.
_na(s::String) = isempty(s) ? "N/A" : s

# A JSON string, as `json.dumps` escapes one: every value here is ASCII text without control characters
# but a quote or backslash.
_json_str(s::AbstractString) = "\"" * replace(s, "\\" => "\\\\", "\"" => "\\\"") * "\""
# A JSON number: the shortest form that reads back, as Python writes a float.
_json_num(x::Float64) = string(x)
_json_list(items) = "[" * join(items, ",") * "]"
_json_obj(pairs) = "{" * join([_json_str(k) * ":" * v for (k, v) in pairs], ",") * "}"
_lonlat(v::Tuple{Float64,Float64}) = _json_list((_json_num(v[1]), _json_num(v[2])))

"""
    item_json(s, fp, id, uri) -> String

The STAC item of product `id`, summarized by `s`, with footprint `fp`, published at `uri`
(`s3://bucket/prefix/<product>.nc`), whose browse image and thumbnail sit beside it: the document
`pystac.Item.to_dict` gives the item `cryoforge`'s `create_stac_item` builds, its keys in that order,
except for a product without a flight direction (an optical one): its item has neither `sat:orbit_state`
nor the sat extension, where `cryoforge` writes `"N/A"`, which the extension's schema rejects
(nasa-jpl/its_live_production#49).
"""
function item_json(s::ItsLiveOffsetProduction.ProductSummary, fp::Footprint4326, id::String, uri::String)
    mission = String(first(split(s.id_img1, '_')))
    f1, f2 = scene_frames(s, mission)
    mid = iso_utc(parse_pair_date(s.date_center))
    created = iso_utc(parse_created_date(s.date_created))
    updated = isempty(s.date_updated) ? created : iso_utc(parse_created_date(s.date_updated))
    orbit = s.flight_direction_img1
    props = Pair{String,String}[
        "mid_datetime" => _json_str(mid), "created" => _json_str(created), "updated" => _json_str(updated),
        "latitude" => _json_num(_py_round(fp.center[2], 4)), "longitude" => _json_num(_py_round(fp.center[1], 4)),
        "date_dt" => _json_num(round(s.date_dt)), "platform" => _json_str(mission),
        "scene_1_id" => _json_str(s.id_img1), "scene_2_id" => _json_str(s.id_img2),
        "scene_1_frame" => _json_str(f1), "scene_2_frame" => _json_str(f2)]
    isempty(orbit) || push!(props, "sat:orbit_state" => _json_str(orbit))
    append!(props, (
        "start_datetime" => _json_str(iso_utc(parse_pair_date(s.acquisition_date_img1))),
        "end_datetime" => _json_str(iso_utc(parse_pair_date(s.acquisition_date_img2))),
        "percent_valid_pixels" => string(round(Int, s.roi_valid_percentage)),
        "proj:code" => _json_str("EPSG:" * string(s.epsg)), "version" => _json_str(ITS_LIVE_DATA_VERSION),
        "datetime" => _json_str(mid)))
    extensions = isempty(orbit) ? STAC_EXTENSIONS : [STAC_EXTENSIONS; STAC_SAT_EXTENSION]
    beside(ext) = replace(uri, ".nc" => ext)
    image(ext, role) = _json_obj(("href" => _json_str(s3_https(beside(ext))), "type" => _json_str("image/png"),
                                  "roles" => _json_list((_json_str(role),))))
    data = _json_obj(("href" => _json_str(s3_https(uri)), "type" => _json_str("application/netcdf"),
                      "alternate" => _json_obj(("s3" => _json_obj(("href" => _json_str(uri),
                                                                   "alternate:name" => _json_str("S3"))),)),
                      "roles" => _json_list((_json_str("data"),))))
    geometry = _json_obj(("type" => _json_str("Polygon"), "coordinates" => _json_list((_json_list(_lonlat.(fp.ring)),))))
    return _json_obj((
        "type" => _json_str("Feature"), "stac_version" => _json_str("1.1.0"),
        "stac_extensions" => _json_list(_json_str.(extensions)), "id" => _json_str(id), "geometry" => geometry,
        "bbox" => _json_list(_json_num.(fp.bbox)), "properties" => _json_obj(props), "links" => "[]",
        "assets" => _json_obj(("data" => data, "overview" => image(".png", "overview"),
                               "thumbnail" => image("_thumb.png", "thumbnail"))),
        "collection" => _json_str("itslive-granules")))
end

# ---------------------------------------------------------------------------
# NSIDC's premet and spatial files
# ---------------------------------------------------------------------------

# Each mission's platform and sensor, as NSIDC names them, by the first field of its scene IDs.
const NSIDC_PLATFORMS = Dict(
    "LC09" => ("LANDSAT-9", "OLI"), "LO09" => ("LANDSAT-9", "OLI"), "LC08" => ("LANDSAT-8", "OLI"),
    "LO08" => ("LANDSAT-8", "OLI"), "LE07" => ("LANDSAT-7", "ETM+"), "LT05" => ("LANDSAT-5", "TM"),
    "LT04" => ("LANDSAT-4", "TM"), "S1A" => ("SENTINEL-1", "Sentinel-1A"), "S1B" => ("SENTINEL-1", "Sentinel-1B"),
    "S1C" => ("SENTINEL-1", "Sentinel-1C"), "S1D" => ("SENTINEL-1", "Sentinel-1D"),
    "S2A" => ("SENTINEL-2", "Sentinel-2A"), "S2B" => ("SENTINEL-2", "Sentinel-2B"),
    "S2C" => ("SENTINEL-2", "Sentinel-2C"), "S2D" => ("SENTINEL-2", "Sentinel-2D"), "NISAR" => ("NISAR", "L-SAR"))

"""
    premet(s, id) -> String

NSIDC's premet file of product `id`: its name, the data version, the two acquisitions' dates and times,
and each scene's platform, instrument and sensor, the scenes' missions read from `id` (`<scene 1>_X_<scene 2>_…`).
"""
function premet(s::ItsLiveOffsetProduction.ProductSummary, id::String)
    halves = split(id, "_X_")
    length(halves) >= 2 || throw(ArgumentError("\"" * id * "\" is not named <scene 1>_X_<scene 2>_…"))
    b, e = parse_pair_date(s.acquisition_date_img1), parse_pair_date(s.acquisition_date_img2)
    lines = ["FileName=" * id, "VersionID_local=" * ITS_LIVE_DATA_VERSION, "Begin_date=" * _ymd(b),
             "End_date=" * _ymd(e), "Begin_time=" * _hms_ms(b), "End_time=" * _hms_ms(e)]
    for h in halves[1:2]
        m = String(first(split(h, '_')))
        p = get(NSIDC_PLATFORMS, m, nothing)
        p === nothing && throw(ArgumentError("NSIDC records no platform for mission " * m * " of " * id))
        append!(lines, ["Container=AssociatedPlatformInstrumentSensor", "AssociatedPlatformShortName=" * p[1],
                        "AssociatedInstrumentShortName=" * p[2], "AssociatedSensorShortName=" * p[2]])
    end
    return join(lines, "\n") * "\n"
end

"""
    spatial(fp) -> String

NSIDC's spatial file: the footprint's four corners from the upper left, clockwise, one
`<longitude>\\t<latitude>` line each to two decimals.
"""
spatial(fp::Footprint4326) =
    join([string(_py_round(c[1], 2)) * "\t" * string(_py_round(c[2], 2)) for c in fp.corners], "\n") * "\n"

"""
    write_catalog_records(product, uri) -> (String, String, String)

The STAC item, premet and spatial file of the product netCDF at `product`, published at `uri`, written
beside it and named as `uri` names it: `<product>.stac.json`, `<product>.nc.premet` and
`<product>.nc.spatial`; their paths.
"""
function write_catalog_records(product::String, uri::String)
    s = ItsLiveOffsetProduction.read_summary(product)
    fp = footprint_4326(s)
    # Named as the product is where it is published, beside the local file.
    id = String(chopsuffix(basename(uri), ".nc"))
    stem = joinpath(dirname(product), id)
    item, pm, sp = stem * ".stac.json", stem * ".nc.premet", stem * ".nc.spatial"
    write(item, item_json(s, fp, id, uri))
    write(pm, premet(s, id))
    write(sp, spatial(fp))
    return item, pm, sp
end
