# NISAR L2 (GSLC) pairs: the crop each geocoded granule is correlated over, its amplitude as the correlator
# reads it, and the identification packaging takes from it.

const GslcProduct = SLCDatasets.GeocodedProduct{SLCDatasets.NisarBackend}

"""
    GslcScene

One GSLC of a pair cropped to the pair's overlap: the opened product, the zero-based `(row, col)` where the
crop starts in the granule's grid, the crop's `(rows, cols)`, and the crop's footprint.
"""
struct GslcScene
    product::GslcProduct
    row0::Int
    col0::Int
    dims::Tuple{Int,Int}
    cropped::Footprint
end

"""
    open_gslc(uri) -> GslcProduct

The GSLC at `uri` — a local path, an `s3://` URI, or an `https://` URL read with an Earthdata Login token
from `EARTHDATA_TOKEN` or `~/.netrc`.
"""
function open_gslc(uri::String)
    src = SLCDatasets.source_for(uri)
    src isa SLCDatasets.LocalFile && return SLCDatasets.open_geocoded(src)::GslcProduct
    src isa SLCDatasets.RemoteHTTP && return SLCDatasets.open_geocoded(src)::GslcProduct
    src isa SLCDatasets.RemoteS3 && return SLCDatasets.open_geocoded(src)::GslcProduct
    throw(ArgumentError("no source for " * uri))
end

# A geocoded product's own map grid as a footprint.
function gslc_footprint(g::GslcProduct)
    m = g.grid
    return ImageFootprint(m.origin, m.spacing, (m.size[2], m.size[1]), GFT.EPSG(m.epsg))
end

"""
    gslc_crop(c::Job, early, late) -> (early, late)

Each GSLC of `c` cropped to the pair's overlap, in acquisition order, as `crop_gslcs` (`nisar_isce3.py`)
does: the window is the envelope of the job's **reference** product's `boundingPolygon`, transformed from
EPSG:4326 (longitude first) into the products' projection; each product is cut to it by its own
geotransform, as `srcwin_for_intersection` does, with offsets and sizes rounded half to even as Python's
`round` does. The polygon's envelope is that of its transformed vertices, which is what OGR transforms.
"""
function gslc_crop(c::Job, early::SceneSource, late::SceneSource)
    pe, pl = open_gslc(early.band), open_gslc(late.band)
    ref = c.reference.name == early.name ? pe : pl
    epsg = ref.grid.epsg
    lon, lat = _polygon_vertices(ref.identification.bounding_polygon)
    xs, ys = copy(lon), copy(lat)
    with_transform(ct -> transform_points!(ct, xs, ys, zeros(length(xs))), 4326, epsg)
    xmin, xmax = minimum(xs), maximum(xs)
    ymin, ymax = minimum(ys), maximum(ys)
    crop(name::String, p::GslcProduct) = begin
        fp = gslc_footprint(p)
        p.grid.epsg == epsg || throw(ArgumentError(
            name * " is in EPSG:" * string(p.grid.epsg) * " and the reference product in EPSG:" * string(epsg) *
            "; the crop needs one projection"))
        px_w, px_h = fp.spacing
        xoff = round(Int, (xmin - fp.origin[1]) / px_w)
        yoff = round(Int, (ymax - fp.origin[2]) / px_h)
        xsize = round(Int, (xmax - xmin) / px_w)
        ysize = round(Int, (ymin - ymax) / px_h)
        (0 <= xoff && 0 <= yoff && xoff + xsize <= fp.size[1] && yoff + ysize <= fp.size[2]) || throw(ArgumentError(
            name * ": the overlap window does not fit the granule's grid"))
        cropped = ImageFootprint((fp.origin[1] + xoff * px_w, fp.origin[2] + yoff * px_h), fp.spacing,
                                 (xsize, ysize), fp.crs)
        GslcScene(p, yoff, xoff, (ysize, xsize), cropped)
    end
    return crop(early.name, pe), crop(late.name, pl)
end

# The `(lon, lat)` vertices of a WKT `POLYGON`'s outer ring, each vertex two or three numbers.
function _polygon_vertices(wkt::String)
    m = match(r"^\s*POLYGON\s*Z?\s*\(\(([^()]*)\)"i, wkt)
    m === nothing && throw(ArgumentError("the boundingPolygon is not a WKT POLYGON"))
    lon, lat = Float64[], Float64[]
    for v in eachsplit(something(m.captures[1]), ',')
        parts = split(strip(v))
        length(parts) >= 2 || throw(ArgumentError("a boundingPolygon vertex has fewer than two coordinates"))
        push!(lon, parse(Float64, parts[1]))
        push!(lat, parse(Float64, parts[2]))
    end
    length(lon) >= 3 || throw(ArgumentError("the boundingPolygon has fewer than three vertices"))
    return lon, lat
end

"""
    gslc_image(g, offset, want; dir, cache_io) -> TileFile

The GSLC's amplitude — the magnitude `convert_slc_to_uint8_amplitude` takes, without its byte rescale, which
the correlator's normalization makes moot — over its crop, cropped again at `offset` to the overlap
`want`, and filled into a tile cache. The crops are index shifts; only the tiles are read, each source
chunk once.
"""
function gslc_image(g::GslcScene, offset::NTuple{2,Int}, want::Tuple{Int,Int}; dir::String, cache_io::Symbol)
    ox, oy = offset
    (g.dims[1] >= oy + want[1] && g.dims[2] >= ox + want[2]) || throw(DimensionMismatch(
        "the overlap does not fit in the GSLC's crop; the crop and the geotransform disagree"))
    px = SLCDatasets.pixels(g.product)
    amp = SLCDatasets.amplitude(px)
    # The cache starts on a boundary of the granule's own chunks, which its 1024-pixel tiles then contain
    # whole, so each chunk is read once however the window falls; the image itself starts `origin` in.
    ch = map(length, first(SLCDatasets.DiskArrays.eachchunk(px)))
    r0, c0 = g.row0 + oy, g.col0 + ox
    ar, ac = r0 - r0 % ch[1], c0 - c0 % ch[2]
    (1024 % ch[1] == 0 && 1024 % ch[2] == 0) || throw(ArgumentError(
        "the GSLC's chunks do not divide the 1024-pixel tiles of its cache"))
    origin = (r0 - ar, c0 - ac)
    w = LazyWindow(amp, ar, ac, (want[1] + origin[1], want[2] + origin[2]))
    return fill_tiles!(TileCache(w; tile = 1024, dir, io = cache_io, origin); inflight = remote_tiles_in_flight())
end

# A NISAR `zeroDopplerStartTime` ("2025-10-28T23:52:01.000000000") as a `DateTime`, its fraction rounded
# to the millisecond.
function nisar_datetime(s::String)
    whole, frac = split(s, '.'; limit = 2)
    return parse_utc(whole) + Dates.Millisecond(round(Int, parse(Float64, "0." * frac) * 1000))
end

"""
    nisar_img_pair_info(early, late, id1, id2, roi_valid, lon, lat) -> ItsLiveOffsetProduction.ImagePairInfo

`IMG_INFO_DICT` for a NISAR pair from the two products' identification. `mission_img`/`satellite_img` are
the constants `"N"`/`1` the reference writes for every NISAR case.
"""
function nisar_img_pair_info(early::String, late::String, id1::SLCDatasets.Identification,
                             id2::SLCDatasets.Identification, roi_valid::Float64, lon::Float64, lat::Float64)
    # Sized for its eight entries, so the attributes iterate in the order an eight-pair `Dict` constructor
    # gives.
    extra = sizehint!(Dict{String,AttrValue}(), 8)
    extra["id_img1"] = early
    extra["id_img2"] = late
    extra["sensor_img1"] = "L"
    extra["sensor_img2"] = "L"
    extra["flight_direction_img1"] = lowercase(id1.pass_direction)
    extra["flight_direction_img2"] = lowercase(id2.pass_direction)
    extra["absolute_orbit_number_img1"] = Int64(id1.absolute_orbit)
    extra["absolute_orbit_number_img2"] = Int64(id2.absolute_orbit)
    return ItsLiveOffsetProduction.ImagePairInfo(nisar_datetime(id1.start_time), nisar_datetime(id2.start_time),
                                                 "N", "N", Int64(1), Int64(1), roi_valid,
                                                 round(lat; digits = 2), round(lon; digits = 2), extra)
end
