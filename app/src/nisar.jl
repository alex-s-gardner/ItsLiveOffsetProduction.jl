# NISAR L2 (GSLC) pairs: the crop each geocoded granule is correlated over, its amplitude as the correlator
# reads it, and the identification packaging takes from it.

const GslcProduct = SLCDatasets.GeocodedProduct{SLCDatasets.NisarBackend}

# The only frequency a NISAR product is read at.
const NISAR_FREQUENCY = "A"

# `p` itself, once it is known to have opened at [`NISAR_FREQUENCY`](@ref): a product opens at the first
# frequency it lists, so a frequency-B-only mode is refused here rather than read at B.
function require_frequency_a(p)
    p.backend.frequency == NISAR_FREQUENCY || throw(ArgumentError(
        p.backend.path * " opens at frequency " * p.backend.frequency * ", not " * NISAR_FREQUENCY))
    return p
end

# The polarizations `b`'s product lists at its frequency.
function nisar_polarizations(b::SLCDatasets.NisarBackend)::Vector{String}
    group = b.product_type in SLCDatasets.GEOCODED_TYPES ? "grids" : "swaths"
    name = string(SLCDatasets.product_path(b), "/", group, "/frequency", b.frequency, "/listOfPolarizations")
    return SLCDatasets._with_file(f -> ChunkManifests.h5_read_strings(f, name), b)
end

"""
    nisar_polarization(reference, secondary) -> String

The polarization a NISAR pair is read at: `"HH"` if `reference` lists it, otherwise `"VV"`
(`nisar_isce3.py:554-570`). Throws if `reference` lists neither, or `secondary` lacks the one chosen.
hyp3-autorift reads a GSLC at HH only, so a VV-only GSLC pair is read here and refused there.
"""
function nisar_polarization(reference::SLCDatasets.NisarBackend, secondary::SLCDatasets.NisarBackend)
    listed = nisar_polarizations(reference)
    pol = "HH" in listed ? "HH" : "VV"
    pol in listed || throw(ArgumentError(
        reference.path * " lists " * join(listed, ", ") * " at frequency " * reference.frequency *
        ", neither HH nor VV"))
    pol in nisar_polarizations(secondary) || throw(ArgumentError(
        secondary.path * " has no " * pol * " at frequency " * secondary.frequency * ", which the reference has"))
    return pol
end

# `p`'s samples at polarization `pol`.
nisar_pixels(p, pol::String) = SLCDatasets.NisarRaster(p.backend, pol)

"""
    GslcScene

One GSLC of a pair cropped to the pair's overlap: the opened product, the polarization the pair is read at
([`nisar_polarization`](@ref)), the zero-based `(row, col)` where the crop starts in the granule's grid, the
crop's `(rows, cols)`, and the crop's footprint.
"""
struct GslcScene
    product::GslcProduct
    polarization::String
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
    src isa SLCDatasets.LocalFile && return open_gslc(src)
    src isa SLCDatasets.RemoteHTTP && return open_gslc(src)
    src isa SLCDatasets.RemoteS3 && return open_gslc(src)
    throw(ArgumentError("no source for " * uri))
end

open_gslc(src::SLCDatasets.AbstractSLCSource) = require_frequency_a(SLCDatasets.open_geocoded(src)::GslcProduct)

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
    ref, sec = c.reference.name == early.name ? (pe, pl) : (pl, pe)
    pol = nisar_polarization(ref.backend, sec.backend)
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
        GslcScene(p, pol, yoff, xoff, (ysize, xsize), cropped)
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
    px = nisar_pixels(g.product, g.polarization)
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
