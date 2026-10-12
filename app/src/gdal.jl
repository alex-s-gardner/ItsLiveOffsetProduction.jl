# The GDAL C API calls the chain makes, through GDAL.jl: datasets opened and closed explicitly, windows read
# into arrays of a named type, OSR identification and point transforms, OGR containment, and GDALWarp.
#
# Every handle is released by the call that took it, with no finalizer, and every value returned has a
# type that does not depend on the file — which is what lets the chain be compiled with `juliac --trim`.

import GDAL

const GDAL_CONFIGURED = Ref(false)

"""
    gdal_network_setup()

Configure GDAL for the `/vsicurl/` and `/vsis3/` reads the chain makes, once per process.

`CURL_CA_BUNDLE` and its two companions go in the process environment, where the curl handle reads them:
GDAL's bundled curl does not find the system trust store, and the failure otherwise surfaces as
`SSL certificate problem` or a NULL dataset. `GDAL_DISABLE_READDIR_ON_OPEN` keeps GDAL from listing a
bucket prefix before every open. The read chunk and cache sizes turn one band from thousands of round
trips into a few hundred. `AWS_REQUEST_PAYER` because the Landsat bucket is requester-pays.
"""
function gdal_network_setup()
    GDAL_CONFIGURED[] && return nothing
    ca = NetworkOptions.ca_roots_path()
    for var in ("CURL_CA_BUNDLE", "SSL_CERT_FILE", "GDAL_HTTP_CAINFO")
        haskey(ENV, var) || (ENV[var] = ca)
    end
    for (key, value) in (("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR"),
                         ("CPL_VSIL_CURL_CHUNK_SIZE", "1048576"),
                         ("CPL_VSIL_CURL_CACHE_SIZE", "268435456"),
                         ("GDAL_HTTP_MULTIPLEX", "YES"),
                         ("GDAL_HTTP_VERSION", "2"),
                         ("VSI_CACHE", "TRUE"),
                         ("AWS_REQUEST_PAYER", "requester"))
        GDAL.cplsetconfigoption(key, value)
    end
    GDAL_CONFIGURED[] = true
    return nothing
end

# GDAL's last error message, for an error naming what failed.
function _gdal_error(what::String)
    m = GDAL.cplgetlasterrormsg()
    return ErrorException(what * (m === nothing || isempty(m) ? "" : ": " * m))
end

"""
    Dataset

An open GDAL raster dataset. Close it with `close`, or use [`with_dataset`](@ref).
"""
struct Dataset
    h::GDAL.GDALDatasetH
end

# Opening a `/vsicurl/` or `/vsis3/` raster and reading a window of one wait on the network inside libgdal,
# so those calls are made GC-safe: a collection another thread starts meanwhile does not wait for them, and
# other tasks — a download among them — keep running.
function open_dataset(path::String)
    gdal_network_setup()
    h = @ccall gc_safe = true GDAL.libgdal.GDALOpenEx(path::Cstring, (GDAL.GDAL_OF_RASTER | GDAL.GDAL_OF_READONLY)::Cuint,
                                                      C_NULL::Ptr{Cstring}, C_NULL::Ptr{Cstring},
                                                      C_NULL::Ptr{Cstring})::GDAL.GDALDatasetH
    h == C_NULL && throw(_gdal_error("GDAL cannot open " * path))
    return Dataset(h)
end

Base.close(d::Dataset) = (GDAL.gdalclose(d.h); nothing)

"""
    with_dataset(f, path)

`f(dataset)` for the raster at `path`, closed when `f` returns or throws.
"""
function with_dataset(f::F, path::String) where {F}
    d = open_dataset(path)
    try
        return f(d)
    finally
        close(d)
    end
end

function geotransform(d::Dataset)
    gt = zeros(Float64, 6)
    GDAL.gdalgetgeotransform(d.h, gt)
    return (gt[1], gt[2], gt[3], gt[4], gt[5], gt[6])
end

# Width and height, in pixels.
rastersize(d::Dataset) = (Int(GDAL.gdalgetrasterxsize(d.h)), Int(GDAL.gdalgetrasterysize(d.h)))

projection_wkt(d::Dataset) = String(GDAL.gdalgetprojectionref(d.h))

# The no-data value of band `b`, or `nothing` where none is set.
function band_nodata(d::Dataset, b::Int = 1)
    ok = Ref{Cint}(0)
    v = GDAL.gdalgetrasternodatavalue(GDAL.gdalgetrasterband(d.h, b), ok)
    return ok[] != 0 ? Float64(v) : nothing
end

_gdt(::Type{UInt8}) = GDAL.GDT_Byte
_gdt(::Type{UInt16}) = GDAL.GDT_UInt16
_gdt(::Type{Int16}) = GDAL.GDT_Int16
_gdt(::Type{Float32}) = GDAL.GDT_Float32
_gdt(::Type{Float64}) = GDAL.GDT_Float64

"""
    read_window(d, T, xoff, yoff, nx, ny; band = 1) -> Matrix{T}

`nx` by `ny` pixels of `band` from the zero-based offset `(xoff, yoff)`, as `T`, indexed `[x, y]` — the
layout GDAL's rows give a column-major array.
"""
function read_window(d::Dataset, ::Type{T}, xoff::Int, yoff::Int, nx::Int, ny::Int; band::Int = 1) where {T}
    out = Matrix{T}(undef, nx, ny)
    b = GDAL.gdalgetrasterband(d.h, band)
    err = GC.@preserve out @ccall gc_safe = true GDAL.libgdal.GDALRasterIO(
        b::GDAL.GDALRasterBandH, GDAL.GF_Read::GDAL.GDALRWFlag, xoff::Cint, yoff::Cint, nx::Cint, ny::Cint,
        pointer(out)::Ptr{Cvoid}, nx::Cint, ny::Cint, _gdt(T)::GDAL.GDALDataType, 0::Cint, 0::Cint)::GDAL.CPLErr
    err == GDAL.CE_None || throw(_gdal_error("GDAL could not read a window"))
    return out
end

# ---------------------------------------------------------------------------
# Spatial references
# ---------------------------------------------------------------------------

# A spatial reference from an EPSG code, in traditional GIS axis order (x/longitude first) when `trad`,
# as `ArchGDAL.importEPSG(code; order = :trad)` builds it.
function _epsg_srs(code::Int; trad::Bool)
    sr = GDAL.osrnewspatialreference(C_NULL)
    GDAL.osrimportfromepsg(sr, Cint(code)) == GDAL.OGRERR_NONE || (GDAL.osrdestroyspatialreference(sr);
        throw(ArgumentError("GDAL does not know EPSG:" * string(code))))
    trad && GDAL.osrsetaxismappingstrategy(sr, GDAL.OAMS_TRADITIONAL_GIS_ORDER)
    return sr
end

"""
    epsg_wkt(code) -> String

EPSG `code`'s definition as WKT1, which a MEM dataset is given as its projection.
"""
function epsg_wkt(code::Int)
    sr = _epsg_srs(code; trad = false)
    try
        p = Ref{Cstring}(C_NULL)
        GDAL.osrexporttowkt(sr, p)
        s = unsafe_string(p[])
        GDAL.vsifree(Ptr{Cvoid}(p[]))
        return s
    finally
        GDAL.osrdestroyspatialreference(sr)
    end
end

"""
    projected_epsg(wkt) -> Int

The EPSG code of a dataset's *projected* CRS, by the reference's own procedure.

`GeogridOptical.getProjectionSystem` (`GeogridOptical.py:93-123`) imports the WKT, calls
`AutoIdentifyEPSG`, requires the result to be projected, and reads the authority code off `PROJCS` —
falling back to a database match if that is empty. A Landsat scene over Antarctica carries a custom
`PROJCRS["PS         WGS84"]` with no authority code of its own; auto-identification resolves it to 3031,
where reading the authority of its base CRS would give 4326.
"""
function projected_epsg(wkt::String)
    sr = GDAL.osrnewspatialreference(C_NULL)
    try
        GDAL.osrimportfromwkt(sr, [wkt]) == GDAL.OGRERR_NONE ||
            throw(ArgumentError("GDAL cannot read a scene's coordinate system"))
        GDAL.osrisprojected(sr) == 1 || throw(ArgumentError(
            "the scene's coordinate system is not projected; geogrid refuses a geographic or local one " *
            "(GeogridOptical.py:110-115) because its pixel arithmetic is in metres"))
        GDAL.osrautoidentifyepsg(sr)
        code = GDAL.osrgetauthoritycode(sr, "PROJCS")
        code === nothing || isempty(code) || return parse(Int, code)
        # The reference's last resort is `gdalsrsinfo -o epsg`, a match against the EPSG database.
        n = Ref{Cint}(0)
        conf = Ref{Ptr{Cint}}(C_NULL)
        matches = GDAL.osrfindmatches(sr, C_NULL, n, conf)
        n[] > 0 || throw(ArgumentError("could not identify an EPSG code for the scene's projected CRS; " *
                                       "the reference raises here too (GeogridOptical.py:123)"))
        best = unsafe_load(Ptr{GDAL.OGRSpatialReferenceH}(matches), 1)
        c = GDAL.osrgetauthoritycode(best, C_NULL)
        say("scene CRS identified by database match rather than by authority code: EPSG:" *
            (c === nothing ? "?" : c) * ", confidence " * string(unsafe_load(conf[], 1)))
        c === nothing && throw(ArgumentError("the database match for the scene's CRS has no EPSG code"))
        return parse(Int, c)
    finally
        GDAL.osrdestroyspatialreference(sr)
    end
end

"""
    with_transform(f, src_epsg, dst_epsg)

`f(ct)` for the OSR coordinate transformation from `src_epsg` to `dst_epsg`, both in traditional GIS axis
order, as `ArchGDAL.crs2transform(EPSG(src), EPSG(dst); order = :trad)` builds it.
"""
function with_transform(f::F, src::Int, dst::Int) where {F}
    s = _epsg_srs(src; trad = true)
    t = _epsg_srs(dst; trad = true)
    ct = GDAL.octnewcoordinatetransformation(s, t)
    try
        ct == C_NULL && throw(_gdal_error("no transformation from EPSG:" * string(src) * " to EPSG:" *
                                          string(dst)))
        return f(ct)
    finally
        ct == C_NULL || GDAL.octdestroycoordinatetransformation(ct)
        GDAL.osrdestroyspatialreference(t)
        GDAL.osrdestroyspatialreference(s)
    end
end

# One point through `ct` as an OGR point geometry, which is how `ArchGDAL.transform!` moves a point
# created with `ArchGDAL.createpoint`: two coordinates, or three when `z` is given.
function transform_point(ct, x::Float64, y::Float64)
    g = GDAL.ogr_g_creategeometry(GDAL.wkbPoint)
    try
        GDAL.ogr_g_addpoint_2d(g, x, y)
        GDAL.ogr_g_transform(g, ct) == GDAL.OGRERR_NONE || throw(_gdal_error("a point did not transform"))
        return (GDAL.ogr_g_getx(g, 0), GDAL.ogr_g_gety(g, 0))
    finally
        GDAL.ogr_g_destroygeometry(g)
    end
end

function transform_point(ct, x::Float64, y::Float64, z::Float64)
    g = GDAL.ogr_g_creategeometry(GDAL.wkbPoint)
    try
        GDAL.ogr_g_addpoint(g, x, y, z)
        GDAL.ogr_g_transform(g, ct) == GDAL.OGRERR_NONE || throw(_gdal_error("a point did not transform"))
        return (GDAL.ogr_g_getx(g, 0), GDAL.ogr_g_gety(g, 0))
    finally
        GDAL.ogr_g_destroygeometry(g)
    end
end

# Many points through `ct` in place, as `ArchGDAL.transform!(xs, ys, zs, ct)` does; `false` if any failed.
function transform_points!(ct, xs::Vector{Float64}, ys::Vector{Float64}, zs::Vector{Float64})
    length(xs) == length(ys) == length(zs) || throw(DimensionMismatch("coordinate vectors differ in length"))
    return GDAL.octtransform(ct, length(xs), xs, ys, zs) != 0
end

# ---------------------------------------------------------------------------
# MEM datasets and the warp
# ---------------------------------------------------------------------------

const VSIMEM_COUNTER = Threads.Atomic{Int}(0)

# A `/vsimem/` path no other call in the process uses: warps run concurrently, one per tile.
vsimem_name(stem::String) = "/vsimem/" * stem * "_" * string(Threads.atomic_add!(VSIMEM_COUNTER, 1))

# A C argument vector over `args`, NULL-terminated; valid while `args` is preserved.
_argv(args::Vector{String}) = push!(Cstring[Base.unsafe_convert(Cstring, a) for a in args], Cstring(C_NULL))

"""
    warp(src, flags, dest) -> Dataset

`gdalwarp` of the one dataset `src` with command-line `flags`, written to `dest` — a `/vsimem/` path for an
`-of VRT` or `-of MEM` result. The caller closes the result.
"""
function warp(src::Dataset, flags::Vector{String}, dest::String)
    argv = _argv(flags)
    opts = GC.@preserve flags GDAL.gdalwarpappoptionsnew(pointer(argv), C_NULL)
    opts == C_NULL && throw(_gdal_error("GDAL rejected the warp options"))
    srcs = [src.h]
    usage = Ref{Cint}(0)
    # GC-safe for the same reason as a network read: a warp runs for a while inside libgdal, and the other
    # threads warping their own tiles meanwhile must be able to collect.
    out = try
        GC.@preserve srcs usage @ccall gc_safe = true GDAL.libgdal.GDALWarp(
            dest::Cstring, C_NULL::GDAL.GDALDatasetH, 1::Cint, pointer(srcs)::Ptr{GDAL.GDALDatasetH},
            opts::Ptr{GDAL.GDALWarpAppOptions}, usage::Ptr{Cint})::GDAL.GDALDatasetH
    finally
        GDAL.gdalwarpappoptionsfree(opts)
    end
    out == C_NULL && throw(_gdal_error("GDAL could not warp to " * dest))
    return Dataset(out)
end

"""
    mem_dataset(data, geotransform, wkt; nodata = nothing) -> Dataset

A one-band MEM dataset holding `data`, indexed `[x, y]`, on `geotransform` in the CRS `wkt` names, with
`nodata` set on the band unless it is `nothing`. The caller closes it.
"""
function mem_dataset(data::Matrix{T}, gt::NTuple{6,Float64}, wkt::String;
                     nodata::Union{Nothing,Float64} = nothing) where {T}
    gdal_network_setup()
    nx, ny = size(data)
    h = GDAL.gdalcreate(GDAL.gdalgetdriverbyname("MEM"), "", nx, ny, 1, _gdt(T), C_NULL)
    h == C_NULL && throw(_gdal_error("GDAL could not create a MEM dataset"))
    d = Dataset(h)
    GDAL.gdalsetgeotransform(h, Float64[gt[1], gt[2], gt[3], gt[4], gt[5], gt[6]])
    GDAL.gdalsetprojection(h, wkt)
    band = GDAL.gdalgetrasterband(h, 1)
    err = GC.@preserve data GDAL.gdalrasterio(band, GDAL.GF_Write, 0, 0, nx, ny, pointer(data), nx, ny,
                                              _gdt(T), 0, 0)
    err == GDAL.CE_None || (close(d); throw(_gdal_error("GDAL could not write a MEM dataset")))
    nodata === nothing || GDAL.gdalsetrasternodatavalue(band, nodata)
    return d
end

# ---------------------------------------------------------------------------
# Vector layers
# ---------------------------------------------------------------------------

"""
    containing_feature(f, path, lon, lat)

`f(feature)` for the first feature of `path`'s first layer whose geometry contains `(lon, lat)`, in feature
order, or `nothing` if none does. The feature is destroyed when `f` returns.
"""
function containing_feature(f::F, path::String, lon::Float64, lat::Float64) where {F}
    gdal_network_setup()
    ds = @ccall gc_safe = true GDAL.libgdal.GDALOpenEx(path::Cstring, (GDAL.GDAL_OF_VECTOR | GDAL.GDAL_OF_READONLY)::Cuint,
                                                       C_NULL::Ptr{Cstring}, C_NULL::Ptr{Cstring},
                                                       C_NULL::Ptr{Cstring})::GDAL.GDALDatasetH
    ds == C_NULL && throw(_gdal_error("GDAL cannot open " * path))
    pt = GDAL.ogr_g_creategeometry(GDAL.wkbPoint)
    try
        GDAL.ogr_g_addpoint_2d(pt, lon, lat)
        layer = GDAL.gdaldatasetgetlayer(ds, 0)
        n = GDAL.ogr_l_getfeaturecount(layer, 1)
        for k in 0:(n - 1)
            feat = GDAL.ogr_l_getfeature(layer, k)
            feat == C_NULL && continue
            try
                GDAL.ogr_g_contains(GDAL.ogr_f_getgeometryref(feat), pt) != 0 && return f(feat)
            finally
                GDAL.ogr_f_destroy(feat)
            end
        end
        return nothing
    finally
        GDAL.ogr_g_destroygeometry(pt)
        GDAL.gdalclose(ds)
    end
end

function field_string(feat, name::String)
    i = GDAL.ogr_f_getfieldindex(feat, name)
    i < 0 && throw(ArgumentError("the feature has no field " * name))
    s = GDAL.ogr_f_getfieldasstring(feat, i)
    return s === nothing ? "" : String(s)
end

function field_int(feat, name::String)
    i = GDAL.ogr_f_getfieldindex(feat, name)
    i < 0 && throw(ArgumentError("the feature has no field " * name))
    return Int(GDAL.ogr_f_getfieldasinteger64(feat, i))
end
