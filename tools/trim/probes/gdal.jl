# Trim probe: the GDAL C API through GDAL.jl, covering every kind of call the production chain makes.
#
#   gdal <cog> <shapefile> <lon> <lat>
#
# Reads a 256² window of `cog`, finds the `shapefile` region containing (lon, lat), exports the
# region's CRS as WKT and reads its semi-major axis, and warps a small in-memory raster to EPSG:4326.

import GDAL

say(s::String) = (write(Core.stdout, s, "\n"); nothing)

# An exception's text without Base's printer, which a trimmed binary lacks. `ErrorException.msg` is
# typed `AbstractString`, so it is narrowed before use.
function message(err)
    m = err isa ErrorException ? err.msg : err isa GDAL.GDALError ? err.msg : nothing
    return m isa String ? m : string(typeof(err))
end

function region(shp::String, lon::Float64, lat::Float64)
    ds = GDAL.gdalopenex(shp, GDAL.GDAL_OF_VECTOR, C_NULL, C_NULL, C_NULL)
    ds == C_NULL && error("cannot open the shapefile")
    pt = GDAL.ogr_g_creategeometry(GDAL.wkbPoint)
    GDAL.ogr_g_setpoint_2d(pt, 0, lon, lat)
    layer = GDAL.gdaldatasetgetlayer(ds, 0)
    GDAL.ogr_l_resetreading(layer)
    name, epsg = "", 0
    while true
        f = GDAL.ogr_l_getnextfeature(layer)
        f == C_NULL && break
        if GDAL.ogr_g_contains(GDAL.ogr_f_getgeometryref(f), pt) != 0
            name = String(GDAL.ogr_f_getfieldasstring(f, GDAL.ogr_f_getfieldindex(f, "name")))
            epsg = Int(GDAL.ogr_f_getfieldasinteger(f, GDAL.ogr_f_getfieldindex(f, "epsg")))
        end
        GDAL.ogr_f_destroy(f)
        isempty(name) || break
    end
    GDAL.ogr_g_destroygeometry(pt)
    GDAL.gdalclose(ds)
    return name, epsg
end

function window(cog::String)
    ds = GDAL.gdalopen(cog, GDAL.GA_ReadOnly)
    ds == C_NULL && error("cannot open the COG")
    gt = zeros(Float64, 6)
    GDAL.gdalgetgeotransform(ds, gt)
    band = GDAL.gdalgetrasterband(ds, 1)
    a = Matrix{Float32}(undef, 256, 256)
    GC.@preserve a GDAL.gdalrasterio(band, GDAL.GF_Read, 32740, 53240, 256, 256, pointer(a), 256, 256,
                                     GDAL.GDT_Float32, 0, 0)
    GDAL.gdalclose(ds)
    return gt, a
end

function crs(epsg::Int)
    sr = GDAL.osrnewspatialreference(C_NULL)
    GDAL.osrimportfromepsg(sr, epsg)
    wkt = Ref{Cstring}(C_NULL)
    GDAL.osrexporttowkt(sr, wkt)
    s = unsafe_string(wkt[])
    GDAL.vsifree(Ptr{Cvoid}(wkt[]))
    err = Ref{GDAL.OGRErr}(GDAL.OGRERR_NONE)
    a = GDAL.osrgetsemimajor(sr, err)
    GDAL.osrdestroyspatialreference(sr)
    return s, a
end

function warp(epsg::Int, wkt::String, x0::Float64, y0::Float64)
    mem = GDAL.gdalgetdriverbyname("MEM")
    src = GDAL.gdalcreate(mem, "", 64, 64, 1, GDAL.GDT_Float32, C_NULL)
    gt = [x0, 120.0, 0.0, y0, 0.0, -120.0]
    GDAL.gdalsetgeotransform(src, gt)
    GDAL.gdalsetprojection(src, wkt)
    data = Float32[i + j for i in 1:64, j in 1:64]
    GC.@preserve data GDAL.gdalrasterio(GDAL.gdalgetrasterband(src, 1), GDAL.GF_Write, 0, 0, 64, 64,
                                        pointer(data), 64, 64, GDAL.GDT_Float32, 0, 0)
    argv = ["-of", "MEM", "-t_srs", "EPSG:4326", "-ts", "32", "32", "-r", "bilinear"]
    cstrs = [Base.unsafe_convert(Cstring, s) for s in argv]
    push!(cstrs, Cstring(C_NULL))
    opts = GC.@preserve argv GDAL.gdalwarpappoptionsnew(pointer(cstrs), C_NULL)
    srcs = [src]
    usage = Ref{Cint}(0)
    dst = GC.@preserve srcs GDAL.gdalwarp("", C_NULL, 1, pointer(srcs), opts, usage)
    GDAL.gdalwarpappoptionsfree(opts)
    dst == C_NULL && error("warp failed")
    out = Matrix{Float32}(undef, 32, 32)
    GC.@preserve out GDAL.gdalrasterio(GDAL.gdalgetrasterband(dst, 1), GDAL.GF_Read, 0, 0, 32, 32,
                                       pointer(out), 32, 32, GDAL.GDT_Float32, 0, 0)
    GDAL.gdalclose(dst)
    GDAL.gdalclose(src)
    return out
end

function (@main)(args::Vector{String})::Cint
    try
        GDAL.cplsetconfigoption("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
        gt, a = window(args[1])
        say(string("cog x0=", gt[1], " dx=", gt[2], " window sum=", sum(Float64, a)))
        name, epsg = region(args[2], parse(Float64, args[3]), parse(Float64, args[4]))
        say(string("region ", name, " epsg ", epsg))
        wkt, semimajor = crs(epsg)
        say(string("wkt chars=", length(wkt), " semi-major=", semimajor))
        out = warp(epsg, wkt, gt[1] + 32740 * gt[2], gt[4] + 53240 * gt[6])
        say(string("warp nonzero=", count(!iszero, out), " sum=", sum(Float64, out)))
        say(string("GDAL_DATA=", GDAL.GDAL_DATA[], " exists=", isdir(GDAL.GDAL_DATA[])))
        say(string("PROJ_LIB=", GDAL.PROJ_LIB[], " exists=", isdir(GDAL.PROJ_LIB[])))
        return Cint(0)
    catch err
        say(string("error: ", message(err)))
        return Cint(1)
    end
end
