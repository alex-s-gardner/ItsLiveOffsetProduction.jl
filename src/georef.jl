# CF grid-mapping attributes from a CRS, through GDAL's OSR C API.
#
# Projection parameters are read through OSR by their WKT names rather than any CF equivalent; the
# reasoning is at `_projparm`.

const GFT = GeoFormatTypes

# An OSR spatial reference, destroyed by `_with_srs` when the call using it returns.
function _with_srs(f::F, init::G) where {F,G}
    sr = GDAL.osrnewspatialreference(C_NULL)
    try
        init(sr)
        return f(sr)
    finally
        GDAL.osrdestroyspatialreference(sr)
    end
end

# A string OSR allocated, copied out and freed.
function _osr_string(export!::F, sr) where {F}
    p = Ref{Cstring}(C_NULL)
    export!(sr, p)
    p[] == C_NULL && throw(ErrorException("GDAL returned no text for a spatial reference export"))
    s = unsafe_string(p[])
    GDAL.vsifree(Ptr{Cvoid}(p[]))
    return s
end

# The WKT GDAL writes for a CRS: an EPSG code's definition from the EPSG database, exported as WKT1
# (`OSRExportToWkt`'s default), or a WKT given as it is.
_crs_wkt(code::Integer) =
    _with_srs(sr -> _osr_string(GDAL.osrexporttowkt, sr),
              sr -> _ogrcheck(GDAL.osrimportfromepsg(sr, Cint(code)), "EPSG:" * string(code)))
_crs_wkt(crs::GFT.EPSG) = _crs_wkt(GFT.val(crs))
_crs_wkt(crs::GFT.WellKnownText) = String(GFT.val(crs))
_crs_wkt(crs) = throw(ArgumentError(
    "cf_grid_mapping needs an integer EPSG code, a GeoFormatTypes.EPSG or a WellKnownText; got a " *
    string(nameof(typeof(crs)))))

function _ogrcheck(err, what::String)
    err == GDAL.OGRERR_NONE || throw(ArgumentError("GDAL cannot read the spatial reference " * what))
    return nothing
end

# One OSR projection parameter, read by its WKT name rather than any CF equivalent — which sidesteps
# a real problem: GDAL's own netCDF driver has changed which CF names it writes for a given
# projection across releases (measured directly, `Polar_Stereographic`'s standard parallel is
# `latitude_of_origin` in the reference container's GDAL and `standard_parallel` in a current one),
# while the WKT parameter name a projection is defined by is stable. `NaN` default: a parameter the
# projection does not use should be absent from the result, not a spurious zero.
_projparm(sr, name::String) = GDAL.osrgetprojparm(sr, name, NaN, C_NULL)

# The EPSG code of a WKT's projected CRS, or else its geographic one; `nothing` when it names no EPSG
# authority, which is a real case (a custom or stripped CRS).
function _wkt_epsg(sr)
    node = "PROJCS"
    auth = GDAL.osrgetauthorityname(sr, node)
    if auth === nothing
        node = "GEOGCS"
        auth = GDAL.osrgetauthorityname(sr, node)
    end
    auth == "EPSG" || return nothing
    code = GDAL.osrgetauthoritycode(sr, node)
    return code === nothing ? nothing : tryparse(Int64, code)
end

"""
    cf_grid_mapping(crs) -> Dict{String,AttrValue}

The CF grid-mapping attributes GDAL's netCDF writer records for `crs` — everything
[`ItsLiveGeoref`](@ref)'s `mapping_attrs` needs except `GeoTransform`, which [`write_product`](@ref)
always recomputes itself. `crs` is an integer EPSG code, a `GeoFormatTypes.EPSG`, or a
`GeoFormatTypes.WellKnownText`.

Confirmed key-for-key against a real captured ITS_LIVE product's `mapping` attributes for
`polar_stereographic` (both `EPSG:3413`/north and, by the same code path, `EPSG:3031`/south — the two
poles differ only in the sign of `latitude_of_origin`, which is what `latitude_of_projection_origin`'s
±90 follows). `transverse_mercator`'s parameter names follow the CF conventions document, which does
not carry the same cross-GDAL-version drift `polar_stereographic`'s did — not yet confirmed against a
captured UTM case.

Reads each projection parameter by its WKT name (`_projparm`) rather than through GDAL's generic
`OSRExportToCF1` CF exporter: that call's own key names for `polar_stereographic` turned out to be the
*current* GDAL's convention, not the reference container's older one (`standard_parallel` where the
reference writes `latitude_of_origin`, plus a `long_name` the reference does not write and a missing
`scale_factor_at_projection_origin`) — a real, measured mismatch, not a hypothetical one.

Supports only the two projections an ITS_LIVE product grid uses, polar stereographic and UTM
(transverse Mercator); any other throws.
"""
function cf_grid_mapping(crs)
    wkt = _crs_wkt(crs)
    return _with_srs(sr -> _cf_grid_mapping(sr, wkt),
                     sr -> _ogrcheck(GDAL.osrimportfromwkt(sr, [wkt]), "in a WKT"))
end

function _cf_grid_mapping(sr, wkt::String)
    proj = GDAL.osrgetattrvalue(sr, "PROJECTION", 0)
    # Entry by entry, sized as `Dict(pairs...)` sizes seven pairs, so the attributes iterate in the
    # order that constructor gives them; see `_img_pair_info_attrs`.
    d = sizehint!(Dict{String,AttrValue}(), 7)
    d["semi_major_axis"] = GDAL.osrgetsemimajor(sr, C_NULL)
    d["inverse_flattening"] = GDAL.osrgetinvflattening(sr, C_NULL)
    d["false_easting"] = _projparm(sr, "false_easting")
    d["false_northing"] = _projparm(sr, "false_northing")
    d["spatial_ref"] = wkt
    d["crs_wkt"] = wkt
    d["proj4text"] = _osr_string(GDAL.osrexporttoproj4, sr)
    # Absent rather than zero when the WKT names no EPSG authority, as GDAL's own netCDF writer
    # leaves it.
    epsg = _wkt_epsg(sr)
    epsg === nothing || epsg == 0 || (d["spatial_epsg"] = epsg)

    if proj == "Polar_Stereographic"
        lat0 = _projparm(sr, "latitude_of_origin")
        d["grid_mapping_name"] = "polar_stereographic"
        d["latitude_of_origin"] = lat0
        d["latitude_of_projection_origin"] = lat0 >= 0 ? 90.0 : -90.0
        d["straight_vertical_longitude_from_pole"] = _projparm(sr, "central_meridian")
        # Not a WKT parameter of this variant (`_projparm` returns `NaN`) — GDAL's netCDF writer adds
        # it anyway, as a fixed `Int64` `1`, and every ITS_LIVE polar grid uses this same
        # standard-parallel-defined variant, so it is a constant here rather than a read.
        d["scale_factor_at_projection_origin"] = Int64(1)
    elseif proj == "Transverse_Mercator"
        d["grid_mapping_name"] = "transverse_mercator"
        d["scale_factor_at_central_meridian"] = _projparm(sr, "scale_factor")
        d["longitude_of_central_meridian"] = _projparm(sr, "central_meridian")
        d["latitude_of_projection_origin"] = _projparm(sr, "latitude_of_origin")
    else
        throw(ArgumentError(
            "cf_grid_mapping supports only Polar_Stereographic and Transverse_Mercator, the two " *
            "projections an ITS_LIVE product grid uses; got " * (proj === nothing ? "none" : proj)))
    end
    return d
end
