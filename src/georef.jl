# CF grid-mapping attributes, and the `ItsLiveGeoref` a `Rasters.Raster` implies.
#
# Projection parameters are read through OSR by their WKT names rather than any CF equivalent; the
# reasoning is at `_projparm`.

using DimensionalData: dims, lookup
using Rasters: AbstractRaster, crs

const GFT = Rasters.GeoFormatTypes
const GDAL = ArchGDAL.GDAL

# ---------------------------------------------------------------------------
# ITS_LIVE grid-mapping attributes from a CRS
# ---------------------------------------------------------------------------

# `Integer` for a bare EPSG code, matching `ImagePairGeometry`'s own `_as_geoformat` — the same
# convenience, written again rather than shared, because that one is private to its package.
_as_geoformat(x::GFT.GeoFormat) = x
_as_geoformat(x::Integer) = GFT.EPSG(Int(x))
_as_geoformat(x) = throw(ArgumentError(
    "cf_grid_mapping needs a GeoFormatTypes.GeoFormat or an integer EPSG code, got a $(typeof(x))"))

# One OSR projection parameter, read by its WKT name rather than any CF equivalent — which sidesteps
# a real problem: GDAL's own netCDF driver has changed which CF names it writes for a given
# projection across releases (measured directly, `Polar_Stereographic`'s standard parallel is
# `latitude_of_origin` in the reference container's GDAL and `standard_parallel` in a current one),
# while the WKT parameter name a projection is defined by is stable. `NaN` default: a parameter the
# projection does not use should be absent from the result, not a spurious zero.
function _projparm(sr, name::AbstractString)
    ptr = Base.unsafe_convert(GDAL.OGRSpatialReferenceH, sr)
    return GDAL.osrgetprojparm(ptr, name, NaN, C_NULL)
end

"""
    cf_grid_mapping(crs) -> Dict{String,Any}

The CF grid-mapping attributes GDAL's netCDF writer records for `crs`.

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
"""
function cf_grid_mapping(crs)
    wkt = convert(GFT.WellKnownText, _as_geoformat(crs)).val
    sr = ArchGDAL.importWKT(wkt)
    proj = ArchGDAL.getattrvalue(sr, "PROJECTION", 0)

    d = Dict{String,Any}(
        "semi_major_axis" => GDAL.osrgetsemimajor(Base.unsafe_convert(GDAL.OGRSpatialReferenceH, sr),
                                                   C_NULL),
        "inverse_flattening" => GDAL.osrgetinvflattening(
            Base.unsafe_convert(GDAL.OGRSpatialReferenceH, sr), C_NULL),
        "false_easting" => _projparm(sr, "false_easting"),
        "false_northing" => _projparm(sr, "false_northing"),
        "spatial_ref" => wkt,
        "crs_wkt" => wkt,
        "proj4text" => ArchGDAL.toPROJ4(sr),
    )
    # `toEPSG` throws — rather than returning `0` or `nothing` — when the WKT carries no authority
    # code at all, a real case (a custom or stripped CRS), not a hypothetical one; confirmed by
    # calling it on one. `spatial_epsg` is then simply absent, matching what GDAL's own netCDF writer
    # does for such a CRS.
    epsg = try
        ArchGDAL.toEPSG(sr)
    catch
        nothing
    end
    epsg === nothing || epsg == 0 || (d["spatial_epsg"] = Int64(epsg))

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
            "projections an ITS_LIVE product grid uses; got $proj"))
    end
    return d
end

"""
    ItsLiveGeoref(r::Rasters.AbstractRaster; pixel_size_x, pixel_size_y) -> ItsLiveGeoref

An [`ItsLiveGeoref`](@ref) built from `r`'s own grid and CRS.

`x`/`y` come from `r`'s own coordinates; `mapping_attrs` from [`cf_grid_mapping`](@ref) on
`Rasters.crs(r)`; `lonlat` from `ArchGDAL.reproject` against `r`'s CRS, `order = :trad` so it returns
`(lon, lat)` regardless of EPSG:4326's own compliant axis order. `pixel_size_x`/`pixel_size_y` default
to `r`'s own lookup step; pass them explicitly for a radar scene, whose ground-range/azimuth pixel
size is not its projected grid's spacing (see [`ItsLiveGeoref`](@ref)).
"""
function ItsLiveGeoref(r::AbstractRaster;
                       pixel_size_x = abs(_grid_step(dims(r, Rasters.X))),
                       pixel_size_y = abs(_grid_step(dims(r, Rasters.Y))))
    x = Float64.(collect(dims(r, Rasters.X)))
    y = Float64.(collect(dims(r, Rasters.Y)))
    georef_crs = crs(r)
    mapping_attrs = cf_grid_mapping(georef_crs)
    lonlat(px, py) = (only(ArchGDAL.reproject([[px, py]], georef_crs, GFT.EPSG(4326);
                                              order = :trad))...,)
    return ItsLiveGeoref(x, y, mapping_attrs, lonlat, Float64(pixel_size_x),
                         Float64(pixel_size_y))
end

_grid_step(d) = step(lookup(d))
