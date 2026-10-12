# The `ItsLiveGeoref` a `Rasters.Raster` implies.

module ItsLiveOffsetProductionRastersExt

import ArchGDAL
using DimensionalData: dims, lookup
using ItsLiveOffsetProduction: ItsLiveOffsetProduction, ItsLiveGeoref, cf_grid_mapping
using Rasters: Rasters, AbstractRaster, crs

const GFT = Rasters.GeoFormatTypes

"""
    ItsLiveGeoref(r::Rasters.AbstractRaster; pixel_size_x, pixel_size_y) -> ItsLiveGeoref

An [`ItsLiveGeoref`](@ref) built from `r`'s own grid and CRS.

`x`/`y` come from `r`'s own coordinates; `mapping_attrs` from [`cf_grid_mapping`](@ref) on
`Rasters.crs(r)`, converted to WKT through ArchGDAL; `lonlat` from `ArchGDAL.reproject` against `r`'s
CRS, `order = :trad` so it returns `(lon, lat)` regardless of EPSG:4326's own compliant axis order.
`pixel_size_x`/`pixel_size_y` default to `r`'s own lookup step; pass them explicitly for a radar scene,
whose ground-range/azimuth pixel size is not its projected grid's spacing (see [`ItsLiveGeoref`](@ref)).
"""
function ItsLiveOffsetProduction.ItsLiveGeoref(r::AbstractRaster;
                                               pixel_size_x = abs(_grid_step(dims(r, Rasters.X))),
                                               pixel_size_y = abs(_grid_step(dims(r, Rasters.Y))))
    x = Float64.(collect(dims(r, Rasters.X)))
    y = Float64.(collect(dims(r, Rasters.Y)))
    georef_crs = crs(r)
    mapping_attrs = cf_grid_mapping(convert(GFT.WellKnownText, georef_crs))
    lonlat(px, py) = (only(ArchGDAL.reproject([[px, py]], georef_crs, GFT.EPSG(4326);
                                              order = :trad))...,)
    return ItsLiveGeoref(x, y, mapping_attrs, lonlat, Float64(pixel_size_x), Float64(pixel_size_y))
end

_grid_step(d) = step(lookup(d))

end # module
