# The ITS_LIVE production parameters, resolved and read the way the reference resolves and reads them.
#
# Every golden product records `autoRIFT_parameter_file` as
# `autorift_parameters/v001/autorift_landice_0120m.shp`, a 122-feature polygon layer covering the UTM zones
# plus the two polar stereographic grids. One feature contains the pair's centroid, and its attributes name
# one cloud-optimized GeoTIFF per parameter and the EPSG code of the output grid.

const PARAMETER_SHAPEFILE =
    "/vsicurl/https://its-live-data.s3.amazonaws.com/autorift_parameters/v001/autorift_landice_0120m.shp"

# `GeometryInputs` field order, and the shapefile attribute naming each raster. The attribute names are
# truncated to ten characters because a DBF field name cannot be longer.
#
# `vx0`/`vy0` rather than `vx`/`vy`: geogrid's `window_offset` is the displacement a *prior* velocity
# implies. `dhdx`/`dhdy` rather than the smoothed `dhdxs`/`dhdys`. Both established by comparison: with
# these every integer band is bitwise identical to the reference container's output.
const PARAMETER_ATTRIBUTES = ("h", "dhdx", "dhdy", "vx0", "vy0", "vxSearchRan", "vySearchRan",
                              "xMinChipSiz", "yMinChipSiz", "xMaxChipSiz", "yMaxChipSiz", "StableSurfa")

"""
    Region

The parameter region covering a pair: its `name`, the output grid's `epsg`, the twelve parameter
rasters' `/vsicurl/` paths in [`PARAMETER_ATTRIBUTES`](@ref) order, and the `/vsicurl/` paths of its
floating-ice and ocean masks (1 on the class, the no-data 0 elsewhere), on the same grid.
"""
struct Region
    name::String
    epsg::Int
    paths::NTuple{12,String}
    floating::String
    ocean::String
end

dem_path(r::Region) = r.paths[1]

_vsicurl(url::String) = startswith(url, "/vsicurl/") ? url : "/vsicurl/" * url

"""
    parameter_region(lon, lat) -> Region

The region whose polygon contains `(lon, lat)`. The UTM-zone and polar polygons do not overlap in the
layer, so there is no priority rule to reproduce. A pair outside the coverage is an error.
"""
function parameter_region(lon::Float64, lat::Float64)
    r = containing_feature(PARAMETER_SHAPEFILE, lon, lat) do feat
        paths = ntuple(i -> _vsicurl(field_string(feat, PARAMETER_ATTRIBUTES[i])), Val(12))
        Region(field_string(feat, "name"), field_int(feat, "epsg"), paths, _vsicurl(field_string(feat, "floatingice")),
               _vsicurl(field_string(feat, "ocean")))
    end
    r === nothing && throw(ArgumentError("no ITS_LIVE parameter region contains (" * string(lon) * ", " *
                                         string(lat) * "); the pair is outside the parameters' coverage"))
    return r
end

const Grid = MapGrid{Float64,GFT.EPSG{1}}

"""
    parameter_grid(r) -> MapGrid

The output grid, from the region's DEM, which is what `geogridOptical` takes its geotransform and size from
(`geogridOptical.cpp:337-339`): every other parameter raster is read at the DEM's window.
"""
function parameter_grid(r::Region)
    return with_dataset(dem_path(r)) do d
        nx, ny = rastersize(d)
        MapGrid(geotransform(d), (nx, ny), r.epsg)
    end
end

"""
    dem_nodata(r) -> Float64

The no-data value the reference applies to five rasters, read from DEM band 1 as it does
(`geogridOptical.cpp:337-339`); an unset value is `0.0`, because the reference passes `pbSuccess = NULL`
and so cannot tell one from zero.
"""
function dem_nodata(r::Region)
    nd = with_dataset(band_nodata, dem_path(r))
    return nd === nothing ? 0.0 : nd
end

"""
    parameter_window(path, window) -> Matrix{Float64}

One parameter raster over `window`, as `Float64` — the type the reference reads every one of them as
(`GDT_Float64` in each `RasterIO` call, `geogridOptical.cpp:700-850`). `window` is one-based; GDAL's
offsets are zero-based.
"""
function parameter_window(path::String, window::GridWindow)
    xs, ys = window.indices
    return with_dataset(path) do d
        read_window(d, Float64, first(xs) - 1, first(ys) - 1, length(xs), length(ys))
    end
end

"""
    mask_window(path, grid, window, margin) -> (Matrix{UInt8}, Matrix{Bool})

A class raster over `window` of `grid` grown by `margin` pixels on every side, indexed `[x, y]`, and where
the raster covers it. The raster may cover any extent on `grid`'s lattice — same pixel size, edges a whole
number of pixels from the grid's — and is located by its own geotransform; a pixel beyond it reads as 0
and is not covered.
"""
function mask_window(path::String, grid::Grid, window::GridWindow, margin::Int)
    xs, ys = window.indices
    g = grid.geotransform
    return with_dataset(path) do d
        e = projected_epsg(projection_wkt(d))
        e == grid.crs.val[1] || throw(ArgumentError(
            path * " is on EPSG:" * string(e) * "; the product grid is on EPSG:" * string(grid.crs.val[1])))
        t = geotransform(d)
        (t[2] == g[2] && t[6] == g[6]) || throw(ArgumentError(
            path * " has pixels of " * string(t[2]) * " m; the product grid's are " * string(g[2]) * " m"))
        ox, oy = (t[1] - g[1]) / g[2], (t[4] - g[4]) / g[6]
        (isinteger(ox) && isinteger(oy)) || throw(ArgumentError(path * " is not on the product grid's lattice"))
        # Columns and rows of the grown window, zero-based on `grid`, then on the raster.
        gx0, gy0 = first(xs) - 1 - margin, first(ys) - 1 - margin
        nx, ny = length(xs) + 2margin, length(ys) + 2margin
        rx0, ry0 = gx0 - Int(ox), gy0 - Int(oy)
        w, h = rastersize(d)
        out = zeros(UInt8, nx, ny)
        inside = falses(nx, ny)
        cx0, cy0 = max(rx0, 0), max(ry0, 0)
        cx1, cy1 = min(rx0 + nx, w) - 1, min(ry0 + ny, h) - 1
        if cx1 >= cx0 && cy1 >= cy0
            cols, rows = cx0 - rx0 + 1:cx1 - rx0 + 1, cy0 - ry0 + 1:cy1 - ry0 + 1
            out[cols, rows] = read_window(d, UInt8, cx0, cy0, cx1 - cx0 + 1, cy1 - cy0 + 1)
            inside[cols, rows] .= true
        end
        (out, Matrix{Bool}(inside))
    end
end

const Inputs = GeometryInputs{Matrix{Float64},Matrix{Float64},Matrix{Float64},Matrix{Float64},
                              Matrix{Float64},Matrix{Float64},Matrix{Float64}}

"""
    geometry_inputs(r, window) -> GeometryInputs

The twelve parameter rasters over `window`, one task each: every one is a `/vsicurl/` read, so run one
after the other they spend most of their time waiting. Each task opens its own dataset, which is how GDAL is
thread-safe.
"""
function geometry_inputs(r::Region, window::GridWindow)
    gdal_network_setup()
    tasks = [Threads.@spawn(parameter_window(r.paths[i], window)) for i in 1:12]
    b = Matrix{Float64}[fetch(t)::Matrix{Float64} for t in tasks]
    return GeometryInputs(; dem = b[1], dhdx = b[2], dhdy = b[3], vx = b[4], vy = b[5], srx = b[6],
                          sry = b[7], csminx = b[8], csminy = b[9], csmaxx = b[10], csmaxy = b[11],
                          ssm = b[12])
end
