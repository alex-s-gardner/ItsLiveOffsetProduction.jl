# Writing a tide model in the band layout `app/src/tidecorrection.jl` reads: shared by `make_cats_cog.jl` and
# `make_otis_cog.jl`.

using ArchGDAL
const AG = ArchGDAL

# The app's constituent order, read from `app/src/tides.jl`'s `TIDE_CONSTITUENTS` without loading the app.
const TIDE_CONSTITUENTS = let src = read(joinpath(@__DIR__, "..", "..", "app", "src", "tides.jl"), String)
    m = match(r"const TIDE_CONSTITUENTS = (\([^)]*\))", src)
    m === nothing && error("app/src/tides.jl defines no TIDE_CONSTITUENTS")
    Tuple(String.(eval(Meta.parse(m.captures[1]))))
end

"""
    write_tide_cog(dst, gt, proj4, bands)

`bands`, pairs of band name and `[x, y]` Float32 array in band order, written to `dst` as one
cloud-optimized GeoTIFF on geotransform `gt` (metres) and projection `proj4`: DEFLATE with the floating-point
predictor, 256-pixel tiles, nearest-neighbour overviews.
"""
function write_tide_cog(dst::String, gt::Vector{Float64}, proj4::String, bands::Vector{Pair{String,Matrix{Float32}}})
    nx, ny = size(last(first(bands)))
    tmp = dst * ".tmp.tif"
    AG.create(tmp; driver = AG.getdriver("GTiff"), width = nx, height = ny, nbands = length(bands),
              dtype = Float32, options = ["TILED=YES", "BIGTIFF=YES"]) do out
        AG.setgeotransform!(out, gt)
        AG.setproj!(out, AG.toWKT(AG.importPROJ4(proj4)))
        for (i, (name, a)) in enumerate(bands)
            AG.write!(AG.getband(out, i), a)
            AG.GDAL.gdalsetdescription(AG.getband(out, i).ptr, name)
        end
    end
    run(`gdal_translate -q -of COG -co COMPRESS=DEFLATE -co PREDICTOR=FLOATING_POINT -co BLOCKSIZE=256
         -co BIGTIFF=YES -co OVERVIEW_RESAMPLING=NEAREST $tmp $dst`)
    rm(tmp)
    return first.(bands)
end
