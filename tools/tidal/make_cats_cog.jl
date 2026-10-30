# CATS2008_v2023's height constants and flexure as one cloud-optimized GeoTIFF, the form the trimmed
# `autorift_itslive` reads them in (`app/src/tidecorrection.jl`).
#
#   julia --project=tools/tidal tools/tidal/make_cats_cog.jl <CATS2008_v2023.nc> <CATS2008_v2023_z.tif>
#
# Bands 1–20 are `hRe` and `hIm` of each constituent in the file's `constituent_order`, interleaved
# (`hRe_m2`, `hIm_m2`, `hRe_s2`, …), in metres, unscaled from the file's Int16 and unchanged in sign.
# Band 21 is `flexure` as a fraction (the file's percent / 100). Each band's description names it.
# The grid is the model's own: 2 km cells on the polar stereographic projection with origin latitude
# −90°, standard parallel −71° and central meridian −70°, written in metres where the file is in km.

using NCDatasets, ArchGDAL
const AG = ArchGDAL

function make_cats_cog(src::String, dst::String)
    NCDataset(src) do ds
        order = split(ds["constituents"].attrib["constituent_order"])
        x = Float64.(ds["x"][:]) .* 1000
        y = Float64.(ds["y"][:]) .* 1000
        dx, dy = x[2] - x[1], y[2] - y[1]
        all(≈(dx), diff(x)) && all(≈(dy), diff(y)) || error("the model grid is not regular")
        dy < 0 || error("the model's y axis ascends; this script writes it north-up as stored")
        nx, ny = length(x), length(y)
        gt = [x[1] - dx / 2, dx, 0.0, y[1] - dy / 2, 0.0, dy]
        proj = replace(ds["mapping"].attrib["spatial_proj4"], "+units=km" => "+units=m")
        names = String[]
        tmp = dst * ".tmp.tif"
        AG.create(tmp; driver = AG.getdriver("GTiff"), width = nx, height = ny, nbands = 2 * length(order) + 1,
                  dtype = Float32, options = ["TILED=YES", "BIGTIFF=YES"]) do out
            AG.setgeotransform!(out, gt)
            AG.setproj!(out, AG.toWKT(AG.importPROJ4(proj)))
            band = 0
            for (k, c) in enumerate(order), v in ("hRe", "hIm")
                band += 1
                # NCDatasets applies `scale_factor`; the variable is (x, y, constituent) in Julia order.
                a = Float32.(coalesce.(ds[v][:, :, k], NaN32))
                AG.write!(AG.getband(out, band), a)
                push!(names, v * "_" * c)
            end
            f = Float32.(coalesce.(ds["flexure"][:, :], 0)) ./ 100
            AG.write!(AG.getband(out, band + 1), f)
            push!(names, "flexure")
            for (i, n) in enumerate(names)
                AG.GDAL.gdalsetdescription(AG.getband(out, i).ptr, n)
            end
        end
        run(`gdal_translate -q -of COG -co COMPRESS=DEFLATE -co PREDICTOR=FLOATING_POINT -co BLOCKSIZE=256
             -co BIGTIFF=YES -co OVERVIEWS=NONE $tmp $dst`)
        rm(tmp)
        return names
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    println(join(make_cats_cog(ARGS[1], ARGS[2]), " "))
end
