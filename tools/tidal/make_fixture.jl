# The test fixture `app/test/tidecorrection.jl` checks a converted tide model against: for each point of a
# `pytmd_reference.py` table, the 4×4-cell window of the model COG around it, with pyTMD's heights at the
# table's times.
#
#   julia --project=tools/tidal tools/tidal/make_fixture.jl <model.tif> <pytmd.csv> <out.json>
#
#   julia --project=tools/tidal tools/tidal/make_fixture.jl ~/data/tides/cogs/Gr1kmTM_v1_z.tif \
#       tools/tidal/tide_pytmd_gr1km.csv app/test/data/tide_pytmd_gr1km.json
#
# The window's first cell centre is `origin` and its cells `spacing` apart, in metres; `bands` is the
# window `[x, y, band]` flattened in Julia's column-major order, `size` its shape.

using ArchGDAL, CSV, JSON3
const AG = ArchGDAL

function make_fixture(model::String, table::String, out::String)
    rows = collect(CSV.File(table; types = Dict(:time => String, :point => String)))
    names = unique(r.point for r in rows)
    points = AG.read(model) do d
        gt = AG.getgeotransform(d)
        nb = AG.nraster(d)
        map(names) do name
            g = filter(r -> r.point == name, rows)
            x, y = 1000 * g[1].x_km, 1000 * g[1].y_km
            # The cell whose centre is just below and left of the point, and one cell to spare each way.
            c0 = floor(Int, (x - gt[1]) / gt[2] - 0.5) - 1
            r0 = floor(Int, (y - gt[4]) / gt[6] - 0.5) - 1
            bands = Array{Float32,3}(undef, 4, 4, nb)
            for b in 1:nb
                # ArchGDAL's `read(band, rows, cols)` returns the window indexed `[x, y]`.
                bands[:, :, b] = AG.read(AG.getband(d, b), r0 + 1:r0 + 4, c0 + 1:c0 + 4)
            end
            (; point = name, lon = g[1].lon, lat = g[1].lat, flexure = g[1].flexure,
               origin = (gt[1] + (c0 + 0.5) * gt[2], gt[4] + (r0 + 0.5) * gt[6]), spacing = (gt[2], gt[6]),
               size = size(bands), bands = vec(bands), times = [r.time for r in g], height = [r.height for r in g])
        end
    end
    write(out, JSON3.write(points))
    return length(points)
end

if abspath(PROGRAM_FILE) == @__FILE__
    println(make_fixture(ARGS[1], ARGS[2], ARGS[3]), " points")
end
