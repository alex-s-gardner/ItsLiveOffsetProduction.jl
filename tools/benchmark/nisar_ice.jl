# NISAR candidates (nisar_pairs.py) over land ice: a pair is kept when any point of a 5x5 grid spanning
# ±1° of latitude (±110 km, about a frame) around its footprint centre falls where the ITS_LIVE parameter
# region's ROI raster is set — the processing mask `roi_valid_percentage` is computed from. The ROI is read
# at its overview 3 (~2 km pixels). Writes nisar_<level>_ice.tsv: the candidate's columns, the region, and
# the fraction of grid points on the ROI.
#   julia nisar_ice.jl SHAPEFILE OUTDIR
using ArchGDAL
const AG = ArchGDAL
shp, out = ARGS
const OVERVIEW = 3

struct RegionROI
    name::String
    geom::AG.IGeometry
    epsg::Int
    roi::String
end

regions = AG.read(shp) do ds
    [RegionROI(AG.getfield(f, "name"), AG.clone(AG.getgeom(f)), AG.getfield(f, "epsg"),
               "/vsicurl/" * AG.getfield(f, "ROI")) for f in AG.getlayer(ds, 0)]
end

# A region's ROI at the overview, with the geotransform of that overview, read on first use.
const CACHE = Dict{String,Tuple{Matrix{UInt8},Vector{Float64}}}()
function roi(r::RegionROI)
    get!(CACHE, r.roi) do
        AG.read(r.roi) do ds
            gt = AG.getgeotransform(ds)
            o = AG.getoverview(AG.getband(ds, 1), OVERVIEW)
            a = AG.read(o)
            sx, sy = AG.width(ds) / AG.width(o), AG.height(ds) / AG.height(o)
            (a, [gt[1], gt[2] * sx, gt[3], gt[4], gt[5], gt[6] * sy])
        end
    end
end

const TRANSFORMS = Dict{Int,Any}()
lonlat_to(epsg) = get!(TRANSFORMS, epsg) do
    AG.unsafe_createcoordtrans(AG.importEPSG(4326; order = :trad), AG.importEPSG(epsg; order = :trad))
end

function region_of(lon, lat)
    p = AG.createpoint(lon, lat)
    i = findfirst(r -> AG.contains(r.geom, p), regions)
    return i === nothing ? nothing : regions[i]
end

function ice_fraction(r::RegionROI, lon0, lat0)
    a, gt = roi(r)
    ct = lonlat_to(r.epsg)
    hits = 0
    for dlat in range(-1.0, 1.0; length = 5), dlon in range(-1.0, 1.0; length = 5)
        lat = clamp(lat0 + dlat, -89.9, 89.9)
        lon = lon0 + dlon / max(cosd(lat), 0.05)
        pt = AG.createpoint(lon, lat)
        AG.transform!(pt, ct)
        x, y = AG.getx(pt, 0), AG.gety(pt, 0)
        col = floor(Int, (x - gt[1]) / gt[2]) + 1
        row = floor(Int, (y - gt[4]) / gt[6]) + 1
        (1 <= col <= size(a, 1) && 1 <= row <= size(a, 2)) || continue
        hits += a[col, row] > 0
    end
    return hits / 25
end

frac = Dict{Tuple{Float64,Float64},Tuple{String,Float64}}()
for level in ("L1", "L2")
    rows = [split(l, '\t') for l in eachline(joinpath(out, "nisar_" * level * "_candidates.tsv"))]
    n = 0
    open(joinpath(out, "nisar_" * level * "_ice.tsv"), "w") do io
        for r in rows
            lon, lat = parse(Float64, r[5]), parse(Float64, r[6])
            name, f = get!(frac, (lon, lat)) do
                reg = region_of(lon, lat)
                reg === nothing ? ("", 0.0) : (reg.name, ice_fraction(reg, lon, lat))
            end
            f > 0 || continue
            println(io, join((r..., name, string(f)), '\t'))
            n += 1
        end
    end
    println(level, ": ", length(rows), " candidates, ", n, " over land ice, ", length(frac), " frame centres tested")
end
