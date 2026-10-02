"""
    ItsLiveOffsetProduction

ITS_LIVE offset products: the inputs a granule pair yields, and the netCDF they are written to.

AutoRIFT.jl correlates an image pair and returns displacements. Turning those into an ITS_LIVE
product is a separate job with its own conventions — plugin and software version strings, scene
naming, swath offset bias, the stable-shift correction, and three netCDF schemas — and this package
owns all of it. [`write_product`](@ref) is the entry point; everything else exists to build the
[`ItsLiveInput`](@ref) it takes.

Every dependency is a hard dependency. The netCDF writer, the geogrid repackaging and the CF grid
mapping are all needed on every run, so there is nothing for a package extension to defer.
"""
module ItsLiveOffsetProduction

using Dates: Dates, DateTime
using Statistics: median, quantile, std

import ArchGDAL
import AutoRIFT
import DimensionalData
import GeoFormatTypes
import ImagePairGeometry
import NCDatasets
import Rasters

export GeogridCoefficients, ReferenceVelocity, ItsLiveGeoref, ImagePairInfo, ItsLiveInput,
       SwathOffsetBias
export write_product, roi_valid_percentage, cf_grid_mapping
# `coefficients`, `reference_velocity` and `image_location` are deliberately not exported: the names are
# generic enough to collide in a caller's namespace, and each is called once per run.

# Types and the function declarations their methods attach to.
include("types.jl")
# `PairGeometry`/`GeometryInputs` to the structs above.
include("geogrid.jl")
# CF grid-mapping attributes, and an `ItsLiveGeoref` from a `Raster`.
include("georef.jl")
# Velocity, stable-shift correction, error estimates, and the netCDF itself.
include("write.jl")

# ---------------------------------------------------------------------------
# Precompilation
# ---------------------------------------------------------------------------

# `write_product` is compiled, not slow. On a 2351x2341 grid it takes 0.78 s warm and 8.5 s on the
# first call in a fresh process; inside the pipeline, where the earlier stages have already compiled
# the shared machinery, the first call still costs 4.86 s of a 6.6 s packaging stage. Production runs
# one process per image pair, so every pair pays that, and a workload here moves it into the package
# image once per install.
#
# The grid is the smallest the writer accepts: `x`/`y` must be spaced at exactly 120 m, and the output
# is padded to a 512-pixel tile regardless, so a 50x30 input compiles the same code a scene-sized one
# would. Both schemas are exercised — the optical one production uses, and the radar one, which adds
# `vr`/`va` and the conversion matrix through branches the optical path never reaches.
using PrecompileTools: @setup_workload, @compile_workload

@setup_workload begin
    nx, ny = 50, 30
    x = collect(0.0:120.0:(120.0 * (nx - 1)))
    y = collect(reverse(0.0:120.0:(120.0 * (ny - 1)))) .+ 1.0e6
    mapping = Dict{String,Any}("grid_mapping_name" => "polar_stereographic",
                               "spatial_epsg" => Int64(3413))
    # A plain function rather than a closure over a projection: the workload must not need `Proj`.
    lonlat(px, py) = (-45.0, 70.0)
    georef = ItsLiveGeoref(x, y, mapping, lonlat, 30.0, 30.0)
    info = ImagePairInfo(DateTime(2020, 1, 1), DateTime(2020, 1, 13), "L", "L", "8", "8",
                         80.0, 70.0, -45.0, Dict{String,Any}())
    refv = ReferenceVelocity(zeros(Float32, ny, nx), zeros(Float32, ny, nx))
    dx = fill(6.0f0, ny, nx)
    dy = fill(-2.0f0, ny, nx)
    chip = fill(UInt16(32), ny, nx)
    interp = falses(ny, nx)
    ssm = trues(ny, nx)

    optical = GeogridCoefficients(fill(0.02, ny, nx), fill(0.001, ny, nx),
                                  fill(0.001, ny, nx), fill(0.02, ny, nx),
                                  ones(ny, nx), ones(ny, nx), nothing, nothing)
    radar = GeogridCoefficients(fill(0.02, ny, nx), fill(0.001, ny, nx),
                                fill(0.001, ny, nx), fill(0.02, ny, nx),
                                ones(ny, nx), ones(ny, nx), fill(0.3, ny, nx), fill(0.3, ny, nx))
    loc = repeat(reshape(Float64.(0:(nx - 1)), 1, nx), ny, 1)
    swath = SwathOffsetBias(loc, copy(loc), 1.0, nx, 20.0, 40.0, false, "S1A")

    @compile_workload begin
        for (kind, coeffs, bias) in ((:optical, optical, nothing), (:radar, radar, swath))
            input = ItsLiveInput(kind, "feature", kind === :radar ? "radar" : "map", "2.1.1",
                                 "params.shp", "precompile", 12.0 * 86400,
                                 georef, info, coeffs, refv, ssm, bias,
                                 dx, dy, chip, 1.0, interp)
            path = tempname() * ".nc"
            try
                write_product(path, input)
            finally
                rm(path; force = true)
            end
        end
    end
end

end # module
