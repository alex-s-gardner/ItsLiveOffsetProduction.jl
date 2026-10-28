"""
    AutoriftItsLiveProduction

ITS_LIVE offset production from two scenes: the sensor's tuned end-to-end workflow, from reading the scenes
at their source to the netCDF product, behind the command-line entry point `autorift_itslive_production`
(see [`julia_main`](@ref)). Built with `juliac --trim=safe` into a self-contained executable; every value
whose type would depend on file contents is resolved by an explicit branch.
"""
module AutoriftItsLiveProduction

using Dates: Dates, Date, DateTime
import Mmap
import NetworkOptions
import AutoRIFT
import ChunkManifests
import FastGeoProjections as FGP
import GeoFormatTypes as GFT
import ImagePairGeometry
using ImagePairGeometry: CoregisteredPair, GeometryInputs, IdentityTransform, ImageFootprint, InMemoryInputs,
                         MapGrid, PairGeometry, ProjectedCoordinate, TransformPair, coregister, footprint_bounds,
                         nested_window, nodata_from, pairgeometry_blocked, transform_pair
using ImagePairGeometry: GridWindow
import ItsLiveOffsetProduction
using ItsLiveOffsetProduction: AttrValue
import OpticalDatasets
import SLCDatasets

include("heap.jl")
include("gdal.jl")
include("io.jl")
include("job.jl")
include("scenes.jl")
include("nisar.jl")
include("mtl.jl")
include("parameters.jl")
include("geometry.jl")
include("tilecache.jl")
include("stage.jl")
include("streamed.jl")
include("imagery.jl")
include("correlate.jl")
include("packaging.jl")
include("radar.jl")
include("s1orbit.jl")
include("s1.jl")
include("s1subset.jl")
include("run.jl")
include("cli.jl")

end
