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

end # module
