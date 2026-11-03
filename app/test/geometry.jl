# The pair geometry (`geometry.jl`): what `coregister` is handed.

using AutoriftItsLiveProduction: same_spacing
using ImagePairGeometry: ImageFootprint
import GeoFormatTypes as GFT

@testset "a pixel spacing rounded in the header still coregisters" begin
    a = ImageFootprint((500_000.0, 7_600_000.0), (15.0, -15.0), (100, 100), GFT.EPSG(32607))
    b = ImageFootprint((500_150.0, 7_600_000.0), (14.999999999999995, -15.0), (100, 100), GFT.EPSG(32607))
    @test same_spacing(a, b).spacing === a.spacing
    @test same_spacing(a, b).origin == b.origin
    c = ImageFootprint((500_150.0, 7_600_000.0), (15.01, -15.0), (100, 100), GFT.EPSG(32607))
    @test same_spacing(a, c) === c
end
