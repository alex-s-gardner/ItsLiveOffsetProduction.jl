# `ImagePairGeometry` geometry repackaged into this package's product structs. Each of these is a
# rename and a transpose, not a computation: the geogrid already derived the values from scene
# geometry, and a product needs them under ITS_LIVE's names and in `[row, col]` order.

using ImagePairGeometry: PairGeometry, GeometryInputs, RadarCoordinate

"""
    coefficients(g::PairGeometry) -> GeogridCoefficients

`g`'s displacement-to-velocity operator, scale factors, and (radar-only) direct axis-velocity band,
repackaged for [`write_product`](@ref) — a rename, not a computation: `g` is what already derives
these from scene geometry.

Bands are **transposed**, the same `[x,y]` to `[row,col]` swap `AutoRIFT.pointset` makes. A sentinel
entry becomes `NaN` rather than passing through as `-32767.0`: `write.jl`'s own helpers
(`_finite_median`, and every velocity-conversion formula) test `isnan`, not the sentinel, so returning
the raw value would leave an invalid point silently included in a median or a velocity sum.

`offset2vr`/`offset2va` come from `g.off2vx_dr`/`g.off2vy_dr` — the radar path's third off2vel band,
which `ImagePairGeometry.velocity_conversion` does not return — and are `nothing` for a projected `g`,
matching `GeogridCoefficients`'s `pair_type`-dependent contract.
"""
function coefficients(g::PairGeometry)
    sentinel = Float64(g.nodata.output)
    band(A) = [v == sentinel ? NaN : v for v in permutedims(A)]

    radar = g.coordinate isa RadarCoordinate
    return GeogridCoefficients(
        band(g.off2vx_dx), band(g.off2vx_dy),
        band(g.off2vy_dx), band(g.off2vy_dy),
        band(g.scale_x), band(g.scale_y),
        radar ? band(g.off2vx_dr) : nothing,
        radar ? band(g.off2vy_dr) : nothing)
end

"""
    reference_velocity(inputs::GeometryInputs) -> (ReferenceVelocity, stable_mask)

`inputs.vx`/`inputs.vy`/`inputs.ssm`, repackaged for [`write_product`](@ref) — the same rasters
`ImagePairGeometry.geometry_inputs` fetched to build the `PairGeometry` [`coefficients`](@ref)
above converts, so a caller who already has `inputs` for one gets the other at no extra fetch.

Transposed the same `[x,y]` to `[row,col]` way `coefficients` is, since `inputs`' arrays cover the same
window in the same layout. Left otherwise as the raw raster values, unlike `coefficients`: a
`GeometryInputs` array is whatever its file's own on-disk sentinel is (not `PairGeometry`'s uniform
`-32767`), and the reference itself reads these rasters the same raw way — masking here would diverge
from what was validated bit-exact against a real captured `VXref`/`VYref`/`SSM`.

Throws rather than returning a placeholder when `inputs` was built without `vx`/`vy`: `ItsLiveAutoRIFT`'s
stable-shift correction has no meaningful bias to remove without a reference velocity.
"""
function reference_velocity(inputs::GeometryInputs)
    inputs.vx === nothing && throw(ArgumentError(
        "reference_velocity needs GeometryInputs built with `vx`/`vy`; got neither"))
    t(A) = permutedims(A)
    stable_mask = inputs.ssm === nothing ? falses(size(t(inputs.vx))) : t(inputs.ssm) .!= 0
    return (ReferenceVelocity(Float32.(t(inputs.vx)), Float32.(t(inputs.vy))), stable_mask)
end

"""
    image_location(g::PairGeometry) -> (location_x, location_y)

`g`'s per-point image position — the pixel index each grid point falls at in the acquisition Geogrid
was run against — transposed the same `[x,y]` to `[row,col]` way [`coefficients`](@ref) is, with a
sentinel entry becoming `NaN`.

Needed only for the Sentinel-1 subswath-offset-bias correction (`SwathOffsetBias`), which masks by
which subswath a point's *range* position falls in — a different thing from `offset_x`/`offset_y`
(the expected displacement) or `location_x`/`location_y`'s role elsewhere, which this package has no
other use for.
"""
function image_location(g::PairGeometry)
    sentinel = Int32(g.nodata.output)
    band(A) = [v == sentinel ? NaN : Float64(v) for v in permutedims(A)]
    return (band(g.location_x), band(g.location_y))
end
