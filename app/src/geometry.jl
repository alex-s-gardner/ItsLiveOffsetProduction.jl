# The pair's geometry: footprints, coregistration, the parameter region, and the geogrid over the overlap.

# The sentinel every geogrid band marks a point outside the image with (`geogridOptical.cpp:1039`).
const SENTINEL = -32767.0

const Pair_ = CoregisteredPair{ProjectedCoordinate{Float64},Nothing}
const Geometry = PairGeometry{Matrix{Int32},Matrix{Float64},GFT.EPSG{1},ProjectedCoordinate{Float64}}

"""
    Setup

Everything the chain derives from a job before reading imagery: the scenes in acquisition order and their
native footprints, whether the pair is warped to one CRS, the coregistered pair, the parameter region,
the output grid and the window on it, the parameter rasters over that window, and the geogrid result. A
NISAR L2 pair also carries its two cropped GSLCs, `nothing` for an optical pair.
"""
struct Setup
    job::Job
    early::SceneSource
    late::SceneSource
    early_fp::Footprint
    late_fp::Footprint
    warp::Bool
    pair::Pair_
    epsg::Int
    region::Region
    grid::Grid
    window::GridWindow
    inputs::Inputs
    geometry::Geometry
    early_gslc::Union{Nothing,GslcScene}
    late_gslc::Union{Nothing,GslcScene}
end

"""
    setup(c::Job) -> Setup

Resolve `c`'s scenes, coregister them, look up the parameter region, and derive the output grid and the
geogrid over it. This is the reference's own order of operations: `coregister` first, because geogrid's
geotransform and size are the *overlap's* (`testGeogridOptical.py:92-100`), then the region lookup on the
overlap's centre, then the grid from the region's DEM.

A pair in two projections is coregistered on the footprints `ensure_same_projection`'s warp would produce
— both scenes onto the earlier one's CRS at its resolution, `-tap`-aligned — and nothing is warped here:
the imagery resamples onto those footprints a window at a time. A NISAR L2 pair is coregistered on its two
GSLCs' crops ([`gslc_crop`](@ref)).
"""
function setup(c::Job)
    early, late = acquisition_order(c)
    if c.platform == "NISAR-L2"
        eg, lg = gslc_crop(c, early, late)
        return _setup(c, early, late, eg.cropped, lg.cropped, false, eg.cropped, lg.cropped, eg, lg)
    end
    efp, lfp = scene_footprint(early.band), scene_footprint(late.band)
    target = footprint_epsg(efp)
    warp = target != footprint_epsg(lfp)
    rfp, sfp = warp ? (warped_footprint(early.band, efp, target), warped_footprint(late.band, efp, target)) :
               (efp, lfp)
    return _setup(c, early, late, efp, lfp, warp, rfp, sfp, nothing, nothing)
end

# Coregistration on the footprints `rfp` and `sfp`, the region lookup on the overlap's centre, and the
# geogrid.
function _setup(c::Job, early::SceneSource, late::SceneSource, efp::Footprint, lfp::Footprint, warp::Bool,
                rfp::Footprint, sfp::Footprint, eg::Union{Nothing,GslcScene}, lg::Union{Nothing,GslcScene})
    target = footprint_epsg(rfp)
    pair = coregister(rfp, sfp; dt = geogrid_seconds(c))
    lon, lat = pair_centroid(pair.coordinate, target)
    region = parameter_region(lon, lat)
    grid = parameter_grid(region)
    nd = dem_nodata(region)
    window, inputs, g = geogrid(grid, pair, region, nd, region.epsg, target)
    return Setup(c, early, late, efp, lfp, warp, pair, target, region, grid, window, inputs, g, eg, lg)
end

# The window the transform's footprint bounds give, the parameters over it, and the geogrid there — the
# only place the transform's type is known. The blocked geogrid is bit-identical to the serial one at any
# block size and task count; each output point depends only on its own inputs.
function geogrid(grid::Grid, pair::Pair_, region::Region, nd::Float64, tf::TransformPair)
    window = grid_window(grid, footprint_bounds(tf, pair.coordinate))
    inputs = geometry_inputs(region, window)
    g = pairgeometry_blocked(grid, pair, InMemoryInputs(inputs, window); transform = () -> tf, window,
                             nodata = nodata_from(nd))
    return (window, inputs, g)
end

function geogrid(grid::Grid, pair::Pair_, region::Region, nd::Float64, g::FGP.NativeCRS, s::FGP.NativeCRS)
    # The inverse from `inv`: the element type and math kernel are the operator's type parameters, so
    # rebuilding from the codes would revert both.
    f = FGP.Transformation(g, s)
    return geogrid(grid, pair, region, nd, TransformPair(f, inv(f)))
end

"""
    geogrid(grid, pair, region, nd, grid_epsg, scene_epsg)

The geogrid with the grid-to-scene transform FastGeoProjections implements for the pair of CRSs — one
concrete transform type per branch. Equal codes are the identity, which is what makes a scene already in
its parameter region's projection bitwise on every band.
"""
function geogrid(grid::Grid, pair::Pair_, region::Region, nd::Float64, gcode::Int, scode::Int)
    gcode == scode && return geogrid(grid, pair, region, nd, transform_pair(IdentityTransform()))
    g, s = FGP.native_crs(gcode), FGP.native_crs(scode)
    if g isa FGP.UTMCRS
        s isa FGP.UTMCRS && return geogrid(grid, pair, region, nd, g, s)
        s isa FGP.PolarStereographicCRS && return geogrid(grid, pair, region, nd, g, s)
    elseif g isa FGP.PolarStereographicCRS
        s isa FGP.UTMCRS && return geogrid(grid, pair, region, nd, g, s)
        s isa FGP.PolarStereographicCRS && return geogrid(grid, pair, region, nd, g, s)
    end
    throw(ArgumentError("no optical transform from the grid's EPSG:" * string(gcode) * " to the scene's EPSG:" *
                        string(scode)))
end

# Every pixel within `r` of a set one: one sweep along each axis carrying the distance since the last set
# pixel. The buffer is `Matrix{Bool}`, not `BitMatrix`: the sweeps write rows and columns from several
# tasks at once, and a `BitMatrix` packs 64 elements into a word, so neighbouring rows would lose each
# other's updates.
function _dilate(src::AbstractMatrix{Bool}, r::Int)
    out = Matrix{Bool}(src)
    nr, nc = size(out)
    Threads.@threads for i in 1:nr
        row = out[i, :]
        _sweep!(row, r)
        out[i, :] = row
    end
    Threads.@threads for j in 1:nc
        col = out[:, j]
        _sweep!(col, r)
        out[:, j] = col
    end
    return out
end

# Set every element within `r` of a set one, forward then backward, judging from the values as they
# stand during the sweep — which is what makes one forward and one backward pass a dilation.
function _sweep!(v::Vector{Bool}, r::Int)
    n = length(v)
    d = r + 1
    for k in 1:n
        d = v[k] ? 0 : d + 1
        d <= r && (v[k] = true)
    end
    d = r + 1
    for k in n:-1:1
        d = v[k] ? 0 : d + 1
        d <= r && (v[k] = true)
    end
    return v
end
