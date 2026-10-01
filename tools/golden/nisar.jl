# The NISAR half of the chain from a granule: where a product's grid sits, and the secondary put on the
# reference's.
#
# Reading the samples is `SLCDatasets`' job and is done there — `pixels(open_slc(rslc))` returns a lazy
# `NisarRaster` and `amplitude` wraps it, `open_geocoded` opens a GSLC and reports its map grid. What is
# here is what the *pipeline* adds: the crop a geocoded pair is correlated over, and the coregistration a
# radar-geometry pair needs.
#
# **Both arrays stay lazy, and `AutoRIFT.ondisk` says so.** An RSLC band is 2.9 billion samples and a GSLC
# band 13 billion, so nothing is materialized and a blocked run sizes its own reads.

# `grid` is deliberately not imported: the harness uses that name for the parameter grid, and
# `GeocodedProduct` carries its own as a field.
using SLCDatasets: NisarRaster, Amplitude, open_geocoded

# Reading one of these costs I/O, so a blocked run windows its reads and an unblocked one is refused
# rather than resolving a pixel at a time.
AutoRIFT.ondisk(::NisarRaster) = true
AutoRIFT.ondisk(a::Amplitude) = AutoRIFT.ondisk(a.parent)

"""
    gslc_footprint(path; frequency = nothing) -> ImageFootprint

A geocoded NISAR product's own map grid, as the footprint [`coregister`](@ref) intersects on.

`SLCDatasets`' `GeocodedProduct` carries it; this is that grid in the harness's own type, which carries the CRS
so a pair in two projections is refused rather than silently intersected.
"""
function gslc_footprint(path::AbstractString; frequency = nothing)
    g = open_geocoded(path; frequency).grid
    return ImageFootprint(origin = g.origin, spacing = g.spacing,
                          size = reverse(g.size), crs = EPSG(g.epsg))
end

"""
    gslc_amplitude(path; frequency = nothing) -> Amplitude

A geocoded product's amplitude, lazily.

What the reference correlates: `convert_slc_to_uint8_amplitude` takes the magnitude and discards the
phase. The byte rescale is not applied — the correlator normalizes every chip, so quantizing here would
only throw resolution away.
"""
gslc_amplitude(path::AbstractString; frequency = nothing) =
    SLCDatasets.amplitude(pixels(open_geocoded(path; frequency)))

"""
    LazyWindow(parent, row0, col0, dims) <: AbstractMatrix{Float32}

A window of `parent` addressed from `(1, 1)`, without materializing it.

A geocoded pair is correlated over the overlap of its two grids, and that overlap is still billions of
samples — so the crop is an index shift and the correlator's blocked reads reach the HDF5 dataset
directly. A `view` would not do: `AutoRIFT._read_window!` indexes with ranges precisely because a view of
a lazy array defers to one read per element.
"""
struct LazyWindow{A<:AbstractMatrix} <: AbstractMatrix{Float32}
    parent::A
    row0::Int
    col0::Int
    dims::Tuple{Int,Int}
end

# The type parameter is computed here and the field declarations coerce, so GDAL's `Int32` width and
# height give a window rather than a `MethodError`.
LazyWindow(parent::AbstractMatrix, row0::Integer, col0::Integer, dims::Tuple{Integer,Integer}) =
    LazyWindow{typeof(parent)}(parent, row0, col0, dims)

Base.size(w::LazyWindow) = w.dims

function Base.getindex(w::LazyWindow, rows::AbstractUnitRange, cols::AbstractUnitRange)
    checkbounds(w, rows, cols)
    return Float32.(w.parent[(w.row0 .+ rows), (w.col0 .+ cols)])
end

Base.getindex(w::LazyWindow, i::Integer, j::Integer) = Float32(w.parent[w.row0 + i, w.col0 + j])
Base.getindex(w::LazyWindow, rows::AbstractUnitRange, j::Integer) =
    Float32.(w.parent[(w.row0 .+ rows), w.col0 + j])
Base.getindex(w::LazyWindow, i::Integer, cols::AbstractUnitRange) =
    Float32.(w.parent[w.row0 + i, (w.col0 .+ cols)])

AutoRIFT.ondisk(w::LazyWindow) = AutoRIFT.ondisk(w.parent)

"""
    gslc_window(fp::ImageFootprint, path) -> (row0, col0)

Where the raster at `path` starts inside the grid `fp` describes, as a zero-based `(row, col)` shift.

**Integer by construction, and checked.** Both grids are the same product's, at the same spacing in the
same projection, so the crop can only be a whole number of pixels from the granule's own origin; a
fractional offset means they are not the grid they claim to be, and is an error rather than a rounded
index.
"""
function gslc_window(fp, path::AbstractString)
    gt = ArchGDAL.getgeotransform(ArchGDAL.read(path))
    (gt[2] ≈ fp.spacing[1] && gt[6] ≈ fp.spacing[2]) || error(
        "$(basename(path)) is at $(gt[2]) x $(gt[6]) m and the granule's grid at " *
        "$(fp.spacing[1]) x $(fp.spacing[2]); the crop is not on the granule's own grid")
    cx = (gt[1] - fp.origin[1]) / fp.spacing[1]
    cy = (gt[4] - fp.origin[2]) / fp.spacing[2]
    for (v, ax) in ((cx, "x"), (cy, "y"))
        isinteger(round(v; digits = 6)) || error(
            "$(basename(path)) sits $v pixels along $ax from the granule's origin; a crop of the same " *
            "grid is a whole number of pixels")
    end
    return (round(Int, cy), round(Int, cx))
end

"""
    ResampledRSLC(reference, secondary, dem) <: AbstractMatrix{Float32}

The secondary RSLC's amplitude on the reference's radar grid, resampled a window at a time.

**An RSLC is not TOPS, and that is the whole of what makes this simpler than the Sentinel-1 path**: one
continuous swath rather than bursts, so there is no seam arithmetic, no per-burst carrier and no
deramp-and-reramp — the azimuth ramp a TOPS burst carries exists because its antenna sweeps, and this
antenna does not. What is left is the part both share: an offset field from `rdr2geo` on the reference and
`geo2rdr` on the secondary, and an eight-tap sinc through it.

The offset field is [`_offset_lattice`](@ref)'s, a node every 64 lines and 512 samples with bilinear
interpolation between — 90,400 nodes over a 57760 x 50511 grid against 2.9 billion per-pixel solves.
Resampling the reference against *itself* through it returns the reference to a correlation of 0.999997
and a mean ratio of 1.000007, which is what says the field and the kernel are right rather than merely
plausible.

Lazy for the reason the Sentinel-1 mosaic is: a materialized `Float32` copy of this grid is 11.7 GB and
the correlator only ever wants a block.
"""
struct ResampledRSLC{S,DL,DS} <: AbstractMatrix{Float32}
    source::S
    dl::DL
    ds::DS
    dims::Tuple{Int,Int}
    extent::Tuple{Int,Int}
    read::Base.RefValue{Int}
end

function ResampledRSLC(reference::AbstractString, secondary::AbstractString, dem)
    rs = open_slc(reference)
    ss = open_slc(secondary)
    src = pixels(ss)
    dims = (nlines(rs), nsamples(rs))
    dl, ds = _offset_lattice(RadarCoordinate(rs), RadarCoordinate(ss), dem, dims[1], dims[2])
    return ResampledRSLC(src, dl, ds, dims, size(src), Ref(0))
end

Base.size(r::ResampledRSLC) = r.dims

function Base.getindex(r::ResampledRSLC, rows::AbstractUnitRange, cols::AbstractUnitRange)
    checkbounds(r, rows, cols)
    # The offset field is addressed in the grid's own zero-based line and sample, which for one continuous
    # swath is the output index less one.
    lines = (first(rows) - 1):(last(rows) - 1)
    samples = (first(cols) - 1):(last(cols) - 1)
    ymin, ymax, xmin, xmax = _support_bounds(r.dl, r.ds, lines, samples)
    ny, nx = r.extent
    brows = clamp(floor(Int, ymin) - 3, 1, ny):clamp(ceil(Int, ymax) + 4, 1, ny)
    bcols = clamp(floor(Int, xmin) - 3, 1, nx):clamp(ceil(Int, xmax) + 4, 1, nx)
    # A window whose support lies wholly outside the secondary reads nothing: `resample_burst`'s bounds
    # test would reject every output pixel of it anyway.
    (isempty(brows) || isempty(bcols)) && return zeros(Float32, length(rows), length(cols))
    src = r.source[brows, bcols]
    r.read[] += length(brows) * length(bcols)
    # No `doppler`: that rejection reproduces `ResampSlc`'s own LUT test, which is a property of the
    # Sentinel-1 burst grid rather than of an interpolation.
    return resample_burst(src, r.dl, r.ds, lines, samples;
                          origin = (first(brows) - 1, first(bcols) - 1), extent = r.extent)
end

Base.getindex(r::ResampledRSLC, i::Integer, j::Integer) = r[i:i, j:j][1, 1]
Base.getindex(r::ResampledRSLC, rows::AbstractUnitRange, j::Integer) = r[rows, j:j][:, 1]
Base.getindex(r::ResampledRSLC, i::Integer, cols::AbstractUnitRange) = r[i:i, cols][1, :]

AutoRIFT.ondisk(::ResampledRSLC) = true
