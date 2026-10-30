# An OTIS-format tide model's height constants as one cloud-optimized GeoTIFF in the band layout
# `app/src/tidecorrection.jl` reads, the one `make_cats_cog.jl` writes for CATS2008_v2023.
#
#   julia --project=tools/tidal tools/tidal/make_otis_cog.jl <h_file> <grid_file> <proj4> <out.tif>
#
# for Gr1kmTM v1 (doi:10.18739/A2B853K18):
#
#   julia --project=tools/tidal tools/tidal/make_otis_cog.jl h_Gr1kmTM_v1 grid_Gr1kmTM_v1 \
#       "+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +datum=WGS84 +units=m" Gr1kmTM_v1_z.tif
#
# Bands are `hRe_<c>` and `hIm_<c>` for each constituent of `TIDE_CONSTITUENTS`, interleaved, then
# `flexure`. A constituent the model does not carry is zero in both bands. `hIm` holds the file's
# imaginary part negated: pyTMD reads an OTIS file's constant as `re + i·im` but a TMD3 file's as
# `hRe − i·hIm`, and the app reads every model the TMD3 way.
#
# An OTIS model has no flexure: band `flexure` is 1 on wet cells. Dry cells within `FILL_CELLS` of a wet
# cell take the nearest wet cell's constants with flexure 1, so a point near the model coastline is
# interpolated from four filled cells. This is pyTMD's nearest-neighbour extrapolation within its default
# 10 km cutoff applied to cells rather than to points; the two agree wherever a point's four surrounding
# cells are wet. Every other cell is NaN with flexure 0.
#
# The file is a Fortran sequential binary in big-endian order: each record framed by 4-byte lengths.
# Heights: header (nx, ny, nc, ylim, xlim, nc 4-character names), then per constituent an (nx, ny) array
# of (re, im) Float32 pairs, x fastest. Grid: header (nx, ny, ylim, xlim, dt, nob), then the open-boundary
# indices, water-column thickness and the 0/1 mask (1 ocean). Limits are cell edges in km.

include("tide_cog.jl")

const FILL_CELLS = 10

_be(io, ::Type{T}) where {T} = ntoh(read(io, T))
_be(io, ::Type{T}, n::Int) where {T} = ntoh.(read!(io, Vector{T}(undef, n)))

function read_otis_heights(path::String)
    open(path) do io
        _be(io, Int32)
        nx, ny, nc = Int.(_be(io, Int32, 3))
        ylim, xlim = Float64.(_be(io, Float32, 2)), Float64.(_be(io, Float32, 2))
        names = [strip(String(read(io, 4))) for _ in 1:nc]
        _be(io, Int32)
        z = Dict{String,Matrix{ComplexF32}}()
        for c in names
            _be(io, Int32)
            v = _be(io, Float32, 2nx * ny)
            _be(io, Int32)
            # (re, im) pairs, x fastest: column-major [x, y].
            z[c] = reshape(complex.(v[1:2:end], v[2:2:end]), nx, ny)
        end
        return (; nx, ny, xlim, ylim, names, z)
    end
end

function read_otis_mask(path::String, nx::Int, ny::Int)
    open(path) do io
        _be(io, Int32)
        gnx, gny = Int.(_be(io, Int32, 2))
        (gnx, gny) == (nx, ny) || error("the grid file is ", gnx, "×", gny, "; the heights are ", nx, "×", ny)
        _be(io, Float32, 5)
        nob = Int(_be(io, Int32))
        _be(io, Int32)
        _be(io, Int32)
        nob == 0 ? _be(io, Int32) : _be(io, Int32, 2nob)
        _be(io, Int32)
        _be(io, Int32)
        h = reshape(_be(io, Float32, nx * ny), nx, ny)
        _be(io, Int32)
        _be(io, Int32)
        m = reshape(_be(io, Int32, nx * ny), nx, ny)
        # pyTMD's rule: a cell is wet where the mask says ocean and the water column is positive.
        return (m .== 1) .& (h .> 0)
    end
end

# The nearest wet cell to every dry cell within `r` cells, as a linear index; 0 where there is none.
function nearest_wet(wet::BitMatrix, r::Int)
    nx, ny = size(wet)
    src = zeros(Int, nx, ny)
    offsets = sort!([(di, dj) for di in -r:r for dj in -r:r if 0 < di^2 + dj^2 <= r^2]; by = o -> o[1]^2 + o[2]^2)
    L = LinearIndices(wet)
    Threads.@threads for j in 1:ny
        for i in 1:nx
            wet[i, j] && continue
            for (di, dj) in offsets
                a, b = i + di, j + dj
                (1 <= a <= nx && 1 <= b <= ny && wet[a, b]) || continue
                src[i, j] = L[a, b]
                break
            end
        end
    end
    return src
end

function make_otis_cog(hfile::String, gfile::String, proj4::String, dst::String)
    m = read_otis_heights(hfile)
    wet = read_otis_mask(gfile, m.nx, m.ny)
    # pyTMD also masks a cell whose constant is NaN or exactly zero.
    for z in values(m.z)
        wet .&= isfinite.(z) .& (z .!= 0)
    end
    extra = setdiff(m.names, TIDE_CONSTITUENTS)
    isempty(extra) || error("constituents ", join(extra, ", "), " are not among ", join(TIDE_CONSTITUENTS, " "))
    src = nearest_wet(wet, FILL_CELLS)
    dx, dy = (m.xlim[2] - m.xlim[1]) / m.nx, (m.ylim[2] - m.ylim[1]) / m.ny
    # North-up raster, metres: row 1 is the file's last y.
    gt = [m.xlim[1] * 1000, dx * 1000, 0.0, m.ylim[2] * 1000, 0.0, -dy * 1000]
    flip(a::Matrix) = a[:, end:-1:1]                  # [x, y ascending] → [x, y descending]
    filled(a) = reshape([wet[k] ? a[k] : src[k] > 0 ? a[src[k]] : NaN32 for k in eachindex(a)], size(a))
    bands = Pair{String,Matrix{Float32}}[]
    for c in TIDE_CONSTITUENTS
        z = get(m.z, c, nothing)
        part(f) = z === nothing ? zeros(Float32, m.nx, m.ny) : flip(filled(f.(z)))
        push!(bands, "hRe_" * c => part(real), "hIm_" * c => part(v -> -imag(v)))
    end
    push!(bands, "flexure" => flip(Float32.(wet .| (src .> 0))))
    names = write_tide_cog(dst, gt, proj4, bands)
    return (; names, wet = count(wet), filled = count(>(0), src))
end

if abspath(PROGRAM_FILE) == @__FILE__
    r = make_otis_cog(ARGS[1], ARGS[2], ARGS[3], ARGS[4])
    println(join(r.names, " "))
    println("wet cells ", r.wet, ", dry cells filled from a wet neighbour ", r.filled)
end
