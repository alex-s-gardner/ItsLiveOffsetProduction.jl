# A product's browse image and thumbnail, drawn as `hyp3-autorift` draws them: the speed `v` through the
# ITS_LIVE color map (`hyp3_autorift.image.make_browse`), and a thumbnail no larger than 100 pixels a side
# (`hyp3lib.image.create_thumbnail`), pixel for pixel. Each step is the reference's own arithmetic:
#
#   1. `v`, its fill read as 0 and clipped to 625 m/yr, mapped onto [0, 1] by a PCHIP curve through the color
#      map's knots (SciPy's `PchipInterpolator`);
#   2. that stretched linearly between its minimum and maximum over the image (matplotlib's `Normalize`, as
#      `imsave` applies it with no limits given);
#   3. the result looked up in a 256-entry table interpolated linearly between the knots' colors, spaced
#      evenly (matplotlib's `LinearSegmentedColormap.from_list` and `Colormap.__call__` with `bytes`);
#   4. the thumbnail resampled with Pillow's bicubic filter over alpha-premultiplied 8-bit channels.
#
# A speed is a 16-bit integer, so steps 1 to 3 are tabulated once per value in the image.

# The ITS_LIVE color map's knots: speed (m/yr), red, green, blue, alpha.
const BROWSE_KNOTS = (
    (0, 255, 255, 255, 0), (2, 166, 238, 255, 255), (4, 97, 195, 219, 255), (9, 84, 169, 254, 255),
    (16, 84, 130, 254, 255), (25, 84, 85, 254, 255), (36, 50, 119, 220, 255), (49, 16, 153, 186, 255),
    (64, 16, 186, 153, 255), (81, 50, 220, 119, 255), (100, 84, 254, 85, 255), (121, 118, 221, 51, 255),
    (144, 153, 186, 16, 255), (169, 187, 152, 17, 255), (196, 221, 118, 51, 255), (225, 255, 85, 85, 255),
    (289, 255, 25, 85, 255), (324, 213, 1, 72, 255), (361, 158, 1, 66, 255), (400, 140, 0, 51, 255),
    (441, 122, 0, 166, 255), (484, 140, 0, 191, 255), (529, 159, 0, 217, 255), (576, 213, 0, 255, 255),
    (625, 255, 0, 138, 255),
)
const THUMBNAIL_SIZE = 100

# A speed of the product, its fill drawn as 0.
_drawn_speed(s::Int16) = s == ItsLiveOffsetProduction.NODATA_I16 ? 0 : Int(s)

# `numpy.linspace(0, 1, n)`: `i * (1 / (n - 1))`, its last value exactly 1.
_linspace01(n::Int) = (step = 1.0 / (n - 1); v = [i * step for i in 0:(n - 1)]; v[end] = 1.0; v)

# SciPy's one-sided three-point end slope for PCHIP (`PchipInterpolator._edge_case`).
function _pchip_edge(h0::Float64, h1::Float64, m0::Float64, m1::Float64)
    d = ((2 * h0 + h1) * m0 - h0 * m1) / (h0 + h1)
    sign(d) != sign(m0) && return 0.0
    (sign(m0) != sign(m1) && abs(d) > 3 * abs(m0)) && return 3 * m0
    return d
end

"""
    pchip_coefficients(x, y) -> Matrix{Float64}

The cubic of each interval of SciPy's `PchipInterpolator(x, y)`, highest power first, one column per
interval: its knot slopes from `_find_derivatives`, its coefficients from `CubicHermiteSpline`.
"""
function pchip_coefficients(x::Vector{Float64}, y::Vector{Float64})
    n = length(x)
    h = [x[k + 1] - x[k] for k in 1:(n - 1)]
    m = [(y[k + 1] - y[k]) / h[k] for k in 1:(n - 1)]
    d = zeros(n)
    for k in 2:(n - 1)
        (sign(m[k]) != sign(m[k - 1]) || m[k] == 0 || m[k - 1] == 0) && continue
        w1 = 2 * h[k] + h[k - 1]
        w2 = h[k] + 2 * h[k - 1]
        d[k] = 1.0 / ((w1 / m[k - 1] + w2 / m[k]) / (w1 + w2))
    end
    d[1] = _pchip_edge(h[1], h[2], m[1], m[2])
    d[n] = _pchip_edge(h[n - 1], h[n - 2], m[n - 1], m[n - 2])
    c = Matrix{Float64}(undef, 4, n - 1)
    for k in 1:(n - 1)
        t = (d[k] + d[k + 1] - 2 * m[k]) / h[k]
        c[1, k] = t / h[k]
        c[2, k] = (m[k] - d[k]) / h[k] - t
        c[3, k] = d[k]
        c[4, k] = y[k]
    end
    return c
end

# The piecewise cubic at `v`, its interval the one holding `v` (the last for the last knot, the end ones
# beyond), summed lowest power first as SciPy's `evaluate_poly1` sums it.
function pchip_value(x::Vector{Float64}, c::Matrix{Float64}, v::Float64)
    k = clamp(searchsortedlast(x, v), 1, length(x) - 1)
    s = v - x[k]
    res = c[4, k]
    z = s
    res += c[3, k] * z
    z *= s
    res += c[2, k] * z
    z *= s
    res += c[1, k] * z
    return res
end

"""
    browse_colormap() -> Matrix{UInt8}

The 256 RGBA colors, one per row, of matplotlib's `LinearSegmentedColormap.from_list` over the knots'
colors spaced evenly on [0, 1], as `Colormap.__call__(…, bytes = True)` reads them: each channel's table
interpolated linearly (`_create_lookup_table`) and truncated to a byte.
"""
function browse_colormap()
    N = 256
    vals = _linspace01(length(BROWSE_KNOTS))
    xk = [v * (N - 1) for v in vals]
    xind = [(N - 1) * v for v in _linspace01(N)]
    byte(c) = unsafe_trunc(UInt8, clamp(c, 0.0, 1.0) * 255)
    lut = Matrix{UInt8}(undef, N, 4)
    for ch in 1:4
        y = [k[ch + 1] / 255 for k in BROWSE_KNOTS]
        lut[1, ch] = byte(y[1])
        lut[N, ch] = byte(y[end])
        for i in 2:(N - 1)
            j = searchsortedfirst(xk, xind[i])
            dist = (xind[i] - xk[j - 1]) / (xk[j] - xk[j - 1])
            lut[i, ch] = byte(dist * (y[j] - y[j - 1]) + y[j - 1])
        end
    end
    return lut
end

"""
    browse_rgba(v) -> Array{UInt8,3}

The browse image of speed `v` (indexed `[x, y]`, `NODATA_I16` where there is no value), as RGBA indexed
`[channel, x, y]`.
"""
function browse_rgba(v::AbstractMatrix{Int16})
    x = Float64[k[1] for k in BROWSE_KNOTS]
    c = pchip_coefficients(x, _linspace01(length(BROWSE_KNOTS)))
    # The curve at each speed the image holds, clipped to the last knot.
    level = [pchip_value(x, c, min(Float64(s), x[end])) for s in 0:maximum(_drawn_speed, v)]
    lo, hi = extrema(level[_drawn_speed(s) + 1] for s in v)
    # The color table entry of each level, as `Normalize` then `Colormap.__call__` index it.
    N = 256
    entry = Vector{Int}(undef, length(level))
    for i in eachindex(level)
        xa = lo == hi ? 0.0 : (level[i] - lo) / (hi - lo) * N
        xa == N && (xa = Float64(N - 1))
        entry[i] = clamp(unsafe_trunc(Int, xa), 0, N - 1) + 1
    end
    lut = browse_colormap()
    out = Array{UInt8,3}(undef, 4, size(v, 1), size(v, 2))
    for j in axes(v, 2), i in axes(v, 1)
        e = entry[_drawn_speed(v[i, j]) + 1]
        for ch in 1:4
            out[ch, i, j] = lut[e, ch]
        end
    end
    return out
end

# ---------------------------------------------------------------------------
# The thumbnail: Pillow's `Image.thumbnail`
# ---------------------------------------------------------------------------

"""
    thumbnail_size(w, h) -> (Int, Int)

The size `Image.thumbnail((100, 100))` gives a `w` by `h` image: unchanged when it already fits, else the
largest that fits with the side it shortens rounded to keep the aspect ratio closest.
"""
function thumbnail_size(w::Int, h::Int)
    n = THUMBNAIL_SIZE
    (n >= w && n >= h) && return (w, h)
    aspect = w / h
    if 1 >= aspect
        x = _round_aspect(n * aspect, k -> abs(aspect - k / n))
        return (x, n)
    end
    y = _round_aspect(n / aspect, k -> k == 0 ? 0.0 : abs(aspect - n / k))
    return (n, y)
end

# Floor or ceiling of `v`, whichever `key` ranks lower (the floor on a tie), at least 1.
function _round_aspect(v::Float64, key::F) where {F}
    a, b = floor(Int, v), ceil(Int, v)
    return max(key(b) < key(a) ? b : a, 1)
end

# Pillow's bicubic kernel, `a = -0.5`.
function _bicubic(x::Float64)
    a = -0.5
    x = abs(x)
    x < 1.0 && return ((a + 2.0) * x - (a + 3.0)) * x * x + 1
    x < 2.0 && return (((x - 5) * x + 8) * x - 4) * a
    return 0.0
end

const PRECISION_BITS = 22

# The taps of each output sample of a resample from `n` to `out` samples (`precompute_coeffs`), as the
# fixed-point weights `normalize_coeffs_8bpc` makes of them, and the first input sample each starts at.
function _resample_taps(n::Int, out::Int)
    scale = n / out
    filterscale = max(scale, 1.0)
    support = 2.0 * filterscale
    taps = Vector{Vector{Int}}(undef, out)
    first_ = Vector{Int}(undef, out)
    for o in 0:(out - 1)
        center = (o + 0.5) * scale
        xmin = max(unsafe_trunc(Int, center - support + 0.5), 0)
        xmax = min(unsafe_trunc(Int, center + support + 0.5), n) - xmin
        ss = 1.0 / filterscale
        w = [_bicubic((k + xmin - center + 0.5) * ss) for k in 0:(xmax - 1)]
        # Summed in order, as Pillow sums them.
        ww = 0.0
        for x in w
            ww += x
        end
        ww != 0.0 && (w ./= ww)
        taps[o + 1] = [p < 0 ? unsafe_trunc(Int, -0.5 + p * (1 << PRECISION_BITS)) :
                       unsafe_trunc(Int, 0.5 + p * (1 << PRECISION_BITS)) for p in w]
        first_[o + 1] = xmin
    end
    return taps, first_
end

_clip8(s::Int) = UInt8(clamp(s >> PRECISION_BITS, 0, 255))

# Pillow's `rgbA2rgba`: each color channel times alpha over 255, rounded.
function _premultiply!(row::Matrix{UInt8}, rgba::Array{UInt8,3}, j::Int)
    for i in axes(row, 2)
        a = UInt32(rgba[4, i, j])
        for ch in 1:3
            t = UInt32(rgba[ch, i, j]) * a + 0x80
            row[ch, i] = UInt8(((t >> 8) + t) >> 8)
        end
        row[4, i] = rgba[4, i, j]
    end
    return row
end

# One output sample of 8-bit channel `ch`: the taps `k` over `get(p)` from input sample `f + 1` on.
function _tap_sum(get::G, k::Vector{Int}, f::Int) where {G}
    s = 1 << (PRECISION_BITS - 1)
    for t in eachindex(k)
        s += Int(get(f + t)) * k[t]
    end
    return _clip8(s)
end

# `rgba` (indexed `[channel, x, y]`), alpha premultiplied a row at a time, resampled along x to `n`.
function _resample_x(rgba::Array{UInt8,3}, n::Int)
    taps, first_ = _resample_taps(size(rgba, 2), n)
    out = Array{UInt8,3}(undef, 4, n, size(rgba, 3))
    row = Matrix{UInt8}(undef, 4, size(rgba, 2))
    for j in axes(rgba, 3)
        _premultiply!(row, rgba, j)
        for i in 1:n, ch in 1:4
            out[ch, i, j] = _tap_sum(p -> row[ch, p], taps[i], first_[i])
        end
    end
    return out
end

# `pre` (premultiplied, indexed `[channel, x, y]`) resampled along y to `n`.
function _resample_y(pre::Array{UInt8,3}, n::Int)
    taps, first_ = _resample_taps(size(pre, 3), n)
    out = Array{UInt8,3}(undef, 4, size(pre, 2), n)
    for j in 1:n, i in axes(pre, 2), ch in 1:4
        out[ch, i, j] = _tap_sum(p -> pre[ch, i, p], taps[j], first_[j])
    end
    return out
end

# `rgba` premultiplied whole, for a thumbnail that keeps its width.
function _premultiplied(rgba::Array{UInt8,3})
    pre = similar(rgba)
    row = Matrix{UInt8}(undef, 4, size(rgba, 2))
    for j in axes(rgba, 3)
        pre[:, :, j] = _premultiply!(row, rgba, j)
    end
    return pre
end

"""
    thumbnail_rgba(rgba) -> Array{UInt8,3}

`rgba` (indexed `[channel, x, y]`) as Pillow's `Image.thumbnail((100, 100))` shrinks an RGBA image: alpha
premultiplied, resampled with the bicubic filter along x and then y — each only when its size changes —
and divided out again.
"""
function thumbnail_rgba(rgba::Array{UInt8,3})
    w, h = size(rgba, 2), size(rgba, 3)
    tw, th = thumbnail_size(w, h)
    (tw, th) == (w, h) && return copy(rgba)
    r = tw == w ? _premultiplied(rgba) : _resample_x(rgba, tw)
    th == h || (r = _resample_y(r, th))
    for j in axes(r, 3), i in axes(r, 2)
        a = Int(r[4, i, j])
        (a == 255 || a == 0) && continue
        for ch in 1:3
            r[ch, i, j] = UInt8(min(255 * Int(r[ch, i, j]) ÷ a, 255))
        end
    end
    return r
end

"""
    write_browse(product) -> (String, String)

The browse image `<product>.png` and thumbnail `<product>_thumb.png` of the product netCDF at `product`,
written beside it, drawn from the speed `v` it holds; their paths.
"""
function write_browse(product::String)
    rgba = browse_rgba(ItsLiveOffsetProduction.read_speed(product))
    stem = endswith(product, ".nc") ? product[1:(end - 3)] : product
    browse, thumb = stem * ".png", stem * "_thumb.png"
    write_png(browse, rgba)
    write_png(thumb, thumbnail_rgba(rgba))
    return browse, thumb
end
