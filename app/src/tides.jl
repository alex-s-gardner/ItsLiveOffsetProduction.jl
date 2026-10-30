# Ocean tide heights from the harmonic constants of CATS2008_v2023 (Howard et al., 2024,
# doi:10.15784/601772), predicted as pyTMD predicts a TMD3 model: `tide_elevations` with
# `corrections = "TMD3"`, delta-T 0 and minor constituents inferred (`pyTMD.predict.ocean_load`,
# `pyTMD.constituents`). Constants are complex amplitudes in metres under pyTMD's convention for TMD3
# files, imaginary part negated from the file's `hIm`, so the height of one constituent is
# `f·(Re·cos θ − Im·sin θ)` with `θ = ω·t + φ₀ + u`.

# The ten constituents CATS2008_v2023 carries, in the order every constant vector here is kept.
const TIDE_CONSTITUENTS = ("m2", "s2", "n2", "k2", "k1", "o1", "p1", "q1", "mf", "mm")

# Angular frequency (rad/s) and phase at the tide epoch (rad), `pyTMD.constituents._constituent_parameters`.
const TIDE_OMEGA = (1.405189e-04, 1.454441e-04, 1.378797e-04, 1.458423e-04, 7.292117e-05,
                    6.759774e-05, 7.252295e-05, 6.495854e-05, 0.053234e-04, 0.026392e-04)
const TIDE_PHASE0 = (1.731557546, 0.0, 6.050721243, 3.487600001, 0.173003674,
                     1.558553872, 6.110181633, 5.877717569, 1.756042456, 1.964021610)

# The tide epoch, 1992-01-01T00:00:00 UTC, as a Modified Julian Day.
const TIDE_EPOCH = DateTime(1992, 1, 1)
const MJD_TIDE_EPOCH = 48622.0

"""
    tide_days(t::DateTime) -> Float64

Days from 1992-01-01T00:00:00 UTC to `t`, the time argument of every function here.
"""
tide_days(t::DateTime) = Dates.value(t - TIDE_EPOCH) / 86_400_000

_normalize_degrees(x::Float64) = mod(x, 360.0)

# Mean longitudes (degrees) of the Moon `s`, Sun `h`, lunar perigee `p` and lunar node `n`, from David
# Cartwright's formulae for 1990–2010 (`pyTMD.astro.mean_longitudes(method = "Cartwright")`), and the
# hour angle `tau` of mean lunar time. The solar perigee is held at 282.8°.
function _mean_longitudes(days::Float64)
    mjd = days + MJD_TIDE_EPOCH
    T = mjd - 51544.4993
    s = _normalize_degrees(218.3164 + 13.17639648 * T)
    h = _normalize_degrees(280.4661 + 0.98564736 * T)
    p = _normalize_degrees(83.3535 + 0.11140353 * T)
    n = _normalize_degrees(125.0445 - 0.05295377 * T)
    tau = 15.0 * 24.0 * mod(mjd, 1.0) - s + h
    return (; tau, s, h, p, n, pp = 282.8)
end

# Nodal factor `f` and angle `u` (rad) of each major constituent, `pyTMD.constituents.nodal_modulation`
# for an OTIS-type model.
function _major_nodal(n::Float64)
    N = deg2rad(n)
    sinn, cosn = sincos(N)
    sin2n, cos2n = sincos(2N)
    sin3n = sin(3N)
    hu(t1, t2) = (hypot(t1, t2), atan(t1, t2))
    m2 = hu(-0.03731 * sinn + 0.00052 * sin2n, 1.0 - 0.03731 * cosn + 0.00052 * cos2n)
    k2 = hu(-0.3108 * sinn - 0.0324 * sin2n, 1.0 + 0.2852 * cosn + 0.0324 * cos2n)
    k1 = hu(-0.1554 * sinn + 0.0029 * sin2n, 1.0 + 0.1158 * cosn - 0.0029 * cos2n)
    o1 = (hypot(0.189 * sinn - 0.0058 * sin2n, 1.0 + 0.189 * cosn - 0.0058 * cos2n),
          deg2rad(10.8 * sinn - 1.3 * sin2n + 0.2 * sin3n))
    q1 = (hypot(1.0 + 0.188 * cosn, 0.188 * sinn), atan(0.189 * sinn / (1.0 + 0.189 * cosn)))
    mf = (1.043 + 0.414 * cosn, deg2rad(-23.7 * sinn + 2.7 * sin2n - 0.4 * sin3n))
    mm = (1.0 - 0.130 * cosn, 0.0)
    one = (1.0, 0.0)
    # m2 s2 n2 k2 k1 o1 p1 q1 mf mm
    return (m2, one, m2, k2, k1, o1, one, q1, mf, mm)
end

# Complex amplitudes of the 18 short-period minor constituents inferred from the majors
# (`pyTMD.predict._infer_short_period`), and their nodal factors, angles and Doodson coefficients
# (`pyTMD.constituents.minor_arguments`). Doodson coefficients multiply (tau, s, h, p, n, pp, 90°).
const MINOR_DOODSON = (
    (1, -3, 0, 2, 0, 0, -1),   # 2q1
    (1, -3, 2, 0, 0, 0, -1),   # sigma1
    (1, -2, 2, -1, 0, 0, -1),  # rho1
    (1, 0, 0, -1, 0, 0, 1),    # m1b
    (1, 0, 0, 1, 0, 0, 1),     # m1a
    (1, 0, 2, -1, 0, 0, 1),    # chi1
    (1, 1, -3, 0, 0, 1, -1),   # pi1
    (1, 1, 2, 0, 0, 0, 1),     # phi1
    (1, 2, -2, 1, 0, 0, 1),    # theta1
    (1, 2, 0, -1, 0, 0, 1),    # j1
    (1, 3, 0, 0, 0, 0, 1),     # oo1
    (2, -2, 0, 2, 0, 0, 0),    # 2n2
    (2, -2, 2, 0, 0, 0, 0),    # mu2
    (2, -1, 2, -1, 0, 0, 0),   # nu2
    (2, 1, -2, 1, 0, 0, 2),    # lambda2
    (2, 1, 0, -1, 0, 0, 2),    # l2a
    (2, 1, 0, 1, 0, 0, 0),     # l2b
    (2, 2, -3, 0, 0, 1, 0),    # t2
)

function _minor_amplitudes(z::NTuple{10,ComplexF64})
    m2, s2, n2, k2, k1, o1, p1, q1 = z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8]
    return (0.263 * q1 - 0.0252 * o1, 0.297 * q1 - 0.0264 * o1, 0.164 * q1 + 0.0048 * o1,
            0.0140 * o1 + 0.0101 * k1, 0.0389 * o1 + 0.0282 * k1, 0.0064 * o1 + 0.0060 * k1,
            0.0030 * o1 + 0.0171 * k1, -0.0015 * o1 + 0.0152 * k1, -0.0065 * o1 + 0.0155 * k1,
            -0.0389 * o1 + 0.0836 * k1, -0.0431 * o1 + 0.0613 * k1, 0.264 * n2 - 0.0253 * m2,
            0.298 * n2 - 0.0264 * m2, 0.165 * n2 + 0.00487 * m2, 0.0040 * m2 + 0.0074 * s2,
            0.0131 * m2 + 0.0326 * s2, 0.0033 * m2 + 0.0082 * s2, 0.0585 * s2)
end

# The u expression for 2q1 subtracts `0.0058 sin 2N` in its denominator where its f uses `cos 2N`;
# that is `pyTMD.constituents.minor_arguments` as written, kept so the inferred heights match it.
function _minor_nodal(n::Float64)
    N = deg2rad(n)
    sinn, cosn = sincos(N)
    sin2n, cos2n = sincos(2N)
    fq = hypot(1.0 + 0.189 * cosn - 0.0058 * cos2n, 0.189 * sinn - 0.0058 * sin2n)
    uq = atan(0.189 * sinn - 0.0058 * sin2n, 1.0 + 0.189 * cosn - 0.0058 * sin2n)
    fn = hypot(1.0 - 0.0373 * cosn, 0.0373 * sinn)
    un = atan(-0.0373 * sinn, 1.0 - 0.0373 * cosn)
    f = (fq, fq, fq, hypot(1.0 + 0.185 * cosn, 0.185 * sinn), hypot(1.0 + 0.201 * cosn, 0.201 * sinn),
         hypot(1.0 + 0.221 * cosn, 0.221 * sinn), 1.0, 1.0, 1.0, hypot(1.0 + 0.198 * cosn, 0.198 * sinn),
         hypot(1.0 + 0.640 * cosn + 0.134 * cos2n, 0.640 * sinn + 0.134 * sin2n), fn, fn, fn, 1.0, fn,
         hypot(1.0 + 0.441 * cosn, 0.441 * sinn), 1.0)
    u = (uq, uq, uq, atan(0.185 * sinn, 1.0 + 0.185 * cosn), atan(-0.201 * sinn, 1.0 + 0.201 * cosn),
         atan(-0.221 * sinn, 1.0 + 0.221 * cosn), 0.0, 0.0, 0.0, atan(-0.198 * sinn, 1.0 + 0.198 * cosn),
         atan(-0.640 * sinn - 0.134 * sin2n, 1.0 + 0.640 * cosn + 0.134 * cos2n), un, un, un, 0.0, un,
         atan(-0.441 * sinn, 1.0 + 0.441 * cosn), 0.0)
    return f, u
end

# Long-period minor constituents inferred by linear interpolation in frequency between the node tide,
# Mm and Mf, normalized by their Cartwright–Edden amplitudes (`pyTMD.predict._admittance_long_period`
# without anelasticity). CATS2008_v2023 has no node constituent, and pyTMD skips this inference when
# one of the three is missing, so it adds nothing here; it is not ported.

"""
    TidePrediction

The time-dependent part of a tide prediction at one instant, shared by every point: per constituent,
`f·cos θ` and `f·sin θ` for the ten majors and the eighteen inferred minors.
"""
struct TidePrediction
    major_cos::NTuple{10,Float64}
    major_sin::NTuple{10,Float64}
    minor_cos::NTuple{18,Float64}
    minor_sin::NTuple{18,Float64}
end

# Tuple lengths are `Val`s so that every tuple's length is part of its type, which `--trim` needs.
function TidePrediction(t::DateTime)
    days = tide_days(t)
    L = _mean_longitudes(days)
    nodal = _major_nodal(L.n)
    θ = ntuple(i -> 86400.0 * TIDE_OMEGA[i] * days + TIDE_PHASE0[i] + nodal[i][2], Val(10))
    major_cos = ntuple(i -> nodal[i][1] * cos(θ[i]), Val(10))
    major_sin = ntuple(i -> nodal[i][1] * sin(θ[i]), Val(10))
    f, u = _minor_nodal(L.n)
    args = (L.tau, L.s, L.h, L.p, L.n, L.pp, 90.0)
    θm = ntuple(Val(18)) do i
        c = MINOR_DOODSON[i]
        G = sum(ntuple(k -> c[k] * args[k], Val(7)))
        deg2rad(G) + u[i]
    end
    minor_cos = ntuple(i -> f[i] * cos(θm[i]), Val(18))
    minor_sin = ntuple(i -> f[i] * sin(θm[i]), Val(18))
    return TidePrediction(major_cos, major_sin, minor_cos, minor_sin)
end

"""
    tide_height(z, p::TidePrediction) -> Float64
    tide_height(z, t::DateTime) -> Float64

Ocean tide height in metres from the ten complex constants `z` (in [`TIDE_CONSTITUENTS`](@ref) order),
with the short-period minor constituents inferred from them.
"""
function tide_height(z::NTuple{10,ComplexF64}, p::TidePrediction)
    h = 0.0
    for i in 1:10
        h += real(z[i]) * p.major_cos[i] - imag(z[i]) * p.major_sin[i]
    end
    zm = _minor_amplitudes(z)
    for i in 1:18
        h += real(zm[i]) * p.minor_cos[i] - imag(zm[i]) * p.minor_sin[i]
    end
    return h
end

tide_height(z::NTuple{10,ComplexF64}, t::DateTime) = tide_height(z, TidePrediction(t))
