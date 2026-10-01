# The Sentinel-1 radar geometry geogrid is given, built the way the pipeline builds it.
#
# `loadMetadataSlc` (`vend/testGeogrid.py:162-212`) is the whole specification, and two things about it
# make the radar path far cheaper than it looks.
#
# **The coregistration does not enter the geometry.** `s1_isce3.process_slc` calls
# `loadMetadataSlc(safe_ref, orbit_ref, ...)` and then sets `meta_s = copy.copy(meta_r)` with only
# `sensingStart` and `sensingStop` replaced. So every number geogrid consumes — and therefore every
# `window_*.tif` band — comes from the *reference* acquisition's burst annotations and its orbit. What
# COMPASS's per-burst resample and hyp3's `merge_swaths` produce is `secondary.tif`'s pixel values, which
# the geometry rungs never read.
#
# **The mosaic's grid is mostly metadata, and one number of it is not.** hyp3 lays all three subswaths
# on the near-range subswath's range origin and zero-fills between them (`s1_isce3.merge_swaths`), and
# `loadMetadataSlc:210` hands geogrid that mosaic's own shape.
#
# **Both extents follow the burst shape `merge_bursts_in_swath` reads, and that is the CSLC's rather than
# the annotation's.** `total_rng_samples` is `last_rng_samples + floor((far.starting_range -
# near.starting_range) / dr)`, and `last_rng_samples` is `num_rng_samples` — taken off the *COMPASS CSLC
# raster* of the far subswath's first burst (`s1_isce3.py:565-566`). The azimuth extent inherits the same
# dependency, since the per-swath merged length is built from that raster's `num_az_samples`.
#
# The two agree exactly when the reference burst was written by `rdr2geo`, and not when it was resampled
# against a cached static topographic layer — `write_yaml`'s `use_static_layer` branch sets
# `bool_reference` false, so COMPASS coregisters the reference burst instead of writing it out, and the
# result lands on its own geocoded grid. Which step ran is readable off the burst directory: `x`/`y`/`z`
# beside a polarization-named SLC for `rdr2geo`, `azimuth.off` beside a stem-named one for the resample.
# On `S1A ... 20170221` the difference is 136 to 195 lines taller and 72 to 97 samples narrower per
# subswath, which is exactly its 181-short, 85-wide disagreement.
#
# [`cslc_grid`](@ref) reads that shape from `product/` when the run kept one, and returns `nothing`
# otherwise — where the annotation is already the right answer. So the derivation stays a derivation: it
# reads COMPASS's own intermediate, which is an *input* to `merge_swaths`, never `reference.tif`, which is
# what the rung compares against.
#
# **`dt` is full precision here, unlike the optical path.** `testGeogrid.py:439` sets
# `repeatTime = (info1.sensingStart - info.sensingStart).total_seconds()` for radar, against `:354`'s
# whole-calendar-day difference for optical. The two conventions are opposite and each is wrong for the
# other path.

using Dates, Printf
using AutoRIFT   # `ondisk`, which the lazy resamplers below answer for themselves
using ImagePairGeometry
using ImagePairGeometry: IdentityTransform, LookRight, incidence_angle
using SLCDatasets
using SLCDatasets: annotation, asf_bursts, burst_raster, bursts, merge_bursts, nbursts,
                   open_slc, orbit, seconds_between, Sentinel1Product

# Sentinel-1 IW covers the swath with three subswaths, and a full-SLC job processes all three
# (`s1_isce3.process_slc` defaults `swaths=(1, 2, 3)`).
const S1_SWATHS = 1:3

"""
    s1_polarization(granule) -> String

The channel the pipeline correlates, from the granule name.

`process.py` asks COMPASS for `co-pol`, which is the like-polarized channel: `HH` for an `SH` or `DH`
product and `VV` for an `SV` or `DV` one. The cross-polarized channel of a dual-pol product is never
the one tracked.
"""
function s1_polarization(granule::AbstractString)
    m = match(r"_1S(S|D)(H|V)_", granule)
    m === nothing && throw(ArgumentError(
        "cannot read a polarization out of \"$granule\"; expected a field like `1SSH` or `1SDV`"))
    return m.captures[2] == "H" ? "HH" : "VV"
end

"""
    s1_orbit(dir, granule) -> String

The orbit file in `dir` whose validity window covers `granule`'s acquisition.

A run directory holds the orbits the container downloaded for that pair, one per granule, and their
names carry the window as `V<start>_<stop>`. Chosen by containment rather than by position, because the
two files are not in a defined order and picking the wrong one puts the state vectors a day away — an
error the geometry would absorb into a plausible-looking footprint rather than reject.
"""
function s1_orbit(dir::AbstractString, granule::AbstractString)
    t = DateTime(split(granule, '_')[6], dateformat"yyyymmddTHHMMSS")
    for path in sort(readdir(dir; join = true))
        endswith(path, ".EOF") || continue
        m = match(r"_V(\d{8}T\d{6})_(\d{8}T\d{6})\.EOF$", basename(path))
        m === nothing && continue
        a = DateTime(m.captures[1], dateformat"yyyymmddTHHMMSS")
        b = DateTime(m.captures[2], dateformat"yyyymmddTHHMMSS")
        a <= t <= b && return path
    end
    throw(ErrorException(
        "no orbit file in $dir covers $t. The container downloads one per granule beside its " *
        "outputs; a pruned or partial run has none, and geogrid cannot place a swath without it."))
end

"""
    AsfSwaths(granule, polarization, orbit_path)

An acquisition whose subswaths are read from ASF's burst extractor.

Two fetches per subswath and no granule transfer: the first burst's metadata file carries the whole
subswath's annotation, including how many bursts it has, and the second call is served from the same
cached file. The rasters are fetched only when read, which the geometry never does.
"""
struct AsfSwaths
    granule::String
    polarization::String
    orbit::String
end

"""
    SafeSwaths(product::Sentinel1Product)

An acquisition whose subswaths are read from an already-parsed local SAFE.

This is the route a burst job takes. Its granule name is `burst2safe`'s own — synthesized from the
requested bursts, with a checksum suffix that is not the one the product name carries — so it names
nothing at ASF and the annotations have to come off the container the container built.
"""
struct SafeSwaths
    product::Sentinel1Product
end

"""
    swath_annotation(src, swath) -> SubswathAnnotation

`swath`'s annotation, which is every number the merged grid is derived from bar the orbit.
"""
swath_annotation(src::AsfSwaths, swath::Integer) =
    only(asf_bursts(src.granule, swath, src.polarization, 1:1; orbit = src.orbit)).backend.annotation
swath_annotation(src::SafeSwaths, swath::Integer) = annotation(src.product, swath)

"""
    merged_swath(src, swath) -> SLC

One subswath of the acquisition, its bursts merged.

Wanted for what an annotation does not carry: the orbit, its epoch and the look side.
"""
function merged_swath(src::AsfSwaths, swath::Integer)
    n = nbursts(swath_annotation(src, swath))
    return merge_bursts(asf_bursts(src.granule, swath, src.polarization, 1:n; orbit = src.orbit))
end
merged_swath(src::SafeSwaths, swath::Integer) = merge_bursts(collect(bursts(src.product, swath)))

"""
    s1_mosaic(src, swaths = S1_SWATHS) -> (RadarCoordinate, UtcTime)

The merged radar grid geogrid is handed for the acquisition `src` reaches.

Built from the near-range subswath and then widened to the mosaic, which is what `loadMetadataSlc`
does: the range origin, sample spacing, PRF, wavelength and orbit are the near subswath's, the sample
count spans to the far subswath's far edge, and the line count comes from the sensing interval across
*all* subswaths rather than from any one of them — subswaths are acquired at slightly different azimuth
times, so the union is wider than each.

The incidence angle is recomputed rather than carried over, because it is the angle at the scene centre
and the centre moves when the range extent widens from one subswath to three.

The absolute sensing start comes back alongside, because the pair's interval is the difference of two
acquisitions' and each acquisition's own `sensing_start` is relative to its own orbit epoch.
"""
function s1_mosaic(src::Union{AsfSwaths,SafeSwaths}, swaths = S1_SWATHS; grid = nothing)
    sws = collect(swaths)
    ann = [swath_annotation(src, sw) for sw in sws]
    # **Both extents follow the burst shape `merge_bursts_in_swath` reads, which is the CSLC's.** It takes
    # `num_az_samples, num_rng_samples` off the first burst's resampled raster (`s1_isce3.py:565-566`), not
    # off the annotation, and the two coincide only when the reference burst was written by `rdr2geo`. A
    # pair whose reference took the resample path instead lands on COMPASS's own grid — 136 to 195 lines
    # taller and 72 to 97 samples narrower — and every extent below inherits that.
    shape = [grid === nothing ? (a.lines_per_burst, a.samples_per_burst) : grid(sw)
             for (a, sw) in zip(ann, sws)]

    # **The mosaic's extents are `merge_swaths`'s, not a subswath's.** `loadMetadataSlc:210` takes both
    # `numberOfLines` and `numberOfSamples` straight from the merged shape when it is given one — the
    # closed-form width at `:197-200` is only the fallback — so the grid geogrid is told about is the
    # mosaic hyp3 laid out, and both extents have to be reproduced from its own arithmetic.
    ilo = argmin(i -> ann[i].starting_range, eachindex(ann))
    ihi = argmax(i -> ann[i].starting_range, eachindex(ann))
    lo, hi = ann[ilo], ann[ihi]
    dr = lo.range_pixel_spacing
    adt = lo.azimuth_time_interval

    # Per subswath, the merged azimuth length: bursts are placed at their valid-line offsets and the
    # merge runs from the first burst's start to the last burst's start plus one burst
    # (`merge_bursts_in_swath:575-579`).
    swath_lines = [1 + round(Int, (seconds_between(first(a.burst_start), last(a.burst_start)) +
                                   (shape[k][1] - 1) * a.azimuth_time_interval) /
                                  a.azimuth_time_interval) for (k, a) in enumerate(ann)]

    # **And then the swath stack adds that whole merged length to the *last* burst's start again**
    # (`merge_swaths:437-439`: `burst_sensing_stop = ref_bursts[-1].sensing_start + burst_length`, where
    # `burst_length` spans `burst_az_samples`, the merged swath rather than one burst). So the mosaic
    # reaches about 1.7 times a subswath's height — 15858 rows against 9145 lines of acquisition on
    # `S1A_IW_SLC__1SSH_20150828`. Reproduced because it is the extent geogrid bounds its azimuth index
    # against: the reference's own `window_location` band 2 runs 0..15857.
    start = minimum(a -> first(a.burst_start), ann)
    span = maximum(zip(ann, swath_lines)) do (a, n)
        seconds_between(start, last(a.burst_start)) + (n - 1) * a.azimuth_time_interval
    end
    nlines = 1 + round(Int, span / adt)

    # `floor` here, against the fallback formula's `round`, and the far subswath's own sample count.
    nsamples = shape[ihi][2] +
               floor(Int, (hi.starting_range - lo.starting_range) / dr)

    prf = 1 / adt

    # The near subswath's merged acquisition supplies what an annotation does not carry: the orbit, its
    # epoch, and the look side.
    base = merged_swath(src, sws[ilo])
    c = RadarCoordinate(base)
    base_start = first(_annotation_of(base).burst_start)

    kwargs = (; orbit = c.orbit, starting_range = lo.starting_range, dr,
              # `c.sensing_start` is seconds against the orbit epoch and corresponds to the near
              # subswath's first burst, so the shift to the earliest subswath is added there.
              sensing_start = c.sensing_start + seconds_between(base_start, start),
              prf, nsamples, nlines, look_side = c.look_side, wavelength = lo.wavelength,
              orbit_epoch_offset = c.orbit_epoch_offset)
    coord = RadarCoordinate(; kwargs..., incidence_angle = incidence_angle(; kwargs...))
    return coord, start
end

_annotation_of(s) = s.backend.annotation

"""
    s1_pair(c::GoldenCase, run) -> CoregisteredPair

`c`'s pair as geogrid receives it: the reference acquisition's geometry, and an interval.

The secondary contributes its sensing start and nothing else, since `process_slc` copies `meta_r` and
replaces only the sensing times. So the secondary's own grid — which is where the coregistration lives —
is deliberately absent from this, and a rung that compares the geogrid bands is comparing the reference
acquisition's geometry alone.

`dt` is the full-precision difference of sensing starts (`testGeogrid.py:439`). Using the optical path's
whole-day convention here would be wrong by up to half a day.
"""
function s1_pair(c::GoldenCase, run::AbstractString)
    # Acquisition order, not the job's: the products report a positive `date_dt` with the earlier scene as
    # `img1` even on the two jobs whose reference is the later acquisition, so the pipeline reorders the
    # pair before geogrid sees it exactly as it does on the optical path.
    rg, sg = acquisition_order(c)
    src(g) = AsfSwaths(g, s1_polarization(g), s1_orbit(run, g))
    return _s1_pair(src(rg), src(sg), S1_SWATHS; grid = cslc_grid(run, S1_SWATHS))
end

"""
    s1_burst_pair(c::GoldenCase, run) -> CoregisteredPair

`c`'s pair for a burst job, which reaches `process_slc` over a synthesized SAFE.

A burst job is the full-SLC path with two substitutions and no third
(`s1_isce3.process_sentinel1_burst_isce3:55-73`). `burst2safe` assembles the requested bursts into a
SAFE, and `process_slc` then runs on it unchanged — so `merge_swaths` mosaics whatever bursts that
container holds, and the merged extents follow from its annotations by the same arithmetic as a full
granule's. The substitutions are:

  * **The container is local.** Its name is `burst2safe`'s, not a granule ASF would serve, so the
    annotations are read from the SAFE the run directory holds rather than fetched.
  * **The subswath set is the bursts'.** `swaths = sorted(set(int(g.split('_')[2][2]) for g in
    reference))` (`:58`), so a job over `IW1` bursts alone mosaics one subswath and the range origin,
    width and incidence angle are that subswath's rather than the three-swath union's.

The two SAFEs are matched to the pair by acquisition time rather than by name: `burst2safe` stamps its
own checksum suffix, which does not agree with the one in the product name.
"""
function s1_burst_pair(c::GoldenCase, run::AbstractString)
    early, late = _burst_safes(run)
    pol(safe) = s1_polarization(basename(safe))
    src(safe) = SafeSwaths(Sentinel1Product(safe; orbit = s1_orbit(run, basename(safe)),
                                            polarization = lowercase(pol(safe)),
                                            swaths = burst_swaths(c)))
    return _s1_pair(src(early), src(late), burst_swaths(c); grid = cslc_grid(run, burst_swaths(c)))
end

# Shared by both routes, because only the annotation source and the subswath set differ: the reference
# acquisition's merged grid, and the interval between the two sensing starts.
function _s1_pair(ref_src, sec_src, swaths; grid = nothing)
    ref, ref_start = s1_mosaic(ref_src, swaths; grid)

    # Only the secondary's sensing start is wanted, so its geometry is built and its grid discarded — and
    # with it the question of which acquisition `product/` holds. The mosaic's shape is one thing rather
    # than one per acquisition, so the grid belongs to the coordinate that is kept.
    _, sec_start = s1_mosaic(sec_src, swaths)

    return CoregisteredPair(ref; dt = seconds_between(ref_start, sec_start))
end

"""
    burst_swaths(c::GoldenCase) -> Vector{Int}

The subswaths a burst job mosaics, from its reference burst list.

`s1_isce3.py:58` reads them out of the burst names — `S1_105602_IW2_...` contributes 2 — and off the
reference list alone, so a pair whose secondary reached a subswath the reference did not still
processes only the reference's.
"""
function burst_swaths(c::GoldenCase)
    # A full-SLC granule names the mode without a subswath — `S1A_IW_SLC__1SSH_...` — and
    # `process_slc` defaults to all three. A burst names exactly one.
    all(g -> occursin("_IW_SLC_", g), c.reference) && return collect(S1_SWATHS)
    sw = Int[]
    for g in c.reference
        m = match(r"_IW(\d)_", g)
        m === nothing && throw(ArgumentError(
            "\"$g\" names neither a full-SLC product nor an IW subswath; a burst is named " *
            "`S1_<id>_IW<n>_<time>_<pol>_<hash>-BURST` and an SLC `S1A_IW_SLC__...`"))
        push!(sw, parse(Int, m.captures[1]))
    end
    return sort!(unique!(sw))
end

# The pair's two SAFEs, earliest first. A burst run holds exactly two, both written by `burst2safe`.
function _burst_safes(run::AbstractString)
    safes = filter(n -> endswith(n, ".SAFE") && isdir(joinpath(run, n)), readdir(run))
    length(safes) == 2 || error(
        "expected the two SAFE products `burst2safe` assembles in $run, found $(length(safes)). " *
        "A burst job's annotations are read from them; a pruned run has none.")
    sort!(safes; by = n -> split(n, '_')[6])
    return joinpath(run, safes[1]), joinpath(run, safes[2])
end

"""
    stage_safe(granule; search = ()) -> String

`granule`'s expanded `.SAFE` tree, downloading it from ASF only if no copy is already on disk.

**Unzipped, because the raster cannot be read in place.** A Sentinel-1 zip stores its measurement TIFFs
deflated, so a line is not addressable without inflating everything before it — `SLCDatasets` says so
rather than reading part of one.

`search` names directories to consult before the granule cache, because the reference's own driver
downloads a full-SLC pair into the run it processes: an acquisition is 2.5 GB zipped and 3.6 GB expanded,
so fetching a second copy of a tree that is already there is 6 GB of disk and a quarter hour of egress for
nothing. A tree found under `search` is returned where it lies and a zip found there is expanded into the
granule cache, both left in place — a run directory holds reference artifacts and is not this function's
to edit.

A zip this function downloads is deleted once expanded: only the tree is read. ASF's download is a URS
redirect chain, so `curl` carries `~/.netrc` through it — the same credential `fetch.jl --check` verifies
and the container mounts — **and a cookie jar**, because URS hands the session back as a cookie and a
`curl` that discards it re-authenticates at every hop until it hits the redirect limit. The symptom is
`curl: (47) Maximum (50) redirects followed` on a URL whose `HEAD` succeeds.
"""
function stage_safe(granule::AbstractString; search = String[])
    dirs = collect(String, search)
    dir = joinpath(CACHE, "granules")
    for d in [dirs; dir]
        isdir(joinpath(d, granule * ".SAFE")) && return joinpath(d, granule * ".SAFE")
    end
    mkpath(dir)
    safe = joinpath(dir, granule * ".SAFE")
    i = findfirst(d -> isfile(joinpath(d, granule * ".zip")), dirs)
    staged = !isnothing(i)
    zip = staged ? joinpath(dirs[i], granule * ".zip") : joinpath(dir, granule * ".zip")
    if !staged && !isfile(zip)
        mission = granule[1:3]                  # `S1A`, `S1B` or `S1C`
        url = "https://sentinel1.asf.alaska.edu/SLC/S$(mission[3])/$(granule).zip"
        tmp = zip * ".partial"
        @info "downloading" granule url
        mktemp() do jar, _
            run(`curl -sS -n -L -f -c $jar -b $jar -o $tmp $url`)
        end
        mv(tmp, zip; force = true)
    end
    @info "expanding" zip = basename(zip)
    run(`unzip -q -o $zip -d $dir`)
    isdir(safe) || error("$(basename(zip)) expanded to no $(basename(safe)); ASF's tree is named " *
                         "`<granule>.SAFE`")
    staged || rm(zip)
    return safe
end

"""
    nisar_l1_pair(c::GoldenCase, run) -> CoregisteredPair

A NISAR L1 RSLC pair, from the two products the run holds.

Much simpler than Sentinel-1 and for one reason: an RSLC is a single acquisition on a single radar grid,
so there is nothing to merge and nothing to mosaic. `loadMetadataRslc` (`testGeogrid.py:240-260`) reads
the zero-Doppler start, the dimensions and the orbit straight off the product, and the orbit travels
*inside* it rather than in a separate `.EOF` — so no orbit file is resolved here.

`dt` is the full-precision difference of the two zero-Doppler starts, which is what `runGeogrid`'s radar
branch takes (`:439`). `CoregisteredPair(::SLC, ::SLC)` computes exactly that.

The products are read from the run directory rather than fetched: they are 11.2 GiB each and the driver
already downloaded them there.
"""
function nisar_l1_pair(c::GoldenCase, run::AbstractString)
    early, late = acquisition_order(c)
    path(name) = begin
        p = joinpath(run, name * ".h5")
        isfile(p) || error("$(name).h5 is not in $run. A NISAR L1 pair is read from the products the " *
                           "driver downloaded there; each is 11.2 GiB.")
        p
    end
    return CoregisteredPair(open_slc(path(early)), open_slc(path(late)))
end

# ---------------------------------------------------------------------------
# Rung 5.2 — the merged radar mosaic the correlator is handed
# ---------------------------------------------------------------------------
#
# `merge_swaths` (`s1_isce3.py:393-530`) is the whole specification, and it is index arithmetic rather
# than signal processing: `read_slc_gdal` takes `np.abs` of each burst raster on the way in, so every
# array downstream is `Float32` amplitude and nothing is resampled. Two nested layouts:
#
#   1. **Bursts into a subswath** (`merge_bursts_in_swath`). Bursts overlap in azimuth, and the seam is
#      put *halfway through* each overlap so no burst contributes its resampling margin.
#   2. **Subswaths into the mosaic** (`merge_swaths`). Each subswath is laid on the near-range one's
#      range origin, and the writer is first-come: a pixel already non-zero is not overwritten.
#
# **The reference burst needs no coregistration.** COMPASS writes the reference burst on its own grid and
# deramping is phase-only, so `abs` of the CSLC is `abs` of the raw burst — measured on
# `S1C_IW_SLC__1SSV_20250416`'s first burst at a median difference of 0.0057 on amplitudes near 200, or
# 3e-5 relative. So the reference side of rung 5.2 is reproducible from the SAFE with no COMPASS run.

"""
    swath_amplitude(p::Sentinel1Product, swath) -> (Matrix{Float32}, Int, Int)

One subswath's bursts merged into a single amplitude raster, with its azimuth and range extents.

Reproduces `merge_bursts_in_swath`. The azimuth seam between two bursts is placed halfway through their
overlap, which is what keeps each burst's resampling margin out of the result, and the range window is
the *first* burst's valid sample range applied to every burst.

A single-burst subswath takes the early return: the burst is written whole, with no valid-region
cropping at all, so its extents are the burst's own.
"""
function swath_amplitude(p::Sentinel1Product, swath::Integer)
    a = annotation(p, swath)
    b = collect(bursts(p, swath))
    n = length(b)
    lpb, spb = a.lines_per_burst, a.samples_per_burst
    dt = a.azimuth_time_interval
    # One handle for the whole subswath: a `.SAFE` stacks every burst of a subswath in one raster, so a
    # burst is a row range of it. `read_pixels` on a single burst cannot be used — it is
    # `first(burst_raster(b))` on a `BurstRaster`, which is not iterable.
    raster = burst_raster(first(b).backend).raster
    amp(i, rows, cols) = abs.(raster[((i - 1) * lpb) .+ rows, cols])

    # The annotation's valid bounds are 1-based inclusive; every index below is the reference's 0-based.
    fvl = a.first_valid_line .- 1
    lvl = a.last_valid_line .- 1
    fvs = a.first_valid_sample .- 1
    lvs = a.last_valid_sample .- 1

    n == 1 && return (Float32.(amp(1, 1:lpb, 1:spb)), lpb, spb)

    # `get_azimuth_reference_offsets`: where each burst's valid region starts and ends in the merged
    # subswath, from its own sensing time and first valid line.
    lims = map(1:n) do i
        s = round(Int, (seconds_between(first(a.burst_start), a.burst_start[i]) + fvl[i] * dt) / dt)
        (s, s + (lvl[i] - fvl[i]) + 1)
    end
    nlines = 1 + round(Int, (seconds_between(first(a.burst_start), last(a.burst_start)) +
                             (lpb - 1) * dt) / dt)
    out = zeros(Float32, nlines, spb)
    for i in 1:n
        # `//` in the reference is floor division and these overlaps are positive, but `fld` says so.
        prev = i > 1 ? fld(lims[i - 1][2] - lims[i][1], 2) : 0
        nxt = i < n ? fld(lims[i][2] - lims[i + 1][1], 2) : 0
        bstart, bend = fvl[i] + prev, 1 + lvl[i] - nxt
        mstart, mend = lims[i][1] + prev, lims[i][2] - nxt
        cols = (fvs[i] + 1):lvs[i]          # `slice(first_valid_sample, last_valid_sample)`
        out[(mstart + 1):mend, cols] = amp(i, (bstart + 1):bend, cols)
    end
    return (out, nlines, spb)
end

"""
    radar_mosaic(p::Sentinel1Product, swaths) -> Matrix{Float32}

The merged amplitude raster `merge_swaths` writes as `reference.tif`.

Each subswath is placed at its own azimuth offset from the earliest sensing start and at its range offset
from the near subswath's origin, and only where the mosaic is still zero — the reference's writer is
first-come, which matters because adjacent subswaths overlap in range.

Two quirks of the extent are reproduced rather than corrected, both already established by
[`s1_mosaic`](@ref): the azimuth span adds a whole merged subswath length to the *last* burst's start, so
the mosaic reaches about 1.7 times a subswath's height; and the width is the **last** subswath's sample
count plus the floored range offset to it, not a union of the three.

A subswath other than the far one is trimmed by 64 samples at its far edge, which is the reference's
`invalid_pixel_buffer` — the resampling margin at a subswath's far range.

`offsets` maps a subswath to a per-burst `(dl, ds)` pair, for a reference that was itself coregistered
and resampled rather than copied; see [`resampled_burst_dirs`](@ref). Each subswath then goes through
[`secondary_swath_amplitude`](@ref) with the product on both sides, since the resampling is the same and
only the source acquisition differs.
"""
function radar_mosaic(p::Sentinel1Product, swaths; offsets = nothing, grid = nothing)
    sws = collect(swaths)
    merged = [offsets === nothing ? swath_amplitude(p, sw) :
              secondary_swath_amplitude(p, p, sw, nothing; offsets = offsets(sw),
                                        grid = grid === nothing ? nothing : grid(sw))
              for sw in sws]
    return _stack_swaths(p, sws, merged)
end

"""
    burst_offsets(dir) -> (dl, ds)

COMPASS's per-pixel azimuth and range offsets for one burst, as callables on a 0-based `(line, sample)`.

`Geo2Rdr` marks a pixel it could not solve with -1e6. Those are moved far outside the burst rather than
passed through, so [`resample_burst`](@ref)'s own bounds guard drops them.
"""
function burst_offsets(dir::AbstractString)
    read_off(name) = permutedims(ArchGDAL.read(ArchGDAL.getband(
                                     ArchGDAL.read(joinpath(dir, name)), 1)))
    lookup(A) = (l, s) -> (v = A[l + 1, s + 1]; v <= -1e5 ? -1e7 : v)
    # An ISCE flat file is described by a sibling `.off.xml`, which GDAL only finds if it may list the
    # directory. `GDAL_DISABLE_READDIR_ON_OPEN` is `EMPTY_DIR` for the S3 reads elsewhere in these tools,
    # so it is lifted for these two opens and put back.
    #
    # **Under a lock, because `CPLSetConfigOption` is process-wide and these reads are concurrent.** A
    # subswath's placements are built one task per burst, so without one a task's restore lands while another
    # is mid-open and GDAL reports `range.off` as "not recognized as being in a supported file format" —
    # measured at 3 failures in 10 concurrent reads where 10 serial reads all succeed.
    #
    # A lock rather than `CPLSetThreadLocalConfigOption`: the thread-local getter returns a pointer into
    # GDAL's own storage, which the following set invalidates, so restoring the previous value means handing
    # back a dangling pointer. Serializing costs little here — offsets are read only to replay the
    # reference's own resample, never on the path a pair is actually processed by.
    return @lock READDIR_OPTION begin
        prev = ArchGDAL.getconfigoption("GDAL_DISABLE_READDIR_ON_OPEN", "")
        ArchGDAL.setconfigoption("GDAL_DISABLE_READDIR_ON_OPEN", "NO")
        try
            (lookup(read_off("azimuth.off")), lookup(read_off("range.off")))
        finally
            ArchGDAL.setconfigoption("GDAL_DISABLE_READDIR_ON_OPEN", prev)
        end
    end
end

# Guards every process-wide flip of `GDAL_DISABLE_READDIR_ON_OPEN`; see [`burst_offsets`](@ref).
const READDIR_OPTION = ReentrantLock()

"""
    cslc_grid(run, swaths) -> (swath -> (lines, samples)) or `nothing`

The shape `merge_bursts_in_swath` reads per subswath: its **first** burst's resampled raster, since that is
the one whose `shape` it takes `num_az_samples` and `num_rng_samples` from and then applies to every burst
of the subswath.

`nothing` when the run kept no `product/`, or when any subswath has no resampled burst in it — a reference
written by `rdr2geo` has no offsets and its CSLC is the annotation's burst, so the annotation is already
right there.
"""
function cslc_grid(run::AbstractString, swaths)
    isdir(joinpath(run, "product")) || return nothing
    dirs = Dict(sw => resampled_burst_dirs(run, sw) for sw in swaths)
    any(isempty, values(dirs)) && return nothing
    shapes = Dict(sw => burst_offset_grid(first(dirs[sw])) for sw in swaths)
    return sw -> shapes[sw]
end

"""
    burst_offset_grid(dir) -> (lines, samples)

The grid COMPASS resampled one burst onto, read from `azimuth.off`'s header alone.

**This is not the annotation's burst, and `merge_bursts_in_swath` uses this one.** It takes
`num_az_samples, num_rng_samples` from the CSLC raster's own shape (`s1_isce3.py:565-566`) and then applies
the *annotation's* `first_valid_line` and `first_valid_sample` as windows into it. On the burst jobs the two
grids coincide, so the distinction is invisible; on the full-SLC pairs that take the resample path they do
not, and every extent in the mosaic depends on this one:

    subswath   annotation           CSLC and offsets
    IW1        1504 x 21530         1640 x 21458
    IW2        1515 x 25376         1696 x 25279
    IW3        1520 x 24454         1715 x 24369

Dimensions only, so this costs a header read rather than the quarter-gigabyte the offsets themselves are.
"""
function burst_offset_grid(dir::AbstractString)
    # `READDIR_OPTION` for the same reason [`burst_offsets`](@ref) holds it: this is called once per burst
    # from the same per-burst tasks, and the option it flips is process-wide.
    return @lock READDIR_OPTION begin
        prev = ArchGDAL.getconfigoption("GDAL_DISABLE_READDIR_ON_OPEN", "")
        ArchGDAL.setconfigoption("GDAL_DISABLE_READDIR_ON_OPEN", "NO")
        try
            ds = ArchGDAL.read(joinpath(dir, "azimuth.off"))
            (ArchGDAL.height(ds), ArchGDAL.width(ds))
        finally
            ArchGDAL.setconfigoption("GDAL_DISABLE_READDIR_ON_OPEN", prev)
        end
    end
end

"""
    resampled_burst_dirs(run, swath; secondary = false) -> Vector{String}

The COMPASS burst directories of one subswath that hold `azimuth.off` and `range.off`, in burst order,
or empty when the acquisition was not resampled.

Whether the *reference* acquisition has them is the whole question for a mosaic. `write_yaml` sets
`bool_reference` false whenever a cached static topographic layer is available, so a reference burst is
then coregistered against that layer and resampled rather than written out; without one it runs `rdr2geo`
instead and its directory holds `x`/`y`/`z` and a polarization-named SLC with no offsets at all. So the
presence of the offsets is both the test for which path ran and the input needed to follow it.

The burst identifier encodes the absolute burst number along the track, so sorting the directory names
puts them in acquisition order.
"""
function resampled_burst_dirs(run::AbstractString, swath::Integer; secondary::Bool = false)
    root = joinpath(run, secondary ? "product_sec" : "product")
    isdir(root) || return String[]
    dirs = String[]
    for id in sort!(filter(endswith("_iw$(swath)"), readdir(root)))
        sub = joinpath(root, id)
        for date in sort!(readdir(sub))
            d = joinpath(sub, date)
            isdir(d) && isfile(joinpath(d, "azimuth.off")) && (push!(dirs, d); break)
        end
    end
    return dirs
end

# `merge_swaths`'s swath stack, over per-subswath rasters a caller has already merged. Shared by the
# reference and the secondary because the layout is the reference's in both cases.
"""
    SwathPlacement

Which window of one merged subswath goes where in the mosaic `merge_swaths` builds.

**The azimuth window is cropped and moved; the range window is cropped and left where it is.** A
subswath's first `first_valid_line` rows are dropped and what remains starts at that subswath's azimuth
offset — zero for the first, which is why a subswath and the mosaic differ by `first_valid_line` rows —
while the range window keeps its own column indices. The asymmetry is `merge_swaths`'s
(`s1_isce3.py:487-500`), not a simplification.
"""
struct SwathPlacement
    swath::Int
    mosaic_rows::UnitRange{Int}
    mosaic_cols::UnitRange{Int}
    swath_rows::UnitRange{Int}
    swath_cols::UnitRange{Int}
end

"""
    _mosaic_layout(p, swaths, shapes) -> (Vector{SwathPlacement}, Tuple{Int,Int})

Where each merged subswath lands in the mosaic, and the mosaic's size.

`shapes[k]` is subswath `k`'s `(rows, columns)` as [`swath_amplitude`](@ref) returned it. Separated from
the filling so a lazily resampled mosaic and a materialized one cannot disagree about an index.
"""
function _mosaic_layout(p::Sentinel1Product, swaths, shapes)
    sws = collect(swaths)
    ann = [annotation(p, sw) for sw in sws]

    dr = first(ann).range_pixel_spacing
    dt = first(ann).azimuth_time_interval
    starts = [first(a.burst_start) for a in ann]
    sensing_start = minimum(starts)
    # `burst_sensing_stop` spans the *merged* subswath rather than one burst, which is the 1.7x quirk.
    span = maximum(zip(ann, shapes)) do (a, sh)
        seconds_between(sensing_start, last(a.burst_start)) + (sh[1] - 1) * a.azimuth_time_interval
    end
    total_az = 1 + round(Int, span / dt)

    # `rng_offsets` are measured from the *first* subswath in the list, and `last_rng_samples` ends as the
    # last subswath's own width — not the widest.
    rng_offsets = [sw == first(sws) ? 0 :
                   floor(Int, (annotation(p, sw).starting_range - first(ann).starting_range) / dr)
                   for sw in sws]
    total_rng = last(shapes)[2] +
                floor(Int, (last(ann).starting_range - first(ann).starting_range) / dr)

    places = map(eachindex(sws)) do k
        a, nrows = ann[k], shapes[k][1]
        fvl, lvl = a.first_valid_line[1] - 1, a.last_valid_line[1] - 1
        fvs, lvs = a.first_valid_sample[1] - 1, a.last_valid_sample[1] - 1
        az_offset = floor(Int, seconds_between(sensing_start, starts[k]) / dt)
        rng_offset = rng_offsets[k] + fvs
        rng_end = rng_offset + (lvs - fvs)
        buffer = sws[k] == maximum(sws) ? 0 : 64
        # A single-burst subswath was written whole, so its azimuth window is the first burst's valid
        # region; a merged one runs to the *last* burst's last valid line counted from the end.
        slc_az_end = nbursts(a) > 1 ? nrows - (lvl_last(a)) : lvl
        merged_az_end = nbursts(a) > 1 ? az_offset + nrows - lvl_last(a) - fvl : az_offset + (lvl - fvl)
        SwathPlacement(sws[k], (az_offset + 1):merged_az_end, (rng_offset + 1):(rng_end - buffer),
                       (fvl + 1):slc_az_end, (fvs + 1):(lvs - buffer))
    end
    return (places, (total_az, total_rng))
end

function _stack_swaths(p::Sentinel1Product, swaths, merged)
    places, dims = _mosaic_layout(p, swaths, [(m[2], m[3]) for m in merged])
    out = zeros(Float32, dims)
    for (k, q) in enumerate(places)
        dst = view(out, q.mosaic_rows, q.mosaic_cols)
        src = view(merged[k][1], q.swath_rows, q.swath_cols)
        for i in eachindex(dst, src)
            # First-come: `cond = merged == 0 & slc != 0`.
            (dst[i] == 0 && src[i] != 0) && (dst[i] = src[i])
        end
    end
    return out
end

"""
    ResampledMosaic(rp, sp, swaths, dem; offsets = nothing, grid = nothing) <: AbstractMatrix{Float32}

The coregistered secondary on the **mosaic** grid, resampled on demand — what a caller hands `autorift`.

[`ResampledSwath`](@ref) is one subswath, which is the reference's `sec_swath_iw<n>.tif`; this is the
`merge_swaths` step on top, so its indices are the ones the geogrid's `window_*` rasters describe. Reading
a window resamples only the bursts that window touches, so nothing holds the mosaic and nothing writes it
— which is the whole of why `secondary.tif` exists.

Equal to `_stack_swaths` of the materialized subswaths, window for window, because both read their indices
from [`_mosaic_layout`](@ref).
"""
struct ResampledMosaic{S} <: AbstractMatrix{Float32}
    swaths::Vector{S}
    places::Vector{SwathPlacement}
    dims::Tuple{Int,Int}
end

function ResampledMosaic(rp::Sentinel1Product, sp::Sentinel1Product, swaths, dem;
                         offsets = nothing, grid = nothing)
    sws = collect(swaths)
    lazy = [ResampledSwath(rp, sp, sw, dem;
                           offsets = offsets === nothing ? nothing : offsets(sw),
                           grid = grid === nothing ? nothing : grid(sw)) for sw in sws]
    places, dims = _mosaic_layout(rp, sws, [size(r) for r in lazy])
    return ResampledMosaic(lazy, places, dims)
end

Base.size(m::ResampledMosaic) = m.dims

function Base.getindex(m::ResampledMosaic, rows::AbstractUnitRange, cols::AbstractUnitRange)
    checkbounds(m, rows, cols)
    out = zeros(Float32, length(rows), length(cols))
    for (k, q) in enumerate(m.places)
        mr = intersect(rows, q.mosaic_rows)
        isempty(mr) && continue
        mc = intersect(cols, q.mosaic_cols)
        isempty(mc) && continue
        # Positional within the placement, as in `_stack_swaths`: the mosaic window and the subswath
        # window have equal lengths by construction.
        r0 = first(q.swath_rows) + (first(mr) - first(q.mosaic_rows))
        c0 = first(q.swath_cols) + (first(mc) - first(q.mosaic_cols))
        got = m.swaths[k][r0:(r0 + length(mr) - 1), c0:(c0 + length(mc) - 1)]
        dst = view(out, (first(mr) - first(rows) + 1):(last(mr) - first(rows) + 1),
                   (first(mc) - first(cols) + 1):(last(mc) - first(cols) + 1))
        for i in eachindex(dst, got)
            # First-come, in subswath order, which is the rule and the order `_stack_swaths` applies.
            (dst[i] == 0 && got[i] != 0) && (dst[i] = got[i])
        end
    end
    return out
end

AutoRIFT.ondisk(::ResampledMosaic) = true

Base.getindex(m::ResampledMosaic, i::Integer, j::Integer) = m[i:i, j:j][1, 1]
Base.getindex(m::ResampledMosaic, rows::AbstractUnitRange, j::Integer) = m[rows, j:j][:, 1]
Base.getindex(m::ResampledMosaic, i::Integer, cols::AbstractUnitRange) = m[i:i, cols][1, :]
Base.getindex(m::ResampledMosaic, ::Colon, ::Colon) = m[axes(m, 1), axes(m, 2)]

"""
    secondary_mosaic(rp, sp, swaths, dem; offsets = nothing, grid = nothing) -> Matrix{Float32}

[`radar_mosaic`](@ref) for the secondary: each subswath resampled onto the reference's, then stacked.

The materialized counterpart of [`ResampledMosaic`](@ref), and what `secondary.tif` holds.
"""
function secondary_mosaic(rp::Sentinel1Product, sp::Sentinel1Product, swaths, dem;
                          offsets = nothing, grid = nothing)
    sws = collect(swaths)
    merged = [secondary_swath_amplitude(rp, sp, sw, dem;
                                        offsets = offsets === nothing ? nothing : offsets(sw),
                                        grid = grid === nothing ? nothing : grid(sw)) for sw in sws]
    return _stack_swaths(rp, sws, merged)
end

# `bursts[-1].last_valid_line` of a subswath, 0-based. Named because the reference reaches it as a
# negative Python index, `slice(first_valid_line, -last_valid_line)`, which counts from the end.
lvl_last(a) = a.last_valid_line[end] - 1

"""
    dem_sampler(path) -> Function

Bilinear terrain height from a geographic DEM, as `(lon_degrees, lat_degrees) -> height`.

The DEM the container downloaded for the pair, which sits in the run directory as `dem.tif` on an
EPSG:4326 grid. Read whole rather than windowed: it is a few thousand cells on a side, and the
coregistration solve asks for scattered points rather than a block.
"""
function dem_sampler(path::AbstractString)
    ds = ArchGDAL.read(path)
    gt = ArchGDAL.getgeotransform(ds)
    z = ArchGDAL.read(ArchGDAL.getband(ds, 1))
    nx, ny = size(z)
    return function (lon_d, lat_d)
        px = (lon_d - gt[1]) / gt[2]
        py = (lat_d - gt[4]) / gt[6]
        i = clamp(floor(Int, px), 0, nx - 2)
        j = clamp(floor(Int, py), 0, ny - 2)
        fx, fy = px - i, py - j
        at(p, q) = Float64(z[p + 1, q + 1])
        return (1 - fx) * (1 - fy) * at(i, j) + fx * (1 - fy) * at(i + 1, j) +
               (1 - fx) * fy * at(i, j + 1) + fx * fy * at(i + 1, j + 1)
    end
end

"""
    coregistration_offset(cr, cs, line, sample, height; iters = 4) -> (dline, dsample, h)

Where the secondary images the ground point the reference images at `(line, sample)`, as an offset in
the reference's own pixels.

The orbit-driven coregistration, and the geometry half of rung 5.2's secondary side: `rdr2geo` on the
reference to reach the ground, `geo2rdr` on the secondary to come back. `line` and `sample` are
zero-based, as the reference's indices are.

**The terrain enters as an outer fixed point.** `rdr2geo` takes a constant height — all its callers in
the geogrid supply one — so the DEM is iterated: solve at the current height, look the DEM up at the
resulting position, solve again. Four passes, which is past convergence for Sentinel-1 geometry.

**The offsets are computed from the solved time and range, not through `azimuth_index`.** Those helpers
round to a whole line and sample for the geogrid's benefit, which is exactly the sub-pixel part a
resampler needs.

# The one line, and where it comes from

**`Geo2Rdr` does not define the azimuth offset as `(t - t0) * prf - line`.** It defines it one line
lower, so that is subtracted here. Measured against ISCE3 itself rather than inferred:
`tools/golden/isce_offsets.py` runs `Rdr2Geo` and `Geo2Rdr` with COMPASS's own arguments and reads the
`azimuth.off` and `range.off` the resampler consumes. Over twenty points spanning a burst of
`S1C_IW_SLC__1SSV_20250416`:

    range:   julia - isce3  mean -0.00000000   sd 8.0e-10
    azimuth: julia - isce3  mean +0.99999904   sd 5.1e-08   before this correction

So the geometry agrees with the reference implementation to eight decimal places on both axes and the
difference is a single exact constant. It is the *resampling position* that matters to a caller — the
input line a given output line reads from — and that is what this returns.

Confirmed independently against the imagery before ISCE3 was consulted, which is what said the residual
was real rather than a bookkeeping artifact: correlating `secondary.tif` against the *raw* secondary
burst locates the offset COMPASS actually used, since those are the same acquisition, and a parabola
through the peak over fifteen points gave +1.0162 +/- 0.0310 lines. The same correlation with the
*reference* on both sides peaks at `(0, 0)` with correlation 1.000, so the mosaic mapping is exact.
"""
function coregistration_offset(cr, cs, line::Integer, sample::Integer, height;
                               iters::Integer = 4)
    el = Ellipsoid()
    az = cr.sensing_start + line / cr.prf
    rg = cr.starting_range + sample * cr.dr
    h = 0.0
    llh = ImagePairGeometry.SVector{3,Float64}(0.0, 0.0, 0.0)
    for _ in 1:iters
        llh = ImagePairGeometry.rdr2geo(cr.orbit, el, az, rg; height = h,
                                       wavelength = cr.wavelength, side = cr.look_side)
        h = height(llh[1] / ImagePairGeometry.DEG2RAD, llh[2] / ImagePairGeometry.DEG2RAD)
    end
    xyz = ImagePairGeometry.lonlat_to_xyz(el,
              ImagePairGeometry.SVector{3,Float64}(llh[1], llh[2], h))
    pm, vm = ImagePairGeometry.interpolate(cs.orbit, ImagePairGeometry.orbit_midtime(cs))
    p = ImagePairGeometry.geo2rdr(cs.orbit, xyz, ImagePairGeometry.midtime(cs),
                                 ImagePairGeometry.orbit_midtime(cs), pm, vm)
    # The `- 1` is `Geo2Rdr`'s convention, measured against it; see above.
    return ((p.aztime - cs.sensing_start) * cs.prf - line - 1,
            (p.range - cs.starting_range) / cs.dr - sample, h)
end

# ---------------------------------------------------------------------------
# The secondary, coregistered onto the reference grid
# ---------------------------------------------------------------------------
#
# `ResampSlc` deramps each chip, interpolates with an eight-tap sinc, and reramps. Only the first two
# reach an amplitude: the reramp is a phase multiply on the finished value.
#
# **Every piece of this was chosen by measurement**, on burst 1 of `S1C_IW_SLC__1SSV_20250416` against
# `secondary.tif` over a 512 x 4096 window, reported as the ratio of means and the correlation:
#
#     no deramp, bicubic                 0.7535   0.73626
#     deramp with the opposite sign      0.7641   0.76199
#     deramp, bicubic                    0.9266   0.98640
#     deramp, eight-tap sinc unwindowed  1.2070   0.98365
#     deramp, eight-tap sinc + Hamming   1.0116   0.99963
#     deramp, eight-tap sinc + Hann      0.9976   0.99957
#
# Hann is also derivable rather than merely best: ISCE's `sinc_coef`, which `Sinc2dInterpolator` is built
# from, weights the sinc by `(1 - pedestal)/2 * cos(pi x / (ns/2)) + (1 + pedestal)/2`, and a pedestal of
# zero makes that exactly `0.5 + 0.5 cos(pi x / 4)` at `ns = 8`.

const CLIGHT = 299792458.0

"""
    range_poly(p::RangePolynomial, range) -> Float64

`p` evaluated at a slant `range` in metres.

**Not `p(range)`.** `RangePolynomial` stores `r0` in metres but keeps the annotation's coefficients in
powers of *slant-range time*, so calling it with a range mixes the two units and returns a number about
sixteen orders of magnitude out. The delta has to be the two-way time: at 838,557 m the FM rate is
−2224.66 Hz/s this way and −3.6e16 Hz/s the other.
"""
range_poly(p, range::Real) = evalpoly(2 * (Float64(range) - p.r0) / CLIGHT, p.coeffs)

"""
    TopsCarrier(slc, coord, lines_per_burst)

The azimuth carrier phase of one Sentinel-1 burst, as `(line, sample) -> radians`.

The TOPS beam sweep puts a quadratic azimuth phase on each burst whose instantaneous frequency reaches a
few kilohertz against a 486 Hz line rate, so the samples of any interpolation chip are aliased with
respect to each other and interpolating the raw complex signal cancels rather than sums — worth 25% of the
amplitude, measured. Removing it first is what `ResampSlc` does.

`s1reader.az_carrier_components` verbatim, and each of its three choices matters:

  * **azimuth time is measured from the middle line index**, `(line - lines_per_burst ÷ 2) * dt` with
    integer division — not from the burst's mid *time*, which is half a line away;
  * **the reference time is a difference**, `dc(r0)/fm(r0) - dc(r)/fm(r)`, so the quadratic's vertex
    carries a constant offset that the beam-centre crossing alone does not;
  * **there is no demodulation term.** The carrier is `pi * kt * (eta - eta_ref)^2` and nothing else.

`kt` is `ks / (1 - ks/ka)` where `ks = 2|v| * steering_rate / wavelength` at the burst mid.
"""
struct TopsCarrier
    fm::Any
    dc::Any
    ks::Float64
    dt::Float64
    r0::Float64
    dr::Float64
    mid_line::Int
    eta_ref0::Float64
end

# Whether the carrier carries the demodulation term as well as the quadratic. `s1reader` has only the
# quadratic; the flag exists because which one reproduces `secondary.tif` better is a measurement.


function TopsCarrier(slc, coord, lines_per_burst::Integer)
    dp = SLCDatasets.deramp_parameters(slc)
    _, v = ImagePairGeometry.interpolate(coord.orbit, ImagePairGeometry.orbit_midtime(coord))
    # The annotation gives the steering rate in degrees per second; `s1reader` holds it in radians.
    ks_rad = dp.azimuth_steering_rate * (abs(dp.azimuth_steering_rate) > 0.1 ? pi / 180 : 1.0)
    ks = ks_rad * 2 * sqrt(sum(abs2, v)) / dp.wavelength
    r0 = dp.starting_range
    eta_ref0 = range_poly(dp.doppler_centroid, r0) / range_poly(dp.azimuth_fm_rate, r0)
    return TopsCarrier(dp.azimuth_fm_rate, dp.doppler_centroid, ks, dp.azimuth_time_interval,
                       r0, dp.range_pixel_spacing, Int(lines_per_burst) ÷ 2, eta_ref0)
end

@inline function (c::TopsCarrier)(line::Integer, sample::Integer)
    r = c.r0 + sample * c.dr
    ka = range_poly(c.fm, r)
    fdc = range_poly(c.dc, r)
    kt = c.ks / (1 - c.ks / ka)
    eta = (line - c.mid_line) * c.dt
    de = eta - (c.eta_ref0 - fdc / ka)
    # The quadratic is `s1reader`'s; the linear term is ISCE3's separate use of the Doppler, which
    # `s1_resample.py` hands `ResampSlc` beside the carrier polynomial. Measured: without it `in_I2`
    # reaches 69.55% of exact bytes against 93.50% with it, and the opposite sign gives 62.78%. Its
    # `eta_ref` subtraction is immaterial — `eta_ref` is ~0.01 s against `eta`'s +/-1.5 — and dropping it
    # changes `in_I2` by 4e-4 of a percent.
    return pi * kt * de * de + 2pi * fdc * de
end

"""
    sinc8(f) -> NTuple{8,Float64}

The eight taps of a Hann-windowed sinc at fractional offset `f`, for taps at −3 through +4.

`Sinc2dInterpolator`'s kernel: ISCE's `sinc_coef` weights `sinc(x)` by
`(1 − pedestal)/2 · cos(πx/(ns/2)) + (1 + pedestal)/2`, which at `pedestal = 0` and `ns = 8` is
`0.5 + 0.5 cos(πx/4)`. Unwindowed the same eight taps overshoot by 21%, and a bicubic undershoots by 7%.
"""
@inline function sinc8(f::Float64)
    # Taps at -3 through +4, which is `Sinc2dInterpolator`'s centring: shifting them to -4 through +3
    # costs `in_I2` 93.50% of exact bytes against 81.76%.
    return ntuple(8) do k
        x = (k - 4) - f
        s = x == 0 ? 1.0 : sinpi(x) / (pi * x)
        s * (0.5 + 0.5 * cospi(x / 4))
    end
end

"""
    SINC_SUBDIVISIONS

Fractional positions per pixel in [`SINC_TAPS`](@ref), the tabulated form of [`sinc8`](@ref).

Evaluating the kernel per pixel costs sixteen transcendentals for every output sample — `sinpi` and
`cospi` eight times each, once per axis — where the taps depend on nothing but the fractional offset. A
table costs two lookups instead, and the quantization it introduces is bounded by `1/2 SINC_SUBDIVISIONS`
of a pixel in interpolation position.

2048, so that bound is 0.00024 px against a resample whose offset field is itself interpolated to about
1e-5 px. Held to the amplitude gate in `tools/golden/coreg_bench.jl`.
"""
const SINC_SUBDIVISIONS = 2048

"""
    SINC_TAPS

[`sinc8`](@ref) on a grid of `SINC_SUBDIVISIONS + 1` fractional offsets, taps down the columns.

**Scaled so each column sums to one**, which is where the per-pixel `sum(ty) * sum(tx)` division goes:
the normalization depends only on the fraction, so it belongs in the table rather than in the inner loop.
"""
const SINC_TAPS = let n = SINC_SUBDIVISIONS
    T = Matrix{Float64}(undef, 8, n + 1)
    for k in 0:n
        t = sinc8(k / n)
        s = sum(t)
        for m in 1:8
            T[m, k + 1] = t[m] / s
        end
    end
    T
end

"""
    sinc8_taps(f) -> NTuple{8,Float64}

[`sinc8`](@ref) at the tabulated offset nearest `f`, already normalized to sum to one.

`f` must lie in `[0, 1]`, which is what `y - floor(y)` gives.
"""
# `floor` to an `Int` without the conversion's own range check.
#
# **The check is 20% of the interpolation kernel** — `Int64(::Float64)` at `float.jl:920` plus the `isinf`
# beside it — and what it guards against cannot reach here: a window's offsets are all finite by the time
# the kernel runs, because [`_support_bounds`](@ref) evaluates the field at the corners of every lattice
# cell the window meets, and `_resample_piece` then calls `floor` on those bounds. A `NaN` anywhere in
# those cells propagates to a corner and throws there, one call earlier and outside the loop.
#
# The values themselves are pixel indices of order 1e3 to 1e5, so nothing here approaches `Int`'s range.
@inline _ifloor(x::Float64) = unsafe_trunc(Int, floor(x))

@inline _tap_bin(f::Float64) = round(Int, f * SINC_SUBDIVISIONS) + 1
@inline sinc8_taps_at(k::Int) = ntuple(m -> @inbounds(SINC_TAPS[m, k]), 8)
@inline sinc8_taps(f::Float64) = sinc8_taps_at(_tap_bin(f))

"""
    deramped_burst(raster, rows, carrier) -> Matrix{ComplexF32}

One burst's samples with its azimuth carrier removed, ready to interpolate.

Done to the whole burst once rather than per interpolation chip: a chip would evaluate the carrier 64
times per output pixel where this evaluates it once per input pixel, which is the difference between
minutes and hours over a subswath.

`rows` and `cols` index `raster`. `origin` is the burst-local zero-based `(line, sample)` of
`raster[first(rows), first(cols)]`, which is what the carrier is a function of: a band of a burst carries
the same phase it does inside the whole burst, so reading one must not restart the carrier at zero.
`(0, 0)` is the whole burst, whose first row *is* its first line.
"""
function deramped_burst(raster, rows::AbstractUnitRange, carrier;
                        cols::AbstractUnitRange = axes(raster, 2),
                        origin::Tuple{Integer,Integer} = (0, 0))
    nr, nc = length(rows), length(cols)
    l0, s0 = Int(origin[1]), Int(origin[2])
    out = Matrix{ComplexF32}(undef, nr, nc)
    # The carrier's range half is hoisted out of each column by `carrier_column`, leaving a degree-5
    # Horner sweep per pixel — see there for what the loop it replaces cost.
    #
    # **The phase stays `Float64` up to the trig.** `cis(-Float32(phi))` rounds a phase of order 1e3 rad
    # to about 5e-4 rad before taking its sine, which is coarser than anything else in this path; the
    # narrowing belongs on the finished phasor, whose parts are in `[-1, 1]`.
    #
    # **The band is read inside the threaded region, a column slab per task.** Reading it once up front
    # was a single-threaded 0.060 s against the 0.170 s of phase work around it — 26% of this function on
    # one core — and it also held a whole extra copy of the burst. A `StripedTiff` is immutable and
    # allocates its own gather buffer per call, so disjoint windows are read concurrently; each slab is
    # then deramped while it is still in cache, and nothing holds the burst twice.
    #
    # Several chunks per thread, so one slow slab cannot set the wall clock, and each is wide enough that
    # a line's copy out of the mapping stays a long contiguous run.
    nchunks = max(1, min(nc, 4 * Threads.nthreads()))
    edges = round.(Int, range(0, nc; length = nchunks + 1))
    Threads.@threads for ci in 1:nchunks
        j0, j1 = edges[ci] + 1, edges[ci + 1]
        j0 > j1 && continue
        slab = raster[rows, cols[j0]:cols[j1]]
        for j in j0:j1
            b = carrier_column(carrier, s0 + j - 1)
            @inbounds for i in 1:nr
                out[i, j] = ComplexF32(slab[i, j - j0 + 1]) *
                            ComplexF32(cis(-carrier_phase(b, carrier, l0 + i - 1)))
            end
        end
    end
    return out
end

# A burst's ramp margins zeroed in place, reproducing the source rectangle `slc_to_vrt_file` exposes.
# Bounds are the annotation's, converted to 0-based by the caller, and the array is 1-based.
#
# `origin` is the burst-local zero-based `(line, sample)` of `A[1, 1]`, so a band of a burst is zeroed
# against the burst's valid rectangle rather than against its own first row.
# **Four rectangle fills rather than a test per sample.** What is zeroed is a frame, so the interior — all
# but a few rows and columns of a burst — is visited by a predicate that is false for every one of its 31
# million samples. Filling the margins instead touches only the margins: 0.044 s per burst becomes the
# cost of the frame, and `fill!` over a view is a store loop rather than a branch per element.
function _zero_outside_valid!(A::AbstractMatrix, fvl::Integer, lvl::Integer, fvs::Integer, lvs::Integer;
                              origin::Tuple{Integer,Integer} = (0, 0))
    nr, nc = size(A)
    oy, ox = Int(origin[1]), Int(origin[2])
    z = zero(eltype(A))
    # The valid rectangle in `A`'s own indices. `r1 < r0` where it misses `A` entirely, which the row
    # fills below then cover between them.
    r0, r1 = clamp(fvl - oy + 1, 1, nr + 1), clamp(lvl - oy + 1, 0, nr)
    c0, c1 = clamp(fvs - ox + 1, 1, nc + 1), clamp(lvs - ox + 1, 0, nc)
    r0 > 1 && fill!(view(A, 1:(r0 - 1), :), z)
    r1 < nr && fill!(view(A, (r1 + 1):nr, :), z)
    if r0 <= r1
        c0 > 1 && fill!(view(A, r0:r1, 1:(c0 - 1)), z)
        c1 < nc && fill!(view(A, r0:r1, (c1 + 1):nc), z)
    end
    return A
end

"""
    resample_burst(deramped, dl, ds, lines, samples) -> Matrix{Float32}

`deramped` read at each reference pixel's own position in the secondary, as amplitude.

`dl` and `ds` are callables giving the offsets at a reference `(line, sample)`; in practice they
interpolate a lattice, since the field moves by 0.0036 lines over 1200 lines and the lattice's own
interpolation error is 1e-5 px.

Amplitude rather than the complex value, so the reramp is omitted: it multiplies the finished sample by a
unit phasor and cannot change a magnitude.

`origin` is the burst-local zero-based `(line, sample)` of `deramped[1, 1]` and `extent` is the burst's
own size, so `deramped` may be a band of the burst rather than all of it.

**The rejection is against `extent` and the read is against the band**, and the two have to be separate
for a windowed result to equal a whole-burst one. A pixel whose eight-tap support leaves the *burst* is
rejected by both, and must be: the reference's resampler has no samples there either. A pixel whose
support leaves the *band* is a band that was sized too small, which is a caller error rather than a pixel
to drop — dropping it silently is exactly how a lazily resampled block would come to differ from the
mosaic — so it throws.
"""
function resample_burst(deramped, dl, ds, lines::AbstractUnitRange,
                        samples::AbstractUnitRange, valid = nothing; doppler = nothing,
                        origin::Tuple{Integer,Integer} = (0, 0),
                        extent::Tuple{Integer,Integer} = size(deramped))
    out = zeros(Float32, length(lines), length(samples))
    oy, ox = Int(origin[1]), Int(origin[2])
    ey, ex = Int(extent[1]), Int(extent[2])
    ny, nx = size(deramped)
    # Over the sample axis's own indices rather than `collect(enumerate(samples))`, which materializes a
    # vector of pairs for `@threads` to index.
    # **Constant-folded, not branched.** `dl` and `ds` are type parameters of this method, so
    # `_cell_linear` is decided at compile time for each specialization and the path not taken is
    # eliminated. A lattice field is advanced per cell — two evaluations every 64 lines rather than a pair
    # per pixel, which is 1.5x of this kernel — and a per-pixel raster is read per pixel, as it must be.
    cellwise = _cell_linear(dl) && _cell_linear(ds)
    Threads.@threads for jj in eachindex(samples)
        s = samples[jj]
        cell = typemin(Int)
        l_cell = 0
        pl = ql = ps = qs = 0.0
        for (ii, l) in enumerate(lines)
            # `_dopplerLUT.contains(az, rng)`, the third of `ResampSlc::_transformTile`'s five rejections
            # and the one the two bounds tests below do not cover.
            #
            # It is evaluated at the **output** pixel — `az = _sensingStart + iRow / _prf` and
            # `rng = _startingRange + iCol * _rangePixelSpacing` — while the LUT is built from the *source
            # burst's* own shape: `doppler_poly1d_to_lut2d` spans `starting_slant_range` to
            # `starting_slant_range + (samples_per_burst - 1) * dr` in range and `0` to
            # `lines_per_burst * dt` in azimuth (`s1_reader.py:131-168`). `ResampSlc` takes its
            # `_sensingStart` and `_startingRange` from the burst's grid, the first constructor argument,
            # rather than from the `ref_rdr_grid` keyword. So the whole test collapses to the output index
            # lying inside the source burst's dimensions, and it bites whenever the output grid is the
            # larger of the two — which is every full-SLC pair, where the CSLC is 136 to 195 lines taller
            # than the burst it came from.
            if doppler !== nothing
                (l > doppler[1] || s > doppler[2] - 1) && continue
            end
            if cellwise
                k = fld(l, _LSTEP)
                if k != cell
                    l_cell = k * _LSTEP
                    pl = dl(l_cell, s)
                    ps = ds(l_cell, s)
                    ql = (dl(l_cell + _LSTEP, s) - pl) / _LSTEP
                    qs = (ds(l_cell + _LSTEP, s) - ps) / _LSTEP
                    cell = k
                end
                d = l - l_cell
                y = l + (pl + ql * d) + 1.0
                x = s + (ps + qs * d) + 1.0
            else
                y = l + dl(l, s) + 1.0
                x = s + ds(l, s) + 1.0
            end
            iy, ix = _ifloor(y), _ifloor(x)
            (iy - 4 < 1 || iy + 4 > ey || ix - 4 < 1 || ix + 4 > ex) && continue
            # The eight-tap support has to sit inside the secondary's own valid region: outside it the
            # reference's resampler has no samples either, and interpolating the zero margin invents a
            # small non-zero value where `secondary.tif` holds none.
            if valid !== nothing
                (iy - 3 < valid[1] || iy + 4 > valid[2] || ix - 3 < valid[3] || ix + 4 > valid[4]) &&
                    continue
            end
            by, bx = iy - oy, ix - ox
            # The taps below read `by - 4 + m` for `m` in `1:8`, so the support is `by - 3` through
            # `by + 4` — which is what the band has to contain and what the message states. The extent
            # test above is deliberately one row wider on the low side, reproducing `ResampSlc`'s own
            # rejection rather than this stencil's true reach.
            (by - 3 < 1 || by + 4 > ny || bx - 3 < 1 || bx + 4 > nx) && throw(ArgumentError(
                "the deramped band covers burst-local lines $(oy + 1):$(oy + ny) and samples " *
                "$(ox + 1):$(ox + nx), but output pixel ($l, $s) reads lines $(iy - 3):$(iy + 4) " *
                "and samples $(ix - 3):$(ix + 4); the band was sized too small for the window"))
            ty = sinc8_taps(y - iy)
            tx = sinc8_taps(x - ix)
            # **`@inbounds` over exactly the range the guard above admits.** The test on the line before
            # this one establishes `by - 3 >= 1`, `by + 4 <= ny`, `bx - 3 >= 1` and `bx + 4 <= nx`, and the
            # rows read below run `by - 3` through `by + 4` — so the accesses are the ones just proved, and
            # no index here comes from anywhere else.
            #
            # Earned rather than assumed: the check is 42% of this kernel, measured, and dropping it is
            # bit-identical on a whole burst — 0 of 30.3 million outputs differ. The indices come from
            # `floor` of an interpolated position, which is why the compiler cannot prove this itself the
            # way it would for `eachindex`. `tools/golden/selftest.jl` pins the guard to the access.
            acc = ComplexF64(0)
            @inbounds for m in 1:8
                row = ComplexF64(0)
                @simd for k in 1:8
                    row += tx[k] * ComplexF64(deramped[by - 4 + m, bx - 4 + k])
                end
                acc += ty[m] * row
            end
            # No tap-sum division: `SINC_TAPS` is normalized, so both axes already carry it.
            out[ii, jj] = Float32(abs(acc))
        end
    end
    return out
end

"""
    secondary_swath_amplitude(rp, sp, swath, dem) -> (Matrix{Float32}, Int, Int)

[`swath_amplitude`](@ref) for the secondary: each burst resampled onto the reference burst's grid first,
then placed by the **reference's** own seam arithmetic.

That the placement is the reference's is not a simplification — `merge_bursts_in_swath` takes
`az_reference_offsets` from `ref_bursts` and the range window from the reference burst's valid samples,
and applies both to the secondary. So the two mosaics share every index and differ only in pixel values.

The offsets come from a lattice of [`coregistration_offset`](@ref) solves, one node every 64 lines and 512
samples, interpolated bilinearly between. That is exact to about 1e-5 px on this field and turns 31
million geometry solves per burst into a few hundred.
"""
function secondary_swath_amplitude(rp::Sentinel1Product, sp::Sentinel1Product, swath::Integer, dem;
                                   offsets = nothing, grid = nothing)
    r = ResampledSwath(rp, sp, swath, dem; offsets, grid)
    out = zeros(Float32, size(r))
    for p in r.places
        out[p.mosaic_rows, p.mosaic_cols] = _resample_piece(r, p, p.lines, p.samples)
    end
    return (out, size(r, 1), size(r, 2))
end

"""
    BurstPlacement

Where one burst's resampled samples go in the merged subswath, and what it takes to produce them.

`lines` and `samples` are the burst's own zero-based output indices; `mosaic_rows` and `mosaic_cols` are
the one-based rows and columns of the merged subswath they land on. The two pairs have equal lengths by
construction, so the mapping is positional: `mosaic_rows[t]` carries `lines[t]`.
"""
struct BurstPlacement{C,DL,DS}
    burst::Int
    source_rows::UnitRange{Int}
    carrier::C
    dl::DL
    ds::DS
    # The secondary's own valid rectangle, zero-based, which bounds what the source offers.
    svalid::NTuple{4,Int}
    lines::UnitRange{Int}
    samples::UnitRange{Int}
    mosaic_rows::UnitRange{Int}
    mosaic_cols::UnitRange{Int}
    extent::Tuple{Int,Int}
end

"""
    ResampledSwath(rp, sp, swath, dem; offsets = nothing, grid = nothing) <: AbstractMatrix{Float32}

The secondary acquisition's amplitude on the reference's merged subswath grid, resampled on demand.

`getindex` over two ranges resamples **only that window**, so nothing holds the mosaic and nothing writes
it to a file. `secondary_swath_amplitude` is this materialized in one pass, and a window read from here
equals the same window of that array — which is a property of sharing one implementation rather than two
that agree.

That the placement is the reference's is not a simplification — `merge_bursts_in_swath` takes
`az_reference_offsets` from `ref_bursts` and the range window from the reference burst's valid samples,
and applies both to the secondary. So the two mosaics share every index and differ only in pixel values.

The offsets come from a lattice of [`coregistration_offset`](@ref) solves, one node every 64 lines and 512
samples, interpolated bilinearly between. That is exact to about 1e-5 px on this field and turns 31
million geometry solves per burst into a few hundred.

**The per-burst setup is done once, at construction.** The carrier fit and the offset lattice cost about a
millisecond each and depend on nothing a window varies, so a caller reading many windows pays for them
once; what a window pays for is the deramp of the band it needs and its own interpolation.

`deramped` counts source pixels deramped so far. Reading the mosaic in blocks deramps the eight-tap halo
of each block twice, so this against the mosaic's area is the read amplification a block size costs.
"""
struct ResampledSwath{S,P} <: AbstractMatrix{Float32}
    source::S
    places::Vector{P}
    dims::Tuple{Int,Int}
    deramped::Base.RefValue{Int}
end

function ResampledSwath(rp::Sentinel1Product, sp::Sentinel1Product, swath::Integer, dem;
                        offsets = nothing, grid = nothing)
    a = annotation(rp, swath)
    sa = annotation(sp, swath)
    n = nbursts(a)
    lpb, spb = a.lines_per_burst, a.samples_per_burst
    dt = a.azimuth_time_interval
    sraster = burst_raster(first(collect(bursts(sp, swath))).backend).raster
    fvl = a.first_valid_line .- 1
    lvl = a.last_valid_line .- 1
    fvs = a.first_valid_sample .- 1
    lvs = a.last_valid_sample .- 1

    lims = map(1:n) do i
        s = round(Int, (seconds_between(first(a.burst_start), a.burst_start[i]) + fvl[i] * dt) / dt)
        (s, s + (lvl[i] - fvl[i]) + 1)
    end
    # **The output grid is the one COMPASS resampled onto, which is the annotation's burst only sometimes.**
    # `merge_bursts_in_swath` reads its `num_az_samples`/`num_rng_samples` off the first burst's CSLC, so
    # the merged height and the merged width both follow that raster and not `lines_per_burst`.
    #
    # `grid` maps a **burst index** to that burst's `(lines, samples)`: the mosaic constructors resolve the
    # per-subswath level and hand this one what is left, matching `offsets`. So a caller supplying it from
    # outside needs two levels — `sw -> (i -> shape)` — which is what `_reference_grid` returns and is why
    # `cslc_grid`'s per-subswath `sw -> shape` cannot be passed here.
    nl, ns = grid === nothing ? (lpb, spb) : grid(1)
    nlines = n == 1 ? nl :
             1 + round(Int, (seconds_between(first(a.burst_start), last(a.burst_start)) +
                             (nl - 1) * dt) / dt)

    # **The annotation is parsed once per acquisition, not once per burst.** `open_slc(path; burst = i)`
    # re-reads and re-parses the whole subswath's annotation XML for each burst it is asked for, and this
    # loop asked fourteen times for a seven-burst subswath. `bursts` walks the bursts of a single parse.
    sb = collect(bursts(sp, swath))
    rb = offsets === nothing ? collect(bursts(rp, swath)) : nothing
    # **A task per burst, because none of this setup is cheap and all of it is independent.** Measured over
    # seven bursts: `RadarCoordinate` 0.242 s, `CarrierFit` 0.191 s and `_offset_lattice` 0.095 s — a fifth
    # of the whole subswath, previously one burst at a time. The lattice threads over its own nodes
    # internally; spawning here lets the scheduler fill the machine with whichever burst is ready rather
    # than draining one lattice before starting the next.
    # `identity.` narrows the element type: `fetch` on a `Task` is `Any`, and `ResampledSwath` holds a
    # vector of a concrete `BurstPlacement`, so the broadcast is what turns twelve `Any`s back into the one
    # type they all share.
    places = identity.(fetch.(map(1:n) do i
        Threads.@spawn begin
        ss = sb[i]
        cs = RadarCoordinate(ss)
        carrier = CarrierFit(TopsCarrier(ss, cs, lpb), lpb, spb)
        dl, ds = if offsets === nothing
            _offset_lattice(RadarCoordinate(rb[i]), cs, dem, lpb, spb)
        else
            offsets(i)
        end

        prev = i > 1 ? fld(lims[i - 1][2] - lims[i][1], 2) : 0
        nxt = i < n ? fld(lims[i][2] - lims[i + 1][1], 2) : 0
        bstart, bend = n == 1 ? (0, nl) : (fvl[i] + prev, 1 + lvl[i] - nxt)
        mstart, mend = n == 1 ? (0, nl) : (lims[i][1] + prev, lims[i][2] - nxt)
        # `min(lvs, ns)`: the valid-sample window is the annotation's and the raster it indexes is the
        # CSLC's, which can be the narrower of the two. NumPy's slicing truncates there in silence.
        cols = n == 1 ? (1:ns) : ((fvs[i] + 1):min(lvs[i], ns))
        BurstPlacement(i, ((i - 1) * lpb + 1):(i * lpb), carrier, dl, ds,
                       (sa.first_valid_line[i] - 1, sa.last_valid_line[i] - 1,
                        sa.first_valid_sample[i] - 1, sa.last_valid_sample[i] - 1),
                       bstart:(bend - 1), (first(cols) - 1):(last(cols) - 1),
                       (mstart + 1):mend, cols, (lpb, spb))
        end
    end))
    # `Int`: a `grid` read from the CSLC's own header carries whatever integer type that header
    # used, and the field is declared `Tuple{Int,Int}`.
    return ResampledSwath(sraster, places, (Int(nlines), Int(ns)), Ref(0))
end

Base.size(r::ResampledSwath) = r.dims

# Where a window's eight-tap support lands in the source, in burst-local one-based coordinates.
#
# **Exact, from the lattice's own cell corners.** Inside one cell of [`_offset_lattice`](@ref)'s grid the
# offsets are bilinear and the index is linear, so `y = l + dl(l, s)` is bilinear there — and a bilinear
# function on a rectangle attains its extremes at that rectangle's corners. The extremes over a window are
# therefore the extremes over the corners of the cells it meets, which are the window's own edges together
# with the lattice lines inside it. Nothing between those points can exceed them.
#
# That makes [`resample_burst`](@ref)'s band assertion unreachable rather than merely unlikely:
# `_resample_piece` reads `floor(ymin) - 3` through `ceil(ymax) + 4`, and a pixel at `y` needs
# `floor(y) - 3` through `floor(y) + 4`, which that range contains for every `y` in `ymin..ymax`.
#
# **A stride through the window is not a bound**, which is what this replaced: the clamp in `lerp`
# extrapolates from the edge cells, so a window reaching past the last node is governed by a cell wider
# than the lattice, and a sampled maximum can sit below the true one. Measured on
# `S1A_IW_SLC__1SSV_20240618` IW2, where the range offset reaches +87 samples, a stride at half the node
# spacing missed it by more than a sample and the assertion fired. Evaluating every pixel instead is also
# exact but costs 0.385 s per burst against this 0.000023 s — 31% of a burst's coregistration, single
# threaded. `tools/golden/selftest.jl` pins the two against each other.
_next_mark(v::Integer, step::Integer) = (fld(v, step) + 1) * step

function _support_bounds(dl, ds, lines::AbstractUnitRange, samples::AbstractUnitRange)
    ymin = xmin = Inf
    ymax = xmax = -Inf
    s = first(samples)
    while true
        l = first(lines)
        while true
            y = l + dl(l, s) + 1.0
            x = s + ds(l, s) + 1.0
            ymin = min(ymin, y); ymax = max(ymax, y)
            xmin = min(xmin, x); xmax = max(xmax, x)
            l == last(lines) && break
            l = min(_next_mark(l, _LSTEP), last(lines))
        end
        s == last(samples) && break
        s = min(_next_mark(s, _SSTEP), last(samples))
    end
    return (ymin, ymax, xmin, xmax)
end

# One burst's contribution to a window, over the burst-local output indices `lines` x `samples`.
#
# The band read from the source is the support those outputs reach, clamped to the burst: a pixel whose
# support leaves the burst is rejected by `resample_burst` and never read, so clamping cannot change an
# answer.
function _resample_piece(r::ResampledSwath, p::BurstPlacement, lines::AbstractUnitRange,
                         samples::AbstractUnitRange)
    lpb, spb = p.extent
    ymin, ymax, xmin, xmax = _support_bounds(p.dl, p.ds, lines, samples)
    brows = clamp(floor(Int, ymin) - 3, 1, lpb):clamp(ceil(Int, ymax) + 4, 1, lpb)
    bcols = clamp(floor(Int, xmin) - 3, 1, spb):clamp(ceil(Int, xmax) + 4, 1, spb)
    origin = (first(brows) - 1, first(bcols) - 1)

    srows = (first(p.source_rows) - 1 + first(brows)):(first(p.source_rows) - 1 + last(brows))
    r.deramped[] += length(brows) * length(bcols)
    deramped = deramped_burst(r.source, srows, p.carrier; cols = bcols, origin)
    # **The resampler's source is the burst zeroed outside its own valid window, not the raw burst.**
    # `slc_to_vrt_file` writes a VRT of the burst's full shape whose `SimpleSource` covers only
    # `first_valid_line:last_valid_line` by `first_valid_sample:last_valid_sample`, with
    # `NoDataValue` 0, so everything outside that rectangle reads as zero
    # (`s1_burst_slc.py:slc_to_vrt_file`). A TOPS burst's ramp-up and ramp-down margins carry small but
    # non-zero amplitudes — a median of 2.24 against a typical 50 on this pair — so reading them fills
    # pixels the reference leaves empty.
    #
    # **Zeroed rather than declined.** Rejecting an output pixel whose interpolation support straddles
    # the boundary is a different rule and a worse one: measured on IW1 burst 1 of
    # `S1A ... 20170221`, guarding declines 161,588 pixels the reference fills. Zeroing the source
    # keeps them, tapered, which is what the reference's own interpolation does with a source rectangle
    # that ends there.
    _zero_outside_valid!(deramped, p.svalid...; origin)
    # No valid-window guard: measured, ISCE3's own criterion is the raster's bounds rather than the
    # annotation's valid region. Guarding on the valid window declines 137,782 pixels `secondary.tif`
    # fills, to avoid filling 2,746 it does not — fifty times the error it removes.
    return resample_burst(deramped, p.dl, p.ds, lines, samples; doppler = p.extent, origin,
                          extent = p.extent)
end

function Base.getindex(r::ResampledSwath, rows::AbstractUnitRange, cols::AbstractUnitRange)
    checkbounds(r, rows, cols)
    out = zeros(Float32, length(rows), length(cols))
    for p in r.places
        mr = intersect(rows, p.mosaic_rows)
        isempty(mr) && continue
        mc = intersect(cols, p.mosaic_cols)
        isempty(mc) && continue
        # Positional, since a placement's mosaic range and its burst range have equal lengths.
        t0 = first(mr) - first(p.mosaic_rows)
        u0 = first(mc) - first(p.mosaic_cols)
        lines = (first(p.lines) + t0):(first(p.lines) + t0 + length(mr) - 1)
        samples = (first(p.samples) + u0):(first(p.samples) + u0 + length(mc) - 1)
        out[(first(mr) - first(rows) + 1):(last(mr) - first(rows) + 1),
            (first(mc) - first(cols) + 1):(last(mc) - first(cols) + 1)] =
            _resample_piece(r, p, lines, samples)
    end
    return out
end

# Scalar and mixed indexing go through the range form, which is the one that does the work. A whole-array
# read is `r[:, :]`, which materializes what `secondary_swath_amplitude` returns.
# **Reading one of these is a resample, not a lookup**, and `AutoRIFT.ondisk` is where an image says that
# a read costs more than memory access. It decides two things for a blocked run: that an unblocked pass is
# refused rather than resolving a pixel at a time, and that a caller may cache what it reads — which
# matters here because a run sweeps every block once per pass and would otherwise resample the mosaic once
# per pass. See `tools/golden/tilecache.jl`.
AutoRIFT.ondisk(::ResampledSwath) = true

Base.getindex(r::ResampledSwath, i::Integer, j::Integer) = r[i:i, j:j][1, 1]
Base.getindex(r::ResampledSwath, rows::AbstractUnitRange, j::Integer) = r[rows, j:j][:, 1]
Base.getindex(r::ResampledSwath, i::Integer, cols::AbstractUnitRange) = r[i:i, cols][1, :]
Base.getindex(r::ResampledSwath, ::Colon, ::Colon) = r[axes(r, 1), axes(r, 2)]

# The offset field on a lattice, with bilinear interpolation between nodes. Two closures rather than two
# matrices so the caller reads a position rather than an index.
const _LSTEP, _SSTEP = 64, 512
function _offset_lattice(cr, cs, dem, lpb::Integer, spb::Integer)
    ls = 0:_LSTEP:(lpb + _LSTEP)
    ss = 0:_SSTEP:(spb + _SSTEP)
    DL = Matrix{Float64}(undef, length(ls), length(ss))
    DS = similar(DL)
    Threads.@threads for jj in eachindex(ss)
        for ii in eachindex(ls)
            DL[ii, jj], DS[ii, jj], _ = coregistration_offset(cr, cs, ls[ii], ss[jj], dem)
        end
    end
    return (LatticeOffsets(DL, length(ls), length(ss)), LatticeOffsets(DS, length(ls), length(ss)))
end

"""
    LatticeOffsets(A, nl, ns)

One offset component interpolated bilinearly on [`_offset_lattice`](@ref)'s grid.

**A type rather than a closure, because the kernel needs a promise the closure could not make.** Inside
one lattice cell this is bilinear, so at a fixed sample it is *linear in the line index* — which lets
[`resample_burst`](@ref) take two values per cell instead of one per pixel. An offset field read straight
from a per-pixel raster, as [`burst_offsets`](@ref) returns, has no such property: interpolating it across
64 lines would invent values it does not hold. [`_cell_linear`](@ref) is where the two are told apart, and
it answers per type so the choice costs nothing at run time.

Indices outside the lattice are clamped to the edge cell, which is also what makes probing a cell's far
end safe at the end of a burst.
"""
struct LatticeOffsets{M<:AbstractMatrix{Float64}}
    A::M
    nl::Int
    ns::Int
end

@inline function (f::LatticeOffsets)(l::Integer, s::Integer)
    fi = l / _LSTEP + 1
    fj = s / _SSTEP + 1
    i = clamp(floor(Int, fi), 1, f.nl - 1)
    j = clamp(floor(Int, fj), 1, f.ns - 1)
    u, v = fi - i, fj - j
    A = f.A
    @inbounds return (1-u)*(1-v)*A[i,j] + u*(1-v)*A[i+1,j] + (1-u)*v*A[i,j+1] + u*v*A[i+1,j+1]
end

"""
    _cell_linear(f) -> Bool

Whether `f` is linear in the line index within one lattice cell, so a caller may advance it per cell.

`false` for anything that has not said otherwise — a per-pixel raster read among them.
"""
_cell_linear(::Any) = false
_cell_linear(::LatticeOffsets) = true

"""
    secondary_mosaic(rp, sp, swaths, dem) -> Matrix{Float32}

The merged amplitude raster `merge_swaths` writes as `secondary.tif`.

The same swath stack as [`radar_mosaic`](@ref) — the reference's extents, offsets and first-come writer —
over [`secondary_swath_amplitude`](@ref)'s coregistered subswaths.
"""
secondary_mosaic(rp::Sentinel1Product, sp::Sentinel1Product, swaths, dem; offsets = nothing,
                 grid = nothing) =
    _stack_swaths(rp, swaths,
                  [secondary_swath_amplitude(rp, sp, sw, dem;
                                             offsets = offsets === nothing ? nothing : offsets(sw),
                                             grid = grid === nothing ? nothing : grid(sw))
                   for sw in collect(swaths)])

"""
    CarrierFit(c::TopsCarrier, lines, samples)

`c` as the **fitted polynomial** ISCE3 actually evaluates, rather than the carrier itself.

`get_az_carrier_poly` samples the carrier every 50 lines and 500 samples and least-squares fits a
polynomial in normalized range and azimuth time, keeping the terms of `x^j y^i` with `i <= 5`, `j <= 3`
and `i + j <= 5` — eighteen of them (`s1reader.polyfit`, `max_order = True`). `ResampSlc` then evaluates
that fit, so reproducing the carrier exactly is *more* accurate than the reference and therefore different
from it.

**The one-sample and one-line shifts fall out of the fit rather than being applied by hand.** The carrier
is evaluated at index `(y, x)` but fitted against the coordinates of `(y + 1, x + 1)`, so the polynomial
returns, at a given pixel, the carrier of the pixel before it on both axes. That is where
`CARRIER_LINE_SHIFT`'s empirical `-1` came from, and with the fit in place the shift belongs at zero.
"""
struct CarrierFit
    coef::Vector{Float64}
    # The same coefficients as a dense `(azimuth order + 1, range order + 1)` table, which is what makes
    # the evaluation below a pair of Horner sweeps rather than a loop over exponent pairs. Absent terms
    # are zero, so the table is read at fixed offsets and the `i + j <= 5` rule costs nothing.
    table::Matrix{Float64}
    xmin::Float64
    xnorm::Float64
    ymin::Float64
    ynorm::Float64
    r0::Float64
    dr::Float64
    dt::Float64
end

# The eighteen exponent pairs, in `polyfit`'s own nesting order.
const _CARRIER_TERMS = [(i, j) for i in 0:5 for j in 0:3 if i + j <= 5]

function CarrierFit(c::TopsCarrier, lines::Integer, samples::Integer;
                   ystep::Integer = 50, xstep::Integer = 500)
    xs = 0:xstep:(samples - 1)
    ys = 0:ystep:(lines - 1)
    n = length(xs) * length(ys)
    # The coordinates the fit is against: one sample and one line past the index the value belongs to.
    rg = [c.r0 + (x + 1) * c.dr for x in xs]
    az = [(y + 1) * c.dt for y in ys]
    xmin, xnorm = minimum(rg), max(maximum(rg) - minimum(rg), 1.0)
    ymin, ynorm = minimum(az), max(maximum(az) - minimum(az), 1.0)
    A = Matrix{Float64}(undef, n, length(_CARRIER_TERMS))
    z = Vector{Float64}(undef, n)
    k = 0
    for (jy, y) in enumerate(ys), (jx, x) in enumerate(xs)
        k += 1
        xn = (rg[jx] - xmin) / xnorm
        yn = (az[jy] - ymin) / ynorm
        for (t, (i, j)) in enumerate(_CARRIER_TERMS)
            A[k, t] = xn^j * yn^i
        end
        z[k] = c(y, x)
    end
    coef = A \ z
    table = zeros(Float64, 6, 4)
    for (t, (i, j)) in enumerate(_CARRIER_TERMS)
        table[i + 1, j + 1] = coef[t]
    end
    return CarrierFit(coef, table, xmin, xnorm, ymin, ynorm, c.r0, c.dr, c.dt)
end

"""
    carrier_column(f::CarrierFit, sample) -> NTuple{6,Float64}

The azimuth polynomial's coefficients at one range sample.

**The fit is separable and the sample index is the outer loop's**, so the range half belongs outside the
column: this evaluates it once per column and leaves a degree-5 polynomial in the line index, which
[`carrier_phase`](@ref) then sweeps per pixel.

What that replaces is a loop over eighteen exponent pairs raising `xn^j * yn^i` with the exponents as
*runtime* integers — each `^` a call to `power_by_squaring`, none of it unrolled. On one burst of
1502 x 20662 the whole deramp goes from 0.355 s to 0.061 s for it, and the phase agrees with the loop it
replaces to 1.3e-11 rad.
"""
@inline function carrier_column(f::CarrierFit, sample::Integer)
    xn = ((f.r0 + sample * f.dr) - f.xmin) / f.xnorm
    T = f.table
    return ntuple(6) do i
        @inbounds T[i, 1] + xn * (T[i, 2] + xn * (T[i, 3] + xn * T[i, 4]))
    end
end

"""
    carrier_phase(b::NTuple{6}, f::CarrierFit, line) -> Float64

The carrier at one line, given [`carrier_column`](@ref)'s coefficients for that column.
"""
@inline function carrier_phase(b::NTuple{6,Float64}, f::CarrierFit, line::Integer)
    yn = (line * f.dt - f.ymin) / f.ynorm
    return ((((b[6] * yn + b[5]) * yn + b[4]) * yn + b[3]) * yn + b[2]) * yn + b[1]
end

@inline (f::CarrierFit)(line::Integer, sample::Integer) =
    carrier_phase(carrier_column(f, sample), f, line)
