# Reading a golden case's scenes, as opposed to establishing that they exist.
#
# `fetch.jl::resolve_inputs` answers reachability and returns no path. This returns a path GDAL can
# open and the footprint behind it, which is what the end-to-end ladder needs in order to start
# where the reference starts: at the granule.
#
# Nothing is downloaded. Both routes read the object in place — `/vsicurl` for the anonymous
# Sentinel-2 mirror, `/vsis3` for requester-pays Landsat — so a case costs the tiles its window
# covers rather than a scene on disk.

using ArchGDAL, Dates, Downloads, GeoFormatTypes, JSON3, ImagePairGeometry
using GeoFormatTypes: EPSG

include("parameters.jl")

"""
    scene_band(platform) -> Symbol

Which band the reference correlates: `:pan` for Landsat 7/8/9, `:green` for Landsat 4/5, `:pan` for
Sentinel-2.

`process.py:77-89`. The STAC asset keys are used rather than band numbers because the number moves
with the sensor — green is B2 on TM and B3 on OLI — while the key does not.
"""
function scene_band(platform::AbstractString)
    platform in ("L4", "L5") && return :green
    platform in ("L7", "L8", "L9", "S2") && return :pan
    throw(ArgumentError("no band rule for platform \"$platform\"; the optical rule is " *
                        "process.py:77-89 and the radar platforms carry no band"))
end

"""
    acquisition_order(c::GoldenCase) -> (early, late)

`c`'s two scene names ordered by acquisition date, which is the order the pipeline correlates them in
whatever order the job listed them.

**The job's `reference`/`secondary` is not that order.** Two of the twenty-two jobs name the later
acquisition as the reference, and their products still report the earlier one as `id_img1` with a
*positive* `date_dt` — on `LE07_L1TP_063018`, `img1` is 20040810 and `date_dt` is +32.000 while the
job's reference is 20040911. So the pipeline sorts the pair before geogrid sees it, and a harness that
follows the job order hands the geogrid a negative interval.

The cost of getting it wrong is not a rejected run: `window_offset` and all four `off2vel` bands come
back negated — a velocity field pointing backwards — and `window_search_range` doubles, because the
short-interval search inflation `max(1, 5 - 4·dt/182)` grows as the interval falls below zero. Every
integer band that does not depend on the interval stays exact, so the case looks two-thirds right.
"""
function acquisition_order(c::GoldenCase)
    a, b = only(c.reference), only(c.secondary)
    da, db = _scene_date(c.platform, a), _scene_date(c.platform, b)
    return da <= db ? (a, b) : (b, a)
end

"""
    STAGED :: Dict{String,String}

Granule name to a local copy of its scene, consulted by [`scene_path`](@ref) before the object store.

Empty unless a caller filled it with [`stage_scene`](@ref). What it exists for is a measurement whose
subject is the computation: a Landsat band is 354 MB behind a requester-pays bucket, so a chain that
reads the scene from there is timing the network as much as the filter.
"""
const STAGED = Dict{String,String}()

"""
    stage_scene(path) -> String

Copy the object at `path` under `<cache>/granules/` and return the local path.

A byte copy rather than a GDAL rewrite, so the staged file *is* the object the reference reads: a
re-encode would turn a JP2 into a GTiff and replace decompression with plain I/O, which changes what a
later read costs. An already-local path is returned untouched.

`/vsis3` goes through `aws s3 cp --request-payer requester`, which pays with whatever `AWS_PROFILE`
names — the same identity [`scene_path`](@ref) reads with.
"""
function stage_scene(path::AbstractString)
    (startswith(path, "/vsis3/") || startswith(path, "/vsicurl/")) || return path
    dir = joinpath(CACHE, "granules")
    dest = joinpath(dir, basename(path))
    isfile(dest) && return dest
    mkpath(dir)
    # Through a temporary, so an interrupted copy cannot leave a truncated file that the check above
    # would then treat as staged.
    tmp = dest * ".partial"
    if startswith(path, "/vsis3/")
        # `--only-show-errors`: the default progress meter writes a carriage-return line per 256 KiB,
        # which is tens of thousands of lines in a redirected harness log.
        run(`aws s3 cp $("s3://" * path[8:end]) $tmp --request-payer requester --only-show-errors`)
    else
        Downloads.download(path[10:end], tmp)
    end
    mv(tmp, dest; force = true)
    return dest
end

"""
    scene_path(c::GoldenCase, which::Symbol) -> String

A GDAL-openable path to `which` scene of `c`, at the band the reference correlates.

`which` is `:reference` or `:secondary` in the pipeline's sense — see [`acquisition_order`](@ref),
which is not the job's sense.

Landsat resolves through the public landsatlook STAC catalogue, whose item names the requester-pays
S3 object; reading it needs `AWS_PROFILE` to point at credentials that can pay, the same identity
`reference.jl` passes to the container.

Sentinel-2 resolves through the granule's `manifest.safe` on the anonymous Google Cloud mirror. The
granule directory name inside a `.SAFE` is not derivable from the product name — it carries its own
datatake identifier — so the manifest is read rather than guessed.
"""
function scene_path(c::GoldenCase, which::Symbol)
    early, late = acquisition_order(c)
    name = which === :reference ? early : late
    haskey(STAGED, name) && return STAGED[name]
    c.platform == "S2" && return _s2_path(name)
    startswith(c.platform, "L") && return _landsat_path(name, scene_band(c.platform))
    throw(ArgumentError("scene_path has no route for platform \"$(c.platform)\""))
end

"""
    gslc_amplitude_paths(c::GoldenCase, run) -> (reference, secondary)

The two byte-amplitude rasters a NISAR L2 GSLC pair is correlated from.

A GSLC is geocoded, so `nisar_isce3.process_gslc` runs the **projected** geogrid over it
(`optical_flag = 1`) and the imagery reaching the correlator is `reference_adjusted.tif` and
`secondary_adjusted.tif` — `convert_slc_to_uint8_amplitude` of the cropped products, already on one map
grid at 2.5 x 5.0 m. So there is no reprojection, no filter and no byte rescale on this path: the
amplitude rasters *are* the correlator's input, which is why rungs 5.2, 5.3 and 5.4 do not apply to it.

Taken from the cached run rather than rebuilt. Producing them means cropping two 11 GB products and
converting 6 GB of complex samples per side; the outputs are 6 GB each and already on disk.
"""
function gslc_amplitude_paths(c::GoldenCase, run::AbstractString)
    paths = (joinpath(run, "reference_adjusted.tif"), joinpath(run, "secondary_adjusted.tif"))
    for p in paths
        isfile(p) || error("$(basename(p)) is not in $run. A NISAR GSLC pair is correlated from the " *
                           "byte-amplitude rasters the driver writes there, and rebuilding them " *
                           "means cropping two 11 GB products.")
    end
    return paths
end

function _landsat_path(name::AbstractString, band::Symbol)
    url = "https://landsatlook.usgs.gov/stac-server/collections/landsat-c2l1/items/$name"
    item = JSON3.read(String(take!(Downloads.download(url, IOBuffer(); timeout = 60))))
    asset = get(item.assets, band, nothing)
    asset === nothing && error("STAC item $name has no `$band` asset")
    href = asset.alternate.s3.href
    startswith(href, "s3://") ||
        error("expected an s3:// href for $name band $band, got $href")
    # Requester pays, and GDAL wants it said explicitly; the account is whatever `AWS_PROFILE`
    # names, matching `run_reference`.
    ArchGDAL.setconfigoption("AWS_REQUEST_PAYER", "requester")
    return "/vsis3/" * href[6:end]
end

const S2_MIRROR = "https://storage.googleapis.com/gcp-public-data-sentinel-2/tiles"

function _s2_path(name::AbstractString)
    # `S2B_MSIL1C_20200612T150759_N0209_R025_T22WEB_20200612T184700`: the tile id is field 6, as
    # `T<zone><band><square>`, and the mirror splits it into three directory levels.
    tile = split(name, '_')[6]
    startswith(tile, "T") && length(tile) == 6 ||
        error("cannot read an MGRS tile out of \"$name\"; expected field 6 like T22WEB")
    safe = "$S2_MIRROR/$(tile[2:3])/$(tile[4:4])/$(tile[5:6])/$name.SAFE"
    manifest = String(take!(Downloads.download("$safe/manifest.safe", IOBuffer(); timeout = 90)))
    # The band the reference correlates is B08 at 10 m. One match is expected; more than one means
    # the manifest layout changed and guessing which is wrong.
    hits = unique(m.match for m in eachmatch(r"GRANULE/[^\"<>]*B08\.jp2", manifest))
    length(hits) == 1 ||
        error("expected one B08 entry in $name's manifest.safe, found $(length(hits))")
    return "/vsicurl/$safe/$(only(hits))"
end

"""
    aligned_scenes(c::GoldenCase) -> (reference, secondary)

Paths to `c`'s two scenes in one projection, reprojecting both if they are not already.

`GeogridOptical.coregister` refuses a pair in two coordinate systems outright
(`GeogridOptical.py:297-298`) — its overlap is index arithmetic in a single system — so the pipeline
brings them together first. `utils.ensure_same_projection` warps **both** scenes, not just the
secondary, to the reference's EPSG at the reference's own resolution, with `lanczos` resampling and
`targetAlignedPixels`. The reference is warped too because `-tap` snaps the output extent to a
multiple of the resolution, which can move its origin.

This must happen before the footprints are read: what geogrid intersects, and what its pixel indices
count in, is the *warped* grid. Four of the twelve optical golden pairs need it — the two cross-path
L7 pairs, the L8×L7 pair and one L5 pair, each straddling UTM 32607 and 32608.

Warped scenes are written to `<cache>/reprojected/<product>/`, mirroring the `reprojected/` directory
the container writes beside its outputs, and are kept: a warp is minutes of compute and, for Landsat,
requester-pays egress. A pair already in one projection is returned untouched, with nothing written.
"""
function aligned_scenes(c::GoldenCase)
    rpath, spath = scene_path(c, :reference), scene_path(c, :secondary)
    rfp, sfp = scene_footprint(rpath), scene_footprint(spath)
    target = footprint_epsg(rfp)
    # The footprints come back with the paths. Reading one is a remote open, and the caller needs exactly
    # these two; recomputing them there cost four opens per run where two will do.
    target == footprint_epsg(sfp) && return (rpath, spath, rfp, sfp)

    dir = joinpath(CACHE, "reprojected", c.product)
    mkpath(dir)
    @info "reprojecting a cross-projection pair" product=c.product target reference=footprint_epsg(rfp) secondary=footprint_epsg(sfp) dir
    # A warped scene has its own geotransform and extent, so its footprint is the warped file's rather
    # than the one measured above.
    wr, ws = _warp_to(rpath, rfp, target, dir), _warp_to(spath, rfp, target, dir)
    return (wr, ws, scene_footprint(wr), scene_footprint(ws))
end

# `reference` supplies the resolution for both scenes, since that is what the reference warps to.
#
# `ArchGDAL.gdalwarp` on the opened dataset rather than `Rasters.warp`, which is the idiomatic call and
# would be preferred but for one thing: it takes a `Raster` and builds a GDAL dataset from it before
# warping (`RastersArchGDALExt/warp.jl:28-45`, whose own TODO is to pass a lazy `FileArray` straight
# through), so the whole scene moves through memory on the way in. Given the dataset, GDAL pulls source
# blocks and writes destination blocks, which is what the command-line tool does.
#
# Measured on one 10980² Sentinel-2 band over `/vsicurl`: the two produce the **identical** output
# geometry — 12099² at the same geotransform — and this path peaks about 200 MiB lower, which is the
# scene at `UInt16`. Switching back is a one-line change once that TODO lands.
function _warp_to(path::AbstractString, reference, target::Integer, dir::AbstractString)
    out = joinpath(dir, replace(basename(path), r"\.(TIF|tif|jp2)$"i => "") * ".tif")
    isfile(out) && (@info "warped scene already cached" out; return out)

    # `abs` on the y resolution: the reference passes the geotransform's own negative value and
    # gdalwarp's `-tr` takes magnitudes. No creation options, because `ensure_same_projection` passes
    # none — a compressed or tiled output would be a different file from the one it wrote.
    #
    # `-of GTiff` is explicit rather than inferred from the extension, because the temporary name below
    # has to carry one GDAL recognizes or it writes nothing and reports no error.
    flags = ["-of", "GTiff", "-t_srs", "EPSG:$target",
             "-tr", string(abs(reference.spacing[1])), string(abs(reference.spacing[2])),
             "-r", "lanczos", "-tap"]
    @info "warping" from=basename(path) to=basename(out)
    # Written through a temporary name and moved, so an interrupted warp cannot leave a truncated file
    # that the `isfile` check above would then treat as cached.
    tmp = out * ".partial.tif"
    ArchGDAL.read(path) do src
        ArchGDAL.gdalwarp([src], flags; dest = tmp) do _
            nothing
        end
    end
    mv(tmp, out; force = true)
    return out
end

"""
    scene_footprint(path) -> ImageFootprint

The image's footprint: origin, signed spacing, size and CRS, read from its geotransform.

`ImageFootprint` carries the CRS so [`coregister`](@ref) can refuse a pair in two projections, which
is what the reference does (`GeogridOptical.py:297-298`) rather than reprojecting on the fly.
"""
function scene_footprint(path::AbstractString)
    gdal_network_setup()
    ds = ArchGDAL.read(path)
    gt = ArchGDAL.getgeotransform(ds)
    (gt[3] == 0 && gt[5] == 0) || error("$path has a rotated geotransform ($(gt[3]), $(gt[5])); " *
                                       "the optical geogrid assumes an axis-aligned grid")
    return ImageFootprint(origin = (gt[1], gt[4]), spacing = (gt[2], gt[6]),
                          size = (ArchGDAL.width(ds), ArchGDAL.height(ds)),
                          crs = EPSG(scene_epsg(ds)))
end

"""
    scene_epsg(ds) -> Int

The EPSG code of a dataset's *projected* CRS, by the reference's own procedure.

`GeogridOptical.getProjectionSystem` (`GeogridOptical.py:93-123`) imports the WKT, calls
`AutoIdentifyEPSG`, requires the result to be projected, and reads the authority code off `PROJCS` —
falling back to a database match if that is empty.

**`ArchGDAL.toEPSG` is not a substitute and fails silently here.** A Landsat scene over Antarctica
carries a custom `PROJCRS["PS         WGS84"]` with no authority code of its own; `toEPSG` walks down
to the `BASEGEOGCRS` and returns **4326**, which is a real code for a different coordinate system.
Nothing then errors: the grid-to-scene transform becomes 4326→4326, a no-op, and the pair's centroid
stays in metres — so the parameter-region lookup is handed a point 2,000 km outside the Earth's
coordinate range. Auto-identification resolves the same scene to **3031** at 100% confidence.
"""
function scene_epsg(ds)
    srs = ArchGDAL.importWKT(ArchGDAL.getproj(ds))
    ArchGDAL.GDAL.osrisprojected(srs.ptr) == 1 || error(
        "the scene's coordinate system is not projected; geogrid refuses a geographic or local " *
        "one (GeogridOptical.py:110-115) because its pixel arithmetic is in metres")
    ArchGDAL.GDAL.osrautoidentifyepsg(srs.ptr)
    code = ArchGDAL.GDAL.osrgetauthoritycode(srs.ptr, "PROJCS")
    isempty(code) || return parse(Int, code)

    # The reference's last resort is `gdalsrsinfo -o epsg`, which is this match against the EPSG
    # database. Done in process rather than by shelling out, and the confidence is reported so a
    # weak match is visible rather than silently adopted.
    n = Ref{Cint}(0)
    conf = Ref{Ptr{Cint}}(C_NULL)
    matches = ArchGDAL.GDAL.osrfindmatches(srs.ptr, C_NULL, n, conf)
    n[] > 0 || error("could not identify an EPSG code for the scene's projected CRS; the " *
                     "reference raises here too (GeogridOptical.py:123)")
    best = unsafe_wrap(Array, matches, n[])[1]
    @info "scene CRS identified by database match rather than by authority code" epsg=ArchGDAL.GDAL.osrgetauthoritycode(best, C_NULL) confidence=unsafe_wrap(Array, conf[], n[])[1]
    return parse(Int, ArchGDAL.GDAL.osrgetauthoritycode(best, C_NULL))
end

"""
    geogrid_seconds(c::GoldenCase) -> Float64

The interval the geogrid is given, in seconds: positive, and a whole number of days.

**Whole calendar days, not the pair's actual separation.** `testGeogridOptical.py:161-165` builds
two `datetime.date` objects from the first eight characters of each scene name's date field and
differences them, so the time of day is discarded. The product's `date_dt` is a different quantity,
computed later from the full timestamps (`netcdf_output.py`), and the two differ by up to half a
day.

Using `date_dt` here is not a rounding nicety: on the golden S2B case it moves 75 points of
`window_offset` and 445 of `window_search_range` by one pixel, because the offset is linear in the
interval and a 6e-5 relative change tips whatever sits nearest a rounding boundary.

Positive because the pair is taken in acquisition order — see [`acquisition_order`](@ref).
"""
function geogrid_seconds(c::GoldenCase)
    early, late = acquisition_order(c)
    d(name) = Date(_scene_date(c.platform, name), dateformat"yyyymmdd")
    return Float64(Dates.value(d(late) - d(early))) * 86400.0
end

# Where the acquisition date sits in a scene name, per platform. Landsat names it as field 4
# (`LC08_L1TP_009011_20200703_...`) and Sentinel-2 as the first eight characters of field 3
# (`S2B_MSIL1C_20200612T150759_...`), which is how the reference reads each.
function _scene_date(platform::AbstractString, name::AbstractString)
    parts = split(name, '_')
    startswith(platform, "L") && return String(parts[4])
    platform == "S2" && return String(parts[3][1:8])
    # Sentinel-1 names the acquisition start as field 6, `20150828T162412`.
    startswith(platform, "S1") && return String(parts[6][1:8])
    # NISAR names it as field 12, `20251028T235201`.
    startswith(platform, "NISAR") && return String(parts[12][1:8])
    throw(ArgumentError("no date rule for platform \"$platform\""))
end

"""
    pair_centroid(coord, epsg) -> (lon, lat)

Centre of `coord` in WGS84 degrees, which is what the parameter-region lookup takes.

The centre of the *coregistered* coordinate rather than of either scene: the region has to be the one
covering the overlap, since that is the only ground the pair reports.

`epsg` is passed rather than read off `coord`, because a `ProjectedCoordinate` carries origin,
spacing and size but no CRS — only the `ImageFootprint` it was built from does. The two scenes share
one CRS by construction, since `coregister` refuses a pair that does not.
"""
function pair_centroid(coord, epsg::Integer)
    cx = coord.origin[1] + (coord.size[1] - 1) * coord.spacing[1] / 2
    cy = coord.origin[2] + (coord.size[2] - 1) * coord.spacing[2] / 2
    lon = lat = 0.0
    # `order = :trad` gives (lon, lat); EPSG:4326's authority order is (lat, lon).
    ArchGDAL.crs2transform(EPSG(epsg), EPSG(4326); order = :trad) do tf
        p = ArchGDAL.createpoint(cx, cy)
        ArchGDAL.transform!(p, tf)
        lon, lat = ArchGDAL.getx(p, 0), ArchGDAL.gety(p, 0)
    end
    return (lon, lat)
end

"""
    footprint_epsg(fp::ImageFootprint) -> Int

The footprint's EPSG code, as an integer.
"""
footprint_epsg(fp) = GeoFormatTypes.val(fp.crs)

"""
    correlator_filter(c::GoldenCase) -> Union{Nothing,PreprocessMethod}

The filter `runAutorift` applies to the cropped pair before quantizing it, or `nothing`.

`testautoRIFT.py:533-540` picks a method **per scene** from the granule name and then
`:293-308` applies **one** of them to both images, most stringent first. So the dispatch is a
reduction over the pair rather than a per-image choice:

  * an `L[EO]07_` scene acquired on or after 2003-05-31 — after the Scan Line Corrector failed —
    contributes `wallis_fill`, and if either scene does, both are filtered with it;
  * an `LT0[45]_` scene contributes `fft`, whose branch is **commented out** and applies no filter
    at all (`:297-306` warns instead, because the band-reject has to run on the native scene before
    geogrid rounds the corners off);
  * anything else contributes `hps`, the plain high-pass.

The width is 5 everywhere except Sentinel-1, which uses 21 (`:526-528`). `StandardDeviationCutoff`
is `0.25` and never overridden (`autoRIFT.py:864`).

Returning `nothing` for an L4/L5 pair is the reference's behaviour, not an omission: those scenes
reach the correlator carrying only the filter [`native_filter`](@ref) already applied.
"""
function correlator_filter(c::GoldenCase)
    width = startswith(c.platform, "S1") ? 21 : 5
    methods = [_scene_method(n) for n in vcat(c.reference, c.secondary)]
    :wallis_fill in methods && return AutoRIFT.WallisGapfill(width, 0.25)
    :fft in methods && return nothing
    return AutoRIFT.Highpass(width)
end

# One scene's contribution to that reduction, by the name tests `testautoRIFT.py:535-540` applies.
function _scene_method(name::AbstractString)
    if occursin(r"^L[EO]07_", name)
        acquired = DateTime(split(name, '_')[4], dateformat"yyyymmdd")
        return acquired >= DateTime(2003, 5, 31) ? :wallis_fill : :hps
    end
    occursin(r"^LT0[45]_", name) && return :fft
    return :hps
end

"""
    native_filter(c::GoldenCase, name) -> Union{Nothing,Symbol}

The filter `process.py` applies to `name`'s **native** scene, writing `filtered/`, or `nothing`.

Gated on the *pair* and then chosen per scene, which is two separate tests and both matter
(`process.py:471-479`, `apply_landsat_filtering`). The pair is filtered only when
`min(reference_platform, secondary_platform)` is `L4`, `L5` or `L7` — so an L8 scene paired with an
L7 one is filtered and the same scene paired with an L9 one is not. Given that, each scene takes the
filter its own platform names: `:fft` for L4 and L5, `:wallis_fill` for L7 and L8.

`:fft` is Wallis at width 5 followed by the band-reject `AutoRIFT.Destripe` reproduces; `:wallis_fill`
is `AutoRIFT.WallisGapfill(5, 0.25)`.
"""
function native_filter(c::GoldenCase, name::AbstractString)
    platforms = [_landsat_platform(n) for n in vcat(c.reference, c.secondary)]
    any(isnothing, platforms) && return nothing
    minimum(platforms) in ("L4", "L5", "L7") || return nothing
    own = _landsat_platform(name)
    own in ("L4", "L5") && return :fft
    own in ("L7", "L8") && return :wallis_fill
    throw(ArgumentError("no native filter rule for \"$name\"; apply_landsat_filtering dispatches " *
                        "on L4, L5, L7 and L8 only"))
end

# `LT05_...` gives "L5": the platform code `get_platform` returns, which drops the leading zero so
# the string comparison `min` relies on orders L4 < L5 < L7 < L8 < L9.
function _landsat_platform(name::AbstractString)
    m = match(r"^L([TEOC])0?(\d)_", name)
    return m === nothing ? nothing : string('L', m.captures[2])
end

"""
    filtered_path(c::GoldenCase, run, name) -> Union{Nothing,String}

The reference's own natively-filtered raster for `name`, or `nothing` if the pair is not filtered.

`process.py` writes `apply_landsat_filtering`'s output into `filtered/` on the scene's **native** grid
and, for a cross-projection pair, `ensure_same_projection` then writes the warped copies into
`reprojected/`. The later directory is preferred when it exists, because that is what geogrid and the
correlator actually read — `reference_path` is reassigned to it (`process.py:488`) before the bounding
box is taken.

Returns the path the reference would have handed downstream, which is what a rung wanting the
*correlator's* input needs. Rung 5.3, which is testing the filter itself, wants `filtered/` on the
native grid instead and asks for it by name.
"""
function filtered_path(c::GoldenCase, run::AbstractString, name::AbstractString)
    native_filter(c, name) === nothing && return nothing
    band = scene_band(c.platform) === :green ? "B2" : "B8"
    for r in cached_runs(c, run), dir in ("reprojected", "filtered")
        p = joinpath(r, dir, "$(name)_$(band).TIF")
        isfile(p) && return p
    end
    error("no filtered raster for $name under any cached run of $(c.product); `filtered/` is what " *
          "`apply_landsat_filtering` writes, and a run pruned of it before " *
          "`prune_run` kept it is restored by another container run")
end

"""
    cached_runs(c::GoldenCase, run) -> Vector{String}

`run` first, then every other cached run of `c`.

**Which run holds an input need not be the run the ladder is comparing against.** `resolve_run` prefers a
run with a capture, and a run made to restore `filtered/` has none — so a rung that looked only in its own
run would ignore the very directory that was just rebuilt, and the advice to re-run the container would
not work. The same holds for a full-SLC pair's `.SAFE` trees, which the driver downloads into whichever
run it processed.
"""
function cached_runs(c::GoldenCase, run::AbstractString)
    root = dirname(run)
    isdir(root) || return [run]
    return [run; [d for d in sort(readdir(root; join = true)) if isdir(d) && d != run]]
end

"""
    native_filtered_paths(c::GoldenCase, run, name) -> (image, zero_mask)

The reference's native-grid filter outputs for `name`: the Float32 raster and its zero mask.

The mask is `nothing` for an L4/L5 pair, which is `apply_fft_filter`'s own return — only the
nodata-infill filter produces one (`process.py:275, 286`).
"""
function native_filtered_paths(c::GoldenCase, run::AbstractString, name::AbstractString)
    band = scene_band(c.platform) === :green ? "B2" : "B8"
    runs = cached_runs(c, run)
    i = findfirst(r -> isfile(joinpath(r, "filtered", "$(name)_$(band).TIF")), runs)
    isnothing(i) && error("no `filtered/$(name)_$(band).TIF` under any cached run of $(c.product); " *
                          "rung 5.3 compares against `filtered/` on the native grid")
    run = runs[i]
    img = joinpath(run, "filtered", "$(name)_$(band).TIF")
    zero = joinpath(run, "filtered", "$(name)_$(band)_zeroMask.TIF")
    return (img, isfile(zero) ? zero : nothing)
end
