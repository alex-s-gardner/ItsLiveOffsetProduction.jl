# The command line.

const USAGE = """
usage: autorift_itslive [<sensor>] <reference> <secondary> [<output.nc>] [options]
       autorift_itslive [<sensor>] [<output.nc>] -r <scene>... -s <scene>... [options]

  sensor       L4 L5 L7 L8 L9 S2 S1 S1BURST NL1 NL2; by default the reference's. A list of SLC zips is
               a burst pair (S1BURST), with --subset; so is a Sentinel-1 SLC given with --subset.
  reference    a scene's band — a local path, an s3:// URI or an https:// URL — or its STAC item
  secondary    (a local .json or a URL). Landsat: the <scene>_B8.TIF (L7, L8, L9) or <scene>_B2.TIF
               (L4, L5) band, with <scene>_MTL.json and <scene>_ANG.txt beside it, or a landsatlook
               STAC item. Sentinel-2: the B08.jp2 inside a .SAFE, or an Earth Search STAC item.
               Sentinel-1 (S1): the IW SLC product's .SAFE directory or its .zip, local or a URL;
               (S1BURST) the SLC zips the bursts are in, comma-separated, with --subset; the bursts'
               extractor URLs, comma-separated
               (https://sentinel1-burst.asf.alaska.edu/<SLC>/IW<n>/<POL>/<index>.tiff); or the
               SAFE burst2safe assembles from them.
               NISAR: the RSLC (NL1) or GSLC (NL2) granule's .h5.
               Or the scene ID alone: a Landsat scene (LC08_L1TP_021248_20210507_20210517_02_T1),
               a Sentinel-2 product, a Sentinel-1 SLC, a NISAR granule, or for S1BURST the bursts
               (S1_342239_IW3_20211202T161351_VV_3876-BURST), comma-separated; each is resolved to
               the input above, NISAR granules and bursts through ASF's catalog (CMR), a burst to
               itself in its SLC's zip.
  -r, --reference SCENE..., -s, --secondary SCENE...
               the two scenes as lists, each item a separate argument (or several in one,
               separated by spaces or commas), as HyP3 passes a burst pair's bursts; in place of
               the positional reference and secondary. A list runs to the next option, so the
               output path, if given, comes before it
  output.nc    where the ITS_LIVE product is written; by default its ITS_LIVE name,
               <ref>_X_<sec>_G0120V02_P<valid>.nc, in the working directory

options:
  --scratch DIR        local scratch for the run's tile caches and a remote Sentinel-2 band (default:
                       \$TMPDIR/autorift_itslive);
                       must be on disk, not a RAM-backed tmpfs such as /tmp on Amazon Linux
  --cache-io MODE      pread (default) or mmap, how the scratch tile caches are accessed
  --filter-tile N      side of the filtered-tile cache's tiles in pixels (default 512)
  --heap-hint-gib N    heap size the collector aims to stay under, GiB (default 12; 0 for none)
  --orbit-dir DIR      S1: the directory holding each product's precise orbit (.EOF) file; without it the
                       files are downloaded from the public s1-orbits bucket
  --dem PATH           S1, NL1: a geographic (EPSG:4326) DEM of heights above the WGS84 ellipsoid for the
                       resample's terrain solve, in place of the ITS_LIVE parameter region's DEM plus
                       the EGM2008 geoid
  --tide-model PATH    S1, S1BURST, NL1, NL2: the tide model the offsets over polar floating ice and ocean
                       are corrected with, for the ocean tide and, from ERA5 sea level pressure, the
                       inverse barometer, in place of the region's own (CATS2008_v2023 in Antarctica,
                       Gr1kmTM in Greenland, both in s3://its-live-data/tidal/)
  --no-tide-correction leave the offsets over floating ice and ocean as measured
  --subset WKT         S1BURST from SLC zips: a POLYGON or MULTIPOLYGON of longitude-latitude vertices,
                       one crossing the antimeridian split at ±180° as OPERA's frames are; the pair is
                       the bursts both scenes hold at least half of whose footprint it covers, read out
                       of the zips
  --bucket NAME        upload the product, its browse image (<product>.png), its thumbnail
                       (<product>_thumb.png) and its catalog records (<product>.stac.json,
                       <product>.nc.premet, <product>.nc.spatial) to this S3 bucket, with the AWS
                       credential chain's credentials
  --bucket-prefix P    the key prefix they are uploaded under (default none)
  --publish-bucket B   publish the six to the ITS_LIVE Open Data bucket B, under
                       velocity_image_pair/<platform>/v02/<region>, with the credentials in
                       PUBLISH_ACCESS_KEY_ID and PUBLISH_SECRET_ACCESS_KEY; with --bucket,
                       publish_info.json, which records where, is uploaded there too
  --quiet              no per-stage timings on stderr

FFTW wisdom: unless AUTORIFT_FFTW_WISDOM_DIR names a directory to read and grow, each run starts from a
copy of the wisdom the bundle ships for this CPU (share/autorift_itslive/fftw-wisdom), so runs
plan their transforms alike, and discards what it measures.

Landsat is read from s3://usgs-landsat, which is requester-pays: the AWS credential chain (environment,
AWS_PROFILE, container or instance role) must name an account that can pay. A NISAR granule over https://
from an Earthdata host is read with the token in EARTHDATA_TOKEN, or one obtained with the
urs.earthdata.nasa.gov credentials in EARTHDATA_USERNAME and EARTHDATA_PASSWORD, or else in ~/.netrc (or
the file NETRC names). Sentinel-1 bursts are fetched from ASF's burst extractor
with those credentials, waiting while the extractor generates them, or inflated out of their SLC
zips, which are read with the same credentials; either way into the run's scratch directory.
"""

"""
    Arguments

A parsed command line: the job's inputs, the output path, the run options, the heap hint, the burst
subset polygon, the bucket and prefix the product is uploaded to, and the Open Data bucket it is published
to. The sensor, the output path, the polygon and the buckets are empty when not given: the sensor is then the reference's
([`infer_sensor`](@ref)), the product named as ITS_LIVE names it in the working directory, and nothing
uploaded.
"""
struct Arguments
    sensor::String
    reference::String
    secondary::String
    output::String
    options::RunOptions
    heap_hint_gib::Float64
    subset::String
    bucket::String
    bucket_prefix::String
    publish_bucket::String
end

# The scenes one value of a `--reference` or `--secondary` list names: separated by white space or commas,
# with `None`, which HyP3 passes for an empty list, naming none.
scene_list_items(v::String) =
    String[String(x) for x in split(v, (' ', '\t', '\n', ','); keepempty = false) if x != "None"]

"""
    parse_args(args) -> Arguments

The command line. Throws an `ArgumentError` naming the problem for anything it cannot use.
"""
function parse_args(args::Vector{String})
    positional = String[]
    scratch = joinpath(tempdir(), "autorift_itslive")
    cache_io = :pread
    filter_tile = 512
    hint = HEAP_HINT_GIB
    verbose = true
    dem = ""
    orbit_dir = ""
    tide = true
    tide_model = ""
    subset = ""
    bucket = ""
    bucket_prefix = ""
    publish_bucket = ""
    lists = (reference = String[], secondary = String[])
    i = 1
    while i <= length(args)
        a = args[i]
        if a in ("-r", "--reference", "-s", "--secondary")
            # Every value up to the next option: the scenes, as HyP3 passes a list.
            list = a in ("-r", "--reference") ? lists.reference : lists.secondary
            i += 1
            while i <= length(args) && !startswith(args[i], "-")
                append!(list, scene_list_items(args[i]))
                i += 1
            end
            isempty(list) && throw(ArgumentError(a * " needs at least one scene"))
        elseif a in ("--scratch", "--cache-io", "--filter-tile", "--heap-hint-gib", "--dem", "--orbit-dir", "--subset",
                 "--bucket", "--bucket-prefix", "--publish-bucket", "--tide-model")
            i < length(args) || throw(ArgumentError(a * " needs a value"))
            v = args[i + 1]
            if a == "--scratch"
                scratch = v
            elseif a == "--dem"
                dem = v
            elseif a == "--orbit-dir"
                orbit_dir = v
            elseif a == "--tide-model"
                # A GDAL virtual path is already what GDAL opens.
                tide_model = startswith(v, "/vsi") ? v : gdal_path(v)
            elseif a == "--subset"
                subset = v
            elseif a == "--bucket"
                bucket = v
            elseif a == "--bucket-prefix"
                bucket_prefix = v
            elseif a == "--publish-bucket"
                # HyP3 passes a null job parameter as `None`.
                publish_bucket = v == "None" ? "" : v
            elseif a == "--cache-io"
                v in ("pread", "mmap") || throw(ArgumentError("--cache-io is pread or mmap, got " * v))
                cache_io = v == "pread" ? :pread : :mmap
            elseif a == "--filter-tile"
                n = tryparse(Int, v)
                (n === nothing || n < 1) && throw(ArgumentError("--filter-tile needs a positive integer, got " * v))
                filter_tile = n
            else
                h = tryparse(Float64, v)
                (h === nothing || h < 0) && throw(ArgumentError("--heap-hint-gib needs a nonnegative number, got " * v))
                hint = h
            end
            i += 2
        elseif a == "--quiet"
            verbose = false
            i += 1
        elseif a == "--no-tide-correction"
            tide = false
            i += 1
        elseif startswith(a, "--")
            throw(ArgumentError("unknown option " * a))
        else
            push!(positional, a)
            i += 1
        end
    end
    (tide || isempty(tide_model)) || throw(ArgumentError(
        "--tide-model names a tide model and --no-tide-correction turns the correction off; give one of them"))
    # The sensor, when given, comes first; no scene input is a sensor code.
    sensor = !isempty(positional) && is_sensor_code(positional[1]) ? uppercase(popfirst!(positional)) : ""
    flagged = !isempty(lists.reference) || !isempty(lists.secondary)
    if flagged
        (isempty(lists.reference) || isempty(lists.secondary)) && throw(ArgumentError(
            "--reference and --secondary are given together, but only one of them is"))
        length(positional) <= 1 || throw(ArgumentError(
            "with --reference and --secondary, expected [<sensor>] [<output.nc>], got " *
            string(length(positional) + !isempty(sensor)) * " positional arguments"))
        reference, secondary = join(lists.reference, ","), join(lists.secondary, ",")
        output = isempty(positional) ? "" : positional[1]
    else
        length(positional) in (2, 3) || throw(ArgumentError(
            "expected [<sensor>] <reference> <secondary> [<output.nc>], got " * string(length(positional) +
            !isempty(sensor)) * " positional arguments"))
        reference, secondary = positional[1], positional[2]
        output = length(positional) == 3 ? positional[3] : ""
    end
    (isempty(bucket) && !isempty(bucket_prefix)) && throw(ArgumentError("--bucket-prefix needs --bucket"))
    return Arguments(sensor, reference, secondary, output,
                     RunOptions(scratch, cache_io, filter_tile, verbose, dem, orbit_dir, tide, tide_model), hint, subset, bucket,
                     bucket_prefix, publish_bucket)
end

# An exception's message without Base's printer, which a trimmed binary does not have for an arbitrary
# exception: the message of the exception types the chain throws, else the exception's type name. A failure
# in a spawned task or a `@threads` loop arrives wrapped, and is reported as the exception it wraps.
# Inlined into each `catch`, and unwrapping in a loop rather than by recursion: the exception has no static
# type, and a call taking it is one `--trim` cannot resolve.
@inline function error_message(wrapped)
    err = wrapped
    for _ in 1:16
        if err isa TaskFailedException
            err = err.task.exception
        elseif err isa CompositeException && !isempty(err.exceptions)
            err = first(err.exceptions)
        else
            break
        end
    end
    m = err isa ArgumentError ? err.msg :
        err isa ErrorException ? err.msg :
        err isa DimensionMismatch ? err.msg :
        err isa SystemError ? err.prefix * ": " * Libc.strerror(err.errnum) :
        err isa Base.IOError ? err.msg :
        err isa GDAL.GDALError ? err.msg : nothing
    return m isa String ? m : string(nameof(typeof(err)))
end

"""
    julia_main(args = ARGS) -> Cint

The executable's entry point: parse the command line, run the job, and return `0`, or print the error and
return `1` (`2` for a usage error).
"""
function julia_main(args::Vector{String} = ARGS)::Cint
    if isempty(args) || args[1] in ("-h", "--help")
        write(Core.stderr, USAGE)
        return isempty(args) ? Cint(2) : Cint(0)
    end
    a = try
        parse_args(args)
    catch err
        say("error: " * error_message(err))
        return Cint(2)
    end
    try
        release_freed_memory!()
        set_heap_hint!(a.heap_hint_gib)
        AutoRIFT.WORKSPACE_POOL_BYTES[] = WORKSPACE_POOL_BYTES
        with_earthdata_netrc() do
            sensor = isempty(a.sensor) ? infer_sensor(a.reference; a.subset) : a.sensor
            job = Job(sensor, a.reference, a.secondary; subset = a.subset)
            path = run_job(job, a.output; options = a.options)
            say("wrote " * path)
            for uri in deliver(path, job.platform; a.bucket, a.bucket_prefix, a.publish_bucket)
                say("uploaded " * uri)
            end
        end
        return Cint(0)
    catch err
        say("error: " * error_message(err))
        return Cint(1)
    end
end

# The path of the running executable.
function _exepath()
    buf = Vector{UInt8}(undef, 4096)
    n = Ref{Csize_t}(length(buf))
    ccall(:uv_exepath, Cint, (Ptr{UInt8}, Ref{Csize_t}), buf, n) == 0 ||
        throw(ErrorException("cannot find the path of the running executable"))
    return String(buf[1:n[]])
end

# Set to the per-run wisdom directory `ensure_environment` made, which the run deletes when it ends.
const TEMP_WISDOM_VAR = "AUTORIFT_ITSLIVE_TEMP_WISDOM"

# The FFTW wisdom the bundle ships: `share/autorift_itslive/fftw-wisdom` beside `bin/`, one file
# per CPU model, named as `AutoRIFT.fftw_wisdom_path` names it.
bundled_wisdom_dir(exe::String) = joinpath(dirname(dirname(exe)), "share", "autorift_itslive", "fftw-wisdom")

"""
    ensure_environment(args)

Two settings a `juliac` executable reads only when its runtime starts, so the executable runs itself again
with them set and the same arguments when either is missing:

  * `JULIA_NUM_THREADS`, without which it runs on one thread: `auto`, every core the process may use. A
    caller limits the threads by setting it.
  * `AUTORIFT_FFTW_WISDOM_DIR`, which `AutoRIFT` reads its FFTW wisdom from when it loads: a new temporary
    directory holding a copy of [`bundled_wisdom_dir`](@ref)'s files, deleted when the run ends. Every
    run therefore starts from the same wisdom, and nothing a run measures outlives it. A caller keeps its
    own wisdom by setting the variable.
"""
function ensure_environment(args::Vector{String})
    need_threads = !haskey(ENV, "JULIA_NUM_THREADS")
    need_wisdom = !haskey(ENV, "AUTORIFT_FFTW_WISDOM_DIR")
    (need_threads || need_wisdom) || return nothing
    exe = _exepath()
    need_threads && (ENV["JULIA_NUM_THREADS"] = "auto")
    if need_wisdom
        dir = mktempdir(; cleanup = false)
        src = bundled_wisdom_dir(exe)
        if isdir(src)
            for f in readdir(src)
                endswith(f, ".wisdom") && cp(joinpath(src, f), joinpath(dir, f))
            end
        end
        ENV["AUTORIFT_FFTW_WISDOM_DIR"] = dir
        ENV[TEMP_WISDOM_VAR] = dir
    end
    argv = String[exe; args]
    ptrs = Ptr{UInt8}[pointer(a) for a in argv]
    push!(ptrs, Ptr{UInt8}(C_NULL))
    GC.@preserve argv ptrs ccall(:execv, Cint, (Cstring, Ptr{Ptr{UInt8}}), exe, ptrs)
    throw(ErrorException("could not run " * exe * " again with its thread count and FFTW wisdom set"))
end

# The per-run wisdom directory, removed once the run is over.
function remove_temp_wisdom()
    dir = get(ENV, TEMP_WISDOM_VAR, "")
    (!isempty(dir) && dir == get(ENV, "AUTORIFT_FFTW_WISDOM_DIR", "")) || return nothing
    rm(dir; recursive = true, force = true)
    return nothing
end

# The entry point `juliac` builds an executable from.
function (@main)(args::Vector{String})::Cint
    try
        ensure_environment(args)
    catch err
        say("error: " * error_message(err))
        return Cint(1)
    end
    try
        return julia_main(args)
    finally
        remove_temp_wisdom()
    end
end
