# The chain on one job: geometry, imagery, grid, correlation, and the product written to disk.

"""
    RunOptions(scratch, cache_io, filter_tile, verbose, dem, orbit_dir, tide_model)

How a run uses local scratch: `scratch` holds the run's tile caches; `cache_io` (`:pread` or `:mmap`) and
`filter_tile` are `TileCache`'s `io` and `AutoRIFT.FilterTileCache`'s tile side. `dem`, when not empty, is
the geographic DEM of ellipsoidal heights a radar resample solves terrain against, in place of
[`parameter_dem`](@ref). `orbit_dir`, when not empty, holds the Sentinel-1 orbit files, in place of those
[`fetch_s1_orbit`](@ref) downloads. `tide_model` is the tide model a radar pair's floating ice is
corrected with ([`TIDE_MODEL`](@ref) by default); empty, no correction is made.
"""
struct RunOptions
    scratch::String
    cache_io::Symbol
    filter_tile::Int
    verbose::Bool
    dem::String
    orbit_dir::String
    tide_model::String
end

RunOptions() = RunOptions(joinpath(tempdir(), "autorift_itslive"), :pread, 512, true, "", "", TIDE_MODEL)

say(s::String) = (write(Core.stderr, s); write(Core.stderr, "\n"); nothing)

_seconds(t0::UInt64) = string(round((time_ns() - t0) / 1e9; digits = 1)) * " s"

_log(o::RunOptions, stage::String, t0::UInt64) = o.verbose && say(rpad_ascii(stage, 11) * _seconds(t0))

# `rpad` for ASCII text: `rpad` and `^` measure `textwidth`, which `juliac --trim` cannot compile.
rpad_ascii(s::String, n::Int) = ncodeunits(s) >= n ? s : s * String(fill(UInt8(' '), n - ncodeunits(s)))

"""
    run_job(c::Job, output; options = RunOptions()) -> String

Run `c` end to end and write its ITS_LIVE product to `output`, returning the path. An empty `output` is the
product's ITS_LIVE name ([`ItsLiveOffsetProduction.product_filename`](@ref)) in the working directory. Every scene is read
from its source; nothing is staged locally except the run's tile caches and a remote Sentinel-2 band, under
a directory of its own, deleted when the run ends, whether it succeeds or not. A scratch directory on a
RAM-backed filesystem is refused before anything is read: see [`ram_backed_mount`](@ref).

A run killed by a signal (SIGTERM from a Spot interruption or AWS Batch, SIGKILL from the kernel's OOM
killer) cannot delete its directory, so each run holds a lock on it ([`lock_run_dir`](@ref)) and the next
run on the same scratch removes every directory whose lock nobody holds ([`remove_stale_runs`](@ref)).
"""
function run_job(c::Job, output::String; options::RunOptions = RunOptions())
    # Assigned inside the `try`, whose scope is its own.
    path = ""
    mkpath(options.scratch)
    fs = ram_backed_mount(options.scratch)
    isempty(fs) || throw(ArgumentError(
        "the scratch directory " * options.scratch * " is on a RAM-backed " * fs * " filesystem, where the " *
        "run's tile caches (gigabytes for a radar pair) would count against memory; pass --scratch a " *
        "directory on disk"))
    remove_stale_runs(options.scratch)
    dir = mktempdir(options.scratch; cleanup = false)
    held = lock_run_dir(dir)
    try
        if startswith(c.platform, "S1")
            t0 = time_ns()
            ss = s1_setup(c, options.orbit_dir, dir)
            _log(options, "geometry", t0)
            t0 = time_ns()
            img = try
                isempty(options.dem) ? s1_imagery(ss, parameter_dem(ss); dir, cache_io = options.cache_io) :
                s1_imagery(ss, geographic_dem(options.dem); dir, cache_io = options.cache_io)
            catch
                # A read that failed because a background fetch did reports that fetch's own error.
                foreach(_await, ss.pending)
                rethrow()
            end
            foreach(_await, ss.pending)
            _log(options, "imagery", t0)
            path = _correlate(c, ss, img, output, options, dir)
        elseif c.platform == "NISAR-L1"
            t0 = time_ns()
            rs = radar_setup(c)
            _log(options, "geometry", t0)
            t0 = time_ns()
            # The reference is read while the DEM is.
            ref = Threads.@spawn radar_reference(rs; dir, cache_io = options.cache_io)
            if isempty(options.dem)
                dem = parameter_dem(rs)
                img = radar_imagery(rs, dem, _await(ref)::TileFile; dir, cache_io = options.cache_io)
            else
                gdem = geographic_dem(options.dem)
                img = radar_imagery(rs, gdem, _await(ref)::TileFile; dir, cache_io = options.cache_io)
            end
            _log(options, "imagery", t0)
            path = _correlate(c, rs, img, output, options, dir)
        else
            t0 = time_ns()
            j = stage_bands(c, dir)
            j === c || _log(options, "download", t0)
            t0 = time_ns()
            s = setup(j)
            _log(options, "geometry", t0)
            t0 = time_ns()
            img = native_imagery(s; dir, cache_io = options.cache_io)
            _log(options, "imagery", t0)
            path = _correlate(j, s, img, output, options, dir)
        end
    finally
        rm(dir; recursive = true, force = true)
        close(held)
    end
    return path
end

const LOCK_EX = Cint(2)
const LOCK_NB = Cint(4)
const RUN_LOCK = ".autorift_run.lock"

_fd(io::IOStream) = reinterpret(Cint, fd(io))
_trylock(io::IOStream) = ccall(:flock, Cint, (Cint, Cint), _fd(io), LOCK_EX | LOCK_NB) == 0

"""
    lock_run_dir(dir) -> IOStream

An exclusive `flock` on `dir`'s lock file, held as long as the returned stream is open. The kernel releases
it when the process ends however it ends, which is what tells [`remove_stale_runs`](@ref) the run is gone.
The lock is taken on a file of another name and the file then renamed, so the lock file exists only while
locked: a run starting meanwhile cannot find it unlocked and remove a directory still in use.
"""
function lock_run_dir(dir::String)
    staging = joinpath(dir, RUN_LOCK * ".new")
    io = open(staging, "w")
    _trylock(io) || (close(io); throw(ErrorException("cannot lock the run's scratch directory " * dir)))
    mv(staging, joinpath(dir, RUN_LOCK))
    return io
end

"""
    remove_stale_runs(scratch)

Delete each run directory under `scratch` whose lock no process holds: what a run killed by a signal leaves.
A directory without a lock file is not a run's and is left alone.
"""
function remove_stale_runs(scratch::String)
    for name in readdir(scratch)
        d = joinpath(scratch, name)
        lockfile = joinpath(d, RUN_LOCK)
        isfile(lockfile) || continue
        io = open(lockfile, "r")
        stale = try
            _trylock(io)
        finally
            close(io)
        end
        stale && rm(d; recursive = true, force = true)
    end
    return nothing
end

"""
    ram_backed_mount(path) -> String

The filesystem type of the mount `path` lies on when it is held in memory (`tmpfs` or `ramfs`), else `""`:
from `/proc/self/mounts`, the deepest mount point containing `path`. `""` where that file does not exist,
as on macOS. `/tmp` is a `tmpfs` of half the memory on Amazon Linux 2023, which makes the default scratch
directory one.
"""
function ram_backed_mount(path::String)
    mounts = "/proc/self/mounts"
    isfile(mounts) || return ""
    p = realpath(path)
    best, fstype = "", ""
    for line in eachline(mounts)
        f = split(line, ' ')
        length(f) >= 3 || continue
        # `/proc/self/mounts` writes a space in a mount point as `\040`.
        mp = replace(String(f[2]), "\\040" => " ")
        inside = mp == "/" || p == mp || startswith(p, mp * "/")
        (inside && length(mp) > length(best)) || continue
        best, fstype = mp, String(f[3])
    end
    return fstype in ("tmpfs", "ramfs") ? fstype : ""
end

# One concrete call chain per correlator filter; the images are closed when it ends. The correlator's
# filter tile cache goes under `dir`, the run's own scratch directory.
function _correlate(c::Job, s, img::Imagery, output::String, options::RunOptions, dir::String)
    try
        f = correlator_filter(c)
        if f === HIGHPASS
            _correlate_and_package(s, img, AutoRIFT.Highpass(filter_width(c)), options.filter_tile, false,
                                   output, options, dir)
        elseif f === WALLIS_FILL
            _correlate_and_package(s, img, AutoRIFT.WallisGapfill(filter_width(c), WALLIS_STD_CUTOFF), nothing,
                                   true, output, options, dir)
        else
            _correlate_and_package(s, img, AutoRIFT.NoPreprocess(), nothing, false, output, options, dir)
        end
    finally
        close(img)
    end
end

function _correlate_and_package(s, img::Imagery, pre::AutoRIFT.PreprocessMethod,
                                filter_tile::Union{Nothing,Int}, masked::Bool, output::String, o::RunOptions,
                                dir::String)
    t0 = time_ns()
    c = correlate(s, img, pre; filter_tile, cache_io = o.cache_io, masked, scratch = dir)
    _log(o, "correlate", t0)
    t0 = time_ns()
    path = package(s, c, output; tide_model = o.tide_model)
    _log(o, "package", t0)
    return path
end
