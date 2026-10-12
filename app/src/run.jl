# The chain on one job: geometry, imagery, grid, correlation, and the product written to disk.

"""
    RunOptions(scratch, cache_io, filter_tile, verbose, dem)

How a run uses local scratch: `scratch` holds the run's tile caches; `cache_io` (`:pread` or `:mmap`) and
`filter_tile` are `TileCache`'s `io` and `AutoRIFT.FilterTileCache`'s tile side. `dem`, when not empty, is
the geographic DEM of ellipsoidal heights a radar resample solves terrain against, in place of
[`parameter_dem`](@ref). `orbit_dir`, when not empty, holds the Sentinel-1 orbit files, in place of those
[`fetch_s1_orbit`](@ref) downloads.
"""
struct RunOptions
    scratch::String
    cache_io::Symbol
    filter_tile::Int
    verbose::Bool
    dem::String
    orbit_dir::String
end

RunOptions() = RunOptions(joinpath(tempdir(), "autorift_itslive_production"), :pread, 512, true, "", "")

say(s::String) = (write(Core.stderr, s); write(Core.stderr, "\n"); nothing)

_seconds(t0::UInt64) = string(round((time_ns() - t0) / 1e9; digits = 1)) * " s"

_log(o::RunOptions, stage::String, t0::UInt64) = o.verbose && say(rpad_ascii(stage, 11) * _seconds(t0))

# `rpad` for ASCII text: `rpad` and `^` measure `textwidth`, which `juliac --trim` cannot compile.
rpad_ascii(s::String, n::Int) = ncodeunits(s) >= n ? s : s * String(fill(UInt8(' '), n - ncodeunits(s)))

"""
    run_job(c::Job, output; options = RunOptions()) -> String

Run `c` end to end and write its ITS_LIVE product to `output`, returning the path. Every scene is read
from its source; nothing is staged locally except the run's tile caches under a directory of its own,
deleted when the run ends, whether it succeeds or not.
"""
function run_job(c::Job, output::String; options::RunOptions = RunOptions())
    mkpath(options.scratch)
    dir = mktempdir(options.scratch; cleanup = false)
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
            _correlate(c, ss, img, output, options)
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
            _correlate(c, rs, img, output, options)
        else
            t0 = time_ns()
            s = setup(c)
            _log(options, "geometry", t0)
            t0 = time_ns()
            img = native_imagery(s; dir, cache_io = options.cache_io)
            _log(options, "imagery", t0)
            _correlate(c, s, img, output, options)
        end
    finally
        rm(dir; recursive = true, force = true)
    end
    return output
end

# One concrete call chain per correlator filter; the images are closed when it ends.
function _correlate(c::Job, s, img::Imagery, output::String, options::RunOptions)
    try
        f = correlator_filter(c)
        if f === HIGHPASS
            _correlate_and_package(s, img, AutoRIFT.Highpass(filter_width(c)), options.filter_tile, false,
                                   output, options)
        elseif f === WALLIS_FILL
            _correlate_and_package(s, img, AutoRIFT.WallisGapfill(filter_width(c), WALLIS_STD_CUTOFF), nothing,
                                   true, output, options)
        else
            _correlate_and_package(s, img, AutoRIFT.NoPreprocess(), nothing, false, output, options)
        end
    finally
        close(img)
    end
    return nothing
end

function _correlate_and_package(s, img::Imagery, pre::AutoRIFT.PreprocessMethod,
                                filter_tile::Union{Nothing,Int}, masked::Bool, output::String, o::RunOptions)
    t0 = time_ns()
    c = correlate(s, img, pre; filter_tile, cache_io = o.cache_io, masked)
    _log(o, "correlate", t0)
    t0 = time_ns()
    package(s, c, output)
    _log(o, "package", t0)
    return nothing
end
