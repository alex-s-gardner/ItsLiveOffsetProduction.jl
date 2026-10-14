# A remote Sentinel-2 band, copied into the run's scratch directory before anything reads it. Read in place,
# a JPEG 2000 band costs a request for each tile-part header and then for each tile.

# The size of each range asked for. The transport merges adjacent ranges into requests of up to 16 MiB and
# keeps 32 in flight; one range larger than that would be a single request.
const STAGE_RANGE_BYTES = 4 << 20

"""
    stage_band(s::SceneSource, dir) -> SceneSource

`s` with its band copied to a file under `dir`, when the band is an `http(s)://` URL; `s` itself otherwise.
"""
function stage_band(s::SceneSource, dir::String)
    is_http(s.band) || return s
    t = ChunkManifests.HTTPTransport()
    total = Int(ChunkManifests.objectsize(t, s.band))
    ranges = ChunkManifests.ByteRange[
        ChunkManifests.ByteRange(UInt64(o), UInt64(min(STAGE_RANGE_BYTES, total - o)))
        for o in 0:STAGE_RANGE_BYTES:(total - 1)]
    path = joinpath(dir, s.name * "_B08.jp2")
    io = open(path, "w")
    try
        for block in ChunkManifests.fetchranges(t, s.band, ranges)
            write(io, block)
        end
    finally
        close(io)
    end
    return SceneSource(s.input, s.name, path, s.mtl, s.ang)
end

# Wait for `t` to stop, whether it succeeded or failed; `_await` reports its error.
function _settle(t::Task)
    try
        wait(t)
    catch
    end
    return nothing
end

"""
    stage_bands(c::Job, dir) -> Job

`c` with the bands of a Sentinel-2 pair copied to `dir` ([`stage_band`](@ref)), both at once; `c` itself for
any other sensor or when both bands are already local.
"""
function stage_bands(c::Job, dir::String)
    c.platform == "S2" || return c
    (is_http(c.reference.band) || is_http(c.secondary.band)) || return c
    a = Threads.@spawn stage_band(c.reference, dir)
    b = Threads.@spawn stage_band(c.secondary, dir)
    # Neither is writing into `dir` any longer when an error is reported.
    _settle(a)
    _settle(b)
    return Job(c.sensor, c.platform, _await(a)::SceneSource, _await(b)::SceneSource, c.product)
end
