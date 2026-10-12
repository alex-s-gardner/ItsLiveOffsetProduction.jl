# Trim probe: windows read through ChunkManifests' typed opens — a local Landsat COG, the same band
# from requester-pays S3, and a Sentinel-2 JPEG 2000 band from a local copy and over HTTPS.
#
#   typedread <landsat.TIF> <s3://bucket/key.TIF> <region> <s2.jp2> <https://.../s2.jp2>

import ChunkManifests as CM
import OpenJpeg_jll
import OpticalDatasets

say(s::String) = (write(Core.stdout, s, "\n"); nothing)

function message(err)
    m = err isa ErrorException ? err.msg : err isa ArgumentError ? err.msg : nothing
    return m isa String ? m : string(typeof(err))
end

function window_sums(a::CM.TypedArray{UInt16})
    s = 0
    for (rx, ry) in ((1:2000, 1:1500), (8000:9123, 4000:5301))
        s += sum(Int, a[rx, ry])
    end
    return s
end

function (@main)(args::Vector{String})::Cint
    try
        local_ls = CM.open_tiff(args[1], UInt16; transport = CM.LocalTransport())
        say(string("landsat local size=", size(local_ls, 1), "x", size(local_ls, 2), " sum=", window_sums(local_ls)))
        s3 = CM.S3Transport("usgs-landsat"; region = args[3], requesterpays = true)
        remote = CM.open_tiff(args[2], UInt16; transport = s3)
        say(string("landsat s3 sum=", window_sums(remote)))
        s2 = CM.open_jp2(args[4], UInt16; transport = CM.LocalTransport())
        say(string("s2 jp2 size=", size(s2, 1), "x", size(s2, 2), " sum=", window_sums(s2)))
        web = CM.open_jp2(args[5], UInt16; transport = CM.HTTPTransport())
        say(string("s2 jp2 https sum=", window_sums(web)))
        # OpticalDatasets' typed entry point: the same band, indexed (row, column).
        r = OpticalDatasets.open_optical(args[1], UInt16; transport = CM.LocalTransport())
        rows = OpticalDatasets.read_window(r, 4000:5301, 8000:9123)
        dest = Matrix{Float32}(undef, 1302, 1124)
        OpticalDatasets.read_window!(dest, r, 4000:5301, 8000:9123)
        say(string("optical rows x cols=", size(r, 1), "x", size(r, 2), " window sum=", sum(Int, rows),
                   " as Float32=", sum(Float64, dest)))
        return Cint(0)
    catch err
        say(string("error: ", message(err)))
        return Cint(1)
    end
end
