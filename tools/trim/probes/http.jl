# Trim probe: ranged HTTPS GETs on libcurl's easy interface, one handle per task.
#
#   http <public-url> <earthdata-url>
#
# Reads four 64 KiB ranges of `public-url` concurrently from four tasks and prints the sum of each
# range's bytes; then reads 100 bytes of `earthdata-url` through the Earthdata Login redirect chain,
# authenticating from `~/.netrc` with libcurl's in-memory cookie engine.

include("curl.jl")

function (@main)(args::Vector{String})::Cint
    try
        # Once, before any task touches libcurl: curl_global_init is not thread-safe.
        ccall((:curl_global_init, libcurl), Cint, (Clong,), Clong(3))
        public, edl = args[1], args[2]
        n = 65536
        tasks = [Threads.@spawn(get_range(public, k * n, (k + 1) * n - 1)) for k in 0:3]
        for (k, t) in enumerate(tasks)
            status, bytes = fetch(t)::Tuple{Int,Vector{UInt8}}
            say(string("range ", k - 1, " status=", status, " length=", length(bytes), " sum=", sum(Int, bytes)))
        end
        status, bytes = get_range(edl, 1000000, 1000099; earthdata = true)
        say(string("earthdata status=", status, " length=", length(bytes), " sum=", sum(Int, bytes)))
        say(string("threads=", Threads.nthreads(), " ca=", NetworkOptions.ca_roots_path()))
        return Cint(0)
    catch err
        say(string("error: ", message(err)))
        return Cint(1)
    end
end
