# A ranged-GET client on libcurl's easy interface, shared by the HTTP and S3 trim probes.

import LibCURL_jll
import NetworkOptions

const libcurl = LibCURL_jll.libcurl
const CURLOPT_URL = Cint(10002)
const CURLOPT_WRITEDATA = Cint(10001)
const CURLOPT_WRITEFUNCTION = Cint(20011)
const CURLOPT_RANGE = Cint(10007)
const CURLOPT_HTTPHEADER = Cint(10023)
const CURLOPT_FOLLOWLOCATION = Cint(52)
const CURLOPT_NETRC = Cint(51)
const CURLOPT_COOKIEFILE = Cint(10031)
const CURLOPT_CAINFO = Cint(10065)
const CURLOPT_NOSIGNAL = Cint(99)
const CURLINFO_RESPONSE_CODE = Cint(0x200002)
const CURL_NETRC_OPTIONAL = Clong(1)

say(s::String) = (write(Core.stdout, s, "\n"); nothing)

function message(err)
    m = err isa ErrorException ? err.msg : nothing
    return m isa String ? m : string(typeof(err))
end

function write_body(ptr::Ptr{UInt8}, size::Csize_t, n::Csize_t, userdata::Ptr{Cvoid})::Csize_t
    buf = unsafe_pointer_to_objref(userdata)::Vector{UInt8}
    len = size * n
    append!(buf, unsafe_wrap(Array, ptr, len))
    return len
end

setopt(h::Ptr{Cvoid}, opt::Cint, v::Clong) =
    ccall((:curl_easy_setopt, libcurl), Cint, (Ptr{Cvoid}, Cint, Clong...), h, opt, v)
setopt(h::Ptr{Cvoid}, opt::Cint, v::String) =
    ccall((:curl_easy_setopt, libcurl), Cint, (Ptr{Cvoid}, Cint, Cstring...), h, opt, v)
setopt(h::Ptr{Cvoid}, opt::Cint, v::Ptr{Cvoid}) =
    ccall((:curl_easy_setopt, libcurl), Cint, (Ptr{Cvoid}, Cint, Ptr{Cvoid}...), h, opt, v)

"""
    get_range(url, first, last; headers = String[], earthdata = false) -> (status, bytes)

Bytes `first:last` (0-based, inclusive, as HTTP spells a range) of `url`. `headers` are full
`"Name: value"` lines; `earthdata` authenticates redirects from `~/.netrc` with in-memory cookies.
"""
function get_range(url::String, first::Int, last::Int; headers::Vector{String} = String[],
                   earthdata::Bool = false)
    h = ccall((:curl_easy_init, libcurl), Ptr{Cvoid}, ())
    h == C_NULL && error("curl_easy_init failed")
    buf = UInt8[]
    range = string(first, "-", last)
    ca = NetworkOptions.ca_roots_path()
    cb = @cfunction(write_body, Csize_t, (Ptr{UInt8}, Csize_t, Csize_t, Ptr{Cvoid}))
    list = C_NULL
    for line in headers
        list = ccall((:curl_slist_append, libcurl), Ptr{Cvoid}, (Ptr{Cvoid}, Cstring), list, line)
    end
    code = Ref{Clong}(0)
    rc = GC.@preserve buf url range ca headers begin
        setopt(h, CURLOPT_URL, url)
        setopt(h, CURLOPT_RANGE, range)
        setopt(h, CURLOPT_NOSIGNAL, Clong(1))
        setopt(h, CURLOPT_FOLLOWLOCATION, Clong(1))
        setopt(h, CURLOPT_CAINFO, ca)
        setopt(h, CURLOPT_WRITEFUNCTION, cb)
        setopt(h, CURLOPT_WRITEDATA, pointer_from_objref(buf))
        list == C_NULL || setopt(h, CURLOPT_HTTPHEADER, list)
        if earthdata
            setopt(h, CURLOPT_NETRC, CURL_NETRC_OPTIONAL)
            setopt(h, CURLOPT_COOKIEFILE, "")
        end
        r = ccall((:curl_easy_perform, libcurl), Cint, (Ptr{Cvoid},), h)
        ccall((:curl_easy_getinfo, libcurl), Cint, (Ptr{Cvoid}, Cint, Ptr{Clong}...), h, CURLINFO_RESPONSE_CODE, code)
        r
    end
    ccall((:curl_easy_cleanup, libcurl), Cvoid, (Ptr{Cvoid},), h)
    list == C_NULL || ccall((:curl_slist_free_all, libcurl), Cvoid, (Ptr{Cvoid},), list)
    rc == 0 || error(string("curl_easy_perform failed with code ", rc, " for ", url))
    return Int(code[]), buf
end
