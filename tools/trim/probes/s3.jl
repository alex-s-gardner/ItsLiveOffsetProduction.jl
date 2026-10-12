# Trim probe: S3 ranged GETs signed with AWS Signature Version 4, without an AWS SDK.
#
#   s3 <region> <bucket> <key> <first> <last> [requester-pays: yes|no]
#
# Credentials come from AWS_ACCESS_KEY_ID / AWS_SECRET_ACCESS_KEY / AWS_SESSION_TOKEN, else from the
# AWS_PROFILE (default "default") section of ~/.aws/credentials. Prints the response status, length
# and the sum of the bytes.

include("curl.jl")
import SHA
using Dates

struct Credentials
    access_key::String
    secret_key::String
    token::String          # empty when the credentials are long-lived
end

function profile_credentials(path::String, profile::String)
    isfile(path) || return nothing
    section = ""
    found = Dict{String,String}()
    for raw in eachline(path)
        line = strip(raw)
        (isempty(line) || startswith(line, '#') || startswith(line, ';')) && continue
        if startswith(line, '[') && endswith(line, ']')
            section = String(strip(line[2:(end - 1)]))
            continue
        end
        section == profile || continue
        i = findfirst('=', line)
        i === nothing && continue
        found[String(strip(line[1:(i - 1)]))] = String(strip(line[(i + 1):end]))
    end
    haskey(found, "aws_access_key_id") && haskey(found, "aws_secret_access_key") || return nothing
    return Credentials(found["aws_access_key_id"], found["aws_secret_access_key"],
                       get(found, "aws_session_token", ""))
end

function credentials()
    akid = get(ENV, "AWS_ACCESS_KEY_ID", "")
    secret = get(ENV, "AWS_SECRET_ACCESS_KEY", "")
    isempty(akid) || isempty(secret) || return Credentials(akid, secret, get(ENV, "AWS_SESSION_TOKEN", ""))
    profile = get(ENV, "AWS_PROFILE", "default")
    file = get(ENV, "AWS_SHARED_CREDENTIALS_FILE", joinpath(homedir(), ".aws", "credentials"))
    c = profile_credentials(file, profile)
    c === nothing && error(string("no AWS credentials: none in the environment and no profile \"",
                                  profile, "\" in ", file))
    return c
end

hmac(key::Vector{UInt8}, data::String) = SHA.hmac_sha256(key, Vector{UInt8}(codeunits(data)))
hexsha(data::String) = bytes2hex(SHA.sha256(data))

# S3 object keys are signed with each path segment percent-encoded, '/' left as is.
function encode_key(key::String)
    io = IOBuffer()
    for b in codeunits(key)
        c = Char(b)
        if ('A' <= c <= 'Z') || ('a' <= c <= 'z') || ('0' <= c <= '9') || c in ('-', '_', '.', '~', '/')
            write(io, b)
        else
            write(io, '%', uppercase(string(b; base = 16, pad = 2)))
        end
    end
    return String(take!(io))
end

"""
    signed_headers(c, region, bucket, key, range, requester_pays, now) -> (url, header lines)

The virtual-hosted URL of `bucket/key` and the request headers of a SigV4-signed GET of `range`.
"""
function signed_headers(c::Credentials, region::String, bucket::String, key::String, range::String,
                        requester_pays::Bool, now::DateTime)
    host = string(bucket, ".s3.", region, ".amazonaws.com")
    path = string("/", encode_key(key))
    amzdate = Dates.format(now, dateformat"yyyymmdd\THHMMSS\Z")
    date = amzdate[1:8]
    payload = "UNSIGNED-PAYLOAD"
    headers = Pair{String,String}["host" => host, "range" => string("bytes=", range),
                                  "x-amz-content-sha256" => payload, "x-amz-date" => amzdate]
    requester_pays && push!(headers, "x-amz-request-payer" => "requester")
    isempty(c.token) || push!(headers, "x-amz-security-token" => c.token)
    sort!(headers; by = first)
    canonical_headers = join([string(k, ":", v, "\n") for (k, v) in headers])
    signed = join([k for (k, _) in headers], ";")
    canonical = string("GET\n", path, "\n\n", canonical_headers, "\n", signed, "\n", payload)
    scope = string(date, "/", region, "/s3/aws4_request")
    to_sign = string("AWS4-HMAC-SHA256\n", amzdate, "\n", scope, "\n", hexsha(canonical))
    k = hmac(Vector{UInt8}(codeunits(string("AWS4", c.secret_key))), date)
    k = hmac(k, region)
    k = hmac(k, "s3")
    k = hmac(k, "aws4_request")
    signature = bytes2hex(hmac(k, to_sign))
    auth = string("AWS4-HMAC-SHA256 Credential=", c.access_key, "/", scope, ", SignedHeaders=", signed,
                  ", Signature=", signature)
    lines = [string(k, ": ", v) for (k, v) in headers if k != "host" && k != "range"]
    push!(lines, string("Authorization: ", auth))
    return string("https://", host, path), lines
end

function (@main)(args::Vector{String})::Cint
    try
        ccall((:curl_global_init, libcurl), Cint, (Clong,), Clong(3))
        region, bucket, key = args[1], args[2], args[3]
        first, last = parse(Int, args[4]), parse(Int, args[5])
        rp = length(args) >= 6 && args[6] == "yes"
        c = credentials()
        range = string(first, "-", last)
        url, headers = signed_headers(c, region, bucket, key, range, rp, Dates.now(Dates.UTC))
        status, bytes = get_range(url, first, last; headers)
        say(string("status=", status, " length=", length(bytes), " sum=", sum(Int, bytes)))
        status == 206 || say(String(copy(bytes)))
        return status == 206 ? Cint(0) : Cint(1)
    catch err
        say(string("error: ", message(err)))
        return Cint(1)
    end
end
