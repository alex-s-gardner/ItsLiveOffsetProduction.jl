# Trim probe: random access into a deflated zip member through a seekable-deflate index (the method of
# zlib's examples/zran.c), with a typed zip central-directory parser. Input is read only through
# `fetch(src, offset, n)`, so the same code serves a remote zip read by byte ranges.
#
#   zipindex <file.zip> <member-substring> <offset> [<offset> ...]
#
# Builds the index of the first member whose name contains `member-substring`, then prints the sum of
# the 4096 uncompressed bytes at each `offset`.

import Zlib_jll

say(s::String) = (write(Core.stdout, s, "\n"); nothing)

function message(err)
    m = err isa ErrorException ? err.msg : nothing
    return m isa String ? m : string(typeof(err))
end

# ---- byte source ---------------------------------------------------------------------------------

struct FileSource
    io::IOStream
end
function fetch(src::FileSource, offset::Int, n::Int)
    buf = Vector{UInt8}(undef, n)
    seek(src.io, offset)
    read!(src.io, buf)
    return buf
end
source_size(src::FileSource) = Int(filesize(src.io))

le16(b::Vector{UInt8}, i::Int) = Int(b[i]) | Int(b[i + 1]) << 8
le32(b::Vector{UInt8}, i::Int) = le16(b, i) | le16(b, i + 2) << 16
le64(b::Vector{UInt8}, i::Int) = le32(b, i) | le32(b, i + 4) << 32

# ---- zip central directory -----------------------------------------------------------------------

struct ZipMember
    name::String
    method::Int              # 0 stored, 8 deflated
    compressed::Int
    uncompressed::Int
    header::Int              # offset of the local file header
end

function zip_members(src::FileSource)
    total = source_size(src)
    tail_len = min(total, 65536 + 22)
    tail = fetch(src, total - tail_len, tail_len)
    e = findlast(i -> le32(tail, i) == 0x06054b50, 1:(length(tail) - 3))
    e === nothing && error("no zip end-of-central-directory record")
    count, cd_size, cd_off = le16(tail, e + 10), le32(tail, e + 12), le32(tail, e + 16)
    if count == 0xffff || cd_size == 0xffffffff || cd_off == 0xffffffff
        e - 20 >= 1 && le32(tail, e - 20) == 0x07064b50 || error("ZIP64 locator missing")
        z64 = fetch(src, le64(tail, e - 20 + 8), 56)
        le32(z64, 1) == 0x06064b50 || error("ZIP64 end record missing")
        count, cd_size, cd_off = le64(z64, 33), le64(z64, 41), le64(z64, 49)
    end
    cd = fetch(src, cd_off, cd_size)
    members = ZipMember[]
    i = 1
    for _ in 1:count
        le32(cd, i) == 0x02014b50 || error("bad central directory entry")
        method = le16(cd, i + 10)
        csize, usize = le32(cd, i + 20), le32(cd, i + 24)
        nlen, xlen, clen = le16(cd, i + 28), le16(cd, i + 30), le16(cd, i + 32)
        hdr = le32(cd, i + 42)
        name = String(cd[(i + 46):(i + 45 + nlen)])
        # ZIP64 extended information: present values replace the 0xffffffff placeholders, in order.
        x = i + 46 + nlen
        xend = x + xlen
        while x + 3 < xend
            id, sz = le16(cd, x), le16(cd, x + 2)
            if id == 0x0001
                j = x + 4
                usize == 0xffffffff && (usize = le64(cd, j); j += 8)
                csize == 0xffffffff && (csize = le64(cd, j); j += 8)
                hdr == 0xffffffff && (hdr = le64(cd, j))
            end
            x += 4 + sz
        end
        push!(members, ZipMember(name, method, csize, usize, hdr))
        i = xend + clen
    end
    return members
end

function data_offset(src::FileSource, m::ZipMember)
    h = fetch(src, m.header, 30)
    le32(h, 1) == 0x04034b50 || error("bad local file header")
    return m.header + 30 + le16(h, 27) + le16(h, 29)
end

# ---- zlib ----------------------------------------------------------------------------------------

mutable struct ZStream
    next_in::Ptr{UInt8}
    avail_in::Cuint
    total_in::Culong
    next_out::Ptr{UInt8}
    avail_out::Cuint
    total_out::Culong
    msg::Ptr{UInt8}
    state::Ptr{Cvoid}
    zalloc::Ptr{Cvoid}
    zfree::Ptr{Cvoid}
    opaque::Ptr{Cvoid}
    data_type::Cint
    adler::Culong
    reserved::Culong
    ZStream() = new(C_NULL, 0, 0, C_NULL, 0, 0, C_NULL, C_NULL, C_NULL, C_NULL, C_NULL, 0, 0, 0)
end

const Z_OK = Cint(0)
const Z_STREAM_END = Cint(1)
const Z_NEED_DICT = Cint(2)
const Z_NO_FLUSH = Cint(0)
const Z_BLOCK = Cint(5)
const WINSIZE = 32768
const CHUNK = 1 << 20

function inflate_init!(z::ZStream)
    v = ccall((:zlibVersion, Zlib_jll.libz), Cstring, ())
    rc = ccall((:inflateInit2_, Zlib_jll.libz), Cint, (Ref{ZStream}, Cint, Cstring, Cint), z, -15, v, sizeof(ZStream))
    rc == Z_OK || error(string("inflateInit2 failed: ", rc))
    return z
end
inflate!(z::ZStream, flush::Cint) = ccall((:inflate, Zlib_jll.libz), Cint, (Ref{ZStream}, Cint), z, flush)
inflate_end!(z::ZStream) = ccall((:inflateEnd, Zlib_jll.libz), Cint, (Ref{ZStream},), z)

# ---- the index -----------------------------------------------------------------------------------

struct AccessPoint
    out::Int                 # uncompressed offset
    in::Int                  # compressed offset of the first full byte after the point
    bits::Int                # bits of the preceding byte that belong to the point, 0–7
    window::Vector{UInt8}    # the 32 KiB of output before `out`
end

struct DeflateIndex
    start::Int               # offset of the compressed data in the source
    compressed::Int
    points::Vector{AccessPoint}
end

function build_index(src::FileSource, start::Int, compressed::Int, span::Int)
    z = inflate_init!(ZStream())
    window = zeros(UInt8, WINSIZE)
    # The start of the stream is an access point that needs neither a primed bit nor a dictionary.
    points = AccessPoint[AccessPoint(0, 0, 0, UInt8[])]
    pos, remaining = start, compressed
    totin, totout, last = 0, 0, 0
    input = UInt8[]
    rc = Z_OK
    GC.@preserve window begin
        while rc != Z_STREAM_END
            if z.avail_in == 0
                remaining > 0 || error("deflate stream ends before its end-of-stream block")
                n = min(CHUNK, remaining)
                input = fetch(src, pos, n)
                pos += n
                remaining -= n
                z.next_in = pointer(input)
                z.avail_in = n
            end
            # Inflate is called again, with no new input, in two cases: a full output buffer may leave
            # output pending after the last input byte is consumed, and `Z_BLOCK` stops at the end of the
            # final block (data_type bits 64 and 128) one call before reporting `Z_STREAM_END`.
            GC.@preserve input while z.avail_in != 0 || z.avail_out == 0 || (z.data_type & 192) == 192
                if z.avail_out == 0
                    z.next_out = pointer(window)
                    z.avail_out = WINSIZE
                end
                totin += Int(z.avail_in)
                totout += Int(z.avail_out)
                rc = inflate!(z, Z_BLOCK)
                totin -= Int(z.avail_in)
                totout -= Int(z.avail_out)
                rc == Z_NEED_DICT && error("raw deflate stream asked for a dictionary")
                rc < 0 && error(string("inflate failed: ", rc))
                rc == Z_STREAM_END && break
                at_block_boundary = (z.data_type & 128) != 0 && (z.data_type & 64) == 0
                if at_block_boundary && totout - last > span
                    left = Int(z.avail_out)
                    w = vcat(window[(WINSIZE - left + 1):WINSIZE], window[1:(WINSIZE - left)])
                    push!(points, AccessPoint(totout, totin, Int(z.data_type & 7), w))
                    last = totout
                end
            end
        end
    end
    inflate_end!(z)
    return DeflateIndex(start, compressed, points)
end

function extract(src::FileSource, idx::DeflateIndex, offset::Int, len::Int)
    k = findlast(p -> p.out <= offset, idx.points)
    k === nothing && error("offset before the first access point")
    p = idx.points[k]
    z = inflate_init!(ZStream())
    pos = idx.start + p.in - (p.bits > 0 ? 1 : 0)
    remaining = idx.start + idx.compressed - pos
    input = fetch(src, pos, min(CHUNK, remaining))
    pos += length(input)
    remaining -= length(input)
    first_in = 1
    if p.bits > 0
        rc = ccall((:inflatePrime, Zlib_jll.libz), Cint, (Ref{ZStream}, Cint, Cint),
                   z, p.bits, Cint(input[1] >> (8 - p.bits)))
        rc == Z_OK || error("inflatePrime failed")
        first_in = 2
    end
    if p.out > 0
        rc = GC.@preserve p ccall((:inflateSetDictionary, Zlib_jll.libz), Cint, (Ref{ZStream}, Ptr{UInt8}, Cuint),
                                  z, pointer(p.window), WINSIZE)
        rc == Z_OK || error("inflateSetDictionary failed")
    end
    skip = offset - p.out
    discard = Vector{UInt8}(undef, WINSIZE)
    out = Vector{UInt8}(undef, len)
    got = 0
    GC.@preserve input discard out begin
        z.next_in = pointer(input, first_in)
        z.avail_in = length(input) - first_in + 1
        while got < len
            if z.avail_in == 0
                remaining > 0 || error("member ends before the requested range")
                input = fetch(src, pos, min(CHUNK, remaining))
                pos += length(input)
                remaining -= length(input)
                z.next_in = pointer(input)
                z.avail_in = length(input)
            end
            if skip > 0
                n = min(skip, WINSIZE)
                z.next_out = pointer(discard)
                z.avail_out = n
                rc = GC.@preserve input inflate!(z, Z_NO_FLUSH)
                skip -= n - Int(z.avail_out)
            else
                z.next_out = pointer(out, got + 1)
                z.avail_out = len - got
                rc = GC.@preserve input inflate!(z, Z_NO_FLUSH)
                got = len - Int(z.avail_out)
            end
            rc < 0 && error(string("inflate failed: ", rc))
            rc == Z_STREAM_END && got < len && skip == 0 && error("member ends before the requested range")
        end
    end
    inflate_end!(z)
    return out
end

function (@main)(args::Vector{String})::Cint
    try
        src = FileSource(open(args[1], "r"))
        members = zip_members(src)
        k = findfirst(m -> occursin(args[2], m.name), members)
        k === nothing && error(string("no member matching ", args[2]))
        m = members[k]
        say(string("members=", length(members), " member=", m.name, " method=", m.method,
                   " compressed=", m.compressed, " uncompressed=", m.uncompressed))
        m.method == 8 || error("member is not deflated")
        t = time()
        idx = build_index(src, data_offset(src, m), m.compressed, 4 << 20)
        say(string("index points=", length(idx.points), " window bytes=", length(idx.points) * WINSIZE,
                   " build seconds=", round(time() - t; digits = 1)))
        for a in args[3:end]
            off = parse(Int, a)
            t = time()
            b = extract(src, idx, off, 4096)
            say(string("offset ", off, " sum=", sum(Int, b), " seconds=", round(time() - t; digits = 3)))
        end
        close(src.io)
        return Cint(0)
    catch err
        say(string("error: ", message(err)))
        return Cint(1)
    end
end
