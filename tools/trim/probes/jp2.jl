# Trim probe: decoding an area of a Sentinel-2 JPEG 2000 band through libopenjp2 directly.
#
#   jp2 <band.jp2> <x0> <y0> <x1> <y1>
#
# Prints the decoded area's size and the sum of its samples.

import OpenJpeg_jll

say(s::String) = (write(Core.stdout, s, "\n"); nothing)

function message(err)
    m = err isa ErrorException ? err.msg : nothing
    return m isa String ? m : string(typeof(err))
end

const OPJ_CODEC_JP2 = Cint(2)
const DPARAMETERS_BYTES = 64 * 1024

# `opj_image_comp_t` and `opj_image_t` from openjpeg.h (2.4 and 2.5 alike).
struct OpjImageComp
    dx::UInt32
    dy::UInt32
    w::UInt32
    h::UInt32
    x0::UInt32
    y0::UInt32
    prec::UInt32
    bpp::UInt32
    sgnd::UInt32
    resno_decoded::UInt32
    factor::UInt32
    data::Ptr{Int32}
    alpha::UInt16
end

struct OpjImage
    x0::UInt32
    y0::UInt32
    x1::UInt32
    y1::UInt32
    numcomps::UInt32
    color_space::Cint
    comps::Ptr{OpjImageComp}
    icc_profile_buf::Ptr{UInt8}
    icc_profile_len::UInt32
end

function decode_area(path::String, x0::Int, y0::Int, x1::Int, y1::Int)
    codec = ccall((:opj_create_decompress, OpenJpeg_jll.libopenjp2), Ptr{Cvoid}, (Cint,), OPJ_CODEC_JP2)
    params = zeros(UInt8, DPARAMETERS_BYTES)
    ccall((:opj_set_default_decoder_parameters, OpenJpeg_jll.libopenjp2), Cvoid, (Ptr{UInt8},), params)
    ccall((:opj_setup_decoder, OpenJpeg_jll.libopenjp2), Cint, (Ptr{Cvoid}, Ptr{UInt8}), codec, params) == 1 ||
        error("opj_setup_decoder failed")
    stream = ccall((:opj_stream_create_default_file_stream, OpenJpeg_jll.libopenjp2), Ptr{Cvoid}, (Cstring, Cint), path, 1)
    stream == C_NULL && error("cannot open the JPEG 2000 file")
    image = Ref{Ptr{OpjImage}}(C_NULL)
    ccall((:opj_read_header, OpenJpeg_jll.libopenjp2), Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Ptr{OpjImage}}),
          stream, codec, image) == 1 || error("opj_read_header failed")
    ccall((:opj_set_decode_area, OpenJpeg_jll.libopenjp2), Cint, (Ptr{Cvoid}, Ptr{OpjImage}, Int32, Int32, Int32, Int32),
          codec, image[], x0, y0, x1, y1) == 1 || error("opj_set_decode_area failed")
    ccall((:opj_decode, OpenJpeg_jll.libopenjp2), Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{OpjImage}),
          codec, stream, image[]) == 1 || error("opj_decode failed")
    img = unsafe_load(image[])
    comp = unsafe_load(img.comps)
    w, h = Int(comp.w), Int(comp.h)
    data = copy(unsafe_wrap(Array, comp.data, w * h))
    ccall((:opj_image_destroy, OpenJpeg_jll.libopenjp2), Cvoid, (Ptr{OpjImage},), image[])
    ccall((:opj_stream_destroy, OpenJpeg_jll.libopenjp2), Cvoid, (Ptr{Cvoid},), stream)
    ccall((:opj_destroy_codec, OpenJpeg_jll.libopenjp2), Cvoid, (Ptr{Cvoid},), codec)
    return w, h, data
end

function (@main)(args::Vector{String})::Cint
    try
        x0, y0, x1, y1 = (parse(Int, a) for a in args[2:5])
        w, h, data = decode_area(args[1], x0, y0, x1, y1)
        say(string("decoded ", w, "x", h, " sum=", sum(Int, data)))
        return Cint(0)
    catch err
        say(string("error: ", message(err)))
        return Cint(1)
    end
end
