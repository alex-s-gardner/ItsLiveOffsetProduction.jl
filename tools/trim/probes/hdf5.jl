# Trim probe: libhdf5 through direct ccalls on HDF5_jll, without HDF5.jl — the calls ChunkManifests'
# HDF5 driver needs to scan and read a NISAR product.
#
#   hdf5 <nisar.h5> <dataset> <row> <col>
#
# Reads the fixed-length `identification/missionId` string, lists the `identification` group's links,
# reports `dataset`'s filter pipeline and chunk count, and reads a 4×4 complex window at (row, col)
# (0-based, HDF5's row-major order).

import HDF5_jll
import Libdl

# Every ccall names `HDF5_jll.libhdf5` itself: the JLL sets that global in its `__init__`, so a `const`
# copy taken at build time would hold the build-time value, not the bundled library.

const hid_t = Int64
const H5P_DEFAULT = hid_t(0)
const H5S_ALL = hid_t(0)
const H5F_ACC_RDONLY = Cuint(0)
const H5_INDEX_NAME = Cint(0)
const H5_ITER_INC = Cint(0)
const H5T_COMPOUND = Cint(6)
const H5S_SELECT_SET = Cint(0)

say(s::String) = (write(Core.stdout, s, "\n"); nothing)

function message(err)
    m = err isa ErrorException ? err.msg : nothing
    return m isa String ? m : string(typeof(err))
end

checked(id::hid_t, what::String) = id < 0 ? error(string("HDF5: ", what, " failed")) : id
checked(status::Cint, what::String) = status < 0 ? error(string("HDF5: ", what, " failed")) : status

function read_fixed_string(f::hid_t, path::String)
    d = checked(ccall((:H5Dopen2, HDF5_jll.libhdf5), hid_t, (hid_t, Cstring, hid_t), f, path, H5P_DEFAULT), "H5Dopen2 $path")
    t = checked(ccall((:H5Dget_type, HDF5_jll.libhdf5), hid_t, (hid_t,), d), "H5Dget_type")
    n = ccall((:H5Tget_size, HDF5_jll.libhdf5), Csize_t, (hid_t,), t)
    buf = zeros(UInt8, n)
    checked(ccall((:H5Dread, HDF5_jll.libhdf5), Cint, (hid_t, hid_t, hid_t, hid_t, hid_t, Ptr{UInt8}),
                  d, t, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf), "H5Dread $path")
    ccall((:H5Tclose, HDF5_jll.libhdf5), Cint, (hid_t,), t)
    ccall((:H5Dclose, HDF5_jll.libhdf5), Cint, (hid_t,), d)
    k = findfirst(iszero, buf)
    return String(buf[1:(k === nothing ? length(buf) : k - 1)])
end

function link_name(group::hid_t, name::Cstring, info::Ptr{Cvoid}, op::Ptr{Cvoid})::Cint
    names = unsafe_pointer_to_objref(op)::Vector{String}
    push!(names, unsafe_string(name))
    return Cint(0)
end

function links(f::hid_t, path::String)
    g = checked(ccall((:H5Gopen2, HDF5_jll.libhdf5), hid_t, (hid_t, Cstring, hid_t), f, path, H5P_DEFAULT), "H5Gopen2")
    names = String[]
    cb = @cfunction(link_name, Cint, (hid_t, Cstring, Ptr{Cvoid}, Ptr{Cvoid}))
    GC.@preserve names checked(ccall((:H5Literate2, HDF5_jll.libhdf5), Cint,
                                     (hid_t, Cint, Cint, Ptr{UInt64}, Ptr{Cvoid}, Ptr{Cvoid}),
                                     g, H5_INDEX_NAME, H5_ITER_INC, C_NULL, cb, pointer_from_objref(names)),
                               "H5Literate2")
    ccall((:H5Gclose, HDF5_jll.libhdf5), Cint, (hid_t,), g)
    return names
end

mutable struct ChunkTally
    n::Int
    bytes::UInt64
end

function tally_chunk(offset::Ptr{UInt64}, mask::Cuint, addr::UInt64, size::UInt64, op::Ptr{Cvoid})::Cint
    t = unsafe_pointer_to_objref(op)::ChunkTally
    t.n += 1
    t.bytes += size
    return Cint(0)
end

function filters_and_chunks(d::hid_t)
    dcpl = checked(ccall((:H5Dget_create_plist, HDF5_jll.libhdf5), hid_t, (hid_t,), d), "H5Dget_create_plist")
    n = ccall((:H5Pget_nfilters, HDF5_jll.libhdf5), Cint, (hid_t,), dcpl)
    ids = Int[]
    for i in 0:(n - 1)
        flags = Ref{Cuint}(0)
        nelmts = Ref{Csize_t}(0)
        config = Ref{Cuint}(0)
        id = ccall((:H5Pget_filter2, HDF5_jll.libhdf5), Cint,
                   (hid_t, Cuint, Ptr{Cuint}, Ptr{Csize_t}, Ptr{Cuint}, Csize_t, Ptr{UInt8}, Ptr{Cuint}),
                   dcpl, i, flags, nelmts, C_NULL, 0, C_NULL, config)
        push!(ids, Int(id))
    end
    ccall((:H5Pclose, HDF5_jll.libhdf5), Cint, (hid_t,), dcpl)
    tally = ChunkTally(0, 0)
    cb = @cfunction(tally_chunk, Cint, (Ptr{UInt64}, Cuint, UInt64, UInt64, Ptr{Cvoid}))
    GC.@preserve tally checked(ccall((:H5Dchunk_iter, HDF5_jll.libhdf5), Cint, (hid_t, hid_t, Ptr{Cvoid}, Ptr{Cvoid}),
                                     d, H5P_DEFAULT, cb, pointer_from_objref(tally)), "H5Dchunk_iter")
    return ids, tally
end

function complex_window(d::hid_t, row::Int, col::Int)
    # `cglobal` with a library chosen at run time fails in a trimmed binary (a MethodError the verifier
    # does not report), so the global is found through the handle the JLL opened.
    native_float = unsafe_load(Ptr{hid_t}(Libdl.dlsym(HDF5_jll.libhdf5_handle, :H5T_NATIVE_FLOAT_g)))
    ct = checked(ccall((:H5Tcreate, HDF5_jll.libhdf5), hid_t, (Cint, Csize_t), H5T_COMPOUND, 8), "H5Tcreate")
    ccall((:H5Tinsert, HDF5_jll.libhdf5), Cint, (hid_t, Cstring, Csize_t, hid_t), ct, "r", 0, native_float)
    ccall((:H5Tinsert, HDF5_jll.libhdf5), Cint, (hid_t, Cstring, Csize_t, hid_t), ct, "i", 4, native_float)
    space = checked(ccall((:H5Dget_space, HDF5_jll.libhdf5), hid_t, (hid_t,), d), "H5Dget_space")
    start = UInt64[row, col]
    count = UInt64[4, 4]
    checked(ccall((:H5Sselect_hyperslab, HDF5_jll.libhdf5), Cint, (hid_t, Cint, Ptr{UInt64}, Ptr{UInt64}, Ptr{UInt64}, Ptr{UInt64}),
                  space, H5S_SELECT_SET, start, C_NULL, count, C_NULL), "H5Sselect_hyperslab")
    mem = checked(ccall((:H5Screate_simple, HDF5_jll.libhdf5), hid_t, (Cint, Ptr{UInt64}, Ptr{UInt64}), 2, count, C_NULL),
                  "H5Screate_simple")
    buf = Matrix{ComplexF32}(undef, 4, 4)   # column-major: buf[j, i] is HDF5 element (row + i, col + j)
    checked(ccall((:H5Dread, HDF5_jll.libhdf5), Cint, (hid_t, hid_t, hid_t, hid_t, hid_t, Ptr{ComplexF32}),
                  d, ct, mem, space, H5P_DEFAULT, buf), "H5Dread window")
    ccall((:H5Sclose, HDF5_jll.libhdf5), Cint, (hid_t,), mem)
    ccall((:H5Sclose, HDF5_jll.libhdf5), Cint, (hid_t,), space)
    ccall((:H5Tclose, HDF5_jll.libhdf5), Cint, (hid_t,), ct)
    return buf
end

function (@main)(args::Vector{String})::Cint
    try
        checked(ccall((:H5open, HDF5_jll.libhdf5), Cint, ()), "H5open")
        ccall((:H5Eset_auto2, HDF5_jll.libhdf5), Cint, (hid_t, Ptr{Cvoid}, Ptr{Cvoid}), hid_t(0), C_NULL, C_NULL)
        path, dset = args[1], args[2]
        row, col = parse(Int, args[3]), parse(Int, args[4])
        f = checked(ccall((:H5Fopen, HDF5_jll.libhdf5), hid_t, (Cstring, Cuint, hid_t), path, H5F_ACC_RDONLY, H5P_DEFAULT), "H5Fopen")
        say(string("missionId=", read_fixed_string(f, "/science/LSAR/identification/missionId")))
        names = links(f, "/science/LSAR/identification")
        say(string("identification links=", length(names), " first=", first(names), " last=", last(names)))
        d = checked(ccall((:H5Dopen2, HDF5_jll.libhdf5), hid_t, (hid_t, Cstring, hid_t), f, dset, H5P_DEFAULT), "H5Dopen2")
        ids, tally = filters_and_chunks(d)
        say(string("filters=", join(string.(ids), ","), " chunks=", tally.n, " stored bytes=", tally.bytes))
        w = complex_window(d, row, col)
        say(string("window sum=", real(sum(w)), " ", imag(sum(w)), "im first=", real(w[1, 1]), " ", imag(w[1, 1]), "im"))
        ccall((:H5Dclose, HDF5_jll.libhdf5), Cint, (hid_t,), d)
        ccall((:H5Fclose, HDF5_jll.libhdf5), Cint, (hid_t,), f)
        return Cint(0)
    catch err
        say(string("error: ", message(err)))
        return Cint(1)
    end
end
