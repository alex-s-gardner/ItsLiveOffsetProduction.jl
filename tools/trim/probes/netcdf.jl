# Trim probe: writing a netCDF-4 file through the NetCDF C library directly, with each feature the
# ITS_LIVE product uses — an unlimited dimension, deflate + shuffle, chunking, fill values, NC_CHAR
# and NC_STRING text attributes, an Int64 attribute, and Int16 data.
#
#   netcdf <output.nc>

import NetCDF_jll

const NC_NETCDF4 = Cint(0x1000)
const NC_CLOBBER = Cint(0)
const NC_UNLIMITED = Csize_t(0)
const NC_GLOBAL = Cint(-1)
const NC_SHORT = Cint(3)
const NC_DOUBLE = Cint(6)
const NC_INT64 = Cint(10)
const NC_CHUNKED = Cint(0)

say(s::String) = (write(Core.stdout, s, "\n"); nothing)

function check(status::Cint, what::String)
    status == 0 && return nothing
    msg = unsafe_string(ccall((:nc_strerror, NetCDF_jll.libnetcdf), Cstring, (Cint,), status))
    error(string(what, ": ", msg))
end

function def_dim(ncid::Cint, name::String, len::Csize_t)
    id = Ref{Cint}(0)
    check(ccall((:nc_def_dim, NetCDF_jll.libnetcdf), Cint, (Cint, Cstring, Csize_t, Ptr{Cint}),
                ncid, name, len, id), "nc_def_dim $name")
    return id[]
end

function def_var(ncid::Cint, name::String, xtype::Cint, dims::Vector{Cint})
    id = Ref{Cint}(0)
    check(ccall((:nc_def_var, NetCDF_jll.libnetcdf), Cint, (Cint, Cstring, Cint, Cint, Ptr{Cint}, Ptr{Cint}),
                ncid, name, xtype, length(dims), dims, id), "nc_def_var $name")
    return id[]
end

put_text!(ncid::Cint, varid::Cint, key::String, val::String) =
    check(ccall((:nc_put_att_text, NetCDF_jll.libnetcdf), Cint, (Cint, Cint, Cstring, Csize_t, Cstring),
                ncid, varid, key, sizeof(val), val), "nc_put_att_text $key")

function put_string!(ncid::Cint, varid::Cint, key::String, val::String)
    p = [Base.unsafe_convert(Cstring, val)]
    GC.@preserve val p check(ccall((:nc_put_att_string, NetCDF_jll.libnetcdf), Cint,
                                   (Cint, Cint, Cstring, Csize_t, Ptr{Cstring}), ncid, varid, key, 1, p),
                             "nc_put_att_string $key")
end

put_int64!(ncid::Cint, varid::Cint, key::String, val::Int64) =
    check(ccall((:nc_put_att_longlong, NetCDF_jll.libnetcdf), Cint, (Cint, Cint, Cstring, Cint, Csize_t, Ref{Int64}),
                ncid, varid, key, NC_INT64, 1, val), "nc_put_att_longlong $key")

put_double!(ncid::Cint, varid::Cint, key::String, val::Float64) =
    check(ccall((:nc_put_att_double, NetCDF_jll.libnetcdf), Cint, (Cint, Cint, Cstring, Cint, Csize_t, Ref{Float64}),
                ncid, varid, key, NC_DOUBLE, 1, val), "nc_put_att_double $key")

function write_probe(path::String)
    ncid = Ref{Cint}(0)
    check(ccall((:nc_create, NetCDF_jll.libnetcdf), Cint, (Cstring, Cint, Ptr{Cint}),
                path, NC_NETCDF4 | NC_CLOBBER, ncid), "nc_create")
    nc = ncid[]
    nx, ny = 300, 200
    # netCDF dimension order is slowest first: (time, y, x).
    t = def_dim(nc, "time", NC_UNLIMITED)
    y = def_dim(nc, "y", Csize_t(ny))
    x = def_dim(nc, "x", Csize_t(nx))
    vx = def_var(nc, "vx", NC_SHORT, Cint[t, y, x])
    check(ccall((:nc_def_var_chunking, NetCDF_jll.libnetcdf), Cint, (Cint, Cint, Cint, Ptr{Csize_t}),
                nc, vx, NC_CHUNKED, Csize_t[1, 128, 128]), "nc_def_var_chunking")
    check(ccall((:nc_def_var_deflate, NetCDF_jll.libnetcdf), Cint, (Cint, Cint, Cint, Cint, Cint),
                nc, vx, 1, 1, 2), "nc_def_var_deflate")
    check(ccall((:nc_def_var_fill, NetCDF_jll.libnetcdf), Cint, (Cint, Cint, Cint, Ref{Int16}),
                nc, vx, 0, Int16(-32767)), "nc_def_var_fill")
    put_text!(nc, vx, "units", "meter/year")
    put_string!(nc, vx, "description", "velocity component in x direction")
    put_int64!(nc, vx, "stable_count", Int64(123456789012))
    put_double!(nc, vx, "error", 1.25)
    put_text!(nc, NC_GLOBAL, "Conventions", "CF-1.8")
    data = Int16[(i - j) % 1000 for i in 1:nx, j in 1:ny]   # Julia (x, y) = netCDF (y, x)
    data[1:10, 1:10] .= Int16(-32767)
    check(ccall((:nc_enddef, NetCDF_jll.libnetcdf), Cint, (Cint,), nc), "nc_enddef")
    check(ccall((:nc_put_vara_short, NetCDF_jll.libnetcdf), Cint, (Cint, Cint, Ptr{Csize_t}, Ptr{Csize_t}, Ptr{Int16}),
                nc, vx, Csize_t[0, 0, 0], Csize_t[1, ny, nx], data), "nc_put_vara_short")
    check(ccall((:nc_close, NetCDF_jll.libnetcdf), Cint, (Cint,), nc), "nc_close")
    return sum(Int, data)
end

function message(err)
    m = err isa ErrorException ? err.msg : nothing
    return m isa String ? m : string(typeof(err))
end

function (@main)(args::Vector{String})::Cint
    try
        s = write_probe(args[1])
        say(string("wrote ", args[1], " data sum=", s))
        return Cint(0)
    catch err
        say(string("error: ", message(err)))
        return Cint(1)
    end
end
