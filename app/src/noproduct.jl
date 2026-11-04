# A pair that yields no velocity, deterministically.

"""
    NoProduct(msg)

Why a pair yields no velocity, where running it again would end the same way: nothing resolved, nothing to
search, no data, or an ASTER pair the geometric correction cannot be trusted on. [`run_job`](@ref) writes such a
pair's product empty, `roi_valid_percentage` 0 and so named `_P000.nc`, which ITS_LIVE production records as
processed; any other exception is an error, and writes nothing.
"""
struct NoProduct <: Exception
    msg::String
end
Base.showerror(io::IO, e::NoProduct) = print(io, e.msg)
