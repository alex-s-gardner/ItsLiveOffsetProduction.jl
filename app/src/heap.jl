# The collector's heap-size target.

# The heap size the collector aims to stay under, GiB. Without one the collector lets the heap grow well
# past what the run holds live. Measured on NISAR L1, the largest golden case, read from ASF: the footprint
# peaks at 14.6 GiB at 12 and 16.4 GiB at 14, against a 16 GiB budget; from local granules, 16.4 GiB with
# no target. 12 costs about 3% in runtime (858 s against 827-849 s from local granules), because the
# collector runs more often. `0` leaves the collector's own default.
const HEAP_HINT_GIB = 12.0

"""
    set_heap_hint!(gib)

Set the collector's heap-size target to `gib` GiB, as `julia --heap-size-hint` does at startup — unless the
process was started with that flag, which takes precedence. `0` sets nothing.
"""
function set_heap_hint!(gib::Real)
    gib > 0 || return nothing
    Base.JLOptions().heap_size_hint == 0 || return nothing
    ccall(:jl_gc_set_max_memory, Cvoid, (UInt64,), round(UInt64, gib * 2^30))
    return nothing
end
