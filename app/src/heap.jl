# The process's memory settings.

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

# Bytes of correlation workspace `AutoRIFT` keeps between uses (`AutoRIFT.WORKSPACE_POOL_BYTES`, 2 GiB by
# default). A NISAR L1 pair visits about a hundred search geometries, most of them 15-20 MiB, and a finished
# level never asks for its own again, so most of a 2 GiB pool is never reused. Measured on NISAR L1 P094 at 4
# threads on Linux: 0.5 GiB lowers the peak from 10.4 to 8.6 GiB, with identical products.
const WORKSPACE_POOL_BYTES = 1 << 29

# glibc's `mallopt` parameters (`malloc.h`).
const M_TRIM_THRESHOLD = Cint(-1)
const M_MMAP_THRESHOLD = Cint(-3)

# glibc otherwise raises its mmap threshold to the size of each mapped block freed, up to 32 MiB, and serves
# later blocks below it from per-thread heaps it rarely shrinks; the correlation's per-block arrays fall in
# that range. Measured on NISAR L1 P094 at 4 threads: without fixed thresholds the resident size climbs from
# 2.4 to 5.5 GiB above the live heap over the correlation and peaks at 15.8 GiB; with them it peaks at 10.4
# GiB, and with the workspace pool above as well at 10.0 GiB in a 16 GiB container, in the same wall time
# (1738 s against 1786 s) and with identical products.
const MALLOC_RELEASE_BYTES = Cint(1 << 20)

"""
    release_freed_memory!()

On Linux, have glibc map every block of `MALLOC_RELEASE_BYTES` or more by itself, so freeing it returns it
to the kernel, and trim a heap's free top above the same size. A no-op elsewhere.
"""
function release_freed_memory!()
    @static if Sys.islinux()
        ccall(:mallopt, Cint, (Cint, Cint), M_MMAP_THRESHOLD, MALLOC_RELEASE_BYTES) == 1 &&
            ccall(:mallopt, Cint, (Cint, Cint), M_TRIM_THRESHOLD, MALLOC_RELEASE_BYTES) == 1 ||
            throw(ErrorException("mallopt refused a malloc threshold"))
    end
    return nothing
end
