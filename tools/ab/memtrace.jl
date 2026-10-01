# Sampling a process's own resident memory while it correlates, and pinning the peak to a stack.
#
# `benchmark/memory.jl` answers "what was the peak" with one `Sys.maxrss()` per subprocess, which is
# the right instrument for regression detection and the wrong one for attribution: a high-water mark
# says nothing about *when* the water was high, so it cannot distinguish a peak set by the imagery
# read from one set by the correlation. This samples the trace instead, and reads the profiler's
# stacks at the sample where the trace peaked.
#
# Included by a child process, which correlates; the parent (`mem_blocks.jl`) renders. Nothing here
# loads AutoRIFT, so the parent can include it too and share the record layout.

using Profile
using Serialization: serialize

# ---------------------------------------------------------------------------
# Reading this process's resident size, without a subprocess
# ---------------------------------------------------------------------------
#
# `Sys.maxrss()` is a high-water mark and `Base.gc_live_bytes()` is the Julia heap, so neither can
# produce a trace of what an instance's memory limit sees. `proc_pid_rusage` reports both current
# figures the OS tracks, and reports them for the calling process, so a sample costs one syscall
# rather than a `ps` fork.
#
# The two differ and both are needed. `resident_size` counts every resident page, including
# clean file-backed ones — so a memory-mapped input inflates it with page cache the kernel would
# evict under pressure rather than OOM-kill for. `phys_footprint` is the figure macOS uses for its
# own memory limits and excludes that clean file-backed residency, so it is what a *lazily read*
# configuration should be judged on. Reporting only `resident_size` would charge a lazy run for the
# page cache its own reads populated.

# `rusage_info_v4`: a 16-byte uuid, then `UInt64` fields. `resident_size` is the 7th and
# `phys_footprint` the 8th, so at 8-byte stride from a `UInt64` buffer they are indices 9 and 10.
const RUSAGE_INFO_V4 = Cint(4)
const RUSAGE_RESIDENT = 9
const RUSAGE_FOOTPRINT = 10
# `ri_user_time` and `ri_system_time`, the first two `UInt64` fields after the uuid.
#
# **In mach absolute ticks, not nanoseconds, despite the field names.** Measured: burning 1.00 s of CPU
# on one thread moves the user field by 23,866,353, which is the same 24 MHz `cntvct_el0` the profiler
# stamps with — so `tick_rate` converts these too. Read as nanoseconds they under-report by ~42x, which
# looks like a process that barely ran rather than like a broken unit.
const RUSAGE_USER_TICKS = 3
const RUSAGE_SYSTEM_TICKS = 4

"""
    rusage!(buf) -> (resident, footprint)

This process's current resident size and physical footprint in bytes.

`buf` is scratch of at least 64 `UInt64`s, reused so a sample allocates nothing — a sampler running
at millisecond cadence must not itself be a source of the pressure it measures.
"""
function rusage!(buf::Vector{UInt64})
    rc = ccall(:proc_pid_rusage, Cint, (Cint, Cint, Ptr{UInt64}),
               getpid(), RUSAGE_INFO_V4, buf)
    rc == 0 || error("proc_pid_rusage failed with $rc")
    return (buf[RUSAGE_RESIDENT], buf[RUSAGE_FOOTPRINT])
end

"""
    cpu_seconds!(buf, hz) -> Float64

This process's user plus system CPU time so far, in seconds, summed over every thread.

`hz` is the tick rate from [`tick_rate`](@ref); the kernel reports these fields in mach absolute ticks
rather than in the nanoseconds their names suggest.

**Reported beside any wall clock that a conclusion rests on**, because the two together separate causes
that wall alone cannot. A run that is slower at the same CPU seconds was *waiting* — on a lock, on an
allocation, on a rebuild; a run that is slower with CPU seconds up in proportion was running at a lower
clock, which is what a power- or thermally-limited machine does. Wall alone leaves those
indistinguishable, and they call for opposite responses: fix the code, or fix the measurement.

Same `buf` and same syscall as [`rusage!`](@ref), so a sampler already holding one pays nothing extra.
"""
function cpu_seconds!(buf::Vector{UInt64}, hz::Real)
    rc = ccall(:proc_pid_rusage, Cint, (Cint, Cint, Ptr{UInt64}),
               getpid(), RUSAGE_INFO_V4, buf)
    rc == 0 || error("proc_pid_rusage failed with $rc")
    return (buf[RUSAGE_USER_TICKS] + buf[RUSAGE_SYSTEM_TICKS]) / hz
end

# ---------------------------------------------------------------------------
# The clock, and why it is this one
# ---------------------------------------------------------------------------
#
# A sample and a profiler stack can only be matched if they carry the same clock. Julia's profiler
# stamps each block with `cpu_cycle_clock`, which on aarch64 is `cntvct_el0` — the architectural
# counter, 24 MHz on Apple silicon — and that is neither `time_ns()` nor `mach_absolute_time()`.
# Verified by reading the register directly and bracketing a profiled call: every stamp falls inside
# the bracket, which neither of the other two clocks does.
#
# So the trace is stamped from the same register. `TICK_HZ` is measured rather than assumed, since
# the rate is a platform property and a hard-coded 24 MHz would silently mislabel every axis on any
# machine that differs.
@static if Sys.ARCH === :aarch64
    cpu_ticks() = Base.llvmcall(
        ("""declare i64 @llvm.read_register.i64(metadata) nounwind
            !0 = !{!"cntvct_el0"}
            define i64 @entry() #0 {
              %v = call i64 @llvm.read_register.i64(metadata !0)
              ret i64 %v
            }
            attributes #0 = { alwaysinline }""", "entry"), UInt64, Tuple{})
else
    # x86-64 stamps the profile with `rdtsc`. Untested here; the alignment check below is what
    # reports it rather than a wrong figure.
    cpu_ticks() = ccall("llvm.readcyclecounter", llvmcall, UInt64, ())
end

"""
    tick_rate(; seconds = 0.2) -> Float64

Ticks per second of the counter [`cpu_ticks`](@ref) reads, measured against `time_ns`.

Measured rather than assumed: the rate is a platform property, and a wrong constant mislabels every
time axis derived from a profile stamp without producing any other symptom.
"""
function tick_rate(; seconds::Real = 0.2)
    t0, c0 = time_ns(), cpu_ticks()
    sleep(seconds)
    t1, c1 = time_ns(), cpu_ticks()
    return (c1 - c0) / ((t1 - t0) / 1e9)
end

# ---------------------------------------------------------------------------
# The sampler
# ---------------------------------------------------------------------------

"""
    MemTrace

A resident-memory trace: one sample per row, stamped on the profiler's clock.

`resident` and `footprint` are bytes (see [`rusage!`](@ref) for the difference and why both are
kept), `live` is `Base.gc_live_bytes` — the Julia heap still reachable, which separates a peak the
process genuinely needs from allocator slack the collector has not returned. `reads` is the running
count of windowed reads the input has served, so a trace can be read against how far through the
blocks the run was.
"""
struct MemTrace
    tick::Vector{UInt64}
    resident::Vector{UInt64}
    footprint::Vector{UInt64}
    live::Vector{Int64}
    reads::Vector{Int32}
    tick_hz::Float64
end

MemTrace(tick_hz::Float64) = MemTrace(UInt64[], UInt64[], UInt64[], Int64[], Int32[], tick_hz)
Base.length(t::MemTrace) = length(t.tick)

"""
    sample!(trace, buf, reads)

Append one sample to `trace`.

Ordered so the memory figures are read as close together as possible: the tick first, then the two
OS figures in one syscall, then the heap counter. `reads` is read last because it is monotonic and a
skew of one block matters far less than a skew between the memory figures.
"""
function sample!(t::MemTrace, buf::Vector{UInt64}, reads::Integer)
    tick = cpu_ticks()
    res, foot = rusage!(buf)
    push!(t.tick, tick)
    push!(t.resident, res)
    push!(t.footprint, foot)
    push!(t.live, Base.gc_live_bytes())
    push!(t.reads, reads % Int32)
    return t
end

"""
    with_trace(f; interval = 0.005, progress = nothing, reads = () -> 0) -> (result, trace, seconds)

Run `f()` while sampling this process's resident memory, and return its value with the trace.

The sampler runs on the **interactive** thread pool, which is what keeps it sampling while every
worker thread is saturated by correlation: a `:default` task would queue behind the run's own tasks
and leave second-long gaps in the trace exactly where the peak is. That requires the process to have
been started with an interactive thread (`-t N,1`), and this checks rather than silently sampling at
whatever cadence it gets.

`progress` is called with the trace and the current read count at roughly ten times `interval`, for
a caller that wants to report while the run proceeds. `reads` supplies the input's read counter.
"""
function with_trace(f; interval::Real = 0.005, progress = nothing, reads = () -> 0)
    Threads.nthreadpools() >= 2 && Threads.nthreads(:interactive) >= 1 || error(
        "the sampler needs an interactive thread so it keeps sampling while the correlation " *
        "saturates the default pool; start Julia with `-t N,1`")
    hz = tick_rate()
    trace = MemTrace(hz)
    buf = zeros(UInt64, 64)
    stop = Threads.Atomic{Bool}(false)
    sampler = Threads.@spawn :interactive begin
        n = 0
        while !stop[]
            sample!(trace, buf, reads())
            n += 1
            isnothing(progress) || n % 10 == 0 && (progress(trace, reads()); nothing)
            sleep(interval)
        end
        # One last sample after the run, so the trace ends at the settled figure rather than at
        # whatever the final `sleep` happened to straddle.
        sample!(trace, buf, reads())
    end
    t0 = time_ns()
    result = try
        f()
    finally
        stop[] = true
        wait(sampler)
    end
    return result, trace, (time_ns() - t0) / 1e9
end

# ---------------------------------------------------------------------------
# Pinning the peak to a stack
# ---------------------------------------------------------------------------

"""
    peak_stacks(data, tick_lo, tick_hi) -> Vector{Pair{String,Int}}

What the profiler recorded running between `tick_lo` and `tick_hi`, most frequent first.

`data` is `Profile.fetch(include_meta = true)`. A profile block is a stack followed by five metadata
words — `[ips..., 0, threadid, taskid, cpu_cycle_clock, thread_sleeping]`, read backwards from the
block-end double NULL as `Profile.is_block_end` lays it out — and the clock word is what makes a
stack addressable by time at all.

**The CPU profiler, not `@profile_walltime`.** The wall-time profiler samples *tasks* rather than
threads, so an idle scheduler task contributes a sample per tick: measured on a two-task saturating
load, 66 of 77 samples were `try_yieldto` and the correlation itself appeared in none of them. Its
sleep field is also the task profiler's fake state, which carries no running/sleeping distinction to
filter on. The CPU profiler's samples are per running thread and flagged, which is what makes the
filtering below possible.

Sleeping samples are dropped, as are threads parked in `__psynch_cvwait` — a thread waiting on a
condition variable is flagged running but is doing nothing, and on a wide machine those outnumber the
working stacks. Measured on this scene: 1399 of 1650 otherwise-unattributable samples were that one
frame.
"""
function peak_stacks(data::Vector{UInt64}, tick_lo::UInt64, tick_hi::UInt64;
                     nframes::Int = 6)
    counts = Dict{String,Int}()
    block_end = 0
    for i in 6:length(data)
        # `Profile.is_block_end`, inlined: two NULL words, preceded by a sleep state in 1:3.
        (data[i] == 0 && data[i - 1] == 0 && data[i - 2] in 1:3) || continue
        state = data[i - 2]
        clock = data[i - 3]
        stack_hi = block_end + 1
        block_end = i
        (tick_lo <= clock <= tick_hi) || continue
        # The field is incremented so it cannot be zero: 1 is running, 2 sleeping, 3 the task
        # profiler's fake state.
        state == 1 || continue
        # The stack is everything from the previous block's end to this block's metadata, innermost
        # frame first as the profiler records it.
        ips = @view data[stack_hi:(i - 6)]
        isempty(ips) && continue
        key = _stack_label(ips, nframes)
        isnothing(key) && continue
        counts[key] = get(counts, key, 0) + 1
    end
    return sort!(collect(counts); by = last, rev = true)
end

# One line naming what a sample was doing, or `nothing` for a sample that was doing nothing.
#
# The innermost frames of a correlation are `mul!`, an FFTW kernel, or a `ccall` at every peak —
# true, and not an answer to "which part of the pipeline is holding this memory". So the label is
# built from the innermost `nframes` frames that are *this package's*, which name a stage.
#
# Three whole-stack classifications come first, because each would otherwise be reported as an
# AutoRIFT frame that is merely the nearest one above it, or dropped as unattributable:
#
#   * A thread parked in `__psynch_cvwait` is idle. Returning `nothing` drops it.
#   * A thread in the collector is attributed to the collector. At a memory peak that is a finding
#     rather than noise — it says the peak is being fought over rather than simply held.
#   * A thread inside FFTW with no Julia frame left on the stack is the transform, which is where a
#     correlation spends much of its time and which unwinding frequently truncates.
function _stack_label(ips, nframes::Int)
    frames = String[]
    gc = false
    fftw = false
    for ip in ips
        for sf in Profile.lookup(ip % UInt)
            fn = string(sf.func)
            if sf.from_c
                fn == "__psynch_cvwait" && return nothing
                startswith(fn, "gc_") && (gc = true)
                # FFTW's generated codelets are named for the transform they perform — `n2fv_14`,
                # `t1fv_12`, `hc2cbdftv_12` — and reached through its own `apply`.
                fn == "apply" && (fftw = true)
                continue
            end
            file = string(sf.file)
            occursin("AutoRIFT", file) || continue
            push!(frames, "$(sf.func) @ $(basename(file)):$(sf.line)")
            length(frames) >= nframes && break
        end
        length(frames) >= nframes && break
    end
    isempty(frames) || return join(frames, " ← ")
    gc && return "(garbage collection)"
    fftw && return "(FFTW transform)"
    return "(no AutoRIFT frame)"
end

"""
    write_trace(path, trace, meta)

Persist a trace and whatever the caller knows about the run that produced it.

Serialized rather than printed: the trace is thousands of samples, the parent renders it, and a
figure re-rendered from a stored trace costs nothing where re-measuring costs minutes.
"""
function write_trace(path::AbstractString, trace::MemTrace, meta::NamedTuple)
    serialize(path, (; trace, meta))
    return path
end
