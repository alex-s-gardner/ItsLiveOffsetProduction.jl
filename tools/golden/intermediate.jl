# Capture the correlator's inputs and outputs from a reference run, and read them into Julia.
#
#   julia --project=tools/golden tools/golden/intermediate.jl S2B_MSIL1C_20200612
#
# This runs the container with `capture.py` in place of the ordinary entry point, so the pipeline is
# the real one but `runAutorift` dumps the arrays it is handed and the arrays it returns. See
# `capture.py` for why the arrays are taken at that boundary rather than reconstructed from the scene
# files: reconstructing them would mean reimplementing the code under test.
#
# What comes back is everything needed to run AutoRIFT.jl on the identical problem — the filtered
# pair, the snapped grid, the per-point search limits and chip bounds, the priors — plus the
# reference's own answer to diff against.

include("manifest.jl")
include("reference.jl")

using JSON3

# `xchg` is the A/B harness's array exchange: element type and both dimensions travel in the file, so
# a reader cannot be handed the wrong shape. Reused rather than reinvented, and the format is
# asserted round-trip by `tools/ab/xchg_test.sh`.
include(joinpath(@__DIR__, "..", "ab", "xchg.jl"))

const AB_DIR = normpath(joinpath(@__DIR__, "..", "ab"))

capture_dir(c::GoldenCase, n::Integer) = joinpath(run_dir(c, n), "capture")

"""
    capture_reference(c::GoldenCase; n = 100, threads = 8, force = false) -> String

Run the reference under `capture.py` and return the directory holding the dumped arrays.

`n` defaults to 100 so a capture run does not collide with the plain runs 1 and 2 that measure
reproducibility: a capture is the same pipeline, but its working directory accumulates extra files
and its stdout differs, and keeping the two kinds separate makes each reproducible on its own.
"""
function capture_reference(c::GoldenCase; n::Integer = 100, threads::Integer = 8, force = false)
    image_present() || error("$IMAGE is not pulled; run `docker pull --platform $PLATFORM $IMAGE`")

    dir = run_dir(c, n)
    cap = capture_dir(c, n)
    if isdir(cap) && !isempty(readdir(cap)) && !force
        # A capture taken before per-level recording existed is complete for the correlator
        # comparison and silently missing the level diagnostic, so say which one is on disk rather
        # than letting a later `length(k.levels) == 0` look like the reference measured no levels.
        haslevels = any(startswith("lvl"), readdir(cap))
        @info "capture already present; pass force = true to redo it" cap levels = haslevels
        haslevels || @warn "this capture has no per-level records; redo it with force = true to " *
                           "attribute a disagreement to a pyramid level" cap
        return cap
    end
    mkpath(dir)

    # `autoRIFT_intermediate.nc` left by a previous run makes the driver load that instead of
    # correlating (`testautoRIFT.py:693-706`), so `runAutorift` is never called and `capture.py`'s
    # patch has nothing to intercept. The run then *succeeds* — a fresh product, a fresh log, no
    # error — while the capture directory keeps whatever the earlier run put there. A forced
    # re-capture that silently returns the arrays it was asked to replace is the worst outcome
    # available, so the file is removed rather than detected.
    stale = joinpath(dir, "autoRIFT_intermediate.nc")
    isfile(stale) && (@info "removing stale intermediate so the correlator runs" stale;
                      rm(stale))

    netrc = joinpath(homedir(), ".netrc")
    awsdir = joinpath(homedir(), ".aws")

    mounts = ["-v", "$dir:/home/ubuntu/work",
              # `capture.py` and the `xchg` module it writes through, read-only.
              "-v", "$(joinpath(@__DIR__, "capture.py")):/opt/capture/capture.py:ro",
              "-v", "$(joinpath(AB_DIR, "xchg.py")):/opt/capture/xchg.py:ro"]
    isfile(netrc) && append!(mounts, ["-v", "$netrc:/home/ubuntu/.netrc:ro"])
    isdir(awsdir) && append!(mounts, ["-v", "$awsdir:/home/ubuntu/.aws:ro"])

    # Landsat inputs are requester-pays on `s3://usgs-landsat`, and `process.py` reaches them through
    # `boto3.client('s3')` and GDAL's `/vsis3/` — both of which read the default profile unless told
    # otherwise. `AWS_PROFILE` names which set of credentials in the mounted `~/.aws` to use, so a
    # machine whose keys live under a named profile needs no edit to that file.
    profile = get(ENV, "AWS_PROFILE", "itslive")

    # `CAPTURE_STAGES` turns on the line-level trace of `runAutorift`'s locals, which is what makes a
    # stage-by-stage comparison possible: the intermediates between the correlator calls are local to
    # the loop body and unreachable any other way. Off by default because the trace is slow.
    stageenv = String[]
    haskey(ENV, "CAPTURE_STAGES") &&
        append!(stageenv, ["-e", "CAPTURE_STAGES=1",
                           "-e", "CAPTURE_STAGE_LEVEL=" * get(ENV, "CAPTURE_STAGE_LEVEL", "1")])
    # `CAPTURE_FLOAT32` runs the same pipeline with `DataType = 1`, so the correlator sees the filtered
    # float field rather than a 256-level quantization of it and the pyramid reaches `arImgDisp_s`. Paired
    # with an ordinary capture of the same case it separates the quantization's contribution to a residual
    # from everything else; see `capture.py`. A separate `--run` is wanted, since the two captures differ
    # in `in_I1`'s element type and nothing else names which is which.
    haskey(ENV, "CAPTURE_FLOAT32") && append!(stageenv, ["-e", "CAPTURE_FLOAT32=1"])

    args = ["--reference", c.reference..., "--secondary", c.secondary...]
    c.frame_id === nothing || append!(args, ["--frame-id", c.frame_id])

    # `pixi run` is how the container's entry point reaches the environment holding autoRIFT and
    # hyp3_autorift; invoking python directly would miss it. It has to run from the pixi project
    # directory, so the working directory is restored explicitly before the pipeline starts —
    # otherwise every relative output path the driver writes (the geogrid rasters, the intermediate,
    # the product itself) lands inside the container and is discarded with it, leaving only the
    # capture directory that `CAPTURE_DIR` names absolutely.
    inner = "cd /hyp3-autorift && pixi run --manifest-path /hyp3-autorift/pyproject.toml " *
            "bash -c 'cd /home/ubuntu/work && exec python /opt/capture/capture.py " *
            join(args, " ") * "'"
    cmd = `docker run --rm --platform $PLATFORM
           $mounts -w /home/ubuntu/work
           -e OMP_NUM_THREADS=$threads -e CAPTURE_DIR=/home/ubuntu/work/capture
           -e AWS_PROFILE=$profile
           -e PYTHONPATH=/opt/capture
           $stageenv
           --entrypoint /bin/bash
           $IMAGE -lc $inner`

    log = joinpath(dir, "capture.log")
    @info "capturing correlator arrays" product=c.product run=n log
    started = time()
    open(log, "w") do io
        try
            run(pipeline(cmd; stdout = io, stderr = io))
        catch
            error("capture run failed; see $log\n$(last_lines(log, 30))")
        end
    end

    # The patch printing "patched" only says it was installed; it has to have *fired*. A run that
    # skipped the correlator writes a product and a log and leaves the capture untouched, so check the
    # manifest is newer than this run rather than merely present.
    manifest = joinpath(cap, "call1.json")
    if !isfile(manifest) || mtime(manifest) < started
        error("the correlator did not run: `capture/call1.json` is " *
              (isfile(manifest) ? "older than this run" : "absent") * ".\n" *
              "`runAutorift` is skipped when an intermediate is present, and the run then " *
              "succeeds while the capture stays stale. See $log")
    end
    isdir(cap) || error("run finished but wrote no capture directory; see $log")
    # The capture is only trustworthy if the run it came from produced the product too. A capture
    # directory alone means the pipeline wrote its outputs somewhere else, and the arrays would then
    # belong to a run whose result cannot be checked against golden.
    try
        run_product(dir)
    catch
        error("capture wrote arrays but the run produced no product in $dir — the pipeline's " *
              "outputs went elsewhere, so this capture cannot be tied to a verified result. " *
              "See $log")
    end
    return cap
end

"""
    LevelRecord

What one pyramid level did, taken from inside the reference's own loop.

`runAutorift`'s inputs and outputs describe the *merged* answer, so which level a disagreement came
from otherwise has to be inferred from the reported `ChipSizeX` — and that inference is weakest
exactly where it matters, at a point the two implementations assign to different levels.

`kind` is `"coarse"` or `"fine"` for a correlator call and `"filtDisp"` for a rejection pass, in the
order the loop made them. A correlator record carries the level's raw `dx`/`dy` before any filtering
or merge, so it separates "this level never measured the point" from "this level measured it and the
merge preferred another".
"""
struct LevelRecord
    seq::Int
    kind::String
    chip_size::Tuple{Float64,Float64}
    oversample::Float64
    grid_shape::Tuple{Int,Int}
    counts::Dict{String,Int}
    # The one filter parameter that is not an integer. `NaN` on a correlator record, which has no
    # filter — comparing a derived `FracValid` against the reference's needs the value, not a rounding
    # of it, since `8/25` and `0.32` differ in the last bits.
    frac_valid::Float64
    arrays::Dict{String,Matrix}
end

"""
    Capture

One `runAutorift` call's arrays and scalars, as Julia values.

`arrays` is keyed as the manifest writes them — `in_I1`, `out_Dx`, and so on — holding each in the
orientation the reference had it. `scalars` carries every attribute that affects the result, so
`Params` can be configured from what the reference used rather than from what the driver is believed
to set. `levels` is the per-level record, empty for a capture taken before it was collected.

`stages` holds the locals `runAutorift` builds *between* those calls, keyed as the trace names them —
`xGrid0_L0`, `SearchLimitX0_rev2_L0`, `MC2_L0`. It is empty unless the capture was taken with
`CAPTURE_STAGES=1`, and it holds one chip-size level per run. This is what makes a stage-by-stage
comparison possible: those arrays are local to the loop body, so the function boundaries are the only
other place two implementations can be compared, and comparing there compares two dozen composed steps
at once.
"""
struct Capture
    call::Int
    arrays::Dict{String,Matrix}
    scalars::Dict{String,Any}
    skipped::Dict{String,String}
    levels::Vector{LevelRecord}
    stages::Dict{String,Matrix}
end


"""
    read_capture(dir; call = 1, mmap = ()) -> Capture

Read one captured call from `dir`.

A pipeline run calls `runAutorift` once, so `call = 1` is the usual case; the argument exists because
a driver that retried would produce more, and silently reading the first of several would compare
against the wrong one.

`mmap` names the arrays to map from the file instead of reading onto the heap, and
`("in_I1", "in_I2")` is the useful value: the imagery is almost all of a capture's bytes — 11.25 GiB
on NISAR L2, against a few hundred MiB for every other array together — and a mapped pair is
file-backed and clean, so the kernel can reclaim it. `xread_mmap` gives the identical array.

**Named rather than a `Bool`, because a mapped array is read-only.** The mapping is opened from a
read-only handle, so writing one is a fault rather than an error; a caller that mutates a captured
array in place must not have that array in this list. The imagery is not mutated by anything here.
"""
function read_capture(dir::AbstractString; call::Integer = 1,
                      mmap::Tuple{Vararg{AbstractString}} = ())
    mpath = joinpath(dir, "call$call.json")
    isfile(mpath) || error("no call$call.json in $dir; captured calls: " *
                           join(filter(f -> startswith(f, "call"), readdir(dir)), ", "))
    m = JSON3.read(read(mpath, String))

    arrays = Dict{String,Matrix}()
    for (name, info) in pairs(m.arrays)
        path = joinpath(dir, String(name))
        arrays[String(name)] = String(name) in mmap ? xread_mmap(path) : xread(path)
    end
    scalars = Dict{String,Any}(String(k) => v for (k, v) in pairs(m.scalars))
    skipped = Dict{String,String}(String(k) => String(v) for (k, v) in pairs(get(m, :skipped, (;))))

    # Absent from a capture taken before per-level recording existed, which is a missing diagnostic
    # rather than a broken capture — the correlator comparison does not read it.
    levels = LevelRecord[]
    for r in get(m, :levels, ())
        la = Dict{String,Matrix}()
        for (name, _) in pairs(get(r, :arrays, (;)))
            la[String(name)] = xread(joinpath(dir, "lvl$(r.seq)_$(name)"))
        end
        counts = Dict{String,Int}()
        for key in (:measured, :in_mask, :kept, :filt_width, :iterations)
            haskey(r, key) && (counts[String(key)] = Int(r[key]))
        end
        # `frac_valid` is the one filter parameter that is not an integer, so it travels beside the
        # counts rather than in them. `NaN` for a correlator record, which has no filter.
        counts_frac = Float64(get(r, :frac_valid, NaN))
        push!(levels, LevelRecord(
            r.seq, String(r.kind),
            (Float64(get(r, :chip_size_x, NaN)), Float64(get(r, :chip_size_y, NaN))),
            Float64(get(r, :oversample, NaN)),
            (Int(r.grid_shape[1]), Int(r.grid_shape[2])),
            counts, counts_frac, la))
    end

    # The stage trace, absent unless the capture asked for it. Keys carry the level, so a directory
    # holding traces of two levels keeps them apart.
    stages = Dict{String,Matrix}()
    for (name, _) in pairs(get(m, :stages, (;)))
        stages[String(name)] = xread(joinpath(dir, "stage_" * String(name)))
    end

    return Capture(m.call, arrays, scalars, skipped, sort!(levels, by = r -> r.seq), stages)
end

read_capture(c::GoldenCase; n::Integer = 100, call::Integer = 1,
             mmap::Tuple{Vararg{AbstractString}} = ()) =
    read_capture(capture_dir(c, n); call, mmap)

# The imagery, which is what a caller maps when it maps anything. Named here rather than spelled at
# each call site so a capture whose imagery keys change is fixed in one place.
const CAPTURE_IMAGERY = ("in_I1", "in_I2")

function main(args)
    isempty(args) && error("usage: intermediate.jl <product-name-fragment> [--run N] [--force]")
    c = only(cases(args[1]))
    n = 100
    i = findfirst(==("--run"), args); i === nothing || (n = parse(Int, args[i + 1]))

    cap = capture_reference(c; n, force = "--force" in args)
    println("capture directory: ", cap)

    k = read_capture(cap)
    println("\ncall ", k.call)
    println("arrays:")
    for name in sort!(collect(keys(k.arrays)))
        a = k.arrays[name]
        println("  ", rpad(name, 18), rpad(string(eltype(a)), 10), size(a))
    end
    isempty(k.skipped) || println("skipped: ", k.skipped)
    println("\nscalars:")
    for name in sort!(collect(keys(k.scalars)))
        println("  ", rpad(name, 24), k.scalars[name])
    end
    return nothing
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS)
