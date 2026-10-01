# Capture the reference's `netCDF_packaging` call for one golden case, and read it into Julia.
#
#   julia --project=tools/golden tools/golden/capture_packaging.jl S2B_MSIL1C_20200612
#
# Mirrors `intermediate.jl`, one level up: that captures `runAutorift`'s arrays, this captures
# `netCDF_packaging`'s — velocity, the geogrid conversion coefficients, the reference velocity, the
# stable-surface mask, the already-determined stable-shift scalars, and `IMG_INFO_DICT`. This is the
# validation input for `write_product`, which replaces this call.

include("manifest.jl")
include("reference.jl")

using JSON3

include(joinpath(@__DIR__, "..", "ab", "xchg.jl"))

const AB_DIR = normpath(joinpath(@__DIR__, "..", "ab"))

packaging_dir(c::GoldenCase, n::Integer) = joinpath(run_dir(c, n), "capture_packaging")

"""
    capture_packaging_call(c::GoldenCase; n = 101, threads = 8, force = false) -> String

Run the reference under `capture_packaging.py` and return the directory holding the dumped
`netCDF_packaging` arguments. The run also produces the ordinary golden-equivalent product (cropped,
in `run_dir(c, n)`), so the same run gives both the packaging inputs and the file to diff the Julia
output against.
"""
function capture_packaging_call(c::GoldenCase; n::Integer = 101, threads::Integer = 8, force = false)
    image_present() || error("$IMAGE is not pulled; run `docker pull --platform $PLATFORM $IMAGE`")

    dir = run_dir(c, n)
    cap = packaging_dir(c, n)
    if isdir(cap) && !isempty(readdir(cap)) && !force
        @info "capture already present; pass force = true to redo it" cap
        return cap
    end
    mkpath(dir)

    stale = joinpath(dir, "autoRIFT_intermediate.nc")
    isfile(stale) && (@info "removing stale intermediate so the correlator runs" stale; rm(stale))

    netrc = joinpath(homedir(), ".netrc")
    awsdir = joinpath(homedir(), ".aws")

    mounts = ["-v", "$dir:/home/ubuntu/work",
              "-v", "$(joinpath(@__DIR__, "capture_packaging.py")):/opt/capture/capture_packaging.py:ro",
              "-v", "$(joinpath(AB_DIR, "xchg.py")):/opt/capture/xchg.py:ro"]
    isfile(netrc) && append!(mounts, ["-v", "$netrc:/home/ubuntu/.netrc:ro"])
    isdir(awsdir) && append!(mounts, ["-v", "$awsdir:/home/ubuntu/.aws:ro"])

    profile = get(ENV, "AWS_PROFILE", "itslive")

    args = ["--reference", c.reference..., "--secondary", c.secondary...]
    c.frame_id === nothing || append!(args, ["--frame-id", c.frame_id])

    inner = "cd /hyp3-autorift && pixi run --manifest-path /hyp3-autorift/pyproject.toml " *
            "bash -c 'cd /home/ubuntu/work && exec python /opt/capture/capture_packaging.py " *
            join(args, " ") * "'"
    cmd = `docker run --rm --platform $PLATFORM
           $mounts -w /home/ubuntu/work
           -e OMP_NUM_THREADS=$threads -e CAPTURE_DIR=/home/ubuntu/work/capture_packaging
           -e AWS_PROFILE=$profile
           -e PYTHONPATH=/opt/capture
           --entrypoint /bin/bash
           $IMAGE -lc $inner`

    log = joinpath(dir, "capture_packaging.log")
    @info "capturing netCDF_packaging arguments" product=c.product run=n log
    started = time()
    open(log, "w") do io
        try
            run(pipeline(cmd; stdout = io, stderr = io))
        catch
            error("capture run failed; see $log\n$(last_lines(log, 40))")
        end
    end

    manifest = joinpath(cap, "call1.json")
    if !isfile(manifest) || mtime(manifest) < started
        error("netCDF_packaging did not run: `capture_packaging/call1.json` is " *
              (isfile(manifest) ? "older than this run" : "absent") * ". See $log")
    end
    return cap
end

function read_packaging_capture(dir::AbstractString; call::Integer = 1)
    mpath = joinpath(dir, "call$call.json")
    isfile(mpath) || error("no call$call.json in $dir")
    m = JSON3.read(read(mpath, String))
    arrays = Dict{String,Matrix}()
    for (name, info) in pairs(m.arrays)
        haskey(info, :skipped_dtype) && continue
        # The manifest's `file` field carries a cosmetic `.abx` suffix that `xchg.write` never adds
        # when given a full path (see `xchg.py::path`); the bare key is the file actually on disk —
        # `intermediate.jl`'s `read_capture` makes the same correction for the same reason.
        arrays[String(name)] = xread(joinpath(dir, String(name)))
    end
    scalars = Dict{String,Any}(String(k) => v for (k, v) in pairs(m.scalars))
    img_pair_info = Dict{String,Any}(String(k) => v for (k, v) in pairs(m.img_pair_info))
    return (; arrays, scalars, img_pair_info)
end

function main(args)
    isempty(args) && error("usage: capture_packaging.jl <product-name-fragment> [--run N] [--force]")
    c = only(cases(args[1]))
    n = 101
    i = findfirst(==("--run"), args); i === nothing || (n = parse(Int, args[i + 1]))

    cap = capture_packaging_call(c; n, force = "--force" in args)
    println("capture directory: ", cap)
    k = read_packaging_capture(cap)
    println("arrays: ", sort(collect(keys(k.arrays))))
    println("scalars: ", k.scalars)
    println("img_pair_info: ", k.img_pair_info)
    dir = run_dir(c, n)
    println("product: ", try run_product(dir) catch e; "(none — $(e.msg))" end)
    return nothing
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS)
