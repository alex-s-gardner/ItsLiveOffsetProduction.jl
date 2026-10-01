# Running the reference container that produced the golden data.
#
#   julia --project=tools/golden tools/golden/reference.jl S2B_MSIL1C_20200612
#   julia --project=tools/golden tools/golden/reference.jl S2B_MSIL1C_20200612 --run 2
#
# `ghcr.io/asfhyp3/hyp3-autorift:0.28.4` is the exact version recorded in every golden product's
# `source` attribute, and an arm64 manifest exists, so the chain that made the golden data runs here
# unmodified — filter, reproject, geogrid, correlate, package, crop.
#
# Two reasons to run it rather than only read its output:
#
#   **Tolerances need a measured basis.** Three fields cannot match between two runs of the reference
#   itself (see `compare.jl`), and asserting that claim is different from believing it. Running the
#   same granule twice and diffing the two products measures the reference's own reproducibility
#   floor, which is the only defensible tolerance for every other field.
#
#   **The intermediate is the diagnostic.** `testautoRIFT.py` writes `autoRIFT_intermediate.nc`
#   holding `Dx`, `Dy`, `InterpMask`, `ChipSizeX`, `SearchLimitX/Y` and `noDataMask` — the correlator
#   output before any conversion to velocity. Keeping it, together with the filtered scenes and the
#   geogrid rasters, is what makes a later product disagreement attributable to the correlator or to
#   the packaging rather than to one of them by elimination.
#
# The working directory is deliberately not cleaned. Everything the container computed stays on disk
# under `runs/<product>/<n>/`, because re-deriving it costs minutes of compute and, for Landsat,
# money.
#
# A run needs credentials for its inputs, which differ by platform: Sentinel-2 needs none, Sentinel-1
# and NISAR need Earthdata (mounted `~/.netrc`), and Landsat needs an AWS identity that can read
# requester-pays `s3://usgs-landsat`. `fetch.jl --check` reports which are reachable.

include("manifest.jl")

const IMAGE = "ghcr.io/asfhyp3/hyp3-autorift:0.28.4"
const PLATFORM = "linux/arm64"

runs_dir(c::GoldenCase) = joinpath(CACHE, "runs", c.product)
run_dir(c::GoldenCase, n::Integer) = joinpath(runs_dir(c), string(n))

"""
    image_present() -> Bool

Whether the reference image is pulled. Checked before a run so a missing image is reported as such
rather than as a container failure.
"""
function image_present()
    try
        # stderr is discarded: "No such image" is the expected answer here, not an error to report.
        out = read(pipeline(`docker image inspect $IMAGE --format "{{.Id}}"`; stderr = devnull), String)
        return !isempty(strip(out))
    catch
        return false
    end
end

"""
    run_reference(c::GoldenCase; n = 1, threads = 8, force = false) -> String

Run the reference container on `c` and return the directory holding its output.

`n` distinguishes repeated runs of the same case, which is how the reference's own reproducibility is
measured — two runs of one granule differ only in the fields that cannot be deterministic.

`threads` sets `OMP_NUM_THREADS` for the correlator's OpenMP loop. It does not affect results: each
grid point writes a distinct output element with no reduction (`REFERENCE.md`), so this is a runtime
choice only.

`prune` deletes the regenerable warped scene copies once the run succeeds — see [`prune_run`](@ref). On
by default because nothing here reads them; `filtered/` is kept either way, since rungs 5.3 and 5.4
compare against it.
"""
function run_reference(c::GoldenCase; n::Integer = 1, threads::Integer = 8, force = false,
                       prune::Bool = true)
    image_present() || error("$IMAGE is not pulled; run `docker pull --platform $PLATFORM $IMAGE`")

    dir = run_dir(c, n)
    if isdir(dir) && !isempty(readdir(dir)) && !force
        @info "run already present; pass force = true to redo it" dir
        return dir
    end
    mkpath(dir)

    netrc = joinpath(homedir(), ".netrc")
    awsdir = joinpath(homedir(), ".aws")

    mounts = ["-v", "$dir:/home/ubuntu/work"]
    isfile(netrc) && append!(mounts, ["-v", "$netrc:/home/ubuntu/.netrc:ro"])
    isdir(awsdir) && append!(mounts, ["-v", "$awsdir:/home/ubuntu/.aws:ro"])

    # Landsat inputs are requester-pays on `s3://usgs-landsat`, and `process.py` reaches them through
    # `boto3.client('s3')` and GDAL's `/vsis3/` — both of which read the default profile unless told
    # otherwise. `AWS_PROFILE` names which set of credentials in the mounted `~/.aws` to use, so a
    # machine whose keys live under a named profile needs no edit to that file.
    profile = get(ENV, "AWS_PROFILE", "itslive")

    args = ["--reference", c.reference..., "--secondary", c.secondary...]
    c.frame_id === nothing || append!(args, ["--frame-id", c.frame_id])

    cmd = `docker run --rm --platform $PLATFORM
           $mounts -w /home/ubuntu/work
           -e OMP_NUM_THREADS=$threads -e AWS_PROFILE=$profile
           $IMAGE ++process hyp3_autorift $args`

    log = joinpath(dir, "container.log")
    @info "running reference" product=c.product run=n log
    open(log, "w") do io
        try
            run(pipeline(cmd; stdout = io, stderr = io))
        catch e
            # The log is the diagnosis, so point at it rather than reproducing a wall of output.
            error("container run failed; see $log\n$(last_lines(log, 25))")
        end
    end
    # Only after the run succeeded: a failed run's scene copies are exactly what a diagnosis needs.
    if prune
        freed = prune_run(dir)
        freed > 0 && @info "pruned regenerable scene copies" freed_GiB=round(freed / 2^30; digits=2)
    end
    return dir
end

function last_lines(path, n)
    isfile(path) || return ""
    ls = readlines(path)
    return join(ls[max(1, end - n + 1):end], "\n")
end

# The warped scene copies the driver writes beside its outputs when a pair straddles two UTM zones.
# Each is an uncompressed Float32 raster of a whole scene, ~2.5 GiB for an L7 pair against ~380 MiB of
# everything else in the run.
#
# Nothing here reads them: `aligned_scenes` builds its own under `<cache>/reprojected/`, and what the
# rungs compare the filter against is `filtered/`, which stays.
const REGENERABLE = ("reprojected",)

# **`filtered/` is not in that list, though it is the larger of the two.** Rungs 5.3 and 5.4 compare the
# Julia filter against the reference's own output of it, and that output is only here: the capture holds
# the bytes the *correlator* was handed, one filter later and quantized. Deleting it costs a container
# run to restore, so the ~4.9 GiB an L7 pair spends on it stays.

# The same, for a driver that writes its scene copies as loose files beside the outputs rather than
# into a directory. `crop_gslcs` (`nisar_isce3.py:508`) writes one uncompressed CFloat32 raster of a
# whole GSLC swath per scene: on the Jakobshavn pair they are **48.3 GiB each**, against 44 GiB for
# everything else in the run put together.
const REGENERABLE_FILES = ("reference_cropped.tif", "secondary_cropped.tif")

"""
    prune_run(dir; dry_run = false) -> Int

Delete the regenerable scene copies under `dir` and return the bytes freed.

Keeps everything a comparison or a diagnosis needs: the product, `autoRIFT_intermediate.nc`, the
geogrid rasters, `offset.tif`/`velocity.tif`, the browse images, the log, and `filtered/` — which rungs
5.3 and 5.4 compare the Julia filter against and no other file holds. What goes is the warped copies and,
on a NISAR L2 pair, the two cropped GSLC rasters, which are 45 GiB each.

The byte images the correlator was handed are in the capture, so pruning cannot cost a correlator
comparison either. What it costs is re-running the driver to inspect a warped scene.
"""
function prune_run(dir::AbstractString; dry_run::Bool = false)
    freed = 0
    for name in REGENERABLE
        path = joinpath(dir, name)
        isdir(path) || continue
        for (root, _, files) in walkdir(path), f in files
            freed += filesize(joinpath(root, f))
        end
        dry_run || rm(path; recursive = true)
    end
    for name in REGENERABLE_FILES
        path = joinpath(dir, name)
        isfile(path) || continue
        freed += filesize(path)
        dry_run || rm(path)
    end
    return freed
end

prune_run(c::GoldenCase, n::Integer; kw...) = prune_run(run_dir(c, n); kw...)

"""
    run_product(dir) -> String

The product `.nc` a container run wrote into `dir`.

The intermediate and the geogrid rasters are also `.nc`/`.tif` files in the same directory, so this
selects on the product naming scheme rather than on the extension alone.
"""
function run_product(dir::AbstractString)
    isdir(dir) || error("no such run directory: $dir")
    hits = filter(readdir(dir)) do f
        endswith(f, ".nc") && occursin("_G0120V02_P", f)
    end
    isempty(hits) && error("no product in $dir; the run may have failed — see container.log")
    length(hits) == 1 || error("$(length(hits)) products in $dir: $(join(hits, ", "))")
    return joinpath(dir, only(hits))
end

"""
    run_intermediate(dir) -> Union{String,Nothing}

The `autoRIFT_intermediate.nc` a run left behind, or `nothing`.

This is the correlator's output before conversion to velocity: `Dx`, `Dy`, `InterpMask`,
`ChipSizeX`, `SearchLimitX/Y`, `noDataMask`. Its presence in the working directory makes
`testautoRIFT.py` skip correlation entirely on a subsequent run, which is worth knowing before
reusing a directory.
"""
function run_intermediate(dir::AbstractString)
    p = joinpath(dir, "autoRIFT_intermediate.nc")
    return isfile(p) ? p : nothing
end

function main(args)
    isempty(args) && error("usage: reference.jl <product-name-fragment> [--run N] [--threads N] [--force]")

    c = only(cases(args[1]))
    n = 1
    i = findfirst(==("--run"), args); i === nothing || (n = parse(Int, args[i + 1]))
    t = 8
    j = findfirst(==("--threads"), args); j === nothing || (t = parse(Int, args[j + 1]))

    dir = run_reference(c; n, threads = t, force = "--force" in args)
    println("run directory: ", dir)
    println("product:       ", try run_product(dir) catch e; "(none — $(e.msg))" end)
    ip = run_intermediate(dir)
    println("intermediate:  ", ip === nothing ? "(none)" : ip)
    return nothing
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main(ARGS)
