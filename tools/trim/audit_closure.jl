# Fails when an environment's dependency closure holds a package whose presence alone breaks a
# trimmed build.
#
#   julia tools/trim/audit_closure.jl <project-dir>
#
# Every package in the closure has its `__init__` compiled and run in the binary whether or not
# anything calls it. The list below is measured: a program that only `import`s the package, built with
# `--trim=safe` on Julia 1.13 (`README.md`, "Loading versus calling"). Packages that load cleanly but
# whose *calls* do not trim — ArchGDAL, NCDatasets, AWS, Downloads, HTTP, JSON, FFTW plans — are not
# listed; the verifier reports those uses where the binary makes them. Weak dependencies whose
# extensions are not triggered are not in the closure and are fine.

import Pkg

const BANNED = Dict(
    "HDF5" => "`@require` blocks in __init__; use HDF5_jll directly",
    "Requires" => "`@require` runs dynamic code in __init__",
    "TiffImages" => "its load reaches Base's task and IO error printers",
    "Rasters" => "its load reaches Base's task and IO error printers",
)

function audit(project::AbstractString)
    Pkg.activate(project; io = devnull)
    deps = Pkg.dependencies()
    bad = sort!([(info.name, BANNED[info.name]) for info in values(deps) if haskey(BANNED, info.name)])
    if isempty(bad)
        println("closure of $project: clean ($(length(deps)) packages)")
        return true
    end
    println("closure of $project holds $(length(bad)) package(s) a trimmed binary cannot load:")
    for (name, why) in bad
        println("  ", rpad(name, 12), why)
    end
    return false
end

audit(ARGS[1]) || exit(1)
