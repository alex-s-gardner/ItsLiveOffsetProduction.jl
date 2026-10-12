#!/usr/bin/env bash
# Build a Julia program with `juliac --trim=safe` and run the result with no Julia depot in sight.
#
#   tools/trim/probe.sh [--project ENV] [--threads N] <program.jl | package-dir> [program args...]
#
# Builds into $TRIM_OUT (default ./trim_out)/<name>/, writes the build log beside it, prints the
# verifier's error count grouped by the package each unresolved call sits in, and — when the build
# succeeds — runs `bin/<name>` with the user depot (`~/.julia`) hidden, so that a binary which still
# reaches into the build machine's depot (artifacts, wisdom, registries) fails here rather than in
# production. Exits with the build's status, or else the run's.
#
# JuliaC is taken from the `@jcdrv` shared environment (`Pkg.add("JuliaC")` there once). TRIM_JULIA names
# the Julia command (default `julia +release`, the juliaup channel); set it where juliaup is absent.
set -uo pipefail
read -r -a jl <<< "${TRIM_JULIA:-julia +release}"

project=""
threads=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --project) project="$2"; shift 2 ;;
        --threads) threads="$2"; shift 2 ;;
        *) break ;;
    esac
done
[[ $# -ge 1 ]] || { echo "usage: $0 [--project ENV] [--threads N] <program.jl|package-dir> [args...]" >&2; exit 2; }

src="$1"; shift
here="$(cd "$(dirname "$0")" && pwd)"
name="$(basename "${src%.jl}")"
out="${TRIM_OUT:-$PWD/trim_out}/$name"
rm -rf "$out" && mkdir -p "$out"
log="$out/build.log"

args=(--experimental --output-exe "$name" --bundle "$out/bundle" --trim=safe)
[[ -n "$project" ]] && args+=(--project "$project")
[[ -n "$threads" ]] && args+=(--jl-option "threads=$threads")

"${jl[@]}" --project=@jcdrv -m JuliaC "${args[@]}" "$src" > "$log" 2>&1
status=$?
"${jl[@]}" --startup-file=no "$here/verifier_report.jl" "$log"
if [[ $status -ne 0 ]]; then
    echo "$name: BUILD FAILED (log: $log)"
    exit $status
fi

empty="$(mktemp -d)"
# The trailing `:` expands to Julia's default depots *except* the user depot, which keeps the
# bundle's own `share/julia` (where `--bundle` puts the artifacts) and hides `~/.julia`.
JULIA_DEPOT_PATH="$empty:" "$out/bundle/bin/$name" "$@" < /dev/null
status=$?
rm -rf "$empty"
echo "$name: run exit $status"
exit $status
