#!/bin/bash
# Run inside the amd64 container `build_linux.sh` starts: copy the repository from `$REPO`, instantiate the
# app with a container-local depot, and build into `/out/build`.
set -eu
apt-get update -qq && apt-get install -y -qq gcc g++ > /dev/null
export JULIA_DEPOT_PATH=/depot JULIA_NUM_THREADS=1 JULIA_PKG_CONCURRENT_DOWNLOADS=4
# Julia's own release multiversioning for x86_64: one binary runs on every x86_64 machine and uses the vector
# width of the one it runs on.
export JULIA_CPU_TARGET="generic;sandybridge,-xsaveopt,clone_all;haswell,-rdrnd,base(1);x86-64-v4,-rdrnd,base(1)"
mkdir -p /work
cp -a "$REPO" /work/
name="$(basename "$REPO")"
rm -rf "/work/$name/.git" "/work/$name/app/build"
julia -e 'import Pkg; Pkg.activate("jcdrv"; shared = true); Pkg.add(name = "JuliaC", version = "0.3.9")'
cd "/work/$name/app"
julia --project=. -e 'import Pkg; Pkg.instantiate(); Pkg.precompile()'
TRIM_JULIA=julia ./build.sh /out/build
