#!/usr/bin/env bash
# Build the `autorift_itslive_production` executable with `juliac --trim=safe`, bundled with the libraries
# and data it loads at run time (GDAL and PROJ data, the CA bundle, the JLL artifacts).
#
#   app/build.sh [outdir]        default: app/build
#
# The executable is `<outdir>/bin/autorift_itslive_production`; the bundle directory is relocatable as a
# whole. `app/wisdom/` holds the FFTW wisdom the bundle ships, measured by running the golden cases on each
# CPU model the executable is deployed on. JuliaC comes from the `@jcdrv` shared environment (`julia +release --project=@jcdrv -e
# 'import Pkg; Pkg.add("JuliaC")'` once). TRIM_JULIA names the Julia command (default `julia +release`).
set -euo pipefail
read -r -a jl <<< "${TRIM_JULIA:-julia +release}"
# The executable runs only on CPUs that one of its compiled targets supports, and juliac compiles for the
# build host's CPU alone unless JULIA_CPU_TARGET says otherwise. The default here is the multiversioning of
# Julia's own release binaries: one executable runs on every CPU of the platform and takes the code path of
# the one it runs on.
cpu_target="${JULIA_CPU_TARGET:-}"
if [ -z "$cpu_target" ]; then
    case "$(uname -sm)" in
        "Linux x86_64")
            cpu_target="generic;sandybridge,-xsaveopt,clone_all;haswell,-rdrnd,base(1);x86-64-v4,-rdrnd,base(1)" ;;
        "Linux aarch64")
            cpu_target="generic;cortex-a57;thunderx2t99;carmel,clone_all;apple-m1,base(3);neoverse-512tvb,-rand,-fpac,base(3)" ;;
        "Darwin arm64")
            cpu_target="generic;apple-m1,clone_all" ;;
        *) echo "no default CPU target for $(uname -sm); set JULIA_CPU_TARGET" >&2; exit 1 ;;
    esac
fi
echo "CPU target: $cpu_target"
here="$(cd "$(dirname "$0")" && pwd)"
out="${1:-$here/build}"
rm -rf "$out"
mkdir -p "$(dirname "$out")"
"${jl[@]}" --project="$here" -e 'import Pkg; Pkg.instantiate()'
JULIA_CPU_TARGET="$cpu_target" "${jl[@]}" --project=@jcdrv -m JuliaC --experimental --trim=safe --output-exe autorift_itslive_production \
    --bundle "$out" "$here"
# Each JLL artifact is bundled whole. Headers, CMake files, static libraries and command-line tools are used
# only to build against or to run by hand; the executable loads shared libraries (`lib/`) and data (`share/`).
find "$out/share/julia/artifacts" -mindepth 2 -maxdepth 2 -type d \( -name include -o -name cmake -o -name bin \) \
    -exec rm -rf {} +
find "$out/share/julia/artifacts" -type f -name '*.a' -delete
"$here/link_libraries.sh" "$out"
# The FFTW wisdom each run starts from (`ensure_environment` in src/cli.jl), one file per CPU model.
mkdir -p "$out/share/autorift_itslive_production/fftw-wisdom"
cp "$here"/wisdom/*.wisdom "$out/share/autorift_itslive_production/fftw-wisdom/"
# What the bundle is and how to run it, at the bundle's root.
cp "$here/RELEASE.md" "$out/README.md"
echo "built $out/bin/autorift_itslive_production"
