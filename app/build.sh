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
here="$(cd "$(dirname "$0")" && pwd)"
out="${1:-$here/build}"
rm -rf "$out"
mkdir -p "$(dirname "$out")"
"${jl[@]}" --project="$here" -e 'import Pkg; Pkg.instantiate()'
"${jl[@]}" --project=@jcdrv -m JuliaC --experimental --trim=safe --output-exe autorift_itslive_production \
    --bundle "$out" "$here"
# The FFTW wisdom each run starts from (`ensure_environment` in src/cli.jl), one file per CPU model.
mkdir -p "$out/share/autorift_itslive_production/fftw-wisdom"
cp "$here"/wisdom/*.wisdom "$out/share/autorift_itslive_production/fftw-wisdom/"
# What the bundle is and how to run it, at the bundle's root.
cp "$here/RELEASE.md" "$out/README.md"
echo "built $out/bin/autorift_itslive_production"
