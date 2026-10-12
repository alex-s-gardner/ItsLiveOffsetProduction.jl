#!/usr/bin/env bash
# Build the linux/x86_64 `autorift_itslive_production` in an amd64 Julia container, from any host Docker
# can run one on (natively, or under emulation on an arm64 Mac).
#
#   app/build_linux.sh <outdir>
#
# The repository is copied into the container, and every other package of `app/Project.toml`'s `[sources]`
# — which must sit beside the repository, in its parent directory — is read from a read-only mount of that
# parent at the same absolute path, so the manifest's paths resolve unchanged. The depot is a Docker volume, `autorift-linux-depot`, kept between
# builds. The bundle is `<outdir>/bin/autorift_itslive_production`.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
repo="$(dirname "$here")"
parent="$(dirname "$repo")"
out="$(mkdir -p "$1" && cd "$1" && pwd)"
docker volume create autorift-linux-depot > /dev/null
docker run --rm --platform linux/amd64 -v "$parent:$parent:ro" -v autorift-linux-depot:/depot -v "$out:/out" \
    -e REPO="$repo" julia:1.13 bash "$repo/app/linux_container_build.sh"
echo "built $out/build/bin/autorift_itslive_production"
