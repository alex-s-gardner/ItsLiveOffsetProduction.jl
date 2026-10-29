#!/bin/bash
# pack.sh BUNDLE OUTDIR: check that BUNDLE's executable starts, and pack BUNDLE as
# OUTDIR/<bundle name>.tar.zst with its .sha256 beside it. The archive holds the bundle directory itself,
# so it unpacks to <bundle name>/bin/autorift_itslive.
set -euo pipefail
bundle=$(cd "$1" && pwd); out=$(mkdir -p "$2" && cd "$2" && pwd)
name=$(basename "$bundle")
# `--help` prints the usage on stderr and exits 0.
"$bundle/bin/autorift_itslive" --help 2>&1 | head -1 | grep -q '^usage:' ||
    { echo "$name: --help did not print the usage"; exit 1; }
tar -C "$(dirname "$bundle")" -cf - "$name" | zstd -T0 -10 -q -f -o "$out/$name.tar.zst"
if command -v sha256sum > /dev/null; then sum=(sha256sum); else sum=(shasum -a 256); fi
(cd "$out" && "${sum[@]}" "$name.tar.zst" > "$name.tar.zst.sha256")
echo "packed $out/$name.tar.zst ($(wc -c < "$out/$name.tar.zst" | tr -d ' ') bytes)"
