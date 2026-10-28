#!/bin/bash
# smoke.sh BUNDLE [STAGE] [GROUP...]
#
# Runs one fast pair of each sensor group (smoke_pairs.tsv: group, sensor, reference, secondary and, for a
# pair that takes one, a --subset polygon) through BUNDLE's executable and checks that each
# reaches STAGE: geometry (the default), imagery or done. A run is stopped once its stage line appears on
# stderr. `geometry` takes under 20 s per group, under three minutes for all eleven, and exercises every
# setup path, which is where a changed signature surfaces; `done` runs the pairs to completion (about 10
# minutes, and 25 more for NL2).
#
#     tools/release/smoke.sh ~/release/autorift_itslive_production-<version>-<platform>
#     tools/release/smoke.sh <bundle> imagery NL1 NL2
#
# The inputs are read as a normal run reads them, so the environment must supply what they need: AWS
# credentials for Landsat (AWS_PROFILE or the instance role), ~/.netrc or EARTHDATA_TOKEN for NISAR and
# Sentinel-1. SMOKE_SCRATCH names the directory the runs' scratch goes under (default: $TMPDIR or /tmp; it
# must be on disk, not tmpfs). Exits nonzero if any group fails.
set -u
bundle=$1; stage=${2:-geometry}; shift $(( $# > 1 ? 2 : 1 )); groups=("$@")
bin=$bundle/bin/autorift_itslive_production
[ -x "$bin" ] || { echo "no executable at $bin"; exit 2; }
case $stage in geometry|imagery) pat="^$stage " ;; done) pat="^wrote " ;; *) echo "STAGE is geometry, imagery or done"; exit 2 ;; esac
here=$(cd "$(dirname "$0")" && pwd)
scratch_root=${SMOKE_SCRATCH:-${TMPDIR:-/tmp}}
limit=${SMOKE_LIMIT:-1800}
fail=0
printf '%-8s %-6s %7s  %s\n' group result seconds detail
while IFS=$'\t' read -r group sensor ref sec subset; do
  [ ${#groups[@]} -eq 0 ] || [[ " ${groups[*]} " == *" $group "* ]] || continue
  work=$(mktemp -d "$scratch_root/smoke.XXXXXX"); log=$work/log; : > "$log"
  t0=$(date +%s); result=FAIL; detail=
  "$bin" "$sensor" "$ref" "$sec" "$work/out.nc" --scratch "$work/scratch" ${subset:+--subset "$subset"} > "$log" 2>&1 &
  pid=$!
  while kill -0 $pid 2>/dev/null; do
    if grep -q -E "$pat" "$log"; then result=PASS; break; fi
    [ $(( $(date +%s) - t0 )) -gt $limit ] && { detail="no $stage stage within ${limit} s"; break; }
    sleep 1
  done
  # The process may have exited between the last check and now.
  [ $result = FAIL ] && grep -q -E "$pat" "$log" && result=PASS
  kill -KILL $pid 2>/dev/null; wait $pid 2>/dev/null
  if [ $result = FAIL ] && [ -z "$detail" ]; then detail=$(grep -m1 '^error' "$log" | cut -c1-110); fi
  [ $result = PASS ] || fail=1
  printf '%-8s %-6s %7s  %s\n' "$group" $result $(( $(date +%s) - t0 )) "$detail"
  rm -rf "$work"
done < "$here/smoke_pairs.tsv"
exit $fail
