#!/bin/bash
# check_inputs.sh DIR: every reference and secondary input of DIR/<group>.tsv requested for its first byte
# (ASF through Earthdata Login with ~/.netrc), its final HTTP status written to DIR/<group>_check.tsv.
cd "$1"
one() {
  jar=$(mktemp); code=$(curl -s -o /dev/null -L --netrc -c $jar -b $jar -r 0-0 -m 120 -w "%{http_code}" "$1"); rm -f $jar
  printf "%s\t%s\n" "$code" "$1"
}
export -f one
for g in L4L5 L7 L8L9 S1 S2 NL1 NL2; do      # not S1BURST: a request makes the extractor generate the burst
  [ -f $g.tsv ] || continue
  cut -f5,6 $g.tsv | tr '\t' '\n' | sort -u | xargs -P 24 -I{} bash -c 'one "$@"' _ {} > ${g}_check.tsv
  echo "$g: $(cut -f1 ${g}_check.tsv | sort | uniq -c | tr '\n' ' ')"
done
