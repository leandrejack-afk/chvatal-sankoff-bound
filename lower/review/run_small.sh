#!/bin/zsh
# Run klow on small (k,h) with a dump, then brute-force check the dumped witness on the full vector.
# usage: review/run_small.sh "k h iters" ...   (run from lower/)
cd "${0:A:h}/.."
for spec in "$@"; do
  set -- ${=spec}
  k=$1 h=$2 it=$3
  out=review/w/k${k}_h${h}
  nice -n 5 ./klow $k $h $it 2 $out > $out.log 2>&1
  cert=$(grep CERT $out.log)
  r=$(echo $cert | sed -E 's/.* r=(-?[0-9]+) .*/\1/')
  res=$(python3 review/bf_low.py $k $h ${out}_w.bin $r 2>&1 | grep -E "PASS|FAIL")
  echo "k=$k h=$h | $cert | $res"
done
