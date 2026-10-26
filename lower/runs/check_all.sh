#!/bin/zsh
# Run the blind checker on every finished witness against its CERT r.
cd "${0:A:h}/.."
for f in runs/k*_h*.log; do
  line=$(grep CERT $f) || continue
  k=$(echo $line | sed -E 's/.*CERT k=([0-9]+) h=.*/\1/'); h=$(echo $line | sed -E 's/.*h=([0-9]+) classes.*/\1/')
  r=$(echo $line | sed -E 's/.* r=([0-9-]+) .*/\1/'); w=runs/k${k}_h${h}_w.bin
  [ -f $w ] || continue
  out=runs/check_k${k}_h${h}.txt
  [ -f $out ] && grep -q PASS $out && continue
  nice -n 5 check/lowcheck $k $h $w 4 $r > $out 2>&1
  echo "k=$k h=$h: $(grep -E 'PASS|FAIL' $out | tail -1)"
done
