#!/bin/zsh
# review5 batch 2 (one core): author's canon_check5 on small cases, author's count_check5 rerun (diff
# against its log), exact finite-n first moment, Monte Carlo.
cd "${0:A:h}"
PY=${PY:-python3}  # needs numpy, scipy, mpmath
echo "== author's canon_check5 binary, small cases"
for a in "2 6 6" "3 4 5" "2 4 9" "4 3 4"; do nice -n 10 ../canon_check5 ${=a}; done
echo "== author's count_check5.py rerun vs count_check5.log"
(cd .. && nice -n 10 $PY -B count_check5.py) > count_check5_rerun.out 2>&1
if cmp -s count_check5_rerun.out ../count_check5.log; then echo "count_check5 rerun: byte-identical to count_check5.log"; else echo "count_check5 rerun DIFFERS:"; diff count_check5_rerun.out ../count_check5.log | head; fi
echo "== exact finite-n first moment"
for a in "2 10 10" "2 12 12" "2 8 12" "3 8 8" "4 6 6" "5 5 5" "8 4 4"; do nice -n 10 $PY -B fm_exact.py ${=a}; done
echo "== Monte Carlo"
for k in 2 3 4 5 6 7 8 16 32; do nice -n 10 ./mc5 $k 4000 30 $((1000 + k)); done
for k in 2 4; do nice -n 10 ./mc5 $k 1000 2000 $((2000 + k)); done
echo "== done"
