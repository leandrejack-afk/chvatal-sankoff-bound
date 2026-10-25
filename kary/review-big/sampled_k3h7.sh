#!/bin/bash
# Reviewer's sampled checks of the k = 3, h = 7 table once the build has finished
# (different sample and block phase from big/checks_k3h7.sh).
KA=kary
cd $KA
PY=python3
until grep -q "accepted transitions=" big/k3h7/build.out 2>/dev/null; do
  pgrep -f "kbig build 3 7" > /dev/null || { echo "build ended without its final line"; exit 2; }
  sleep 10
done
cat big/k3h7/build.out
tail -1 big/k3h7/build.err
date
RVK_THREADS=4 RVCMP_STRIDE=101 RVCMP_PHASE=50 nice -n 10 ./rvcmp_big 3 7 big/k3h7/kbig_k3_h7.bin; echo "rvcmp_big exit $?"
date
nice -n 10 $PY check_table_big.py 3 7 big/k3h7 600 2027; echo "check_table_big exit $?"
date
