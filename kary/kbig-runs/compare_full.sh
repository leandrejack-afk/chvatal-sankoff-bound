#!/bin/bash
# kbig-runs/compare_full.sh k h (from kary/): every state of kbig_k<k>_h<h>.bin against the earlier review's own
# transition function (rvcmp_big), as four quarter runs (blocks of 1,024 states, stride 4, phases 0..3)
cd "$(dirname "$0")/.."
mkdir -p big/logs
K=$1; H=$2
for P in 0 1 2 3; do
  RVCMP_STRIDE=4 RVCMP_PHASE=$P nice -n 10 ./rvcmp_big $K $H big/k${K}h${H}/kbig_k${K}_h${H}.bin > big/logs/rvcmp_k${K}_h${H}_q$P.log 2>&1
  echo "phase $P exit $?" >> big/logs/rvcmp_k${K}_h${H}_q$P.log
done
