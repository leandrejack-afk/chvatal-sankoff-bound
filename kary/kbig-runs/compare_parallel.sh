#!/bin/bash
# kbig-runs/compare_parallel.sh k h P1 P2 ... (from kary/): like compare_full.sh, but the given quarters run side by side
cd "$(dirname "$0")/.."
mkdir -p big/logs
K=$1; H=$2; shift 2
for P in "$@"; do
  ( RVCMP_STRIDE=4 RVCMP_PHASE=$P nice -n 10 ./rvcmp_big $K $H big/k${K}h${H}/kbig_k${K}_h${H}.bin > big/logs/rvcmp_k${K}_h${H}_q$P.log 2>&1
    echo "phase $P exit $?" >> big/logs/rvcmp_k${K}_h${H}_q$P.log ) &
done
wait
