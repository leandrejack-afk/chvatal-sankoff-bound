#!/bin/bash
# Rerun of every exact certificate with certify_exact.py (directed decimals), 2026-10-25.
P=${PYTHON:-python3}   # needs numpy and mpmath
cd "$(dirname "$0")"
(
  for kz in "3 0.1442" "4 0.1147" "5 0.0966" "6 0.0841" "7 0.0748" "8 0.0675" "9 0.0616" \
            "10 0.0568" "11 0.0527" "12 0.0492" "13 0.0462" "14 0.0436" "15 0.0412"; do
    set -- $kz; nice -n 10 $P certify_exact.py $1 2 $2
  done
  for kz in "3 0.1652" "4 0.1285" "5 0.1063" "6 0.0912" "7 0.0803" "8 0.0719" "9 0.0652" \
            "10 0.0597" "11 0.0552"; do
    set -- $kz; nice -n 10 $P certify_exact.py $1 3 $2
  done
) > certify_small.log 2>&1 &
for kz in "15 0.0426" "14 0.0452" "13 0.0480" "12 0.0514"; do
  set -- $kz; /usr/bin/time -l nice -n 10 $P certify_exact.py $1 3 $2
done > certify_big.log 2>&1
wait
echo ALL-DONE >> certify_big.log
