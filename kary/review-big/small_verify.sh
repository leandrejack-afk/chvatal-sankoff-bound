#!/bin/bash
# For each (k, h): a folder of read-only symlinks (my kbig rebuild + data/ vector, weights, wtypes),
# then kverify's recorded numbers (data/certify_*.out), kbverify, and rvb_verify.
KA=kary
D=$KA/data
RB=$KA/review-big
for kh in "$@"; do
  set -- $kh; k=$1; h=$2; T=$RB/t/k${k}h${h}; mkdir -p $T
  ln -sf $RB/repro/kbig_k${k}_h${h}.bin $T/
  for f in kvec_k${k}_h${h}.bin weights_k${k}_h${h}.txt wtypes_k${k}_h${h}.txt; do ln -sf $D/$f $T/; done
  echo "== k=$k h=$h"
  grep -E "lambda_cert|rho_lower" $D/certify_k${k}_h${h}.out | sed 's/^/  kverify (data\/certify): /'
  nice -n 10 $KA/kbverify $k $h $T 2>&1 | grep -E "lambda_cert|rho_lower|streamed|bad|error" | sed 's/^/  kbverify: /'
  nice -n 10 $RB/bin/rvb_verify $k $h $T/kbig_k${k}_h${h}.bin $T/kvec_k${k}_h${h}.bin $T/weights_k${k}_h${h}.txt $T/wtypes_k${k}_h${h}.txt 4 | grep -E "lambda_emu|rho_emu|nearest|OK|NOT|bad|start" | sed 's/^/  rvb_verify: /'
done
