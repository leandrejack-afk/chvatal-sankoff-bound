#!/bin/bash
# Reviewer's check of a kbig certificate, without writing into big/:
#   review-big/check_cert.sh k h  (reads big/k<k>h<h>/, writes review-big/c<k><h>/)
# 1. kbcertify.py re-run on symlinks (weights are regenerated into review-big/c<k><h>/);
# 2. the weights must equal the ones beside the certificate;
# 3. rvb_verify (emulated directed rounding + a-priori bound) and rvb_bound.py;
# 4. rv_theirweights.py (exact Fractions) for the referee's weights;
# 5. with FULL=1: rvcmp_big's certificate mode (every transition recomputed independently),
#    and rv_certmode_bound.py: its 6-decimal bound must equal the printed one.
set -u
K=$1; H=$2
KA=kary
PY=python3
S=$KA/big/k${K}h${H}
C=$KA/review-big/c${K}${H}
mkdir -p $C
for f in kbig_k${K}_h${H}.bin kvec_k${K}_h${H}.bin wtypes_k${K}_h${H}.txt; do ln -sf $S/$f $C/; done
cd $KA
echo "== certificate beside the table"; cat $S/certify_k${K}_h${H}.out
echo "== kbcertify.py re-run"; nice -n 10 $PY kbcertify.py $K $H $C | tee $C/certify_rerun.out
echo "== weights equal: $(cmp -s $C/weights_k${K}_h${H}.txt $S/weights_k${K}_h${H}.txt && echo yes || echo NO)"
echo "== rvb_verify"
nice -n 10 $KA/review-big/bin/rvb_verify $K $H $C/kbig_k${K}_h${H}.bin $C/kvec_k${K}_h${H}.bin $C/weights_k${K}_h${H}.txt $C/wtypes_k${K}_h${H}.txt 6 | tee $C/rvb.out
nice -n 10 $PY $KA/review-big/rvb_bound.py $K $C/rvb.out $S/certify_k${K}_h${H}.out
echo "== rv_theirweights"
(cd $KA/review && nice -n 10 $PY rv_theirweights.py $K $H $C $C/refweights.txt)
if [ "${FULL:-0}" = 1 ]; then
  echo "== rvcmp_big certificate mode (all transitions recomputed)"
  RVK_THREADS=${THREADS:-6} nice -n 10 $KA/rvcmp_big $K $H $C/kbig_k${K}_h${H}.bin $C/kvec_k${K}_h${H}.bin $C/refweights.txt | tee $C/rvcmp_cert.txt
  echo "== bound from the independent transitions against the printed certificate"
  $PY $KA/review-big/rv_certmode_bound.py $K $C/rvcmp_cert.txt $S/certify_k${K}_h${H}.out
fi
echo "== sha256 (sha256nc: F_NOCACHE reads, same digest as shasum -a 256)"
$KA/review-big/bin/sha256nc $S/kbig_k${K}_h${H}.bin $S/kvec_k${K}_h${H}.bin $S/weights_k${K}_h${H}.txt $S/wtypes_k${K}_h${H}.txt
echo "== against $S/SHA256_k${K}_h${H}.txt"; cat $S/SHA256_k${K}_h${H}.txt 2>/dev/null
