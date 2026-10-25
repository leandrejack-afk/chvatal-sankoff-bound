#!/bin/bash
# Eigen phase and certificate for one kbig table, run from kary/ (the table in big/k<k>h<h>/,
# kary.c's (h-1) table and vector in data/):
#   kbig-runs/pipeline.sh k h zcenter        (TOL1, TOL2 optional, see 2 and 4)
# 1. start vector: the (h-1) vector of data/ projected onto the new states;
# 2. stage 1: 8 values of z, 0.0004 apart, centred on zcenter, to max/min ratio 1 + TOL1 (default 1e-8);
# 3. parabola through the three best bounds (zfit.py); its vertex rounded to 4 decimals = zc;
# 4. stage 2: the 8 values zc - 0.0003 .. zc + 0.0004, from the stage-1 vectors, to 1 + TOL2 (default 1e-10);
# 5. the z with the smallest bound is certified (kweights.py, kbverify, interval arithmetic).
set -euo pipefail
K=$1; H=$2; C=$3
D=big/k${K}h${H}
PY=${PY:-python3}                 # needs mpmath
NICE="nice -n 10"
export KBIG_QOS=user
log() { echo "$(date '+%H:%M:%S') $*" | tee -a $D/pipeline.log; }

log "project k=$K h=$H"
$NICE ./kbig project $K $H $D data | tee -a $D/pipeline.log
Z1=$($PY -c "c=$C; print(' '.join(f'{c + 0.0004 * (i - 3.5):.4f}' for i in range(8)))")
log "stage 1: z = $Z1"
$NICE ./kbig eigen $K $H $D $D/start_k${K}_h${H}.bin ${TOL1:-1e-8} 120 $Z1 > $D/eigen1.out 2> $D/eigen1.err
cat $D/eigen1.out >> $D/pipeline.log
mv $D/kvecs_k${K}_h${H}.bin $D/kvecs_stage1.bin
$PY zfit.py $D/eigen1.out | tee $D/zfit1.out | tee -a $D/pipeline.log
ZC=$(sed -n 's/^4-decimal z nearest the vertex: //p' $D/zfit1.out)
if [ -z "$ZC" ]; then log "no vertex inside the stage-1 range; stopping"; exit 3; fi
Z2=$($PY -c "c=$ZC; print(' '.join(f'{c + 0.0001 * (i - 3):.4f}' for i in range(8)))")
log "stage 2: z = $Z2"
$NICE ./kbig eigen $K $H $D $D/kvecs_stage1.bin ${TOL2:-1e-10} 120 $Z2 > $D/eigen2.out 2> $D/eigen2.err
cat $D/eigen2.out >> $D/pipeline.log
$PY zfit.py $D/eigen2.out | tee $D/zfit2.out | tee -a $D/pipeline.log || true
ZB=$($PY - "$D/eigen2.out" <<'EOF'
import re, sys
best = min((float(m.group(2)), m.group(1)) for m in
           (re.match(r"z=([0-9.]+) lambda in \[[0-9.]+, [0-9.]+\] bound=([0-9.]+)", l) for l in open(sys.argv[1])) if m)
print(best[1])
EOF
)
log "best z = $ZB; certifying"
cp $D/kvec_k${K}_h${H}_z${ZB}.bin $D/kvec_k${K}_h${H}.bin
$NICE $PY kbcertify.py $K $H $D | tee $D/certify_k${K}_h${H}.out | tee -a $D/pipeline.log
log "done"
