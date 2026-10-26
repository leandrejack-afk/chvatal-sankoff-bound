#!/bin/zsh
# After the h = 7 batch: k = 5 at h = 8, k = 4 at h = 9 (each warm-started from h - 1), k = 3 at h = 10.
cd "${0:A:h}/.."
until [ -f runs/h7.done ]; do sleep 20; done
run() { local k=$1 h=$2; shift 2; nice -n 10 env "$@" ./klow $k $h 150 8 runs/k${k}_h${h} > runs/k${k}_h${h}.log 2>&1; }
run 5 7
run 5 8 KLOW_INIT=runs/k5_h7_w.bin
run 4 8
run 4 9 KLOW_INIT=runs/k4_h8_w.bin
run 3 9
run 3 10 KLOW_INIT=runs/k3_h9_w.bin
echo done > runs/big.done
