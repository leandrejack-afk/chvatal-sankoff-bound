#!/bin/zsh
# k = 6..15 at h = 7 (Rosenfeld 2024 stopped at h = 6 for these alphabets)
cd "${0:A:h}/.."
for k in 6 7 8 9 10 11 12 13 14 15; do
  nice -n 10 ./klow $k 7 150 8 runs/k${k}_h7 > runs/k${k}_h7.log 2>&1
done
echo done > runs/h7.done
