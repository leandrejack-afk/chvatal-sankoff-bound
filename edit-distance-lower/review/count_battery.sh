#!/bin/sh
# one job: rules k n m ; compares C brute-force totals with the Python automaton totals
r=$1; k=$2; n=$3; m=$4
P=${PYTHON:-python3}  # a Python with numpy, scipy, mpmath
c=$(nice -n 10 ./count_canon $r $k $n $m)
a=$(nice -n 10 $P automaton_count.py $r $k $n $m)
ct=$(echo "$c" | grep '^totals'); head=$(echo "$c" | head -1); fm=$(echo "$c" | grep first-moment)
if [ "$ct" = "$a" ]; then v=IDENTICAL; else v=DIFFERENT; fi
echo "$head | $fm | brute vs automaton: $v | $ct"
