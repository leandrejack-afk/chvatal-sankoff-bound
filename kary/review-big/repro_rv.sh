#!/bin/bash
# Reviewer copy of big/repro.sh: same steps, output under review-big/repro, binary review-big/bin/kbig
# (text section identical to kary/kbig). Converted ktab copies are deleted after the comparison.
T=review-big/repro
KB=review-big/bin/kbig
mkdir -p $T
for kh in "$@"; do
  set -- $kh
  $KB build $1 $2 $T > $T/build_k$1_h$2.out 2> $T/build_k$1_h$2.err || { echo "k=$1 h=$2 build failed"; continue; }
  $KB klue $1 $2 $T > /dev/null || { echo "k=$1 h=$2 klue failed"; continue; }
  if cmp -s $T/ktab_k$1_h$2.bin data/ktab_k$1_h$2.bin && cmp -s $T/wtypes_k$1_h$2.txt data/wtypes_k$1_h$2.txt; then
    echo "k=$1 h=$2 identical: $(shasum -a 256 $T/ktab_k$1_h$2.bin | cut -c1-16) $(sed 's/ build .*//' $T/build_k$1_h$2.out) | $(tail -1 $T/build_k$1_h$2.err)"
  else
    echo "k=$1 h=$2 DIFFERENT"
  fi
  rm -f "review-big/repro/ktab_k$1_h$2.bin"
done
