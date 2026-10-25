#!/bin/bash
# Rebuild with kbig, rewrite in kary.c's format, compare byte for byte with kary.c's tables in data/.
T=big/repro
mkdir -p $T
for kh in "$@"; do
  set -- $kh
  ./kbig build $1 $2 $T > $T/build_k$1_h$2.out 2> $T/build_k$1_h$2.err || { echo "k=$1 h=$2 build failed"; continue; }
  ./kbig klue $1 $2 $T > /dev/null || { echo "k=$1 h=$2 klue failed"; continue; }
  if cmp -s $T/ktab_k$1_h$2.bin data/ktab_k$1_h$2.bin && cmp -s $T/wtypes_k$1_h$2.txt data/wtypes_k$1_h$2.txt; then
    echo "k=$1 h=$2 identical: $(shasum -a 256 $T/ktab_k$1_h$2.bin | cut -c1-16) $(cat $T/build_k$1_h$2.out | sed 's/ build .*//')"
  else
    echo "k=$1 h=$2 DIFFERENT"
  fi
done
