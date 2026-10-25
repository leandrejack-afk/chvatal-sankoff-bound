#!/bin/bash
# Logs for publication, with the final binaries: (1) kbig rebuilds of all 16 earlier tables,
# rewritten by "kbig klue", compared byte for byte with kary.c's tables in data/;
# (2) kbverify on the kbig files against kverify's recorded certificates.
cd "$(dirname "$0")/.."
mkdir -p big/logs
{
  echo "kbig binary: $(shasum -a 256 kbig | cut -c1-64)  (built from kbig.c $(shasum -a 256 kbig.c | cut -c1-64))"
  echo "$(date '+%Y-%m-%d %H:%M')"
  nice -n 10 ./kbig-runs/repro.sh "2 2" "2 3" "2 4" "2 5" "2 6" "2 7" "2 8" "3 2" "3 3" "3 4" "3 5" "3 6" "4 2" "4 3" "4 4" "4 5"
  for kh in "3 6" "4 5"; do set -- $kh; echo "sha256 of the rebuilt k=$1 h=$2 table: $(./sha256nc big/repro/ktab_k$1_h$2.bin | cut -c1-64); recorded in data/SHA256_kary.txt: $(grep ktab_k$1_h$2 data/SHA256_kary.txt | cut -c1-64)"; done
} > big/logs/repro_identity.log 2>&1
{
  echo "kbverify binary: $(shasum -a 256 kbverify | cut -c1-64)  (built from kbverify.c $(shasum -a 256 kbverify.c | cut -c1-64))"
  for kh in "2 8" "3 2" "3 3" "3 4" "3 5" "3 6" "4 2" "4 3" "4 4" "4 5"; do
    set -- $kh
    cp data/kvec_k$1_h$2.bin data/weights_k$1_h$2.txt big/repro/
    a=$(nice -n 10 ./kbverify $1 $2 big/repro 2>/dev/null | grep -E "^(lambda_cert|rho_lower)=" | tr '\n' ' ')
    b=$(grep -E "^(lambda_cert|rho_lower)=" data/certify_k$1_h$2.out | tr '\n' ' ')
    if [ "$a" = "$b" ]; then echo "k=$1 h=$2 identical: $a"; else echo "k=$1 h=$2 DIFFERENT: kbverify $a / kverify $b"; fi
  done
} > big/logs/kbverify_vs_kverify.log 2>&1
