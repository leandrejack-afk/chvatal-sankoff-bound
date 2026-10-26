#!/bin/zsh
# Runs every lowcheck test from lower/check/ and prints each command before its output.
cd "${0:A:h}/.." || exit 1  # the check/ folder

run() {
    print -r -- "\$ $1"
    eval "$1"'; print -r -- "[exit status ${pipestatus[*]}]"' 2>&1
    print
}

TMP=$(mktemp -d)
print "[exit status] lists every command of a pipeline in order, so 'lowcheck ... | grep' shows lowcheck's status first."
EDIT='import struct,sys; d=bytearray(open(sys.argv[1],"rb").read()); i=int(sys.argv[3]); v=struct.unpack_from("<I",d,4*i)[0]+int(sys.argv[4]); struct.pack_into("<I",d,4*i,v); open(sys.argv[2],"wb").write(bytes(d)); print("entry",i,"is now",v)'
CASES=(2:3 2:6 3:3 3:5 4:4 5:4 6:5 7:4 10:3 3:7 4:6)

print "=== BUILD"
run "clang -O3 -march=native -lpthread lowcheck.c -o lowcheck"
run "clang -O3 -march=native -Wall -Wextra -Wconversion -Wno-sign-conversion -lpthread lowcheck.c -o lowcheck && echo no warnings"

print "=== TEST 3: every witness against claims.txt (also the class counts of TEST 2 and the bounds of TEST 5)"
while read -r _ kk hh _ _ rr _; do
    run "./lowcheck ${kk#k=} ${hh#h=} witness/k${kk#k=}_h${hh#h=}_w.bin 6 ${rr#r=}"
done < witness/claims.txt

print "=== TEST 1: brute force over all k^(2h) pairs for every witness"
for kh in $CASES; do
    k=${kh%%:*}; h=${kh##*:}
    run "./lowcheck --brute $k $h witness/k${k}_h${h}_w.bin 6 | sed -n '/^r = /p;/--brute:/,\$p'"
done

print "=== TEST 2: class counts without a witness (k=6 h=7 claim, then the production sizes)"
for kh in 6:7 15:7 5:8 4:9 3:11; do
    k=${kh%%:*}; h=${kh##*:}
    run "./lowcheck --count $k $h 6"
done

print "=== TEST 4: tamper tests on copies (TMP=$TMP)"
run "python3 -I -c '$EDIT' witness/k4_h6_w.bin $TMP/a.bin 2514 1048576 && ./lowcheck 4 6 $TMP/a.bin 6 58317659 | grep -E '^r =|^min|^bound|^claimed'"
run "python3 -I -c '$EDIT' witness/k4_h6_w.bin $TMP/b.bin 9684 -1048576 && ./lowcheck 4 6 $TMP/b.bin 6 58317659 | grep -E '^r =|^min|^bound|^claimed'"
run "python3 -I -c '$EDIT' witness/k4_h6_w.bin $TMP/c.bin 2514 67108864 && ./lowcheck 4 6 $TMP/c.bin 6 58317659 | grep -E '^r =|^min|^bound|^claimed'"
run "python3 -I -c '$EDIT' witness/k2_h6_w.bin $TMP/d.bin 63 1048576 && ./lowcheck 2 6 $TMP/d.bin 6 84071083 | grep -E '^r =|^min|^bound|^claimed'"
run "python3 -I -c '$EDIT' witness/k3_h7_w.bin $TMP/e.bin 82 1048576 && ./lowcheck 3 7 $TMP/e.bin 6 68534490 | grep -E '^r =|^min|^bound|^claimed'"
run "head -c \$((4*350890)) witness/k4_h6_w.bin > $TMP/short.bin && ./lowcheck 4 6 $TMP/short.bin 6 58317659"
run "cat witness/k4_h6_w.bin witness/k2_h3_w.bin | head -c \$((4*350892)) > $TMP/long.bin && ./lowcheck 4 6 $TMP/long.bin 6 58317659"
run "./lowcheck 4 5 witness/k4_h6_w.bin 6"
run "./lowcheck 4 6 witness/k4_h6_w.bin 6 58317660 | tail -1"

print "=== TEST 5 extra: hand-derived h=1 witness w(match)=2^27, w(mismatch)=0; expect r = floor(2^27/k), bound 2/(k+1)"
run "python3 -I -c 'import struct,sys; open(sys.argv[1],\"wb\").write(struct.pack(\"<2I\", 1<<27, 0))' $TMP/h1.bin"
for k in 2 3 4; do
    run "./lowcheck $k 1 $TMP/h1.bin 6 | grep -E '^classes|^r =|^bound'"
done
run "python3 -I -c 'from fractions import Fraction as F
for k in (2, 3, 4):
    r = (1 << 27) // k; print(k, r, (F(2 * r, (1 << 27) + r) * 10**9).__floor__())'"

print "=== default thread count"
run "./lowcheck 2 3 witness/k2_h3_w.bin | head -1"

SP=${TMPDIR:-/tmp}/lowcheck-tests; mkdir -p $SP
print "=== EXTRA: naive Python reference (dictionary ranks, no bitmap, no shared code), all witnesses"
for kh in $CASES; do
    k=${kh%%:*}; h=${kh##*:}
    run "python3 -I $SP/impl-tamper/naive.py $k $h witness/k${k}_h${h}_w.bin"
done

print "=== EXTRA: witness facts (value range, the zero entry, first-letter identity of match pairs)"
for kh in $CASES; do
    k=${kh%%:*}; h=${kh##*:}
    run "python3 -I $SP/impl-tests/witness_facts.py $k $h witness/k${k}_h${h}_w.bin"
done

rm -rf "$TMP"
