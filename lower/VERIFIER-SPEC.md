# Spec: independent checker for k-ary Chvátal–Sankoff lower-bound certificates

Write a checker, in C (one file, `check/lowcheck.c`, built with
`clang -O3 -march=native -lpthread`), that verifies certificates produced by another
program. **Do not read any other `.c` or `.py` file under `c31-lcs/` (in particular
`lower/klow.c`, `lower/klow_plain.c`).** The point is a second implementation written
only from this spec and the paper. You may read the papers in `c31-lcs/papers/`
(`lueker2009.txt` Section 2 is the method; `2407.18113.pdf`, Section 3, Facts 9 and 10,
states it for any alphabet size).

## The mathematics

gamma_k is the limit of E[LCS(X, Y)]/n for independent uniform random strings X, Y of
length n over an alphabet of k letters. Fix h. For strings s, t of length h, Lueker's
recurrences give, for a vector w indexed by all pairs (s, t),

- if s[0] == t[0]:  T(w)(s, t) = 1 + (1/k^2) * sum over letters c, c' of w(s[1..]c, t[1..]c')
- if s[0] != t[0]:  T(w)(s, t) = max( (1/k) * sum_c w(s[1..]c, t),  (1/k) * sum_c w(s, t[1..]c) )

where s[1..]c is s with its first letter removed and c appended. Lemma 2.3 of Lueker
(stated there for k = 2; the proof carries over to any k unchanged): if T(w) >= w + r
componentwise, then gamma_k >= 2r/(1+r).

## What a certificate is

A witness file of little-endian `uint32` values. Each value is a fixed-point number with
unit 2^-27 (value v means v / 2^27). There is one value per **class** of pairs:

1. For a pair (s, t), let z = s.t, the concatenation, of length 2h.
2. RGS(z) relabels the letters of z by order of first appearance: the first letter
   becomes 0, the next letter not seen before becomes 1, and so on. The result is a
   sequence with at most k distinct labels.
3. Consider ALL sequences of length 2h that are such relabelings (restricted growth
   strings: first entry 0, each entry at most 1 + the maximum of the earlier entries)
   with at most k distinct labels, sorted lexicographically. The rank of a sequence is
   its position in this order, starting from 0.
4. The two orientations of the pair give rank(RGS(s.t)) and rank(RGS(t.s)). The class of
   (s, t) is identified by the smaller of the two. Two pairs are in the same class iff
   they have the same smaller rank (this is the same as: one is obtained from the other by
   a permutation of the alphabet, possibly combined with swapping s and t).
5. The witness lists one value per class, in increasing order of that smaller rank. So
   the number of values equals the number of RGS sequences z of length 2h (at most k
   labels) with rank(z) <= rank(RGS(swap(z))), where swap(z) = z[h..2h-1] . z[0..h-1].

w(s, t) for an arbitrary pair means the value stored for its class.

## What the checker computes

For every class, take its representative (s, t) (the orientation with the smaller rank;
any member of the class gives the same result, see below), and compute with exact
integer arithmetic (64-bit or 128-bit sums, never floating point):

- if s[0] == t[0]:  T = 2^27 + floor( S / k^2 ),  S = sum over all k^2 letter pairs (c, c')
  of w(s[1..]c, t[1..]c')
- if s[0] != t[0]:  T = floor( max(A, B) / k ),  A = sum over the k letters c of
  w(s[1..]c, t),  B = sum over c of w(s, t[1..]c)

Every floor rounds down, so T is at most the true value of Lueker's T at that pair.
Then r = min over classes of (T - w(class)), an exact integer (it may be negative for a
bad witness). Print r, the class where the minimum occurs (its rank and the strings),
the maximum of T - w, and the bound 2r / (2^27 + r) rounded DOWN to 9 decimals using
exact integer arithmetic (for r <= 0 print that no bound follows). Also print the number
of classes the checker computed and check the file holds exactly that many values (fail
loudly otherwise).

Why a class representative suffices: w is constant on classes, and T commutes with
permuting the alphabet and with swapping s and t, so T(w) is constant on classes too.

## Command line

`lowcheck k h witness.bin [threads] [claimed_r]` : print the result; if claimed_r is
given, exit 0 iff the computed r >= claimed_r, and print PASS or FAIL.

## Tests you must run and report

1. **Brute force for small (k, h):** for k^(2h) up to about 10^8, also expand the witness
   to the full vector over all pairs (via the class map), compute T at EVERY pair (s, t)
   directly from the definition above, and check that the minimum of T - w over all pairs
   equals the class-based r, and that every pair's value of T equals its class
   representative's T. Also check, for small k and h, that every alphabet permutation and
   the swap map each pair into the same class (exhaustively over all permutations for
   k <= 5).
2. **Class counts:** your computed class count must equal the witness length. Known counts
   (from the other program, so treat them as claims to test, not facts): k=2 h=3: 20;
   k=2 h=6: 1056; k=3 h=3: 70; k=3 h=5: 5002; k=4 h=4: 1451; k=5 h=4: 1993; k=6 h=5:
   55091; k=7 h=4: 2151; k=10 h=3: 117; k=3 h=7: 399310; k=4 h=6: 350891;
   k=6 h=7: 57365191.
3. **The witnesses in `lower/check/witness/`** (files `k<k>_h<h>_w.bin`; the claimed r
   values are in `claims.txt`). Report your r for each. Your r must be >= the claimed r
   (your T rounds less than the producing program does, so it may be larger by a few
   units; report the difference).
4. **Tamper test:** flip one value in a copy of a witness (for example subtract 2^20 from
   one entry, or add 2^20) and show the checker reports a smaller r or a negative one.
5. **Reproduce known bounds:** with the witness files above, the bounds for k = 2 must
   agree with Lueker's Table I (h = 3: 0.747922, h = 6: 0.770273) to 6 decimals.

## Performance

The real witnesses will be large: k=6..15 at h=7 (up to about 95 million classes,
0.4 GB), k=5 h=8 (about 640 million classes, 2.6 GB), k=4 h=9 (about 1.43 billion classes,
5.7 GB), maybe k=3 h=11 (about 2.6 billion, 10.5 GB). The Mac has 24 GB of memory with
other programs running: keep the checker's memory to the witness (read it with `mmap`
read-only, or into memory) plus at most about 1.5 GB of index structures, and use
pthreads (8 threads by default). A class lookup must cost O(h) arithmetic, not a search
through a table of all classes: compute the lexicographic rank of an RGS from a table of
completion counts, and turn the smaller rank into the dense position with a bitmap over
all ranks plus per-word prefix counts (or any other O(1) or O(log) structure you design).
Target: the k=4 h=9 check in under an hour on 8 threads.

## Deliverables

- `lower/check/lowcheck.c` and a short `lower/check/TESTS.txt` that lists each test you
  ran with its command and output (the harness may refuse a file named like a report; if
  it does, put the test log in your final reply instead).
- Your final reply: what you built, how the rank and class map work, the results of every
  test above, anything in this spec you found ambiguous or wrong, and anything suspicious
  in the witnesses.
