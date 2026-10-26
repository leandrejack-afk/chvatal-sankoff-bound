# Lower bounds on the Chvátal–Sankoff constants for 3 to 15 letters

For two independent uniformly random strings of length n over an alphabet of k letters,
the expected length of a longest common subsequence is γ_k·n + o(n). This folder raises the
best known **lower** bounds on γ_k for every k from 3 to 15.

| k | previous best | its h | this work | h | classes |
|---|---|---|---|---|---|
| 3 | 0.68422 [R24] | 10 | **0.686509** | 11 | 2,615,147,350 |
| 4 | 0.614333 [H24] | 8 | **0.617356** | 9 | 1,431,743,147 |
| 5 | 0.56206 [R24] | 7 | **0.566171** | 8 | 637,624,313 |
| 6 | 0.51850 [R24] | 6 | **0.524151** | 7 | 57,365,191 |
| 7 | 0.48712 [R24] | 6 | **0.493089** | 7 | 82,035,085 |
| 8 | 0.46074 [R24] | 6 | **0.466991** | 7 | 92,494,661 |
| 9 | 0.43806 [R24] | 6 | **0.444595** | 7 | 95,063,689 |
| 10 | 0.41826 [R24] | 6 | **0.425055** | 7 | 95,440,625 |
| 11 | 0.40072 [R24] | 6 | **0.407778** | 7 | 95,473,791 |
| 12 | 0.38504 [R24] | 6 | **0.392336** | 7 | 95,475,506 |
| 13 | 0.37088 [R24] | 6 | **0.378406** | 7 | 95,475,555 |
| 14 | 0.35798 [R24] | 6 | **0.365743** | 7 | 95,475,556 |
| 15 | 0.34616 [R24] | 6 | **0.354154** | 7 | 95,475,556 |

The previous bounds are from Rosenfeld [R24, Table 1b] and, for k = 4, from Heineman et al.
[H24, Table 3] (also Soiffer et al., ISIT 2025). h is the length of the strings in the states of
Lueker's recurrence. Each bound is rounded down from a certificate computed with every division
rounded down (see below). The full certified values are the CERT lines in `runs/` and
`runs/server/`.

## Method

Nothing in the method is new: it is Lueker's lower-bound technique [L09, Section 2] in the
k-letter form of Kiwi and Soto [KS09] and Rosenfeld [R24, Section 3]. For strings s, t of
length h and a vector w indexed by all pairs (s, t), let

- T(w)(as′, at′) = 1 + (1/k²) Σ_{c,c′} w(s′c, t′c′)
- T(w)(as′, bt′) = max( (1/k) Σ_c w(s′c, bt′), (1/k) Σ_c w(as′, t′c) ) for a ≠ b.

Lemma 2.3 of [L09] gives γ_k ≥ 2r/(1+r) whenever T(w) ≥ w + r componentwise. Lueker states
Section 2 for k = 2 only, but the proof uses nothing else about the alphabet: the mismatched
operator is a maximum of averages, so it commutes with adding a constant and with scaling by
α > 0; the matched operator satisfies T(αv) = αT(v) + 1 − α; and his equation (3) holds for every
k. So Observation 2.1 and Lemmas 2.2 and 2.3 hold as stated for any k. (Rosenfeld's Lemma 11
[R24] is the same statement in the form T(w, w − r) ≥ w + r.) A good w is found by iterating T
from 0; the witness is the last vector, and one more step, in which every division rounds down
(fixed point, unit 2⁻²⁷), gives r. So the computed T is never above the true one.

What is new is the storage, which makes one more step of h fit in a laptop's memory.

- **Classes.** w(s, t) does not change when the alphabet is permuted or s and t are swapped.
  The programs store one value per class: the letters of s·t are relabelled in order of first
  appearance (a restricted growth string with at most k labels), the string is ranked
  lexicographically from a table of completion counts, and the smaller of the ranks of s·t
  and t·s names the class. A bitmap over all ranks, with prefix counts, turns a rank into a
  position. For large k the number of classes stops growing with k (it is about half the Bell
  number B_{2h}), so every k ≥ 2h costs the same.
- **One step in two passes.** With F(x, y) = ⌊Σ_c w(xc, y) / k⌋ over pairs of lengths h − 1
  and h (alphabet symmetry only), the mismatched case is max(F(s′, bt′), F(t′, as′)) and the
  matched case is 1 + ⌊Σ_c F(t′, s′c) / k⌋, taking the larger of the two rounding orders so that
  the computed T keeps the swap symmetry. A letter that occurs in neither string stands for all
  k − u such letters at once.

For example k = 6, h = 7 needs 0.3 GB of memory (57,365,191 classes); Rosenfeld reports 3 to
4 GB for his h = 6 runs. The k = 3, h = 11 run needed about 19 GB and was made on a rented
32 GB server; all others ran on a 24 GB laptop.

## Checks

1. **Known values.** `klow_plain.c` (plain layout, no symmetry) reproduces Lueker's Table I
   [L09] for h = 1, …, 8 (0.666666 to 0.776860) and the values of [H24, Table 3] for
   k = 5, h = 5 (0.549817) and k = 6, …, 10, h = 4 (0.499229, 0.466481, 0.438799, 0.414876,
   0.393811), to all six digits printed there. `klow.c` also gives the Kiwi–Soto values that
   [H24] lists as the previous best (k = 3, h = 6: 0.671697; k = 4, h = 5: 0.599248;
   k = 5, h = 4: 0.539129; k = 6, h = 3: 0.479452) and the values of [H24] at k = 3, h = 9
   (0.682218) and k = 4, h = 7 and 8 (0.610590, 0.614333).
2. **Symmetric against plain.** With the same rounding route (`KLOW_FROUTE=1`), `klow_plain`
   and `klow` print identical lines (every step's minimum, maximum, spread and bound) on 17
   cases with k = 2, …, 10 and h = 3, …, 8.
3. **Rosenfeld's h.** At the history lengths of [R24] the program gives, for k = 6, …, 15 and
   h = 6, 0.518617, 0.487246, 0.460856, 0.438183, 0.418378, 0.400847, 0.385159, 0.370991,
   0.358097, 0.346283, each 0.00010 to 0.00013 above the published value (Rosenfeld explains
   the gap by the fixed precision p = 100000 he used); for k = 3, h = 10 it gives 0.684537
   against 0.68422, and for k = 5, h = 7, 0.562175 against 0.56206. So the gains in the table
   come from the larger h, apart from about 10⁻⁴. Logs: `runs/h6/`.
4. **Independent checker.** `check/lowcheck.c` was written from `VERIFIER-SPEC.md` and the
   papers only, without access to the producing programs. It recomputes T from the definition
   (sums over all k² letter pairs, no intermediate F) with exact integer arithmetic and its own
   class map, and checks its class counts against a Burnside orbit count. Its own tests
   (`check/TESTS.txt`) include a brute force over every pair for small cases, all alphabet
   permutations for k ≤ 5, tampered witnesses, and Lueker's Table I. On every witness in the
   table it reproduces r exactly (`runs/check_*.txt`, `runs/server/check_*.txt`).
5. **Review.** A separate fresh-context review (`review/`, summary in `review/REVIEW.md`) wrote
   its own exact-rational checker (`review/rcheck.c`) and found the certificates sound; on the
   witnesses for k = 3 (h = 11), 4 (h = 9), 5 (h = 8) and 6 to 15 (h = 7) its exact r is equal to
   or slightly above the claimed r (`runs/rcheck_*.txt`, `runs/server/rcheck_k3_h11.out`).
6. **Two machines.** The iteration is integer arithmetic only. Cold-start runs of three small
   cases on ARM macOS (clang) and x86 Linux (gcc) give byte-identical witnesses, and a rerun of
   k = 6, h = 7 gives the same file. `witness-sha256.txt` lists the SHA-256 of every witness
   with the command that regenerates it.

## Reproducing

```
make                           # klow, klow_plain, check/lowcheck
./klow 6 7 150 8 k6_h7         # iterate; writes the witness k6_h7_w.bin and prints the CERT line
check/lowcheck 6 7 k6_h7_w.bin 8 47667784    # independent check of that witness
KLOW_FROUTE=1 ./klow_plain 3 5 100 8         # plain layout, same rounding route, small k and h
```

Arguments of `klow`: k, h, maximum number of steps, threads, output prefix. The run stops when
T(w) − w is constant to within 2 units, then writes the witness and prints its r and bound.
`KLOW_INIT=<witness for h − 1>` starts from the previous witness instead of 0.

## References

- [H24] G. T. Heineman et al., Improved lower bounds on the expected length of longest common
  subsequences, arXiv:2407.10925 (2024); D. Soiffer, A. Salls, C. Miller, D. Reichman,
  G. Sárközy, G. Heineman, ISIT 2025.
- [KS09] M. Kiwi, J. Soto, On a speculated relation between Chvátal–Sankoff constants of
  several sequences, Combin. Probab. Comput. 18 (2009), doi:10.1017/S0963548309009900.
- [L09] G. S. Lueker, Improved bounds on the average length of longest common subsequences,
  J. ACM 56(3) (2009) 17.
- [R24] M. Rosenfeld, Upper bounds on the average edit distance between two random strings,
  arXiv:2407.18113 (2024).

These programs were written with the assistance of a large language model; the results rest
on the checks listed above.
