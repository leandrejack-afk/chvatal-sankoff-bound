# Fresh-context review, 26 October 2026

A reviewer with no part in writing `klow.c` and `kedit.c` was asked to find any way the
certificates could be wrong. The files in this folder are its tools; this note summarizes its
verdict, written up from its report.

**Verdict.** `klow.c` certificates: sound. `kedit.c` certificates: sound, with one required
edit outside the program: `results.py` printed the edit-distance upper bounds truncated, which
rounds an upper bound down; it now rounds them up.

**What it checked.**

- The mathematics: the computed operator is never above Lueker's T at any class, including the
  matched case (both rounding orders are lower bounds), the multiplicity k − u of letters absent
  from both strings, the use of swap symmetry inside F, and the shift by the minimum. Lemma 2.3
  of Lueker's paper carries over to every k (the mismatched operator is translation invariant and
  positively homogeneous, the matched one satisfies T(αv) = αT(v) + 1 − α).
- The class machinery: Burnside orbit counts equal every logged class count; `klow.c`'s own rank,
  unrank and successor functions, compiled into a harness (`klow_harness.c`), agree exhaustively
  at 23 small (k, h), including k > 2h and k = 30, and on random tests up to k = 3, h = 12.
- The printed bounds: 40 CERT lines checked with exact fractions; the dumped witness is the input
  of the step whose r is printed.
- Overflow: the spread of values grows about 1.2h (about 14 at h = 11, against a ceiling of 32);
  sums are 64-bit and an overflow stops the run.
- An exact-rational checker (`rcheck.c`) on the witnesses of k = 5 to 15 at h = 7: exact r equal
  to the claimed r or above it by under one unit; 100,000 to 300,000 random raw pairs per witness
  agreed with their class representative.
- Mutants (`mut/`): ceiling instead of floor, a wrong multiplicity, a wrong lookup; the exact
  checkers catch over-claims as small as 0.19 units.

**Minor findings, all addressed.** A guard against overflow of the class counts was added (it
could only trigger at h ≥ 13 with k ≥ 9, as a crash rather than a wrong result); the bound
printers refuse impossible values; `klow_plain.c` checks its dump; the README now explains why
Lueker's Lemma 2.3 holds for every k, has a Makefile and a list of witness hashes.
