# Three and four letters

Lueker's automaton M_h (J. ACM 56(3), 2009, Section 3.2) gives upper bounds on the
Chvátal–Sankoff constant γ_k for every alphabet size k, but his computations are for
k = 2 only. Run with k = 3 and k = 4, it gives

| k | h | states | classes | z | λ_cert (rounded up) | bound |
|---|---|---|---|---|---|---|
| 3 | 2 | 10 | 21 | 0.1442 | 0.111829497 | 0.764972 |
| 3 | 3 | 249 | 63 | 0.1652 | 0.155666008 | 0.754807 |
| 3 | 4 | 5,116 | 141 | 0.1821 | 0.198535252 | 0.747887 |
| 3 | 5 | 98,041 | 279 | 0.1948 | 0.235651554 | 0.743282 |
| 3 | 6 | 1,711,682 | 525 | 0.2048 | 0.268056913 | **0.739975** |
| 4 | 2 | 10 | 52 | 0.1147 | 0.110575331 | 0.707648 |
| 4 | 3 | 309 | 292 | 0.1285 | 0.148063383 | 0.696851 |
| 4 | 4 | 9,427 | 1,108 | 0.1402 | 0.186602820 | 0.689064 |
| 4 | 5 | 290,442 | 3,460 | 0.1496 | 0.222621807 | **0.683622** |

so **γ_3 < 0.739975** and **γ_4 < 0.683622**. The best upper bounds before were
γ_3 < 0.765803 and γ_4 < 0.708236 (V. Dančík, *Expected length of longest common
subsequences*, PhD thesis, University of Warwick, 1994, Figure 4.6), and the automata with
h = 2, which have ten states, already improve on them. States are counted after the
reduction by renaming letters; z is the double nearest to the value shown; the full
certificate output for each row is in `certificates/certify_k<k>_h<h>.out`.

## What changes for k > 2

The automaton itself (Figure 7 of the paper) does not depend on k. What does is the set of
inputs and their weights. The minimal match pairs over k letters are the pairs (a'c, b'c)
in which a' and b' are strings over the other k − 1 letters with no letter in common.
After truncation to the last h characters of each string, a class (α'c, β'c) has weight
z^(|α'|+|β'|+2) times

- 1 if both strings are shorter than h,
- 1/(1 − (k − 1 − |B|)z) if only the top string has length h,
- 1/(1 − (k − 1 − |A|)z) if only the bottom string has length h,
- Σ_{s+t≤f} f!/(s! t! (f−s−t)!) g(|A|, s) g(|B|, t) if both have length h,

where A and B are the letter sets of α' and β', f = k − 1 − |A| − |B|, and
g(a, s) = Σ_j (−1)^(s−j) C(s, j)/(1 − (a + j)z) counts strings over a + s letters that use
each of s given letters. For k = 2 this is Lueker's z^(h+1)/(1 − z). `kweights.py`
evaluates these at the exact double z to 80 digits and rounds up (and down, for the lower
bound on the spectral radius).

## Checks

1. Run with k = 2, `kary` gives Lueker's Table II state counts for h = 2..8, his λ to all
   printed digits, the 91,508,397 accepted transitions at h = 8, and re-certifies 0.826280
   (`certificates/certify_k2_h8.out`).
2. `kref.py` finds the input classes from the definition of a minimal match pair (lcs
   computations, not the closed form); they agree with `kary.c` for every run.
3. `kweights.py check`: for every class, the number of minimal match pairs of each size up
   to a limit between 11 and 15 equals the coefficient of the closed form.
4. `check_table.py` recomputes transitions with `kref.py` (Python, no shared code): every
   state for k = 3, h ≤ 5 and k = 4, h ≤ 4; 5,000 random states at (3, 6) and 2,000 at
   (4, 5) (`logs/check_table_*.log`). No mismatch.
5. `indep_check.py` builds the unreduced automaton with `kref.py` for k = 3, h ≤ 4 and
   k = 4, h ≤ 3; spectral radii agree to 12 or more digits and the orbit counts equal our
   state counts.
6. `kcover.py`: canonical decompositions of random string pairs (uniform, and biased toward
   long runs and few letters so that truncation occurs) are never rejected; at (3, 6) and
   (4, 5) alone, 529,139 inputs, 42,811 of them truncated (`logs/cover_*.log`).
   `redcontrols_cover.py`: two deliberately broken automata fail within 17 pairs.
7. `kverify` sets the rounding mode upward (then downward) in every worker block, with a
   self-test, checks that every entry of v is positive, and reports λ_cert and a lower bound
   on the spectral radius; they agree to at least 12 digits in every row.
8. `review/` is a separate implementation written in a fresh review pass from Lueker's
   paper, before reading the code above: its own state encoding, transition function,
   class grouping and weights (`rvk.c`, `rv_*.py`). It reproduces every state count,
   accepted-transition count and bound in the table, matches both headline tables state
   by state, re-runs the certificates, and gives exact rational certificates for h = 2
   (`review/log_h2cert.txt`: rational z, λ and v with A_z v ≤ λ v checked in exact
   arithmetic, giving γ_3 ≤ 0.764972 and γ_4 ≤ 0.707648).

Known limits, closed by the review checks: `kweights.py` rounds an 80-digit value rather
than the exact rational (review/log_theirweights.txt checks every rounded weight against
the exact rational), and `kverify` does not cross-check the weight types stored in the
table header against `wtypes_*.txt` (they agree where checked).

These computations have not yet been repeated by anyone else.

## Reproducing

```
cd kary
clang -O3 -march=native -o kary kary.c -lm
clang -O2 -ffp-model=strict -march=native -o kverify kverify.c -lm
mkdir -p data
./kary build 3 6 data                 # 2.8 GB table, under a minute on 10 cores
cp vectors/kvec_k3_h6.bin data/       # or: ./kary scan 3 6 data 0.19 0.23
python3 kcertify.py 3 6 data          # needs mpmath; CERTIFIED: gamma_3 <= 0.739975
./kary build 4 5 data                 # 2.2 GB table
cp vectors/kvec_k4_h5.bin data/       # or: ./kary scan 4 5 data 0.14 0.165
python3 kcertify.py 4 5 data          # CERTIFIED: gamma_4 <= 0.683622
```

The tables are deterministic; `certificates/SHA256_kary.txt` lists the checksums of the
two headline tables, vectors and weight files. Any positive vector is a valid certificate,
so the stored vectors only save the power iteration.

I use a large language model to help write and check this code.
