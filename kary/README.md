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
| 3 | 6 | 1,711,682 | 525 | 0.2048 | 0.268056913 | 0.739975 |
| 3 | 7 | 28,966,771 | 975 | 0.2129 | 0.296474627 | **0.737495** |
| 4 | 2 | 10 | 52 | 0.1147 | 0.110575331 | 0.707648 |
| 4 | 3 | 309 | 292 | 0.1285 | 0.148063383 | 0.696851 |
| 4 | 4 | 9,427 | 1,108 | 0.1402 | 0.186602820 | 0.689064 |
| 4 | 5 | 290,442 | 3,460 | 0.1496 | 0.222621807 | 0.683622 |
| 4 | 6 | 8,502,177 | 9,748 | 0.1570 | 0.254391027 | **0.679701** |

so **γ_3 < 0.737495** and **γ_4 < 0.679701**. The rows k = 3, h = 7 and k = 4, h = 6 (added
later on 25 October 2026) need tables of 76 GB and 169 GB, which are built and used from
disk; see [Larger automata](#larger-automata-tables-on-disk). The best upper bounds before
were γ_3 < 0.765803 and γ_4 < 0.708236 (V. Dančík, *Expected length of longest common
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
   paper and the summary of results (the table above and the weight formula), before
   reading the code above: its own state encoding, transition function, class grouping
   and weights (`rvk.c`, `rv_*.py`). It reproduces every state count,
   accepted-transition count and bound in the table, matches both headline tables state
   by state, re-runs the certificates, and gives exact rational certificates for h = 2
   (`review/log_h2cert.txt`: rational z, λ and v with A_z v ≤ λ v checked in exact
   arithmetic, giving γ_3 ≤ 0.764972 and γ_4 ≤ 0.707648).

Known limits, closed by the review checks: `kweights.py` rounds an 80-digit value rather
than the exact rational (review/log_theirweights.txt checks every rounded weight against
the exact rational), and `kverify` does not cross-check the weight types stored in the
table header against `wtypes_*.txt` (they agree where checked).

All of this code, the review implementation included, was written with the help of the same
large language model, and these computations have not yet been repeated by anyone else.
Logs of the checks rerun on 25 October 2026 are in `logs/rerun-2026-10-25/`.

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

## Larger automata: tables on disk

For three letters with h = 7 and four letters with h = 6 the tables are 75.9 GB and 169.3 GB,
far more than the 24 GB of memory of the laptop used, so they are built and used from disk.

- `kbig.c` builds the same table as `kary.c`. Its transition function and hash table are
  copied from kary.c unchanged; each worker sorts and merges the transitions of its state as
  kary.c does and looks the targets up read-only, and the keys not yet known are then inserted
  one at a time in kary.c's order, so every state number is kary.c's. The records go to disk
  as they are produced (`kbig_k<k>_h<h>.bin`: header, records, offsets, keys); `kbig klue`
  rewrites a table in kary.c's format.
- `kbig eigen` streams the records from disk and iterates eight values of z in each pass.
  `kbig-runs/pipeline.sh` runs the eigen phase: a start vector projected from the h − 1 vector,
  eight values of z 0.0004 apart, a parabola through the three best bounds, eight values
  0.0001 apart around its vertex (started from the interpolated coarse vectors), and the
  certificate for the best one. The vector only has to be good: any positive vector gives a
  valid certificate.
- `kbverify.c` is kverify.c reading the table in chunks of whole states; it also checks that
  the weight types in the table header are those the weights were computed for.
  `kbcertify.py` is kcertify.py with kbverify, and compares the printed six-decimal bound with
  the exact upper end of the interval in rational arithmetic.

| k | h | states | accepted transitions | records | z | λ_cert | ρ_lower | exact value | bound |
|---|---|---|---|---|---|---|---|---|---|
| 3 | 7 | 28,966,771 | 9,520,016,086 | 9,427,252,448 | 0.2129 | 0x1.2f970b6938e14p-2 | 0x1.2f970b68f84efp-2 | 0.7374942650 | **0.737495** |
| 4 | 6 | 8,502,177 | 27,720,529,578 | 21,142,620,884 | 0.1570 | 0x1.047f14d17bdedp-2 | 0x1.047f14cf05242p-2 | 0.6797004676 | **0.679701** |

ρ_lower is a lower bound on the spectral radius, so no vector can certify less than
0.73749426499 and 0.679700467298 at these values of z. The certificate outputs are
`certificates/certify_k3_h7.out` and `certificates/certify_k4_h6.out`.

Checks (logs in `logs/kbig-2026-10-25/`):

1. kbig rebuilds all 16 earlier tables (k = 2, h ≤ 8; k = 3, h ≤ 6; k = 4, h ≤ 5) byte for
   byte, with the SHA-256 sums recorded for (3, 6) and (4, 5) (`repro_identity.log`), and
   kbverify gives kverify's λ_cert and ρ_lower bit for bit on the ten earlier certificates
   (`kbverify_vs_kverify.log`).
2. Every state of both new tables was recomputed with the transition function of the earlier
   review implementation (`review/rvcore.h`, through `rvcmp_big.c`, in four quarter runs per
   table) and compared with the stored records: all 28,966,771 states and 9,520,016,086
   transitions at (3, 7), all 8,502,177 states and 27,720,529,578 transitions at (4, 6),
   no difference (`rvcmp_k*_q*.log`).
3. `check_table_big.py` recomputes sampled states with `kref.py`: 5,400 states at (3, 7)
   and 2,400 at (4, 6), 9,273,303 transitions, no mismatch.
4. `kcover.py`: 100,000 random string pairs at each of (3, 7) and (4, 6), 2,191,009
   canonical steps, 142,788 of them truncated inputs: never rejected.
5. `review/rv_theirweights.py`: every rounded weight brackets the exact rational weight at
   the exact double z, as adjacent doubles (24 of 24 and 28 of 28).
6. A second review pass with fresh context (`review-big/REVIEW.md`) read the new code against
   kary.c and the paper, attacked the streaming and the rounding (corrupted tables and
   vectors, thousands of chunk boundaries), and checked both certificates with its own
   verifier (`review-big/rvb_verify.c`: rounding emulated with error-free transformations,
   different streaming), which gives the same λ_cert bit for bit, and with an a priori error
   bound that needs no rounding mode at all; both give the same six decimals. Like the first
   review pass, it is part of my own checking with the same language model, not an independent
   replay.

Reproducing, from `kary/` (on a 10-core laptop shared with other work: about 3 hours of
building and 2 of eigen phase; 260 GB of free disk):

```
clang -O3 -march=native -o kary kary.c -lm
clang -O3 -march=native -o kbig kbig.c -lm
clang -O2 -ffp-model=strict -march=native -o kbverify kbverify.c -lm
mkdir -p data big/k3h7 big/k4h6
./kary build 3 6 data && cp vectors/kvec_k3_h6.bin data/   # the h - 1 table and vector (start vector)
./kbig build 3 7 big/k3h7
kbig-runs/pipeline.sh 3 7 0.2126                         # CERTIFIED: gamma_3 <= 0.737495
./kary build 4 5 data && cp vectors/kvec_k4_h5.bin data/
./kbig build 4 6 big/k4h6
TOL1=1e-7 TOL2=1e-9 kbig-runs/pipeline.sh 4 6 0.1570     # CERTIFIED: gamma_4 <= 0.679701
clang -O3 -ffp-model=strict -o rvcmp_big rvcmp_big.c -lpthread
kbig-runs/compare_full.sh 3 7; kbig-runs/compare_full.sh 4 6   # check 2, logs in big/logs/
```

The certified vectors `kvec_k3_h7.bin` and `kvec_k4_h6.bin` (232 MB and 68 MB) are too large for
the repository and are attached to the release `ipl-2026`; `certificates/SHA256_kbig.txt`
has the SHA-256 sums of the two tables, the two vectors and the weight files. The eigen phase
is deterministic for a given binary, and the exact values lie more than 5 × 10⁻⁷ below the
printed bounds, so a vector computed with another compiler gives the same six decimals.

I use a large language model to help write and check this code.
