# Lower bounds on the average edit distance between two random strings

For independent uniform random strings X, Y of length n over k letters, the expected edit distance
(unit-cost insertions, deletions, substitutions) is α_k·n + o(n). Certified lower bounds from this
folder, against the best previously known:

| k | previous best | source | new lower bound | estimate of α_k (BS26) |
|---|---|---|---|---|
| 2 | 0.178037 | 1 − γ₂, γ₂ ≤ 0.821963 (Tomar 2026, preprint) | **0.213213** | 0.27496 |
| 3 | 0.283660 | Bilardi–Schimd 2026 | **0.335199** | 0.41532 |
| 4 | 0.359783 | Bilardi–Schimd 2026 | **0.417664** | 0.50419 |
| 5 | 0.415173 | Bilardi–Schimd 2026 | **0.478381** | 0.56678 |
| 6 | 0.457766 | Bilardi–Schimd 2026 | **0.525535** | 0.61390 |
| 7 | 0.491836 | Bilardi–Schimd 2026 | **0.563528** | 0.65089 |
| 8 | 0.519901 | Bilardi–Schimd 2026 | **0.594981** | 0.68082 |
| 16 | 0.644758 | Bilardi–Schimd 2026 | **0.735355** | 0.80586 |
| 32 | 0.738677 | Bilardi–Schimd 2026 | **0.837606** | 0.88619 |

Bilardi–Schimd 2026: G. Bilardi, M. Schimd, Computable bounds and Monte Carlo estimates of the
expected edit distance, Information 17 (2026) 685; the previous bounds are from their Table 1 and
the estimates from their Table 6. For k = 2 the best peer-reviewed bound is 1 − 0.826280 = 0.17372
from Lueker's γ₂ ≤ 0.826280; the 0.178037 uses T. S. Tomar's preprint bound γ₂ ≤ 0.821963
(doi:10.5281/zenodo.23610945).

## The argument

**Canonical optimal alignments.** An alignment of X and Y is a sequence of moves: M (a match of two
equal letters, cost 0), S (a substitution of two unequal letters, cost 1), I (an insertion, cost 1),
D (a deletion, cost 1). Every pair (X, Y) has an optimal alignment in which

1. no I is next to a D (a pair I D or D I can be replaced by one diagonal move, saving at least 1, so
   this holds in every optimal alignment);
2. no S comes directly after an I or a D (I S can be replaced by S I and D S by S D at the same cost;
   if the swapped diagonal pair had equal letters it would be a match and the alignment was not
   optimal);
3. the match that ends a maximal run of r insertions (or deletions) uses a letter different from all
   r letters of the run (if it equals the letter at offset s of the run, match it there instead:
   I^r M becomes I^s M I^(r−s) at the same cost);
4. no single match sits directly between a deletion and an insertion (D M I on x_i, (x_{i+1}, y_j),
   y_{j+1} can be replaced by the diagonal moves (x_i, y_j), (x_{i+1}, y_{j+1}), which cost at most 2),
   and likewise no I M D.

Rewrites 2 and 3 keep the number of insertions and deletions and move a diagonal move (M or S)
earlier past one of them; rewrite 4 removes two of them. So the pair (number of insertions and
deletions, sum over diagonal moves of the number of insertions and deletions before it) decreases
lexicographically with every rewrite, and rewriting stops. Costs never increase, and an optimal
alignment cannot get cheaper, so the final alignment is optimal and satisfies 1 to 4.
`canon_check3.c` (rules 1 to 4) and `canon_check2.c` (rules 1 to 3) confirm exhaustively that the
minimum cost over alignments obeying the rules equals the edit distance for every pair with
(k, n, m) = (2, 10, 10), (2, 8, 12), (3, 7, 7), (3, 6, 8), (4, 5, 6), (5, 5, 5), (6, 4, 5), about
36 million pairs; for rules 1 to 4 also (2, 12, 12) and (3, 8, 8), another 60 million pairs.

**First moment.** P(ed(X, Y) ≤ an) ≤ E[number of alignments obeying 1 to 4, valid for (X, Y), with
cost ≤ an]. An alignment is valid if its M moves join equal letters and its S moves unequal ones;
since every letter is used by at most one move, the probability is a product over moves: 1/k for M,
(k − 1)/k for S, and for the match ending a run of r indels, (1/k)((k − 1)/k)^r (it must also avoid
the r letters of the run, which no other move constrains). With x marking a letter of either
string and t a unit of cost, for any x > 0 and 0 < t ≤ 1,

  E[...] ≤ x^(−2n) t^(−an) G(x, t),

where G(x, t) sums the weights M = x²/k (times ((k − 1)/k)^r after a run of r), S = x²t(k − 1)/k,
I = D = xt over all move sequences accepted by the automaton with states diag, MD and MI (a match
right after a deletion run or an insertion run, after which the opposite indel is forbidden by
rule 4), I_1..I_R and D_1..D_R (the current run length, capped at R; ignoring constraints beyond R
only increases the count). `count_check2.py` (and `count_check.py` for rules 1 to 3) confirm that
these weights give exactly the brute-force expected counts, as polynomials in the cost, for small n. G is finite when the spectral radius of the automaton's
weight matrix A(x, t) is below 1. So if ρ(A) < 1 and x²t^a > 1, then P(ed ≤ an) → 0 exponentially,
and α_k ≥ a.

**Certificate.** For each k (`certify_runs4.py`, R = 32): rational x and t near the optimum, and a
rational positive vector v with (Av)_i < v_i for every state i, checked in exact arithmetic (this
proves ρ(A) < 1 for A ≥ 0); then a = −2 ln x / ln t evaluated in interval arithmetic and rounded
down. Output: `certify_runs4_R32.log` (rules 1 to 3: `certify_runs.py`, `certify_runs_R32.log`, with
0.209789, 0.328355, 0.408774, 0.468384, 0.515005, 0.552816, 0.584300, 0.726970, 0.832820).

With rules 1 and 2 only (a 3-state automaton, `certify_simple.py`) the bounds are
already 0.295355 (k = 3), 0.380300 (k = 4), 0.567660 (k = 8), 0.829531 (k = 32), above Bilardi and
Schimd's for every k ≥ 3.

**Earlier work.** Counting canonical edit scripts is not new for this constant. Bilardi and Schimd
(2026, Section 6.2) count "canonical simple scripts", with all deletions first, then substitutions,
then insertions, and bound the number of positions with binomial coefficients; Ganguly, Mossel and
Rácz (arXiv:1601.07086, Appendix A) used a similar volume argument. The idea of counting only
alignments in a canonical form also underlies Dančík and Paterson's upper bounds on the
Chvátal–Sankoff constants. What is new here is the choice of local rules (adjacency of insertions,
deletions and substitutions, and the letters at the end of a run) and their exact count by a
small transfer matrix.

**A check that needs no matrix.** Eliminating the run states reduces ρ(A) < 1 to one scalar
inequality: with a = x²/k + x²t(k − 1)/k, u = xt(k − 1)/k and b = (x²/k)·u/(1 − u), it is a + 2b < 1
for rules 1 to 3 and a + b + ab < 1 for rules 1 to 4 (with R = ∞). `review/verify_rule4.py`
(`r3` and `r4` modes) checks every certified (x, t) in the logs this way, in exact fractions, and
recomputes a* to 18 digits.

**Review.** A fresh-context review re-derived the argument and the certificates independently
(`review/`): random rewriting with adversarial tie-breaks (about 70 million rewrites), its own
exact counts in 34 configurations, Monte Carlo at n = 2000 (no sample at or below the bounds), and
a literature search. Verdict: sound, with edits to this README, which are made.
