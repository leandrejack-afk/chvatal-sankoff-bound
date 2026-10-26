# Lower bounds on the average edit distance between two random strings

For independent uniform random strings X, Y of length n over k letters, the expected edit distance
(unit-cost insertions, deletions, substitutions) is α_k·n + o(n). Certified lower bounds from this
folder, against the best previously known:

| k | previous best | source | new lower bound | estimate of α_k (BS26) |
|---|---|---|---|---|
| 2 | 0.178037 | 1 − γ₂, γ₂ ≤ 0.821963 (Tomar 2026, preprint) | **0.216276** | 0.27496 |
| 3 | 0.283660 | Bilardi–Schimd 2026 | **0.340971** | 0.41532 |
| 4 | 0.359783 | Bilardi–Schimd 2026 | **0.424710** | 0.50419 |
| 5 | 0.415173 | Bilardi–Schimd 2026 | **0.485876** | 0.56678 |
| 6 | 0.457766 | Bilardi–Schimd 2026 | **0.533057** | 0.61390 |
| 7 | 0.491836 | Bilardi–Schimd 2026 | **0.570869** | 0.65089 |
| 8 | 0.519901 | Bilardi–Schimd 2026 | **0.602039** | 0.68082 |
| 16 | 0.644758 | Bilardi–Schimd 2026 | **0.740007** | 0.80586 |
| 32 | 0.738677 | Bilardi–Schimd 2026 | **0.839930** | 0.88619 |

Bilardi–Schimd 2026: G. Bilardi, M. Schimd, Computable bounds and Monte Carlo estimates of the
expected edit distance, Information 17 (2026) 685; the previous bounds are from their Table 1 and
the estimates from their Table 6. For k = 2 the best peer-reviewed bound is 1 − 0.826280 = 0.17372
from Lueker's γ₂ ≤ 0.826280; the 0.178037 uses T. S. Tomar's preprint bound γ₂ ≤ 0.821963
(doi:10.5281/zenodo.23610945).

The first version of this folder (commit f834334, rules 1 to 4 below) gave 0.213213, 0.335199,
0.417664, 0.478381, 0.525535, 0.563528, 0.594981, 0.735355 and 0.837606. Rule 5, added the same
day, raises each of these by 0.0023 to 0.0075.

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
   and likewise no I M D;
5. the run of insertions that directly follows a substitution on (x_i, y_j) uses no letter equal to
   x_i, and the run of deletions that directly follows it uses no letter equal to y_j (if
   y_{j+s} = x_i, then S I^r can be replaced by I^s M I^(r−s), which costs one less; likewise for
   deletions, so this holds in every optimal alignment).

Rewrites 2 and 3 keep the number of insertions and deletions and move a diagonal move (M or S)
earlier past one of them; rewrite 4 removes two of them. So the pair (number of insertions and
deletions, sum over diagonal moves of the number of insertions and deletions before it) decreases
lexicographically with every rewrite, and rewriting stops. Costs never increase, and an optimal
alignment cannot get cheaper, so the final alignment is optimal and satisfies 1 to 4, and, being
optimal, also 5. `canon_check5.c` (rules 1 to 5), `canon_check3.c` (rules 1 to 4) and
`canon_check2.c` (rules 1 to 3) confirm exhaustively that the minimum cost over alignments obeying
the rules equals the edit distance for every pair with (k, n, m) = (2, 10, 10), (2, 8, 12), (3, 7, 7),
(3, 6, 8), (4, 5, 6), (5, 5, 5), (6, 4, 5), about 36 million pairs; for rules 1 to 4 and 1 to 5 also
(2, 12, 12) and (3, 8, 8), another 60 million pairs (`canon_check5.log`).

**First moment.** P(ed(X, Y) ≤ an) ≤ E[number of alignments obeying 1 to 5, valid for (X, Y), with
cost ≤ an]. An alignment is valid if its M moves join equal letters and its S moves unequal ones.
Each letter is used by one move; besides its own move, a run letter is constrained only by the
match that ends its run (rule 3) and, if the run directly follows a substitution, by the substituted
letter (rule 5). So the probability is a product over moves and runs: 1/k for M, (k − 1)/k for S,
and

- for a run of r indels after a match (or at the start) together with the match that ends it,
  (1/k)((k − 1)/k)^r;
- for a run of r indels directly after a substitution together with the match that ends it,
  ((k − 1)/k)^r · (1/k)·h(r) with h(r) = 1/k + ((k − 1)/k)((k − 2)/(k − 1))^r: each run letter avoids
  the substituted letter, and the match must also avoid the run letters; if it equals the
  substituted letter (probability 1/k) that holds already, and otherwise (probability (k − 1)/k)
  each run letter has to avoid one more letter, (k − 2)/(k − 1) each;
- for a run at the end of the strings, 1, or ((k − 1)/k)^r directly after a substitution.

With x marking a letter of either string and t a unit of cost, for any x > 0 and 0 < t ≤ 1,

  E[...] ≤ x^(−2n) t^(−an) G(x, t),

where G(x, t) sums these weights, with M = x²/k, S = x²t(k − 1)/k and I = D = xt, over all move
sequences accepted by an automaton with states for a diagonal move (one after a match or at the
start, one after a substitution), MD and MI (a match right after a deletion run or an insertion run,
after which the opposite indel is forbidden by rule 4), and the run states, by the current run
length capped at R, kept separately for runs directly after a substitution. Capping drops some
constraints, since ((k − 1)/k)^r and h(r) never increase with r, which only increases the count.
`count_check5.py` (rules 1 to 5), `count_check2.py` (rules 1 to 4) and `count_check.py` (rules 1 to
3) confirm that these weights give exactly the brute-force expected counts, as polynomials in the
cost, for small n. G is finite when the spectral radius of the automaton's weight matrix A(x, t) is
below 1. So if ρ(A) < 1 and x²t^a > 1, then P(ed ≤ an) → 0 exponentially, and α_k ≥ a.

**Certificate.** For each k (`certify_rule5.py`, R = 32): rational x and t near the optimum, and a
rational positive vector v with (Av)_i < v_i for every state i, checked in exact arithmetic (this
proves ρ(A) < 1 for A ≥ 0); then a = −2 ln x / ln t evaluated in interval arithmetic and rounded
down. Output: `certify_rule5_R32.log`. Earlier stages: rules 1 to 4, `certify_runs4.py` and
`certify_runs4_R32.log` (the first version's values above); rules 1 to 3, `certify_runs.py` and
`certify_runs_R32.log` (0.209789, 0.328355, 0.408774, 0.468384, 0.515005, 0.552816, 0.584300,
0.726970, 0.832820).

With rules 1 and 2 only (a 3-state automaton, `certify_simple.py`) the bounds are
already 0.295355 (k = 3), 0.380300 (k = 4), 0.567660 (k = 8), 0.829531 (k = 32), above Bilardi and
Schimd's for every k ≥ 3.

**Earlier work.** Counting canonical edit scripts is not new for this constant. Bilardi and Schimd
(2026, Section 6.2) count "canonical simple scripts", with all deletions first, then substitutions,
then insertions, and bound the number of positions with binomial coefficients; Ganguly, Mossel and
Rácz (arXiv:1601.07086, Appendix A) used a similar volume argument. The idea of counting only
alignments in a canonical form also underlies Dančík and Paterson's upper bounds on the
Chvátal–Sankoff constants. What is new here is the choice of local rules (adjacency of insertions,
deletions and substitutions, and the letters of a run, against the match that ends it and the
substitution before it) and their exact count by a small transfer matrix.

**A check that needs no matrix.** Eliminating the run states reduces ρ(A) < 1 to one scalar
inequality. With R = ∞, M = x²/k, S = x²t(k − 1)/k, a = M + S, u = xt(k − 1)/k, w = xt(k − 2)/k,
b = M·u/(1 − u) and c = M·((1/k)·u/(1 − u) + ((k − 1)/k)·w/(1 − w)), it is a + 2b < 1 for rules
1 to 3, a + b + ab < 1 for rules 1 to 4, and a + b + Mb + S(2c − b) < 1 together with xt < 1 and
a + b < 2 for rules 1 to 5 (with c = b this is the rules 1 to 4 test). `review/verify_rule4.py`
(`r3` and `r4` modes) and `review5/cert_check5.py` check every certified (x, t) in the logs this
way, in exact fractions (the latter with R = 32, like the certificate), and recompute a* to 18 and
90 digits.

**Review.** A fresh-context review re-derived the argument and the certificates independently
(`review/`): random rewriting with adversarial tie-breaks (about 70 million rewrites), its own
exact counts in 34 configurations, Monte Carlo at n = 2000 (no sample at or below the bounds), and
a literature search. Verdict: sound, with edits to this README, which are made. A second
fresh-context review checked rule 5 (`review5/`): it tested the rule in every optimal alignment
directly (exhaustively at 13 sizes and on about 175,000 random pairs up to n = 150), found a
rule-obeying optimal alignment for every pair with its own search, enumerated the probabilities of
the new blocks and the independence between blocks on 1,834 move sequences, rebuilt the matrix from
its own transition table (equal entry by entry), certified it with its own vectors and the scalar
test, and ran Monte Carlo (at n = 4000, mean ed/n 0.29087 for k = 2 and 0.52067 for k = 4).
Verdict: all nine bounds hold, with edits to this README, which are made.
