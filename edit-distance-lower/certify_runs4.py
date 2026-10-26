"""Rigorous lower bounds on the edit-distance constant alpha_k (run-aware canonical alignments,
rules 1-4).

alpha_k = lim E[ed(X, Y)]/n for independent uniform X, Y in [k]^n (unit-cost Levenshtein).

Canonical form. Every pair (X, Y) has an optimal alignment (moves M: match of equal letters,
cost 0; S: substitution of unequal letters; I: insertion; D: deletion; cost 1 each) such that
  (1) no I is adjacent to a D (one diagonal move would be cheaper);
  (2) no S directly follows an I or a D (I S -> S I and D S -> S D keep the cost; if the
      swapped diagonal pair had equal letters the alignment was not optimal);
  (3) the match after a maximal run of r insertions (deletions) uses a letter different from all r
      letters of the run (if it equals the letter at offset s, match it there instead:
      I^r M -> I^s M I^(r-s) keeps the cost);
  (4) no D M I and no I M D (a single match between a deletion and an insertion): D M I on
      x_i, (x_{i+1}, y_j), y_{j+1} can be replaced by the two diagonal moves (x_i, y_j),
      (x_{i+1}, y_{j+1}) at cost at most 2.
Rewrites (2) and (3) keep the number of indels and move a diagonal move earlier; rewrite (4) removes
two indels; so the pair (number of indels, sum over diagonal moves of the indels before it) drops
lexicographically, and the rewriting stops.
Rule (1) holds in every optimal alignment; rewrites (2)-(4) never raise the cost, so the result is
still optimal and satisfies (1)-(4). (Exhaustively confirmed on about 96 million pairs by
canon_check3.c.)

First moment. P(ed <= a n) <= E[#canonical alignments of cost <= a n valid for (X, Y)]
<= x^(-2n) t^(-a n) G(x, t) for x > 0, 0 < t <= 1, where G sums the weights of all move sequences
accepted by the automaton with states diag, MD, MI (a match right after a deletion run or an
insertion run), I_1..I_R, D_1..D_R (run length, capped at R):
  M = x^2/k (times ((k-1)/k)^r when it ends a run of length r, r capped at R), S = x^2 t (k-1)/k,
  I = D = x t,
with transitions diag -> {M, S, I, D}, MD -> {M, S, D}, MI -> {M, S, I}, I_r -> {M -> MI, I},
D_r -> {M -> MD, D}. (The weights reproduce the brute-force expected counts exactly for small n:
count_check2.py.) Capping at R drops some of the
constraints, which only increases the count. G is finite when rho(A) < 1. So if rho(A) < 1 and
x^2 t^a > 1 then P(ed <= a n) -> 0, hence alpha_k >= a.

Certificate: rational x, t and a rational vector v > 0 with A v < v componentwise (exact
arithmetic; this implies rho(A) < 1 for A >= 0), then a = -2 ln x / ln t in interval arithmetic,
rounded down to 6 decimals.
usage: certify_runs4.py R k [k ...]
"""
import sys
from fractions import Fraction as Fr

import mpmath
import numpy as np
from scipy.optimize import brentq, minimize_scalar


def build(k, x, t, R, exact=False):
    """States: 0 diag, 1 MD, 2 MI, 3..2+R the I-run states I_1..I_R, 3+R..2+2R the D-run states."""
    if exact:
        M, S, I, q = x * x / k, x * x * t * Fr(k - 1, k), x * t, Fr(k - 1, k)
        A = [[Fr(0)] * (3 + 2 * R) for _ in range(3 + 2 * R)]
    else:
        M, S, I, q = x * x / k, x * x * t * (k - 1) / k, x * t, (k - 1) / k
        A = np.zeros((3 + 2 * R, 3 + 2 * R))
    DG, MD, MI, i0, d0 = 0, 1, 2, 3, 3 + R
    for s in (DG, MD, MI):
        A[s][DG] = M + S
        if s != MD:
            A[s][i0] = I  # no insertion right after (deletion run, match)
        if s != MI:
            A[s][d0] = I  # no deletion right after (insertion run, match)
    for r in range(1, R + 1):
        A[i0 + r - 1][MI] = M * q ** r
        A[i0 + r - 1][i0 + min(r + 1, R) - 1] = I
        A[d0 + r - 1][MD] = M * q ** r
        A[d0 + r - 1][d0 + min(r + 1, R) - 1] = I
    return A


def rho(k, x, t, R):
    return max(abs(np.linalg.eigvals(build(k, x, t, R))))


def optimum(k, R):
    def x_of_t(t):
        return brentq(lambda x: rho(k, x, t, R) - 1, 1e-9, 64.0)

    def neg_a(t):
        x = x_of_t(t)
        return 2 * np.log(x) / np.log(t) if x > 1 else 0.0

    grid = np.linspace(0.005, 0.995, 199)
    vals = [neg_a(t) for t in grid]
    i = int(np.argmin(vals))
    r = minimize_scalar(neg_a, bounds=(grid[max(i - 1, 0)], grid[min(i + 1, 198)]),
                        method="bounded", options={"xatol": 1e-13})
    return x_of_t(r.x), r.x, -r.fun


def certify(k, R):
    xf, tf, af = optimum(k, R)
    x = Fr(xf * (1 - 1e-9)).limit_denominator(10 ** 12)
    t = Fr(tf).limit_denominator(10 ** 12)
    Af = build(k, float(x), float(t), R)
    w, V = np.linalg.eig(Af)
    v = np.abs(V[:, int(np.argmax(np.abs(w)))].real)
    v = v / v.max()
    vq = [Fr(float(c)).limit_denominator(10 ** 15) + Fr(1, 10 ** 18) for c in v]
    A = build(k, x, t, R, exact=True)
    assert all(a >= 0 for row in A for a in row) and all(c > 0 for c in vq)
    Av = [sum(A[i][j] * vq[j] for j in range(len(vq))) for i in range(len(vq))]
    worst = max(Av[i] / vq[i] for i in range(len(vq)))
    if not worst < 1:
        raise SystemExit(f"k={k}: A v < v not certified (max ratio {float(worst)})")
    iv = mpmath.iv
    iv.dps = 50
    X = iv.mpf(x.numerator) / x.denominator
    T = iv.mpf(t.numerator) / t.denominator
    a = -2 * iv.log(X) / iv.log(T)
    mpmath.mp.dps = 60
    lo = mpmath.mpf(a.a)  # lower end of the interval, at full precision
    return x, t, af, float(worst), lo, mpmath.floor(lo * 10 ** 6) / 10 ** 6


if __name__ == "__main__":
    R = int(sys.argv[1])
    for k in map(int, sys.argv[2:]):
        x, t, af, worst, lo, a6 = certify(k, R)
        print(f"k={k} R={R}: x={x} t={t} max (Av)_i/v_i = {worst:.12f} < 1; "
              f"a >= {mpmath.nstr(mpmath.floor(lo * 10 ** 12) / 10 ** 12, 12)}; alpha_{k} >= {mpmath.nstr(a6, 6)}  (float {af:.9f})",
              flush=True)
