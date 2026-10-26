"""Rigorous lower bounds on the edit-distance constant alpha_k (run-aware canonical alignments,
rules 1-5).

alpha_k = lim E[ed(X, Y)]/n for independent uniform X, Y in [k]^n (unit-cost Levenshtein).

Canonical form: rules (1)-(4) as in certify_runs4.py (every pair has an optimal alignment obeying
them), plus
  (5) after a substitution S on (a, b), the following run of insertions avoids the letter a and the
      following run of deletions avoids b (if an inserted y_{j+s} equals a, then I^s M I^(r-s) in
      place of S I^r costs one less; likewise for deletions). This holds in every optimal alignment,
      so the optimal alignment obeying (1)-(4) also obeys (5). (Exhaustively confirmed by
      canon_check5.c.)

First moment: as in certify_runs4.py, with the automaton states G0 (diagonal after M, or start),
G1 (diagonal after S), MD, MI, IR_r and DR_r (runs after G0, MD or MI) and IS_r and DS_r (runs right
after S), r capped at R, and weights
  M = x^2/k, times q^r when it ends a run IR_r or DR_r, and times h(r) = 1/k + q (p/q)^r when it ends a
      run IS_r or DS_r;
  S = x^2 t q;  I = D = x t, times q for every letter of a run IS or DS;
with q = (k-1)/k, p = (k-2)/k. A letter of a run after S avoids the S letter (probability q); the
match ending the run must also avoid the run letters: when it equals the S letter (probability 1/k)
that holds already, otherwise (probability q) each run letter must avoid one more letter (p/q each).
These are exact (count_check5.py compares them with brute-force counts); capping r at R uses
q^R >= q^r and h(R) >= h(r) for r >= R, which only increases the count.

Certificate: rational x, t and a rational vector v > 0 with A v < v componentwise (exact
arithmetic; this implies rho(A) < 1 for A >= 0), then a = -2 ln x / ln t in interval arithmetic,
rounded down to 6 decimals.
usage: certify_rule5.py R k [k ...]
"""
import sys
from fractions import Fraction as Fr

import mpmath
import numpy as np
from scipy.optimize import brentq, minimize_scalar


def build(k, x, t, R, exact=False):
    """States: 0 G0 (diagonal after M, or start), 1 G1 (diagonal after S), 2 MD, 3 MI, then the run
    states IR_1..IR_R and DR_1..DR_R (runs after G0, MD or MI) and IS_1..IS_R and DS_1..DS_R (runs
    right after S), indexed by run length capped at R."""
    if exact:
        M, S, I, q, p = x * x / k, x * x * t * Fr(k - 1, k), x * t, Fr(k - 1, k), Fr(k - 2, k)
        n = 4 + 4 * R
        A = [[Fr(0)] * n for _ in range(n)]
    else:
        M, S, I, q, p = x * x / k, x * x * t * (k - 1) / k, x * t, (k - 1) / k, (k - 2) / k
        A = np.zeros((4 + 4 * R, 4 + 4 * R))
    h = lambda r: 1 / k + q * (p / q) ** r if not exact else Fr(1, k) + q * (p / q) ** r
    G0, G1, MD, MI = 0, 1, 2, 3
    IR, DR, IS, DS = 4, 4 + R, 4 + 2 * R, 4 + 3 * R  # state of run length r is base + r - 1
    for s in (G0, G1, MD, MI):
        A[s][G0] += M
        A[s][G1] += S
        if s != MD:  # no insertion right after (deletion run, match)
            A[s][IS if s == G1 else IR] += I * q if s == G1 else I
        if s != MI:  # no deletion right after (insertion run, match)
            A[s][DS if s == G1 else DR] += I * q if s == G1 else I
    for r in range(1, R + 1):
        nxt = min(r + 1, R) - 1
        A[IR + r - 1][MI] += M * q ** r
        A[IR + r - 1][IR + nxt] += I
        A[DR + r - 1][MD] += M * q ** r
        A[DR + r - 1][DR + nxt] += I
        A[IS + r - 1][MI] += M * h(r)
        A[IS + r - 1][IS + nxt] += I * q
        A[DS + r - 1][MD] += M * h(r)
        A[DS + r - 1][DS + nxt] += I * q
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
              f"a >= {mpmath.nstr(mpmath.floor(lo * 10 ** 12) / 10 ** 12, 12, strip_zeros=False)}; "
              f"alpha_{k} >= {mpmath.nstr(a6, 6, strip_zeros=False)}  (float {af:.9f})",
              flush=True)
