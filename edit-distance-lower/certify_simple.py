"""Rigorous lower bounds on the edit-distance constant alpha_k from canonical alignments.

alpha_k = lim E[ed(X, Y)]/n for independent uniform X, Y in [k]^n, ed = unit-cost Levenshtein.

Canonical alignments. Every pair (X, Y) has an optimal alignment (moves M: match of equal letters,
cost 0; S: substitution of unequal letters, cost 1; I: insertion; D: deletion, cost 1 each) with no
adjacent I-D or D-I pair (replacing the pair by one diagonal move saves at least 1) and no S directly
after an I or a D (I S -> S I and D S -> S D keep the cost: if the swapped diagonal pair had equal
letters it would be a match and the alignment would not be optimal). Repeating the swaps moves every
S earlier, so it terminates.

First moment. P(ed <= a n) <= E[#canonical alignments of cost <= a n valid for (X, Y)]
<= x^(-2n) t^(-a n) G(x, t) for x > 0, 0 < t <= 1, where G sums over all move sequences accepted
by the 3-state automaton (last move: diagonal, I, D) the weights M = x^2/k, S = x^2 t (k-1)/k,
I = D = x t. G is finite when the spectral radius of
    A = [[M+S, I, D], [M, I, 0], [M, 0, D]]
is below 1. So if rho(A) < 1 and x^2 t^a > 1 then P(ed <= a n) -> 0 and alpha_k >= a.

Certificate for each k: rational x, t; rho(A) < 1 is checked exactly (I - A is a nonsingular
M-matrix iff its leading principal minors are positive); then the bound a = -2 ln x / ln t is
evaluated in interval arithmetic and rounded down.
usage: certify_simple.py k [k ...]
"""
import sys
from fractions import Fraction as Fr

import mpmath
import numpy as np
from scipy.optimize import brentq, minimize_scalar


def matrix(k, x, t):
    M, S, I = x * x / k, x * x * t * (k - 1) / k, x * t
    return [[M + S, I, I], [M, I, 0 * I], [M, 0 * I, I]]


def rho_float(k, x, t):
    return max(abs(np.linalg.eigvals(np.array(matrix(k, x, t), dtype=float))))


def optimum(k):
    """Float search: maximize a(t) = -2 ln x(t) / ln t where rho(x(t), t) = 1."""
    def x_of_t(t):
        f = lambda x: rho_float(k, x, t) - 1
        return brentq(f, 1e-9, 64.0)
    def neg_a(t):
        x = x_of_t(t)
        return 2 * np.log(x) / np.log(t) if x > 1 else 0.0
    grid = np.linspace(0.005, 0.995, 199)
    vals = [neg_a(t) for t in grid]
    i = int(np.argmin(vals))
    lo_t, hi_t = grid[max(i - 1, 0)], grid[min(i + 1, len(grid) - 1)]
    r = minimize_scalar(neg_a, bounds=(lo_t, hi_t), method="bounded", options={"xatol": 1e-13})
    t = r.x
    return x_of_t(t), t, -r.fun


def m_matrix_ok(A):
    """I - A nonsingular M-matrix (A >= 0) iff all leading principal minors of I - A are > 0."""
    B = [[(1 if i == j else 0) - A[i][j] for j in range(3)] for i in range(3)]
    m1 = B[0][0]
    m2 = B[0][0] * B[1][1] - B[0][1] * B[1][0]
    m3 = (B[0][0] * (B[1][1] * B[2][2] - B[1][2] * B[2][1])
          - B[0][1] * (B[1][0] * B[2][2] - B[1][2] * B[2][0])
          + B[0][2] * (B[1][0] * B[2][1] - B[1][1] * B[2][0]))
    return m1 > 0 and m2 > 0 and m3 > 0


def certify(k):
    xf, tf, af = optimum(k)
    # back off x slightly so rho < 1 strictly, as an exact rational
    x = Fr(xf * (1 - 1e-9)).limit_denominator(10 ** 12)
    t = Fr(tf).limit_denominator(10 ** 12)
    A = matrix(k, x, t)
    assert all(a >= 0 for row in A for a in row)
    if not m_matrix_ok(A):
        raise SystemExit(f"k={k}: rho(A) < 1 not certified")
    iv = mpmath.iv
    iv.dps = 50
    X = iv.mpf(x.numerator) / x.denominator
    T = iv.mpf(t.numerator) / t.denominator
    a = -2 * iv.log(X) / iv.log(T)  # an interval containing the exact value
    mpmath.mp.dps = 60
    lo = mpmath.mpf(a.a)  # its lower end, kept at full precision
    lo12 = mpmath.floor(lo * 10 ** 12) / 10 ** 12  # printed lower bound, rounded down
    a6 = mpmath.floor(lo * 10 ** 6) / 10 ** 6
    return x, t, af, lo, a6


if __name__ == "__main__":
    for k in map(int, sys.argv[1:]):
        x, t, af, lo, a6 = certify(k)
        print(f"k={k}: x={x} t={t} rho(A)<1 certified; a >= {mpmath.nstr(mpmath.floor(lo * 10 ** 12) / 10 ** 12, 12)}; "
              f"alpha_{k} >= {mpmath.nstr(a6, 6)}  (float optimum {af:.8f})")
