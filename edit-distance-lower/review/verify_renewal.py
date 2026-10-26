"""Independent re-verification of certify_runs.py / certify_simple.py certificates (review code).

Method (different from the implementer's A v < v test): in the run automaton every cycle passes
through the diag state except the self-loops at I_R and D_R (weight I = x t).  For an irreducible
nonnegative matrix, rho(A) < 1 iff sum_L (A^L)_{00} < oo, and by the renewal decomposition
sum_L (A^L)_{00} = 1/(1 - F) whenever I < 1 and the first-return series F < 1, where

  F_R = M + S + 2 M [ sum_{r=1}^{R} (I q)^r + q^R I^{R+1} / (1 - I) ]      (run automaton, cap R)
  F_3 = M + S + 2 M I / (1 - I)                                          (3-state automaton)
  F_oo = M + S + 2 M (I q) / (1 - I q)                                    (uncapped, R = oo)

with M = x^2/k, S = x^2 t (k-1)/k, I = x t, q = (k-1)/k.  All checks are exact Fractions.
Then a* = -2 ln x / ln t in decimal arithmetic (correctly rounded ln) at 90 digits, compared with
the published 6-decimal value as an exact decimal string.
usage: verify_renewal.py <file with lines 'k=.. x=p/q t=p/q ... alpha_k >= v'> {runs R | simple}
"""
import re
import sys
from decimal import Decimal, getcontext
from fractions import Fraction as Fr

getcontext().prec = 90


def weights(k, x, t):
    q = Fr(k - 1, k)
    return x * x / k, x * x * t * q, x * t, q


def F_runs(k, x, t, R):
    M, S, I, q = weights(k, x, t)
    assert I < 1
    geo = sum((I * q) ** r for r in range(1, R + 1))
    return M + S + 2 * M * (geo + q ** R * I ** (R + 1) / (1 - I))


def F_simple(k, x, t):
    M, S, I, q = weights(k, x, t)
    assert I < 1
    return M + S + 2 * M * I / (1 - I)


def F_inf(k, x, t):
    M, S, I, q = weights(k, x, t)
    assert I * q < 1
    return M + S + 2 * M * (I * q) / (1 - I * q)


def dec(fr):
    return Decimal(fr.numerator) / Decimal(fr.denominator)


def main():
    path, mode = sys.argv[1], sys.argv[2]
    R = int(sys.argv[3]) if mode == "runs" else None
    for line in open(path):
        m = re.search(r"k=(\d+).*?x=(\d+)/(\d+) t=(\d+)/(\d+).*alpha_\d+ >= ([0-9.]+)", line)
        if not m:
            continue
        k = int(m.group(1))
        x, t = Fr(int(m.group(2)), int(m.group(3))), Fr(int(m.group(4)), int(m.group(5)))
        claimed = Fr(m.group(6))  # exact decimal string as printed
        assert x > 1 and 0 < t < 1
        F = F_runs(k, x, t, R) if mode == "runs" else F_simple(k, x, t)
        Fi = F_inf(k, x, t)
        a = -2 * (dec(x).ln()) / (dec(t).ln())
        a_lo = a - Decimal("1e-80")  # ln correctly rounded at 90 digits: error far below this
        ok = F < 1 and Fr(a_lo) >= claimed
        print(f"k={k}: x*t = {float(x*t):.9f} < 1; 1 - F = {float(1 - F):.3e} (F<1: {F < 1}); "
              f"F_oo <= F: {Fi <= F}; a* = {str(a)[:22]}; claimed {m.group(6)}; "
              f"slack a*-claimed = {float(Fr(a_lo) - claimed):.3e}; {'OK' if ok else 'FAIL'}")


if __name__ == "__main__":
    main()
