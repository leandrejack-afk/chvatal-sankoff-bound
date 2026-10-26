"""Independent certificate check for the rules 1-3 and rules 1-4 automata (review code).

Censoring the run states (Schur complement; rho(A22) = I = x t < 1 is needed) leaves a matrix K on
the diagonal-type states.  With a = M + S, b = M W_R, W_R = sum_{r=1}^R (Iq)^r + q^R I^(R+1)/(1-I):
  rules 1-3: K = [a + 2b]                         -> rho(A) < 1  iff  a + 2b < 1
  rules 1-4: K = [[a, b, b], [a, b, 0], [a, 0, b]] (diag, MI, MD)
             eigenvalues b and roots of l^2 - (a+b) l - ab = 0  ->  rho < 1  iff  a + b + ab < 1
(R = oo: W = Iq/(1 - Iq).)  Exact Fractions; a* = -2 ln x / ln t by decimal at 90 digits.
Also: the optimum of a over (x, t) for R = 32 and R = oo from these closed forms (mpmath, 40 digits).
usage: verify_rule4.py r3|r4 <rerun file with x=, t=>   |   verify_rule4.py opt
"""
import re
import sys
from decimal import Decimal, getcontext
from fractions import Fraction as Fr

import mpmath

getcontext().prec = 90


def W(I, q, R):
    if R is None:
        return I * q / (1 - I * q)
    return sum((I * q) ** r for r in range(1, R + 1)) + q ** R * I ** (R + 1) / (1 - I)


def ab(k, x, t, R, one=1):
    q = one * (k - 1) / k
    M, S, I = x * x / k, x * x * t * q, x * t
    return M + S, M * W(I, q, R), I


def margin(rules, k, x, t, R, one=1):
    a, b, I = ab(k, x, t, R, one)
    assert I < 1
    return 1 - (a + 2 * b) if rules == "r3" else 1 - (a + b + a * b)


def dec(fr):
    return Decimal(fr.numerator) / Decimal(fr.denominator)


def check(rules, path):
    for line in open(path):
        m = re.search(r"k=(\d+).*?x=(\d+)/(\d+) t=(\d+)/(\d+).*alpha_\d+ >= ([0-9.]+)", line)
        if not m:
            continue
        k = int(m.group(1))
        x, t = Fr(int(m.group(2)), int(m.group(3))), Fr(int(m.group(4)), int(m.group(5)))
        claimed = Fr(m.group(6))
        mg = margin(rules, k, x, t, 32, Fr(1))
        a = -2 * dec(x).ln() / dec(t).ln()
        ok = mg > 0 and Fr(a - Decimal("1e-80")) >= claimed
        print(f"k={k}: 1 - (rho condition) = {float(mg):.3e} > 0: {mg > 0}; a* = {str(a)[:20]}; "
              f"claimed {m.group(6)}; {'OK' if ok else 'FAIL'}")


def optimum(rules, k, R):
    mpmath.mp.dps = 40

    def x_of_t(t):
        f = lambda x: margin(rules, k, x, t, R, mpmath.mpf(1))
        lo, hi = mpmath.mpf(1), mpmath.mpf(1) / t * 0.999999  # need I = x t < 1
        if R is None:
            hi = min(hi, k / ((k - 1) * t) * 0.999999)
        return mpmath.findroot(f, (lo, hi), solver="anderson")

    def a_of_t(t):
        return -2 * mpmath.log(x_of_t(t)) / mpmath.log(t)

    # golden-section on t in (0.02, 0.6)
    lo, hi = mpmath.mpf("0.02"), mpmath.mpf("0.6")
    g = (mpmath.sqrt(5) - 1) / 2
    c, d = hi - g * (hi - lo), lo + g * (hi - lo)
    fc, fd = a_of_t(c), a_of_t(d)
    for _ in range(160):
        if fc > fd:
            hi, d, fd = d, c, fc
            c = hi - g * (hi - lo)
            fc = a_of_t(c)
        else:
            lo, c, fc = c, d, fd
            d = lo + g * (hi - lo)
            fd = a_of_t(d)
    return (c + d) / 2, a_of_t((c + d) / 2)


if __name__ == "__main__":
    if sys.argv[1] == "opt":
        for rules in ("r3", "r4"):
            for k in (2, 3, 4, 5, 6, 7, 8, 16, 32):
                t32, a32 = optimum(rules, k, 32)
                tinf, ainf = optimum(rules, k, None)
                print(f"{rules} k={k}: R=32 optimum a = {mpmath.nstr(a32, 12)} (t={mpmath.nstr(t32, 8)}); "
                      f"R=oo optimum a = {mpmath.nstr(ainf, 12)}; gain from R=oo {mpmath.nstr(ainf - a32, 3)}",
                      flush=True)
    else:
        check(sys.argv[1], sys.argv[2])
