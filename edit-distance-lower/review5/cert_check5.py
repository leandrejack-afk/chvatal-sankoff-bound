"""Review of rule 5 (independent): re-check every certificate in certify_rule5_R32.log.

For each logged k, x, t:
 (1) the author's matrix certify_rule5.build(exact=True) equals the reviewer's matrix from the named
     transition table (table5.matrix), entry by entry, for R in {1, 2, 3, 32} (also at other x, t);
 (2) own vector: v = (Id - A)^(-1) 1 solved at 60 digits (mpmath LU) on the reviewer's matrix,
     converted to exact rationals, then A v < v checked in exact Fractions (implies rho(A) < 1);
 (3) scalar test with no matrix: censoring the run states leaves K on (G0, G1, MD, MI),
       K = [[M, S, b, b], [M, S, c, c], [M, S, b, 0], [M, S, 0, b]],
     b = M W_R (plain runs), c = M V_R (runs right after S), with
       W_R = sum_{r<=R} (Iq)^r + q^R I^(R+1)/(1-I),
       V_R = sum_{r<=R} (Iq)^r h(r) + h(R) (Iq)^(R+1)/(1-Iq),
     eigenvalues 0, b and the roots of l^2 - (a+b) l + (ab - 2Mb - 2Sc), a = M + S; rho(A) < 1 iff
     I = x t < 1, a + b < 2 and f(1) = 1 - a - b + ab - 2Mb - 2Sc > 0.  Checked exactly.
 (4) a* = -2 ln x / ln t at 90 digits (decimal); printed 12- and 6-digit values must be <= a*.
 (5) discrimination: the rules 1-4 scalar test (a + b + ab < 1) at the rule-5 (x, t) must FAIL, and
     the rule-5 scalar test at x (1 + 1e-8) must FAIL.
usage: cert_check5.py <log>
"""
import importlib.util
import re
import sys
from decimal import Decimal, getcontext
from fractions import Fraction as Fr

sys.dont_write_bytecode = True
import mpmath  # noqa: E402

from table5 import h, matrix, qp  # noqa: E402

HERE = __file__.rsplit("/", 1)[0]
spec = importlib.util.spec_from_file_location("certify_rule5", HERE + "/../certify_rule5.py")
author = importlib.util.module_from_spec(spec)
spec.loader.exec_module(author)

getcontext().prec = 90


def scalar(k, x, t, R):
    q, p = qp(k)
    M, S, I = x * x / k, x * x * t * q, x * t
    a = M + S
    W = sum((I * q) ** r for r in range(1, R + 1)) + q ** R * I ** (R + 1) / (1 - I)
    V = sum((I * q) ** r * h(k, r) for r in range(1, R + 1)) + h(k, R) * (I * q) ** (R + 1) / (1 - I * q)
    b, c = M * W, M * V
    f1 = 1 - a - b + a * b - 2 * M * b - 2 * S * c
    return I, a, b, c, M, S, f1


def scalar_ok(k, x, t, R):
    I, a, b, c, M, S, f1 = scalar(k, x, t, R)
    return I < 1 and a + b < 2 and f1 > 0, f1


def rules14_ok(k, x, t, R):
    q, _ = qp(k)
    M, S, I = x * x / k, x * x * t * q, x * t
    W = sum((I * q) ** r for r in range(1, R + 1)) + q ** R * I ** (R + 1) / (1 - I)
    a, b = M + S, M * W
    return 1 - (a + b + a * b) > 0


def compare_matrices(k, x, t, R):
    A1 = author.build(k, x, t, R, exact=True)
    st, A2 = matrix(k, x, t, R)
    n = len(st)
    assert len(A1) == n and all(len(row) == n for row in A1)
    diffs = [(st[i], st[j], A1[i][j], A2[i][j]) for i in range(n) for j in range(n) if A1[i][j] != A2[i][j]]
    return diffs


def own_vector_check(k, x, t, R):
    st, A = matrix(k, x, t, R)
    n = len(st)
    mpmath.mp.dps = 60
    B = mpmath.matrix(n, n)
    for i in range(n):
        for j in range(n):
            B[i, j] = (1 if i == j else 0) - mpmath.mpf(A[i][j].numerator) / A[i][j].denominator
    v = mpmath.lu_solve(B, mpmath.matrix([1] * n))
    vq = []
    for i in range(n):
        man, exp = v[i].man_exp
        vq.append(Fr(int(man)) * (Fr(2) ** int(exp)))
    if not all(c > 0 for c in vq):
        return False, None
    worst = max(sum(A[i][j] * vq[j] for j in range(n) if A[i][j]) / vq[i] for i in range(n))
    return worst < 1, worst


if __name__ == "__main__":
    allok = True
    for line in open(sys.argv[1]):
        m = re.search(r"k=(\d+) R=(\d+): x=(\d+)/(\d+) t=(\d+)/(\d+).*a >= ([0-9.]+); alpha_\d+ >= ([0-9.]+)", line)
        if not m:
            continue
        k, R = int(m.group(1)), int(m.group(2))
        x, t = Fr(int(m.group(3)), int(m.group(4))), Fr(int(m.group(5)), int(m.group(6)))
        a12, a6 = m.group(7), m.group(8)
        # (1)
        nd = 0
        for RR in (1, 2, 3, R):
            for (xx, tt) in ((x, t), (Fr(3, 2), Fr(1, 7)), (Fr(11, 10), Fr(2, 5))):
                d = compare_matrices(k, xx, tt, RR)
                nd += len(d)
                for e in d[:3]:
                    print("   MATRIX DIFF", k, RR, e)
        # (2)
        ok2, worst = own_vector_check(k, x, t, R)
        # (3)
        ok3, f1 = scalar_ok(k, x, t, R)
        ok3inf = None
        # (4)
        astar = -2 * (Decimal(x.numerator) / Decimal(x.denominator)).ln() / (Decimal(t.numerator) / Decimal(t.denominator)).ln()
        ok4 = Decimal(a12) <= astar and Decimal(a6) <= astar
        # (5)
        disc1 = not rules14_ok(k, x, t, R)
        disc2 = not scalar_ok(k, x * Fr(100000001, 100000000), t, R)[0]
        ok = nd == 0 and ok2 and ok3 and ok4 and disc1 and disc2
        allok &= ok
        print(f"k={k:2d}: matrix == own table (R=1,2,3,{R}; 3 points): {nd == 0}; own v: max (Av)_i/v_i = "
              f"{float(worst):.12f} < 1: {ok2}; scalar f(1) = {float(f1):.4e} > 0: {ok3}; a* = {str(astar)[:22]} "
              f">= printed {a12} and {a6}: {ok4}; rules1-4 test fails here: {disc1}; x(1+1e-8) fails: {disc2}"
              f" -> {'OK' if ok else 'FAIL'}", flush=True)
    print("ALL OK" if allok else "FAILURES")
