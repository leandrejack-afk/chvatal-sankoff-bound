"""Exact rational certificate for h = 2 (k = 3, 4), from rv_machine.py only.

  python rv_h2cert.py K Znum Zden SYM
Builds M_2 (orbit machine if SYM=1, raw machine if SYM=0), exact weights at z = Znum/Zden,
a rational v > 0 (power iteration rounded to 1e-12 grid), a rational lambda (ceil to 1e-15),
verifies A_z v <= lambda v in Fractions, then bounds 2 ln(kz)/ln(lambda) with mpmath
interval arithmetic and, independently, with the decimal module (correctly rounded ln).
"""
import sys
from decimal import Decimal, getcontext
from fractions import Fraction
from math import ceil

from mpmath import iv, mp

import rv_machine as M


def main():
    k, zn, zd, sym = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4] == "1"
    h = 2
    z = Fraction(zn, zd)
    assert 0 < z < Fraction(1, k)
    states, rows, cl = M.build(k, h, sym)
    n = len(states)
    w = [M.weight(k, M.wkey(k, h, c), z) for c in cl]
    A = [dict() for _ in range(n)]
    for s, row in enumerate(rows):
        for t, ci in row:
            A[s][t] = A[s].get(t, 0) + w[ci]
    # float power iteration for a good v
    vf = [1.0] * n
    for _ in range(400):
        nv = [sum(float(a) * vf[t] for t, a in A[s].items()) for s in range(n)]
        mx = max(nv)
        vf = [x / mx for x in nv]
    D = 10 ** 12
    v = [Fraction(round(x * D), D) for x in vf]
    assert all(x > 0 for x in v)
    Av = [sum(a * v[t] for t, a in A[s].items()) for s in range(n)]
    ratios = [Av[s] / v[s] for s in range(n)]
    lam_max = max(ratios)
    lam = Fraction(ceil(lam_max * 10 ** 15), 10 ** 15)
    ok = all(Av[s] <= lam * v[s] for s in range(n))
    assert ok and 0 < lam < 1
    print(f"k={k} h=2 sym={int(sym)} z={z} states={n} classes={len(cl)} "
          f"accepted={sum(len(r) for r in rows)}")
    print(f"  exact max ratio = {float(lam_max)!r}, min ratio = {float(min(ratios))!r}")
    print(f"  lambda = {lam}   A_z v <= lambda v exactly: {ok}")
    # rigorous bound, mpmath intervals
    iv.prec = 400
    kz = iv.mpf(k * z.numerator) / iv.mpf(z.denominator)
    li = iv.mpf(lam.numerator) / iv.mpf(lam.denominator)
    B = 2 * iv.log(kz) / iv.log(li)
    mp.dps = 30
    print(f"  bound interval (mpmath iv, 400 bits): [{mp.nstr(B.a, 20)}, {mp.nstr(B.b, 20)}]")
    # independent: decimal, 60 digits, ln correctly rounded
    getcontext().prec = 60
    d = 2 * (Decimal(k * z.numerator) / Decimal(z.denominator)).ln() / (Decimal(lam.numerator) / Decimal(lam.denominator)).ln()
    print(f"  bound (decimal 60 digits) = {d}")
    sg, man, ex, _ = B._mpi_[1]  # exact upper endpoint
    hi = (-1) ** sg * Fraction(man) * (Fraction(2) ** ex)
    up6 = Fraction(ceil(hi * 10 ** 6), 10 ** 6)
    print(f"  certified (rounded up to 6 decimals): {float(up6):.6f}")
    if sym:
        print("  states (x, y, dx, dy) and v:")
        for s in range(n):
            print("   ", states[s], v[s])


if __name__ == "__main__":
    main()
