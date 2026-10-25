"""Rigorous 2 ln(kz)/ln(lambda) for (k, z hex double, lambda hex double) triples.

  python rv_bound.py K ZHEX LAMHEX [LABEL]
mpmath interval arithmetic at 300 bits; cross-checked with decimal (60 digits, ln correctly
rounded). Prints the interval and the bound rounded up to 6 decimals, with the exact margin.
"""
import sys
from decimal import Decimal, getcontext
from fractions import Fraction
from math import ceil

from mpmath import iv, mp


def exact(x):
    return Fraction(float.fromhex(x))


def main():
    k = int(sys.argv[1])
    z, lam = exact(sys.argv[2]), exact(sys.argv[3])
    label = sys.argv[4] if len(sys.argv) > 4 else ""
    assert 0 < z < Fraction(1, k) and 0 < lam < 1
    iv.prec = 300
    kz = iv.mpf(k * z.numerator) / iv.mpf(z.denominator)
    li = iv.mpf(lam.numerator) / iv.mpf(lam.denominator)
    B = 2 * iv.log(kz) / iv.log(li)
    lo_t, hi_t = B._mpi_
    lo = (-1) ** lo_t[0] * Fraction(lo_t[1]) * Fraction(2) ** lo_t[2]
    hi = (-1) ** hi_t[0] * Fraction(hi_t[1]) * Fraction(2) ** hi_t[2]
    getcontext().prec = 60
    d = 2 * (Decimal(k * z.numerator) / Decimal(z.denominator)).ln() / (Decimal(lam.numerator) / Decimal(lam.denominator)).ln()
    up6 = Fraction(ceil(hi * 10 ** 6), 10 ** 6)
    assert Fraction(d) <= up6 and lo <= hi and hi - lo < Fraction(1, 10 ** 60)
    mp.dps = 20
    print(f"{label} k={k} z={float(z)!r} lambda={float(lam)!r}: bound in [{mp.nstr(mp.mpf(lo.numerator) / lo.denominator, 16)}, "
          f"{mp.nstr(mp.mpf(hi.numerator) / hi.denominator, 16)}], decimal {str(d)[:18]}, "
          f"certified {float(up6):.6f} (margin {float(up6 - hi):.2e})")


if __name__ == "__main__":
    main()
