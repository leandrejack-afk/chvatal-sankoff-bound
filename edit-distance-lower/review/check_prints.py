"""Review code: compare every printed number in the certify logs with the exact a* = -2 ln x/ln t
(decimal, 90 digits): the 12-digit 'a >= ...' print and the 6-decimal 'alpha_k >= ...' claim.
Also report the side on which mpmath.mpf(a.a) (53-bit nearest) lands relative to the interval end."""
import re, sys
from decimal import Decimal, getcontext
from fractions import Fraction as Fr
import mpmath
from mpmath.libmp import to_rational
getcontext().prec = 90
for path in sys.argv[1:]:
    print("==", path)
    for line in open(path):
        m = re.search(r"k=(\d+).*?x=(\d+)/(\d+) t=(\d+)/(\d+).*?a (?:>=|in \[) ?([0-9.]+).*alpha_\d+ >= ([0-9.]+)", line)
        if not m: continue
        k = int(m.group(1)); x = Fr(int(m.group(2)), int(m.group(3))); t = Fr(int(m.group(4)), int(m.group(5)))
        a = -2 * (Decimal(x.numerator) / Decimal(x.denominator)).ln() / (Decimal(t.numerator) / Decimal(t.denominator)).ln()
        p12, p6 = Decimal(m.group(6)), Decimal(m.group(7))
        iv = mpmath.iv; iv.dps = 50
        A = -2 * iv.log(iv.mpf(x.numerator) / x.denominator) / iv.log(iv.mpf(t.numerator) / t.denominator)
        lo53 = mpmath.mpf(A.a)
        pq = to_rational(lo53._mpf_); end = to_rational(A._mpi_[0])
        side = Fr(*pq) - Fr(*end)
        print(f"k={k}: a*={str(a)[:19]}  printed12={m.group(6):<15} {'OK' if p12 <= a else 'ABOVE a* by %.1e' % float(p12 - a):<20} "
              f"claim6={m.group(7):<9} {'OK' if p6 <= a else 'FAIL'}  mpf(a.a)-end={float(side):+.1e}")
