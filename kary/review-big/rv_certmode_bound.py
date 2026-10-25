"""The bound from rvcmp_big's certificate mode, compared with kbcertify's (reviewer, 2026-10-25).

  python rv_certmode_bound.py k rvcmp_output.txt certify_k<k>_h<h>.out

rvcmp_big's certificate mode recomputes every transition with the earlier referee's own delta and
sums with upward rounding, so its "max ratio (up)" is a lambda for an independent certificate of the
same vector. This script turns it into the 6-decimal bound (120-digit interval arithmetic, exact
decimal round-up) and requires it to equal kbcertify's printed bound. A table with missing or wrong
rows shows up here as a larger lambda and a larger bound (review-big/REVIEW.md, finding 1).
"""
import re
import sys
from fractions import Fraction

import mpmath
from mpmath import iv

k = int(sys.argv[1])
txt = open(sys.argv[2]).read()
assert "SELF-TEST OR LOOKUP FAILURE" not in txt, "rvcmp_big reported a self-test or lookup failure"
orb = re.search(r"orbit map: (\d+) stored states share an orbit with another, (\d+) not in kary's own canonical form; "
                r"start state index (-?\d+)", txt)
assert orb and orb.groups() == ("0", "0", "0"), f"orbit map line: {orb.group(0) if orb else 'missing'}"
m = re.search(r"z=(0x[0-9a-fp.+-]+): max ratio \(up\) (0x[0-9a-fp.+-]+)", txt)
z, lam = float.fromhex(m.group(1)), float.fromhex(m.group(2))
ct = open(sys.argv[3]).read()
cz = float.fromhex(re.search(r"z=(0x\S+)", ct).group(1))
lc = float.fromhex(re.search(r"lambda_cert=(\S+)", ct).group(1))
cb = Fraction(re.search(r"CERTIFIED: gamma_\d <= ([0-9.]+)", ct).group(1))
assert cz == z, "different z"
assert 0 < z < 1 / k and 0 < lam < 1

iv.dps = 120
mpmath.mp.dps = 120


def frac(x):
    mm, e = mpmath.mpf(x).man_exp
    return Fraction(int(mm)) * Fraction(2) ** e


zq, lq = Fraction(z), Fraction(lam)
B = 2 * iv.log(k * iv.mpf(zq.numerator) / zq.denominator) / iv.log(iv.mpf(lq.numerator) / lq.denominator)
up = frac(B.b)
six = Fraction((up * 10**6).__ceil__(), 10**6)
print(f"referee lambda (independent transitions, upward) {lam.hex()} = {lam!r}; kbverify lambda_cert {lc.hex()}; "
      f"ratio - 1 = {lam / lc - 1:.3e}")
print(f"referee bound {iv.nstr(B, 15)} -> {float(six):.6f}; kbcertify printed {float(cb):.6f}; "
      f"{'SAME' if six == cb else 'DIFFERENT'}")
sys.exit(0 if six == cb else 1)
