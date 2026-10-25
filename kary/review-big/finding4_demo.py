"""Finding 4: the old and the new rounding-up step of kbcertify.py on a crafted interval end.

The interval's upper end is put 1e-30 above a 6-decimal value d whose nearest double lies above d.
Old code (kcertify.py and the 12:42 kbcertify.py): prints d, which is below the interval end.
New code (kbcertify.py, 13:49): prints d + 1e-6.
"""
import math
from fractions import Fraction

import mpmath
from mpmath import iv

iv.dps = 80
mpmath.mp.dps = 80
# a 6-decimal d near the certified values whose nearest double is above d
d = None
for n in range(739975, 740975):
    dec = Fraction(n, 10**6)
    if Fraction(float(dec)) > dec:
        d = dec
        break
print(f"d = {float(d):.6f}: nearest double {float(d).hex()} is above d by {float(Fraction(float(d)) - d):.3e}")
y = mpmath.mpf(d.numerator) / d.denominator + mpmath.mpf("1e-30")
B = iv.mpf([y - mpmath.mpf("1e-20"), y])
my, ey = mpmath.mpf(y).man_exp
yq = Fraction(int(my)) * Fraction(2) ** ey
assert Fraction(float(d)) > yq > d

# old step
six = math.ceil(float(B.b) * 1e6) / 1e6
if mpmath.mpf(six) < B.b:
    six += 1e-6
assert mpmath.mpf(six) >= B.b
old = f"{six:.6f}"

# new step (copied from kbcertify.py, 13:49)
ub = mpmath.mpf(B._mpi_[1])
assert ub == B.b
m, e = ub.man_exp
m = int(m)
upper = Fraction(m) * Fraction(2) ** e
six = Fraction(math.ceil(upper * 10**6), 10**6)
assert six >= upper and six - upper < Fraction(1, 10**6)
six = float(six)
assert Fraction(f"{six:.6f}") >= upper
new = f"{six:.6f}"
print(f"interval upper end exceeds d by {float(upper - d):.3e}")
print(f"old step prints {old}: exact decimal >= upper end? {Fraction(old) >= upper}")
print(f"new step prints {new}: exact decimal >= upper end? {Fraction(new) >= upper}")
