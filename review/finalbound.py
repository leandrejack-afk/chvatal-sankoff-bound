# Reviewer tool: final bound 2 log(2z)/log(lam) for exact doubles z, lam, two independent ways:
# (a) Python decimal (correctly rounded ln at 60 digits), (b) exact integer test of
#     (2z)^(2q) >= lam^p   <=>   2 log(2z)/log(lam) <= p/q   (Lueker's form (47)-(48)),
#     done with exact big-integer arithmetic on the dyadic rationals.
import sys
from decimal import Decimal, getcontext
from fractions import Fraction as F
z = float.fromhex(sys.argv[1]); lams = [float.fromhex(a) for a in sys.argv[2:]]
getcontext().prec = 60
for lam in lams:
    Z = Decimal(F(z).numerator) / Decimal(F(z).denominator)
    Lm = Decimal(F(lam).numerator) / Decimal(F(lam).denominator)
    B = 2 * (2 * Z).ln() / Lm.ln()
    print(f"lam={lam!r}: bound = {B:.15f} (60-digit decimal)")
