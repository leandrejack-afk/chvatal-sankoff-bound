"""kcertify.py for kbig tables: certified bound gamma_k <= 2 log(kz) / log(lambda_cert).

1. kweights.py write: exact weights at the vector's z (an exact double), 80 digits,
   rounded up and down to doubles;
2. kbverify (directed rounding, records streamed from kbig_k<k>_h<h>.bin): lambda_cert >= max_s (A_z v)_s / v_s and
   rho_lower <= min_s (A_z v)_s / v_s;
3. the bound in interval arithmetic (mpmath.iv, 80 digits), rounded up to 6 decimals.
   Also the value at rho_lower: no certificate for this machine and this z
   can give a bound below it.

Usage: python kbcertify.py <k> <h> <dir>
"""
import math
from fractions import Fraction
import os
import struct
import subprocess
import sys

import mpmath
from mpmath import iv

k, h, d = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
here = os.path.dirname(os.path.abspath(__file__))
py = sys.executable
print(subprocess.run([py, f"{here}/kweights.py", "write", str(k), str(h), d],
                     capture_output=True, text=True, check=True).stdout.strip())
out = subprocess.run([f"{here}/kbverify", str(k), str(h), d], capture_output=True, text=True, check=True).stdout
print(out.strip())
lam = float.fromhex(out.split("lambda_cert=")[1].split()[0])
rho = float.fromhex(out.split("rho_lower=")[1].split()[0])
with open(f"{d}/kvec_k{k}_h{h}.bin", "rb") as f:
    f.read(8 + 8)
    (z,) = struct.unpack("<d", f.read(8))
assert 0 < z < 1.0 / k and 0 <= rho <= lam < 1
iv.dps = 80
mpmath.mp.dps = 80
Z = iv.mpf(mpmath.mpf(z))
B = 2 * iv.log(k * Z) / iv.log(iv.mpf(mpmath.mpf(lam)))
Bl = 2 * iv.log(k * Z) / iv.log(iv.mpf(mpmath.mpf(rho))) if rho > 0 else iv.mpf(0)
# the printed 6-decimal value, as an exact decimal, must be >= the exact upper end of B
# (review-big/REVIEW.md, finding 4: comparing the double nearest to it is not enough)
ub = mpmath.mpf(B._mpi_[1])          # the upper end of B as an exact binary number
assert ub == B.b
m, e = ub.man_exp
m = int(m)
upper = Fraction(m) * Fraction(2) ** e
six = Fraction(math.ceil(upper * 10**6), 10**6)
assert six >= upper and six - upper < Fraction(1, 10**6)
six = float(six)
assert Fraction(f"{six:.6f}") >= upper
print(f"bound interval: {iv.nstr(B, 15)}")
print(f"best possible at this z (from rho_lower): {iv.nstr(Bl, 12)}")
print(f"CERTIFIED: gamma_{k} <= {six:.6f}   (k={k}, h={h}, z={z!r}, lambda_cert={lam!r})")
