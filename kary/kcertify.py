"""Certified bound gamma_k <= 2 log(kz) / log(lambda_cert) from a stored vector.

1. kweights.py write: exact weights at the vector's z (an exact double), 80 digits,
   rounded up and down to doubles;
2. kverify (directed rounding): lambda_cert >= max_s (A_z v)_s / v_s and
   rho_lower <= min_s (A_z v)_s / v_s;
3. the bound in interval arithmetic (mpmath.iv, 80 digits), rounded up to 6 decimals.
   Also the value at rho_lower: no certificate for this machine and this z
   can give a bound below it.

Usage: python kcertify.py <k> <h> <dir>
"""
import math
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
out = subprocess.run([f"{here}/kverify", str(k), str(h), d], capture_output=True, text=True, check=True).stdout
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
six = math.ceil(float(B.b) * 1e6) / 1e6
if mpmath.mpf(six) < B.b:
    six += 1e-6
assert mpmath.mpf(six) >= B.b
print(f"bound interval: {iv.nstr(B, 15)}")
print(f"best possible at this z (from rho_lower): {iv.nstr(Bl, 12)}")
print(f"CERTIFIED: gamma_{k} <= {six:.6f}   (k={k}, h={h}, z={z!r}, lambda_cert={lam!r})")
