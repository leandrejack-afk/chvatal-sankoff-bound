"""Turn a vector file into a rigorous bound on gamma_2.

1. read z (an exact double) from eigen_h<h>.bin;
2. compute each weight w_{mu,z}(T_{c,l}) = z^{|l|+2} (|l| < h-1) or
   z^{h+1}/(1-z) (|l| = h-1) to 80 digits and round it UP to a double;
3. run verify_cert (directed rounding) to get lambda_cert with A_z v <= lambda_cert v;
4. evaluate 2 log(2z)/log(lambda_cert) in interval arithmetic (mpmath.iv, 80 digits)
   and print a 6-decimal upper bound that is certainly >= it.
Usage: python certify.py <h> <data dir>
"""
import math
import struct
import subprocess
import sys

import mpmath
from mpmath import iv

h = int(sys.argv[1])
d = sys.argv[2]
with open(f"{d}/eigen_h{h}.bin", "rb") as f:
    eh, n, z = struct.unpack("<iQd", f.read(4 + 8 + 8))
assert eh == h
mpmath.mp.dps = 80
iv.dps = 80
zr = mpmath.mpf(z)          # exact value of the stored double


def up(x):
    fx = float(x)
    if mpmath.mpf(fx) < x:
        fx = math.nextafter(fx, math.inf)
    assert mpmath.mpf(fx) >= x
    return fx


lines = []
for l in range(-(h - 1), h):
    a = abs(l)
    wexact = zr ** (a + 2) if a < h - 1 else zr ** (h + 1) / (1 - zr)
    lines.append(float.hex(up(wexact)))
wfile = f"{d}/weights_up_h{h}.txt"
with open(wfile, "w") as f:
    f.write("\n".join(lines) + "\n")
out = subprocess.run([__import__("os").environ.get("VERIFIER", "./verify_cert"), f"{d}/trans_h{h}.bin", f"{d}/eigen_h{h}.bin", wfile],
                     capture_output=True, text=True, check=True).stdout
print(out.strip())
lam = float.fromhex(out.split("lambda_cert=")[1].split()[0])
assert 0 < lam < 1 and 0 < z < 0.5
Z = iv.mpf(mpmath.mpf(z))
LAM = iv.mpf(mpmath.mpf(lam))
B = 2 * iv.log(2 * Z) / iv.log(LAM)
print(f"bound interval: [{mpmath.nstr(B.a, 15)}, {mpmath.nstr(B.b, 15)}]")
six = math.ceil(float(B.b) * 1e6) / 1e6
if mpmath.mpf(six) < B.b:
    six += 1e-6
assert mpmath.mpf(six) >= B.b
print(f"CERTIFIED: gamma_2 <= {six:.6f}   (h={h}, z={z!r}, lambda_cert={lam!r})")
