"""Bound from rvb_verify's output (reviewer, 2026-10-25).

  python rvb_bound.py k rvb_verify.out [kbcertify.out]

1. From lambda_emu (directed rounding, emulated): 2 log(kz)/log(lambda) in mpmath interval
   arithmetic (120 digits); the smallest 6-decimal number >= its upper end, compared as an exact
   decimal (not as a double).
2. From the round-to-nearest maximum alone: lambda_apriori = qn_max * (1 - 2^-53)^-(nmax + 2),
   evaluated exactly (Fraction), then the same bound. This needs no directed rounding at all.
3. With a kbcertify output: its lambda_cert must equal lambda_emu bit for bit, and its
   printed bound must equal the 6-decimal value of step 1.
"""
import re
import sys
from fractions import Fraction

import mpmath
from mpmath import iv

k = int(sys.argv[1])
txt = open(sys.argv[2]).read()
assert "\nOK" in txt, "rvb_verify did not print OK"
z = float.fromhex(re.search(r"z=(\S+)", txt).group(1))
lam = float.fromhex(re.search(r"lambda_emu=(\S+)", txt).group(1))
rho = float.fromhex(re.search(r"rho_emu=(\S+)", txt).group(1))
qn = float.fromhex(re.search(r"nearest_max_up_weights=(\S+)", txt).group(1))
nmax = int(re.search(r"nmax=(\d+)", txt).group(1))
assert 0 < z < 1 / k and 0 < rho <= lam < 1

iv.dps = 120
mpmath.mp.dps = 120


def frac(x):
    """exact value of an mpmath number"""
    m, e = mpmath.mpf(x).man_exp
    return Fraction(m) * Fraction(2) ** e


def bound(lam_exact):
    Z = iv.mpf(Fraction(z).numerator) / iv.mpf(Fraction(z).denominator)
    L = iv.mpf(lam_exact.numerator) / iv.mpf(lam_exact.denominator)
    B = 2 * iv.log(k * Z) / iv.log(L)
    hi = frac(B.b)                          # exact binary value of the interval's upper end
    six = Fraction((hi * 10**6).__ceil__(), 10**6)
    assert six >= hi
    return B, six


u = Fraction(1, 2**53)
lam_ap = Fraction(qn) / (1 - u) ** (nmax + 2)
B1, s1 = bound(Fraction(lam))
B2, s2 = bound(lam_ap)
B3, _ = bound(Fraction(rho))
print(f"z = {z!r} ({z.hex()}), k = {k}")
print(f"lambda_emu      = {lam.hex()} = {lam!r}")
print(f"lambda_apriori  = qn_max (1-u)^-(nmax+2), nmax = {nmax}: {float(lam_ap)!r} "
      f"(lambda_apriori / lambda_emu - 1 = {float(lam_ap / Fraction(lam) - 1):.3e})")
print(f"bound from lambda_emu:     {iv.nstr(B1, 15)} -> certified {float(s1):.6f}")
print(f"bound from lambda_apriori: {iv.nstr(B2, 15)} -> certified {float(s2):.6f}")
print(f"bound at rho_emu (no certificate at this z can go below): {iv.nstr(B3, 12)}")
print(f"same 6-decimal bound by both routes: {s1 == s2}")
if len(sys.argv) > 3:
    ct = open(sys.argv[3]).read()
    lc = re.search(r"lambda_cert=(\S+)", ct).group(1)
    cb = re.search(r"CERTIFIED: gamma_\d <= ([0-9.]+)", ct).group(1)
    print(f"kbcertify lambda_cert {lc} equals lambda_emu: {float.fromhex(lc) == lam}")
    print(f"kbcertify printed {cb}; exact decimal >= interval upper end: {Fraction(cb) >= frac(B1.b)}; "
          f"equals this script's value: {Fraction(cb) == s1}")
