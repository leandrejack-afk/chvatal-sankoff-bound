# Reviewer tool: exact (Fraction) check that each weight in a weights file is >= the true
# w_{mu,z}(T_{c,l}) for the exact double z stored in the eigen file, and how tight it is.
import struct, sys
from fractions import Fraction as F
h = int(sys.argv[1]); eig = sys.argv[2]; wf = sys.argv[3]
with open(eig, "rb") as f:
    eh, n, z = struct.unpack("<iQd", f.read(20))
assert eh == h
Z = F(z)  # exact value of the double
ws = [float.fromhex(s.strip()) for s in open(wf) if s.strip()]
assert len(ws) == 2*h - 1, len(ws)
worst = None; ok = True
for k, l in enumerate(range(-(h-1), h)):
    a = abs(l)
    exact = Z**(a+2) if a < h-1 else Z**(h+1) / (1 - Z)
    W = F(ws[k])
    rel = (W - exact) / exact
    if W < exact: ok = False
    # is it the smallest double >= exact?  (next double below W must be < exact)
    import math
    below = F(math.nextafter(ws[k], 0.0))
    tight = below < exact
    print(f"l={l:+d} |l|={a} exact~{float(exact):.17g} w_up={ws[k]!r} (w_up-exact)/exact={float(rel):.3e} {'>=OK' if W>=exact else 'BELOW!'} {'tight' if tight else 'NOT-TIGHT'}")
print(f"h={h} z={z!r} ({z.hex()}) all weights >= exact: {ok}")
