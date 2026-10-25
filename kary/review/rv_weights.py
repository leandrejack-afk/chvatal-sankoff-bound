"""Exact class weights at the double nearest Z, rounded down / nearest / up, one line per class.

  ./rvk classes K H | python rv_weights.py K H Z > weights.txt
Weights are exact Fractions (rv_machine.weight with z = the exact double), then rounded.
"""
import math
import sys
from fractions import Fraction

import rv_machine as M


def round3(x):
    f = float(x)  # int/int true division: correctly rounded
    up = f if Fraction(f) >= x else math.nextafter(f, math.inf)
    dn = f if Fraction(f) <= x else math.nextafter(f, -math.inf)
    assert Fraction(dn) <= x <= Fraction(up) and math.nextafter(dn, math.inf) >= up
    return dn, f, up


def main():
    k, h = int(sys.argv[1]), int(sys.argv[2])
    z = Fraction(float(sys.argv[3]))
    assert 0 < z < Fraction(1, k)
    cache = {}
    for line in sys.stdin:
        la, lb, nA, nB = map(int, line.split())
        key = (la, lb, nA, nB, la == h, lb == h)
        if key not in cache:
            cache[key] = round3(M.weight(k, key, z))
        print(" ".join(x.hex() for x in cache[key]))
    print(f"z = {float(z).hex()} ({float(z)!r}), {len(cache)} weight types", file=sys.stderr)


if __name__ == "__main__":
    main()
