"""Check kary's rounded weights against the referee's exact weights, per kary weight type.

  python rv_theirweights.py K H DATADIR OUTFILE
Reads DATADIR/wtypes_kK_hH.txt, DATADIR/weights_kK_hH.txt and z from DATADIR/kvec_kK_hH.bin.
For each type: exact weight (Fraction, rv_machine.weight at the exact double z); kary's up must
be >= exact, kary's down <= exact, and they must be the two neighbouring doubles (or equal).
Writes OUTFILE with the referee's own "down up" hex pair per kary weight type (for rvcmp).
"""
import math
import struct
import sys
from fractions import Fraction

import rv_machine as M


def main():
    k, h, d, outp = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4]
    with open(f"{d}/kvec_k{k}_h{h}.bin", "rb") as f:
        f.read(16)
        (z,) = struct.unpack("<d", f.read(8))
    zq = Fraction(z)
    types = [tuple(map(int, l.split())) for l in open(f"{d}/wtypes_k{k}_h{h}.txt")]
    theirs = [tuple(float.fromhex(x) for x in l.split()) for l in open(f"{d}/weights_k{k}_h{h}.txt")]
    assert len(types) == len(theirs)
    bad = 0
    out = []
    for (i, size, ta, tb, na, nb), (up, dn) in zip(types, theirs):
        la = h if ta else (size - h if tb else min(h - 1, size - 1))
        lb = h if tb else size - la
        key = (la, lb, na, nb, bool(ta), bool(tb))
        assert la + lb == size and 1 <= la <= h and 1 <= lb <= h
        w = M.weight(k, key, zq)
        ok = Fraction(dn) <= w <= Fraction(up) and (up == dn or math.nextafter(dn, math.inf) == up)
        bad += not ok
        f = float(w)
        mu = f if Fraction(f) >= w else math.nextafter(f, math.inf)
        md = f if Fraction(f) <= w else math.nextafter(f, -math.inf)
        out.append(f"{md.hex()} {mu.hex()}")
        if not ok:
            print("BAD", i, size, ta, tb, na, nb, up.hex(), dn.hex(), float(w))
    open(outp, "w").write("\n".join(out) + "\n")
    print(f"k={k} h={h} z={z.hex()}: {len(types)} weight types, kary up >= exact >= kary down and adjacent: "
          f"{len(types) - bad}/{len(types)}")


if __name__ == "__main__":
    main()
