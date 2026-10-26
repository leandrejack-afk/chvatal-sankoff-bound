"""Check every CERT line: the printed bound must equal floor(1e9 * 2r/(2^27 + r)) / 1e9
(LCS) or ceil(1e9 * 2r / 2^27) / 1e9 (edit), computed here with exact fractions, and must be
<= (resp. >=) the exact value. Usage: python3 boundcheck.py <log files or witness claims>"""
import re, sys
from fractions import Fraction
from math import floor, ceil
pat = re.compile(r"CERT k=(\d+) h=(\d+) classes=(\d+) shift=(\d+) r=(-?\d+) (gamma|alpha)_(\d+) (>=|<=) ([0-9.<=]+)")
bad = 0
for f in sys.argv[1:]:
    for line in open(f):
        m = pat.search(line)
        if not m: continue
        k, h, cl, sh, r = map(int, m.groups()[:5]); kind, k2, rel, printed = m[6], int(m[7]), m[8], m[9]
        u = 2 ** sh
        if kind == "gamma":
            exact = Fraction(2 * r, u + r) if r > 0 else None
            want = f"0.{floor(exact * 10**9):09d}" if exact is not None else "<=0"
            ok = (printed == want) and (exact is None or Fraction(printed) <= exact)
        else:
            exact = Fraction(2 * r, u)
            q = ceil(exact * 10**9)
            want = f"{q // 10**9}.{q % 10**9:09d}"
            ok = printed == want and Fraction(printed) >= exact
        bad += not ok
        print(f"{f.split('/')[-1]}: k={k} h={h} r={r} printed {printed} exact {float(exact) if exact else None:.15f} want {want} {'OK' if ok else 'BAD'}")
print("bad:", bad)
