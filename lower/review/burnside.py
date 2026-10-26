"""Independent class counts for klow.c / kedit.c.

Orbits of pairs (s, t) of length-h strings over k letters under S_k x Z_2 (alphabet
permutations and swapping s, t), by Burnside's lemma over cycle types:
  #orbits = (1/(2 k!)) * sum_pi [ fix(pi)^(2h) + fix(pi^2)^h ]
fix(pi) = number of fixed letters of pi; a pair fixed by pi.swap has t fixed by pi^2, s = pi(t).
Also the number of RGS of length L with at most k labels = sum_{j<=k} S(L, j).
Usage: python3 burnside.py   (compares with every 'k=.. h=.. rgs=.. classes=.. F=..' line in runs/*.log)
"""
import glob, re, sys
from fractions import Fraction
from math import factorial
from functools import lru_cache

def partitions(n, maxpart=None):
    if maxpart is None: maxpart = n
    if n == 0: yield []; return
    for p in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - p, p):
            yield [p] + rest

def class_size(part, k):
    # number of permutations of S_k with cycle type part
    from collections import Counter
    c = Counter(part)
    den = 1
    for length, mult in c.items():
        den *= length ** mult * factorial(mult)
    return factorial(k) // den

def orbits(k, h):
    tot = 0
    for part in partitions(k):
        n = class_size(part, k)
        f1 = part.count(1)                      # fixed points of pi
        f2 = sum(p for p in part if p in (1, 2))  # fixed points of pi^2
        tot += n * (f1 ** (2 * h) + f2 ** h)
    assert tot % (2 * factorial(k)) == 0
    return tot // (2 * factorial(k))

@lru_cache(None)
def S2(n, j):
    if n == 0 and j == 0: return 1
    if n == 0 or j == 0: return 0
    return j * S2(n - 1, j) + S2(n - 1, j - 1)

def nrgs(L, k):
    return sum(S2(L, j) for j in range(1, k + 1))

if __name__ == "__main__":
    root = sys.argv[1] if len(sys.argv) > 1 else "."
    pat = re.compile(r"k=(\d+) h=(\d+) rgs=(\d+) classes=(\d+) F=(\d+)")
    seen = set()
    bad = 0
    for f in sorted(glob.glob(root + "/runs/*.log") + glob.glob(root + "/runs/*/*.log")):
        for line in open(f):
            m = pat.search(line)
            if not m: continue
            k, h, rg, cl, F = map(int, m.groups())
            key = (k, h, rg, cl, F)
            if key in seen: continue
            seen.add(key)
            o, a, b = orbits(k, h), nrgs(2 * h, k), nrgs(2 * h - 1, k)
            ok = (o == cl and a == rg and b == F)
            bad += not ok
            print(f"{f.split('/')[-1]:>14}: k={k:2d} h={h:2d} classes {cl} burnside {o} | rgs {rg} stirling {a} | F {F} stirling {b}  {'OK' if ok else 'MISMATCH'}")
    # the spec's list of known counts
    spec = {(2,3):20,(2,6):1056,(3,3):70,(3,5):5002,(4,4):1451,(5,4):1993,(6,5):55091,(7,4):2151,(10,3):117,(3,7):399310,(4,6):350891,(6,7):57365191}
    for (k, h), c in spec.items():
        o = orbits(k, h); bad += (o != c)
        print(f"spec k={k} h={h}: claimed {c} burnside {o} {'OK' if o == c else 'MISMATCH'}")
    # projected big runs
    for k, h in [(5,8),(4,8),(4,9),(3,9),(3,10),(3,11),(13,7),(14,7),(15,7),(30,2),(30,20)]:
        print(f"projection k={k} h={h}: classes {orbits(k,h)} rgs {nrgs(2*h,k)} F {nrgs(2*h-1,k)}")
    print("mismatches:", bad)
