"""Weights w_{mu,z}(class) for the k-letter machine, exactly, and a brute-force check.

A class is mu(p) = (alpha, beta) = (a'c, b'c) with |alpha|, |beta| <= h.  If
|alpha| < h the top string of p is alpha itself; if |alpha| = h, p's top
padding is u a'' with alpha = a''c and u any string such that p is still a
minimal match pair.  With A = letters(a''), B = letters(b'') and
f = k - 1 - |A| - |B| free letters:

  top truncated only:     z^(|alpha|+|beta|) / (1 - (k-1-|B|) z)
  bottom truncated only:  z^(|alpha|+|beta|) / (1 - (k-1-|A|) z)
  both truncated:         z^(|alpha|+|beta|) * sum_{s+t<=f} f!/(s! t! (f-s-t)!) g(|A|,s) g(|B|,t)
  where g(a, s) = sum_j (-1)^(s-j) C(s,j) / (1 - (a+j) z)
  is the generating function of strings over a+s letters using all s.

Usage:
  python kweights.py check <k> <h> <max size>     brute force against the closed form
  python kweights.py write <k> <h> <dir>          weights at the vector's z, rounded up and down
"""
import itertools
import math
import struct
import sys
from collections import defaultdict

import mpmath

import kref

mpmath.mp.dps = 80


def closed_mp(k, desc, z):
    size, ta, tb, na, nb = desc
    v = z ** size
    if ta and not tb:
        return v / (1 - (k - 1 - nb) * z)
    if tb and not ta:
        return v / (1 - (k - 1 - na) * z)
    if ta and tb:
        f = k - 1 - na - nb

        def g(a, s):
            return sum((-1) ** (s - j) * math.comb(s, j) / (1 - (a + j) * z) for j in range(s + 1))
        tot = mpmath.mpf(0)
        for s in range(f + 1):
            for t in range(f - s + 1):
                mult = math.factorial(f) // (math.factorial(s) * math.factorial(t) * math.factorial(f - s - t))
                tot += mult * g(na, s) * g(nb, t)
        return v * tot
    return v


def closed_series(k, desc, S):
    """Integer coefficients of the closed form up to z^S."""
    size, ta, tb, na, nb = desc
    out = [0] * (S + 1)
    if size > S:
        return out

    def geo(r):
        return [r ** m for m in range(S + 1)]

    def g(a, s):
        return [sum((-1) ** (s - j) * math.comb(s, j) * (a + j) ** m for j in range(s + 1)) for m in range(S + 1)]

    def mul(p, q):
        r = [0] * (S + 1)
        for i, x in enumerate(p):
            if x:
                for j in range(S + 1 - i):
                    r[i + j] += x * q[j]
        return r
    if ta and not tb:
        tail = geo(k - 1 - nb)
    elif tb and not ta:
        tail = geo(k - 1 - na)
    elif ta and tb:
        f = k - 1 - na - nb
        tail = [0] * (S + 1)
        for s in range(f + 1):
            for t in range(f - s + 1):
                mult = math.factorial(f) // (math.factorial(s) * math.factorial(t) * math.factorial(f - s - t))
                pr = mul(g(na, s), g(nb, t))
                tail = [x + mult * y for x, y in zip(tail, pr)]
    else:
        tail = [1] + [0] * S
    for m in range(S + 1 - size):
        out[size + m] = tail[m]
    return out


def mmps_constructive(k, S):
    """All minimal match pairs of total size <= S: (a'c, b'c), c not in a'b', letters disjoint."""
    for c in range(k):
        others = [t for t in range(k) if t != c]
        for la in range(S - 1):
            for a in itertools.product(others, repeat=la):
                sa = set(a)
                rest = [t for t in others if t not in sa]
                for lb in range(S - 1 - la):
                    for b in itertools.product(rest, repeat=lb):
                        yield a + (c,), b + (c,)


def check(k, h, S):
    # 1. constructive enumeration equals the definition, for small sizes
    S0 = min(S, 9 if k <= 3 else 8)
    cons = set(p for p in mmps_constructive(k, S0))
    defin = set()
    strs = [s for m in range(1, S0) for s in itertools.product(range(k), repeat=m)]
    for a in strs:
        for b in strs:
            if len(a) + len(b) <= S0 and kref.is_mmp(a, b):
                defin.add((a, b))
    assert cons == defin, (len(cons), len(defin))
    print(f"k={k}: constructive minimal match pairs = definition up to size {S0} ({len(cons)} pairs)")
    # 2. per class, counts by size against the closed-form series
    counts = defaultdict(lambda: [0] * (S + 1))
    n = 0
    for a, b in mmps_constructive(k, S):
        cls = (a[-h:], b[-h:])
        counts[cls][len(a) + len(b)] += 1
        n += 1
    classes = kref.inputs_by_definition(k, h)
    assert set(counts) == set(classes), "class sets differ"
    bad = 0
    for cls in classes:
        desc = kref.weight_descriptor(k, h, *cls)
        if counts[cls] != closed_series(k, desc, S):
            bad += 1
            if bad <= 3:
                print("MISMATCH", cls, desc, counts[cls], closed_series(k, desc, S))
    print(f"k={k} h={h}: {len(classes)} classes, {n} minimal match pairs of size <= {S}: {bad} mismatches")
    # 3. the mpmath closed form against its own series at a small z
    z = mpmath.mpf(1) / (4 * k)
    for cls in classes[:: max(1, len(classes) // 50)]:
        desc = kref.weight_descriptor(k, h, *cls)
        ser = closed_series(k, desc, 120)
        approx = sum(c * z ** i for i, c in enumerate(ser))
        assert abs(approx - closed_mp(k, desc, z)) < mpmath.mpf(10) ** -40 * closed_mp(k, desc, z)
    return bad == 0


def up(x):
    fx = float(x)
    if mpmath.mpf(fx) < x:
        fx = math.nextafter(fx, math.inf)
    assert mpmath.mpf(fx) >= x
    return fx


def down(x):
    fx = float(x)
    if mpmath.mpf(fx) > x:
        fx = math.nextafter(fx, -math.inf)
    assert mpmath.mpf(fx) <= x
    return fx


def write(k, h, d):
    with open(f"{d}/kvec_k{k}_h{h}.bin", "rb") as f:
        kk, hh = struct.unpack("<ii", f.read(8))
        (n,) = struct.unpack("<Q", f.read(8))
        (z,) = struct.unpack("<d", f.read(8))
    assert (kk, hh) == (k, h) and 0 < z < 1.0 / k
    zr = mpmath.mpf(z)
    lines = []
    for line in open(f"{d}/wtypes_k{k}_h{h}.txt"):
        i, size, ta, tb, na, nb = map(int, line.split())
        w = closed_mp(k, (size, ta, tb, na, nb), zr)
        assert w > 0
        lines.append(f"{float.hex(up(w))} {float.hex(down(w))}")
    path = f"{d}/weights_k{k}_h{h}.txt"
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"wrote {path}: {len(lines)} weight types at z = {z!r}")


if __name__ == "__main__":
    if sys.argv[1] == "check":
        ok = check(int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
        sys.exit(0 if ok else 1)
    write(int(sys.argv[2]), int(sys.argv[3]), sys.argv[4])
