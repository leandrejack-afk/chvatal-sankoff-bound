"""Referee checks on rv_machine.py.

  python rv_checks.py mmp K N        minimal match pairs: definition vs closed form, all pairs of size <= N
  python rv_checks.py weights K H N  class weights: brute force (definition) vs closed form, sizes <= N
  python rv_checks.py build K H SYM Z [Z...]   state count, accepted transitions, spectral radius at each z
"""
import itertools
import sys
import time
from fractions import Fraction

import numpy as np
import scipy.sparse as sp

import rv_machine as M


def lcs(a, b):
    prev = [0] * (len(b) + 1)
    for i in range(1, len(a) + 1):
        cur = [0] * (len(b) + 1)
        for j in range(1, len(b) + 1):
            if a[i - 1] == b[j - 1]:
                cur[j] = prev[j - 1] + 1
            else:
                cur[j] = max(cur[j - 1], prev[j])
        prev = cur
    return prev[-1]


def is_mmp_def(a, b):
    """Lueker's definition: match pair, ends on an essential match, lcs = 1."""
    if not a or not b or a[-1] != b[-1]:
        return False
    l = lcs(a, b)
    return l == 1 and lcs(a[:-1], b) < l and lcs(a, b[:-1]) < l


def is_mmp_closed(a, b):
    if not a or not b or a[-1] != b[-1]:
        return False
    c = a[-1]
    A, B = set(a[:-1]), set(b[:-1])
    return c not in A and c not in B and not (A & B)


def all_strings(k, maxlen):
    for L in range(maxlen + 1):
        yield from itertools.product(range(k), repeat=L)


def check_mmp(k, N):
    n = bad = good = 0
    for na in range(N + 1):
        for a in itertools.product(range(k), repeat=na):
            for nb in range(N - na + 1):
                for b in itertools.product(range(k), repeat=nb):
                    n += 1
                    d = is_mmp_def(a, b)
                    good += d
                    if d != is_mmp_closed(a, b):
                        bad += 1
    print(f"mmp k={k} size<={N}: {n} pairs, {good} minimal match pairs, {bad} disagreements")


class Ser:
    """Truncated power series with Fraction coefficients, order < N."""
    N = 10

    def __init__(self, c):
        self.c = (list(c) + [Fraction(0)] * Ser.N)[:Ser.N]

    @staticmethod
    def lift(o):
        return o if isinstance(o, Ser) else Ser([Fraction(o)])

    def __add__(self, o):
        o = Ser.lift(o)
        return Ser([p + q for p, q in zip(self.c, o.c)])

    __radd__ = __add__

    def __neg__(self):
        return Ser([-p for p in self.c])

    def __sub__(self, o):
        return self + (-Ser.lift(o))

    def __rsub__(self, o):
        return Ser.lift(o) + (-self)

    def __mul__(self, o):
        o = Ser.lift(o)
        r = [Fraction(0)] * Ser.N
        for i, p in enumerate(self.c):
            if p:
                for j in range(Ser.N - i):
                    r[i + j] += p * o.c[j]
        return Ser(r)

    __rmul__ = __mul__

    def inv(self):
        assert self.c[0] == 1
        r = [Fraction(0)] * Ser.N
        r[0] = Fraction(1)
        for n in range(1, Ser.N):
            r[n] = -sum(self.c[i] * r[n - i] for i in range(1, n + 1))
        return Ser(r)

    def __truediv__(self, o):
        return self * Ser.lift(o).inv()

    def __rtruediv__(self, o):
        return Ser.lift(o) * self.inv()

    def __pow__(self, e):
        r = Ser([Fraction(1)])
        for _ in range(e):
            r = r * self
        return r


def check_weights(k, h, N):
    Ser.N = N + 1
    z = Ser([0, 1])
    cl = M.classes(k, h)
    idx = {c: i for i, c in enumerate(cl)}
    tally = [[0] * (N + 1) for _ in cl]
    npairs = 0
    for na in range(1, N):
        for a in itertools.product(range(k), repeat=na):
            for nb in range(1, N - na + 1):
                for b in itertools.product(range(k), repeat=nb):
                    if is_mmp_def(a, b):
                        npairs += 1
                        tally[idx[(a[-h:], b[-h:])]][na + nb] += 1
    bad = 0
    for i, c in enumerate(cl):
        w = M.weight(k, M.wkey(k, h, c), z)
        for n in range(N + 1):
            if w.c[n] != tally[i][n]:
                bad += 1
                if bad < 10:
                    print("MISMATCH", c, n, w.c[n], tally[i][n])
    print(f"weights k={k} h={h} size<={N}: {len(cl)} classes, {npairs} minimal match pairs, {bad} mismatches")


def spectral(states, rows, cl, k, h, z):
    n = len(states)
    wk = [M.wkey(k, h, c) for c in cl]
    wcache = {}
    r, cidx, val = [], [], []
    for s, row in enumerate(rows):
        for t, ci in row:
            key = wk[ci]
            if key not in wcache:
                wcache[key] = float(M.weight(k, key, Fraction(z)))
            r.append(s)
            cidx.append(t)
            val.append(wcache[key])
    A = sp.csr_matrix((val, (r, cidx)), shape=(n, n))
    v = np.ones(n)
    for it in range(20000):
        w = A @ v
        ratio = w / v
        hi, lo = ratio.max(), ratio.min()
        v = w / w.max()
        if hi / lo - 1 < 1e-13:
            break
    return lo, hi, it


def main():
    cmd = sys.argv[1]
    if cmd == "mmp":
        check_mmp(int(sys.argv[2]), int(sys.argv[3]))
    elif cmd == "weights":
        check_weights(int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
    elif cmd == "build":
        k, h, sym = int(sys.argv[2]), int(sys.argv[3]), sys.argv[4] == "1"
        t0 = time.time()
        states, rows, cl = M.build(k, h, sym)
        acc = sum(len(r) for r in rows)
        rec = sum(len({(t, M.wkey(k, h, cl[ci])) for t, ci in r}) for r in rows)
        print(f"k={k} h={h} sym={int(sym)}: states {len(states)}, classes {len(cl)}, accepted {acc}, "
              f"records {rec}  ({time.time() - t0:.1f}s)", flush=True)
        for zs in sys.argv[5:]:
            lo, hi, it = spectral(states, rows, cl, k, h, float(zs))
            print(f"   z={zs}: rho in [{lo:.17g}, {hi:.17g}] after {it} iterations", flush=True)


if __name__ == "__main__":
    main()
