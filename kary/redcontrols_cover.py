"""Red controls: the covering test must catch a wrong machine.

(a) backtracking prefers left over up (wrong canonical order);
(b) the column differences dx are decoded in reverse order.
Each mutated machine must reject some canonical decomposition.
"""
import random
import sys

import kcover
import kref

orig = kref.delta


def delta_left_first(h, q, alpha, beta):
    x, y, dx, dy = q
    lcs = {(h, h): h}
    for i in range(h, 0, -1):
        lcs[(i - 1, h)] = lcs[(i, h)] - dx[i - 1]
    for j in range(h, 0, -1):
        lcs[(h, j - 1)] = lcs[(h, j)] - dy[j - 1]
    a, b = x + tuple(alpha), y + tuple(beta)
    n, m = len(a), len(b)
    for i in range(n + 1):
        for j in range(m + 1):
            if i > h or j > h:
                if i == 0:
                    lcs[(i, j)] = lcs[(i, j - 1)]
                elif j == 0:
                    lcs[(i, j)] = lcs[(i - 1, j)]
                elif a[i - 1] == b[j - 1]:
                    lcs[(i, j)] = lcs[(i - 1, j - 1)] + 1
                else:
                    lcs[(i, j)] = max(lcs[(i, j - 1)], lcs[(i - 1, j)])
    i, j = n, m
    while i >= h and j >= h and not (i == h and j == h):
        if lcs[(i, j)] == lcs[(i, j - 1)]:
            j -= 1
        elif lcs[(i, j)] == lcs[(i - 1, j)]:
            i -= 1
        else:
            i -= 1
            j -= 1
    if (i, j) != (h, h):
        return None
    ndx = tuple(lcs[(t, m)] - lcs[(t - 1, m)] for t in range(n - h + 1, n + 1))
    ndy = tuple(lcs[(n, t)] - lcs[(n, t - 1)] for t in range(m - h + 1, m + 1))
    return (a[-h:], b[-h:], ndx, ndy)


def delta_reversed_dx(h, q, alpha, beta):
    x, y, dx, dy = q
    return orig(h, (x, y, tuple(reversed(dx)), dy), alpha, beta)


def run(k, h, mutant, pairs=2000):
    kref.delta = mutant
    rng = random.Random(5)
    for t in range(pairs):
        a = kcover.rand_string(k, rng.randint(1, 60), rng, t % 4)
        b = kcover.rand_string(k, rng.randint(1, 60), rng, t % 4)
        pieces, tail, L = kcover.canonical_pieces(a, b)
        q = kref.start(h)
        for (pa, pb) in pieces:
            q = kref.delta(h, q, pa[-h:], pb[-h:])
            if q is None:
                kref.delta = orig
                return t + 1
    kref.delta = orig
    return None


fired = 0
for name, mut in (("left before up", delta_left_first), ("dx reversed", delta_reversed_dx)):
    for k, h in ((3, 4), (4, 3)):
        r = run(k, h, mut)
        print(f"{name}, k={k} h={h}: " + (f"caught after {r} pairs" if r else "NOT CAUGHT"))
        fired += r is not None
print(f"{fired} of 4 red controls fired")
sys.exit(0 if fired == 4 else 1)
