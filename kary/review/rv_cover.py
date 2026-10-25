"""Covering test and deep weight test for rv_machine.py.

  python rv_cover.py cover K H NPAIRS SEED
     random string pairs -> canonical decomposition from the full DP table (Lueker 3.2 rules)
     -> every piece a minimal match pair (by definition), tail a null pair -> raw M_h must accept.
  python rv_cover.py wdeep K H E
     for one class per weight key: enumerate every extension (u alpha, v beta) with |u|+|v| <= E
     (u only if |alpha| = h, v only if |beta| = h), keep those that are minimal match pairs by
     definition, tally by size, compare with the closed-form power series.
"""
import itertools
import random
import sys

import rv_machine as M
from rv_checks import Ser, is_mmp_def, lcs


def canonical_decomposition(a, b):
    n, m = len(a), len(b)
    T = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            T[i][j] = T[i - 1][j - 1] + 1 if a[i - 1] == b[j - 1] else max(T[i - 1][j], T[i][j - 1])
    i, j = n, m
    cuts = []
    while i > 0 and j > 0:
        if T[i][j] == T[i - 1][j]:
            i -= 1
        elif T[i][j] == T[i][j - 1]:
            j -= 1
        else:
            cuts.append((i, j))
            i -= 1
            j -= 1
    cuts.reverse()
    pieces, pi, pj = [], 0, 0
    for (ci, cj) in cuts:
        pieces.append((a[pi:ci], b[pj:cj]))
        pi, pj = ci, cj
    return pieces, (a[pi:], b[pj:]), T[n][m]


def rand_string(k, rng):
    mode = rng.random()
    L = rng.randint(0, 40)
    if mode < 0.4:
        return tuple(rng.randrange(k) for _ in range(L))
    letters = rng.sample(range(k), rng.randint(1, k))
    out = []
    while len(out) < L:  # long runs so that mu truncates
        out += [rng.choice(letters)] * rng.randint(1, 9)
    return tuple(out[:L])


def cover(k, h, npairs, seed):
    rng = random.Random(seed)
    steps = trunc = 0
    for _ in range(npairs):
        a, b = rand_string(k, rng), rand_string(k, rng)
        pieces, tail, l = canonical_decomposition(a, b)
        assert len(pieces) == l and lcs(a, b) == l
        assert all(is_mmp_def(x, y) for x, y in pieces), (a, b)
        assert lcs(*tail) == 0
        q = M.start_state(h)
        for (x, y) in pieces:
            p = (x[-h:], y[-h:])
            trunc += len(x) > h or len(y) > h
            q = M.delta(h, q, p)
            steps += 1
            if q is None:
                print("REJECTED canonical list", a, b)
                return
    print(f"cover k={k} h={h}: {npairs} pairs, {steps} steps, {trunc} truncated inputs, 0 rejections")


def wdeep(k, h, E):
    Ser.N = 2 * h + E + 1
    z = Ser([0, 1])
    reps = {}
    for c in M.classes(k, h):
        reps.setdefault(M.wkey(k, h, c), c)
    bad = tested = 0
    for key, (al, be) in reps.items():
        la, lb = len(al), len(be)
        if la < h and lb < h:
            continue
        tally = [0] * Ser.N
        for e in range(E + 1):
            for eu in range(e + 1):
                ev = e - eu
                if (eu and la < h) or (ev and lb < h):
                    continue
                for u in itertools.product(range(k), repeat=eu):
                    for v in itertools.product(range(k), repeat=ev):
                        if is_mmp_def(u + al, v + be):
                            tally[la + lb + e] += 1
        w = M.weight(k, key, z)
        for n in range(la + lb + E + 1):
            tested += 1
            if w.c[n] != tally[n]:
                bad += 1
                print("MISMATCH", key, n, w.c[n], tally[n])
    print(f"wdeep k={k} h={h} E={E}: {len(reps)} weight keys, {tested} coefficients, {bad} mismatches")


if __name__ == "__main__":
    if sys.argv[1] == "cover":
        cover(*map(int, sys.argv[2:6]))
    else:
        wdeep(*map(int, sys.argv[2:5]))
