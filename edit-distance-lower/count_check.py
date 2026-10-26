"""Validate the probability weights: for small k, n, compare (a) the brute-force average over all
pairs (X, Y) in [k]^n x [k]^n of the number of canonical alignments, as a polynomial in t (cost),
with (b) the automaton count: sum over move sequences consuming n letters of each string of the
product of weights M = 1/k (times ((k-1)/k)^r after a run of r indels), S = t (k-1)/k, I = D = t."""
import itertools
import sys
from fractions import Fraction as Fr
from collections import defaultdict

def brute(k, n):
    total = defaultdict(int)  # cost -> count summed over pairs
    for X in itertools.product(range(k), repeat=n):
        for Y in itertools.product(range(k), repeat=n):
            # DP: state (i, j, type, runset) -> {cost: number of canonical partial alignments}
            dp = defaultdict(lambda: defaultdict(int))
            dp[(0, 0, 0, frozenset())][0] = 1
            for i in range(n + 1):
                for j in range(n + 1):
                    for (ii, jj, s, rs) in [key for key in list(dp.keys()) if key[0] == i and key[1] == j]:
                        for c, w in list(dp[(ii, jj, s, rs)].items()):
                            if i < n and j < n:
                                if X[i] == Y[j]:
                                    if not (s != 0 and X[i] in rs):
                                        dp[(i + 1, j + 1, 0, frozenset())][c] += w
                                elif s == 0:
                                    dp[(i + 1, j + 1, 0, frozenset())][c + 1] += w
                            if j < n and s != 2:
                                dp[(i, j + 1, 1, (rs if s == 1 else frozenset()) | {Y[j]})][c + 1] += w
                            if i < n and s != 1:
                                dp[(i + 1, j, 2, (rs if s == 2 else frozenset()) | {X[i]})][c + 1] += w
            for (ii, jj, s, rs), d in dp.items():
                if ii == n and jj == n:
                    for c, w in d.items():
                        total[c] += w
    npairs = k ** (2 * n)
    return {c: Fr(w, npairs) for c, w in total.items()}

def automaton(k, n):
    q = Fr(k - 1, k)
    # state (i, j, type, run length) -> {cost: weight}
    dp = defaultdict(lambda: defaultdict(Fr))
    dp[(0, 0, 0, 0)][0] = Fr(1)
    for i in range(n + 1):
        for j in range(n + 1):
            for key in [kk for kk in list(dp.keys()) if kk[0] == i and kk[1] == j]:
                _, _, s, r = key
                for c, w in list(dp[key].items()):
                    if i < n and j < n:
                        dp[(i + 1, j + 1, 0, 0)][c] += w * Fr(1, k) * (q ** r if s != 0 else 1)
                        if s == 0:
                            dp[(i + 1, j + 1, 0, 0)][c + 1] += w * q
                    if j < n and s != 2:
                        dp[(i, j + 1, 1, r + 1 if s == 1 else 1)][c + 1] += w
                    if i < n and s != 1:
                        dp[(i + 1, j, 2, r + 1 if s == 2 else 1)][c + 1] += w
    out = defaultdict(Fr)
    for (ii, jj, s, r), d in dp.items():
        if ii == n and jj == n:
            for c, w in d.items():
                out[c] += w
    return dict(out)

for k, n in [(2, 1), (2, 2), (2, 3), (2, 4), (3, 2), (3, 3), (4, 2), (4, 3)]:
    b, a = brute(k, n), automaton(k, n)
    same = set(b) == set(a) and all(b[c] == a[c] for c in b)
    print(f"k={k} n={n}: {'IDENTICAL' if same else 'DIFFERENT'}  E[#canonical alignments] by cost: {dict(sorted(a.items()))}" if same else f"k={k} n={n}: DIFFERENT\n brute {dict(sorted(b.items()))}\n auto  {dict(sorted(a.items()))}", flush=True)
