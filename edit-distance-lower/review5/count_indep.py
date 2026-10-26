"""Review of rule 5 (independent): brute force over ALL alignments of ALL pairs, rules checked by a
direct predicate on the move sequence (weights_enum.obeys), against the reviewer's own two-counter
automaton (table5.two_counter).  Also, per pair: min cost over rule-obeying alignments == ed
(existence of a canonical optimal alignment), and P(ed <= c) <= E[N_<=c] for every c.
usage: count_indep.py [k,n,m ...]
"""
import itertools
import sys
from collections import defaultdict
from fractions import Fraction as Fr

sys.dont_write_bytecode = True
from table5 import two_counter  # noqa: E402
from weights_enum import obeys  # noqa: E402


def paths(n, m):
    """All move sequences over G (diagonal), I, D from (0,0) to (n,m) with no I next to a D."""
    out = []

    def rec(i, j, seq):
        if i == n and j == m:
            out.append("".join(seq))
            return
        last = seq[-1] if seq else None
        if i < n and j < m:
            seq.append("G"); rec(i + 1, j + 1, seq); seq.pop()
        if j < m and last != "D":
            seq.append("I"); rec(i, j + 1, seq); seq.pop()
        if i < n and last != "I":
            seq.append("D"); rec(i + 1, j, seq); seq.pop()

    rec(0, 0, [])
    return out


def label(p, X, Y):
    i = j = 0
    w = []
    cost = 0
    for c in p:
        if c == "G":
            if X[i] == Y[j]:
                w.append("M")
            else:
                w.append("S"); cost += 1
            i += 1; j += 1
        elif c == "I":
            w.append("I"); j += 1; cost += 1
        else:
            w.append("D"); i += 1; cost += 1
    return "".join(w), cost


def ed(X, Y):
    n, m = len(X), len(Y)
    prev = list(range(m + 1))
    for i in range(1, n + 1):
        cur = [i] + [0] * m
        for j in range(1, m + 1):
            cur[j] = min(prev[j - 1] + (X[i - 1] != Y[j - 1]), prev[j] + 1, cur[j - 1] + 1)
        prev = cur
    return prev[m]


def run(k, n, m):
    P = paths(n, m)
    total = defaultdict(int)
    edcount = defaultdict(int)
    nobad = 0
    for X in itertools.product(range(k), repeat=n):
        for Y in itertools.product(range(k), repeat=m):
            best = None
            for p in P:
                w, c = label(p, X, Y)
                if obeys(w, X, Y):
                    total[c] += 1
                    best = c if best is None or c < best else best
            e = ed(X, Y)
            edcount[e] += 1
            if best != e:
                nobad += 1
                if nobad <= 3:
                    print(f"  NO CANONICAL OPTIMUM: X={X} Y={Y} ed={e} min rule-obeying={best}")
    K = k ** (n + m)
    brute = {c: Fr(v, K) for c, v in total.items()}
    auto = two_counter(k, n, m)
    same = brute == auto
    # first moment at every cost level
    fm_ok = True
    cum_p, cum_e = Fr(0), Fr(0)
    worst = None
    for c in range(0, n + m + 1):
        cum_p += Fr(edcount.get(c, 0), K)
        cum_e += auto.get(c, Fr(0))
        if cum_p > cum_e:
            fm_ok = False
        if cum_p > 0:
            r = cum_e / cum_p
            worst = r if worst is None or r < worst else worst
    print(f"k={k} n={n} m={m}: {len(P)} paths x {K} pairs; brute == own automaton: {same}; "
          f"pairs without canonical optimum: {nobad}; P(ed<=c) <= E[N<=c] for all c: {fm_ok} "
          f"(min ratio E/P = {float(worst):.4f})", flush=True)
    if not same:
        print("  brute", dict(sorted(brute.items())))
        print("  auto ", dict(sorted(auto.items())))
    return same and nobad == 0 and fm_ok


if __name__ == "__main__":
    cases = [(2, 1, 1), (2, 2, 2), (2, 3, 3), (2, 4, 4), (2, 5, 5), (2, 2, 6), (2, 3, 6), (2, 1, 7),
             (2, 4, 6), (3, 2, 2), (3, 3, 3), (3, 4, 4), (3, 2, 5), (3, 1, 5), (4, 2, 2), (4, 3, 3),
             (4, 2, 4), (5, 2, 3), (5, 3, 3), (6, 2, 3)]
    if len(sys.argv) > 1:
        cases = [tuple(map(int, a.split(","))) for a in sys.argv[1:]]
    allok = True
    for k, n, m in cases:
        allok &= run(k, n, m)
    print("ALL OK" if allok else "FAILURES")
