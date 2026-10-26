"""Reviewer's own named transition table for the rules 1-5 automaton (independent of certify_rule5.py).

States: 'G0' (diagonal-type after M, or start), 'G1' (after S), 'MD' (M right after a deletion run),
'MI' (M right after an insertion run), and run states ('IR', r), ('DR', r) (runs after G0/MD/MI) and
('IS', r), ('DS', r) (runs right after S), r = current run length.
A transition is (dst, dx, dy, cost, prob): dx/dy letters of X/Y consumed, prob the exact factor.
Derivation (see weights_enum.py for the enumeration check):
  P(S I^r M) = q * q^r * (1/k) h(r), h(r) = 1/k + q (p/q)^r  (S letter a; each run letter avoids a;
  the match letter c = a w.p. 1/k (run letters then avoid one letter), else they avoid two);
  P(S I^r <end>) = q * q^r;  P(I^r M) = (1/k) q^r;  M = 1/k, S = q.
"""
from fractions import Fraction as Fr

DIAG = ("G0", "G1", "MD", "MI")


def qp(k):
    return Fr(k - 1, k), Fr(k - 2, k)


def h(k, r):
    q, p = qp(k)
    return Fr(1, k) + q * (p / q) ** r


def transitions(k, s, R=None):
    """Outgoing transitions of state s. With R, run lengths are capped at R (the run state R keeps
    the factor of length R: q^R, h(R)), which over-counts."""
    q, p = qp(k)
    cap = (lambda r: r) if R is None else (lambda r: min(r, R))
    out = []
    if s in DIAG:
        out.append(("G0", 1, 1, 0, Fr(1, k)))
        out.append(("G1", 1, 1, 1, q))
        if s != "MD":
            out.append((("IS", 1) if s == "G1" else ("IR", 1), 0, 1, 1, q if s == "G1" else Fr(1)))
        if s != "MI":
            out.append((("DS", 1) if s == "G1" else ("DR", 1), 1, 0, 1, q if s == "G1" else Fr(1)))
        return out
    kind, r = s
    after_s = kind in ("IS", "DS")
    ins = kind in ("IR", "IS")
    end = Fr(1, k) * (h(k, r) if after_s else q ** r)
    out.append(("MI" if ins else "MD", 1, 1, 0, end))
    out.append(((kind, cap(r + 1)), 0 if ins else 1, 1 if ins else 0, 1, q if after_s else Fr(1)))
    return out


def states(R):
    return list(DIAG) + [(kd, r) for kd in ("IR", "DR", "IS", "DS") for r in range(1, R + 1)]


def two_counter(k, n, m, R=None):
    """Sum over accepted move sequences consuming exactly n letters of X and m of Y of the product
    of the factors, by cost: {cost: Fraction}."""
    from collections import defaultdict
    layer = {(0, 0): {"G0": {0: Fr(1)}}}
    total = defaultdict(Fr)
    # process cells in order of i + j
    cells = defaultdict(lambda: defaultdict(lambda: defaultdict(Fr)))
    cells[(0, 0)]["G0"][0] = Fr(1)
    for d in range(n + m + 1):
        for i in range(max(0, d - m), min(n, d) + 1):
            j = d - i
            if (i, j) not in cells:
                continue
            for s, poly in list(cells[(i, j)].items()):
                for dst, dx, dy, c, pr in transitions(k, s, R):
                    ii, jj = i + dx, j + dy
                    if ii > n or jj > m:
                        continue
                    for cost, w in poly.items():
                        cells[(ii, jj)][dst][cost + c] += w * pr
    for s, poly in cells[(n, m)].items():
        for cost, w in poly.items():
            total[cost] += w
    return {c: w for c, w in total.items() if w != 0}


def matrix(k, x, t, R):
    """Weighted adjacency matrix with x per letter (either string), t per unit cost, run cap R."""
    st = states(R)
    idx = {s: i for i, s in enumerate(st)}
    A = [[Fr(0)] * len(st) for _ in st]
    for s in st:
        for dst, dx, dy, c, pr in transitions(k, s, R):
            A[idx[s]][idx[dst]] += x ** (dx + dy) * t ** c * pr
    return st, A
