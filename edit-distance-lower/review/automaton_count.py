"""Review code: expected number of canonical alignments by cost, from the move automaton with two
separate letter counters (exact Fractions), times k^(n+m), to compare with count_canon.c totals.

Independent of count_check.py / count_check2.py: the state is a named tuple-free encoding
(kind, r) with kind in {'diag', 'MD', 'MI', 'I', 'D'}; weights are probabilities per move:
M 1/k (times ((k-1)/k)^r if it ends a run of r), S (k-1)/k, I = D = 1 (no letter condition).
Rules 1-3: MD and MI behave exactly like diag.  Rules 1-4: no I from MD, no D from MI.
usage: automaton_count.py rules k n m
"""
import sys
from collections import defaultdict
from fractions import Fraction as Fr


def count(rules, k, n, m):
    q = Fr(k - 1, k)
    # layer by (i, j); value: dict state -> dict cost -> weight
    table = defaultdict(lambda: defaultdict(lambda: defaultdict(Fr)))
    table[(0, 0)][("diag", 0)][0] = Fr(1)
    for i in range(n + 1):
        for j in range(m + 1):
            cell = table.get((i, j))
            if not cell:
                continue
            for (kind, r), costs in list(cell.items()):
                for c, w in list(costs.items()):
                    if i < n and j < m:
                        if kind in ("I", "D"):
                            nk = ("MI" if kind == "I" else "MD") if rules == 4 else "diag"
                            table[(i + 1, j + 1)][(nk, 0)][c] += w * Fr(1, k) * q ** r
                        else:
                            table[(i + 1, j + 1)][("diag", 0)][c] += w * Fr(1, k)
                            table[(i + 1, j + 1)][("diag", 0)][c + 1] += w * q
                    can_I = kind != "D" and not (rules == 4 and kind == "MD")
                    can_D = kind != "I" and not (rules == 4 and kind == "MI")
                    if j < m and can_I:
                        table[(i, j + 1)][("I", r + 1 if kind == "I" else 1)][c + 1] += w
                    if i < n and can_D:
                        table[(i + 1, j)][("D", r + 1 if kind == "D" else 1)][c + 1] += w
    out = defaultdict(Fr)
    for (kind, r), costs in table[(n, m)].items():
        for c, w in costs.items():
            out[c] += w
    scale = k ** (n + m)
    res = []
    for c in range(n + m + 1):
        v = out.get(c, Fr(0)) * scale
        assert v.denominator == 1, (c, v)
        res.append(int(v))
    return res


if __name__ == "__main__":
    rules, k, n, m = map(int, sys.argv[1:5])
    print("totals " + " ".join(map(str, count(rules, k, n, m))))
