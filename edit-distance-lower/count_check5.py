"""(Rules 1-5.) Validate the probability weights: for small k, n, compare (a) the brute-force average
over all pairs (X, Y) in [k]^n x [k]^n of the number of alignments obeying rules 1 to 5, as a
polynomial in t (cost), with (b) the automaton count: the sum over move sequences consuming n letters
of each string of the product of the weights
  M = 1/k, times q^r when it ends a run of r indels that follows M (or the start), and times
      h(r) = 1/k + q (p/q)^r when it ends a run of r indels that follows S;
  S = t q;  I = D = t, times q inside a run that follows S;
with q = (k-1)/k and p = (k-2)/k.
Rule 5: after S on (a, b), the following insertion run avoids a and the following deletion run avoids
b. Each letter of such a run avoids the S letter (probability q); the match ending the run must also
avoid the run's letters: if it equals the S letter (probability 1/k) that is already true, and
otherwise (probability q) each run letter avoids one more letter, (p/q) each.
usage: count_check5.py"""
import itertools
from collections import defaultdict
from fractions import Fraction as Fr


def brute(k, n):
    """States as in canon_check5.c: 0 diagonal after M (or start), 1 I-run, 2 D-run, 3 M after a
    D-run, 4 M after an I-run, 5 diagonal after S, 6 I-run after S, 7 D-run after S."""
    total = defaultdict(int)  # cost -> number of rule-obeying alignments, summed over pairs
    for X in itertools.product(range(k), repeat=n):
        for Y in itertools.product(range(k), repeat=n):
            dp = defaultdict(lambda: defaultdict(int))
            dp[(0, 0, 0, frozenset())][0] = 1
            for i in range(n + 1):
                for j in range(n + 1):
                    for key in [kk for kk in list(dp) if kk[0] == i and kk[1] == j]:
                        _, _, s, rs = key
                        irun, drun = s in (1, 6), s in (2, 7)
                        for c, w in list(dp[key].items()):
                            if i < n and j < n:
                                if X[i] == Y[j]:
                                    if not ((irun or drun) and X[i] in rs):
                                        ns = 3 if drun else 4 if irun else 0
                                        dp[(i + 1, j + 1, ns, frozenset())][c] += w
                                elif not (irun or drun):
                                    dp[(i + 1, j + 1, 5, frozenset())][c + 1] += w
                            if j < n and not drun and s != 3:
                                after_s = s in (5, 6)
                                if not (after_s and Y[j] == X[i - 1]):
                                    ns = 6 if after_s else 1
                                    dp[(i, j + 1, ns, (rs if irun else frozenset()) | {Y[j]})][c + 1] += w
                            if i < n and not irun and s != 4:
                                after_s = s in (5, 7)
                                if not (after_s and X[i] == Y[j - 1]):
                                    ns = 7 if after_s else 2
                                    dp[(i + 1, j, ns, (rs if drun else frozenset()) | {X[i]})][c + 1] += w
            for (ii, jj, s, rs), d in dp.items():
                if ii == n and jj == n:
                    for c, w in d.items():
                        total[c] += w
    return {c: Fr(w, k ** (2 * n)) for c, w in total.items()}


def automaton(k, n):
    q, p = Fr(k - 1, k), Fr(k - 2, k)
    h = lambda r: Fr(1, k) + q * (p / q) ** r
    dp = defaultdict(lambda: defaultdict(Fr))  # (i, j, state, run length) -> {cost: weight}
    dp[(0, 0, 0, 0)][0] = Fr(1)
    for i in range(n + 1):
        for j in range(n + 1):
            for key in [kk for kk in list(dp) if kk[0] == i and kk[1] == j]:
                _, _, s, r = key
                for c, w in list(dp[key].items()):
                    if s in (0, 3, 4, 5):  # diagonal-type states
                        if i < n and j < n:
                            dp[(i + 1, j + 1, 0, 0)][c] += w / k
                            dp[(i + 1, j + 1, 5, 0)][c + 1] += w * q
                        if j < n and s != 3:
                            dp[(i, j + 1, 6 if s == 5 else 1, 1)][c + 1] += w * (q if s == 5 else 1)
                        if i < n and s != 4:
                            dp[(i + 1, j, 7 if s == 5 else 2, 1)][c + 1] += w * (q if s == 5 else 1)
                    else:  # inside a run
                        ends = q ** r if s in (1, 2) else h(r)
                        if i < n and j < n:
                            dp[(i + 1, j + 1, 4 if s in (1, 6) else 3, 0)][c] += w / k * ends
                        per = q if s in (6, 7) else 1
                        if s in (1, 6) and j < n:
                            dp[(i, j + 1, s, r + 1)][c + 1] += w * per
                        if s in (2, 7) and i < n:
                            dp[(i + 1, j, s, r + 1)][c + 1] += w * per
    out = defaultdict(Fr)
    for (ii, jj, s, r), d in dp.items():
        if ii == n and jj == n:
            for c, w in d.items():
                out[c] += w
    return dict(out)


if __name__ == "__main__":
    for k, n in [(2, 1), (2, 2), (2, 3), (2, 4), (2, 5), (2, 6), (2, 7), (3, 2), (3, 3), (3, 4), (3, 5),
                 (4, 2), (4, 3), (4, 4), (5, 3), (6, 3)]:
        b, a = brute(k, n), automaton(k, n)
        same = set(b) == set(a) and all(b[c] == a[c] for c in b)
        if same:
            print(f"k={k} n={n}: IDENTICAL  E[#alignments obeying rules 1-5] by cost: {dict(sorted(a.items()))}", flush=True)
        else:
            print(f"k={k} n={n}: DIFFERENT\n brute {dict(sorted(b.items()))}\n auto  {dict(sorted(a.items()))}", flush=True)
