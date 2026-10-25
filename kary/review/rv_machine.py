"""Referee's own implementation of Lueker's machine M_h for a k-letter alphabet.

Written from Lueker, J. ACM 56(3) 2009, Section 3 and Figure 7 only (no code from
kary/ or chvatal-sankoff-bound/ was read before this file was written).

State: (x, y, dx, dy), x, y tuples of h letters in 0..k-1, dx, dy tuples of h bits.
dx[t-1] = lcs[t, h] - lcs[t-1, h]   (last column, rows 0..h, local indices)
dy[t-1] = lcs[h, t] - lcs[h, t-1]   (last row, columns 0..h)
Input class: (alpha, beta) = mu(p), a minimal match pair with both sides of length <= h.
Orbit form ("canon"): letters renamed by first occurrence in x1 y1 x2 y2 ... xh yh
(interleaved; deliberately not the x-then-y order used in kary/).
"""
import itertools
from fractions import Fraction
from math import comb


def classes(k, h):
    """All (a'c, b'c) with |a'|,|b'| <= h-1, c not in a'b', letters(a') and letters(b') disjoint."""
    out = []
    for c in range(k):
        others = [l for l in range(k) if l != c]
        strs = [()]
        for L in range(1, h):
            strs += list(itertools.product(others, repeat=L))
        for a1 in strs:
            sa = set(a1)
            for b1 in strs:
                if sa & set(b1):
                    continue
                out.append((a1 + (c,), b1 + (c,)))
    return out


def wkey(k, h, cls):
    """Everything the weight depends on."""
    a, b = cls
    return (len(a), len(b), len(set(a[:-1])), len(set(b[:-1])), len(a) == h, len(b) == h)


def surj_gf(nfix, s, z):
    """GF of strings over (nfix fixed letters) + (s extra letters) that use every extra letter."""
    return sum((-1) ** (s - j) * comb(s, j) / (1 - (nfix + j) * z) for j in range(s + 1))


def weight(k, key, z):
    """w(alpha, beta) = sum over minimal match pairs p with mu(p) = (alpha, beta) of z^|p|.

    Derivation (referee's): a preimage is (u alpha, v beta). u may be nonempty only if
    |alpha| = h, v only if |beta| = h. The pair is minimal iff c is not in u v and
    letters(u alpha') and letters(v beta') are disjoint.
    - top long only: u is any string avoiding c and B: k-1-|B| letters.
    - both long: free letters F = Sigma - A - B - {c}, f = |F|. Sum over U = the set of
      F-letters that u uses: u ranges over strings over A+U using all of U, and then v is
      any string over B + (F - U). (Different grouping from the one in RESULTS.md.)
    """
    la, lb, nA, nB, tl, bl = key
    base = z ** (la + lb)
    if not tl and not bl:
        return base
    if tl and not bl:
        return base / (1 - (k - 1 - nB) * z)
    if bl and not tl:
        return base / (1 - (k - 1 - nA) * z)
    f = k - 1 - nA - nB
    W = 0
    for s in range(f + 1):
        W += comb(f, s) * surj_gf(nA, s, z) / (1 - (nB + f - s) * z)
    return base * W


def start_state(h):
    return ((0,) * h, (0,) * h, (1,) * h, (1,) * h)


def canon(q):
    x, y, dx, dy = q
    mp = {}
    for t in range(len(x)):
        for ch in (x[t], y[t]):
            if ch not in mp:
                mp[ch] = len(mp)
    return (tuple(mp[c] for c in x), tuple(mp[c] for c in y), dx, dy)


def delta(h, q, p):
    """Figure 7, line by line. Returns the new state or None (the rejecting state)."""
    x, y, dx, dy = q
    al, be = p
    a = x + al
    b = y + be
    n, m = len(a), len(b)
    L = [[None] * (m + 1) for _ in range(n + 1)]
    L[h][h] = h                                   # line 3
    for i in range(h, 0, -1):                     # line 4, decodeVertDiffs
        L[i - 1][h] = L[i][h] - dx[i - 1]
    for j in range(h, 0, -1):                     # decodeHorzDiffs
        L[h][j - 1] = L[h][j] - dy[j - 1]
    for i in range(n + 1):                        # lines 7-13
        for j in range(m + 1):
            if i > h or j > h:
                if i == 0:
                    L[i][j] = L[i][j - 1]
                elif j == 0:
                    L[i][j] = L[i - 1][j]
                elif a[i - 1] == b[j - 1]:
                    L[i][j] = L[i - 1][j - 1] + 1
                else:
                    L[i][j] = max(L[i][j - 1], L[i - 1][j])
    i, j = n, m                                   # lines 14-22
    while i >= h and j >= h:
        if L[i][j] == L[i - 1][j]:
            i -= 1
        elif L[i][j] == L[i][j - 1]:
            j -= 1
        else:
            i -= 1
            j -= 1
        if i == h and j == h:
            ndx = tuple(L[n - h + t][m] - L[n - h + t - 1][m] for t in range(1, h + 1))
            ndy = tuple(L[n][m - h + t] - L[n][m - h + t - 1] for t in range(1, h + 1))
            return (a[n - h:], b[m - h:], ndx, ndy)
    return None


def build(k, h, sym=True):
    """BFS from the start state. Returns (states, rows): rows[s] = list of (target, class index)."""
    cl = classes(k, h)
    s0 = start_state(h)
    if sym:
        s0 = canon(s0)
    index = {s0: 0}
    states = [s0]
    rows = []
    qi = 0
    while qi < len(states):
        q = states[qi]
        row = []
        for ci, p in enumerate(cl):
            t = delta(h, q, p)
            if t is None:
                continue
            if sym:
                t = canon(t)
            ti = index.get(t)
            if ti is None:
                ti = len(states)
                index[t] = ti
                states.append(t)
            row.append((ti, ci))
        rows.append(row)
        qi += 1
    return states, rows, cl
