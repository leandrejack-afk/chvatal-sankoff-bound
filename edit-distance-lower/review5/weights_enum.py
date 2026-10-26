"""Review of rule 5 (independent): block probabilities by exhaustive enumeration of letters.

For a move sequence w (moves G = diagonal, I, D; a G is M or S according to the letters), the
probability that w is a valid alignment obeying rules 1-5 is computed here by enumerating EVERY
assignment of letters to the letters w consumes, checking the rules directly on the sequence
(no automaton), and dividing by k^(#letters).  It is compared with the product of the per-move
automaton weights (this reviewer's own derivation):
  M: 1/k; S: q;  plain run letter (after M, start, MD, MI): 1;  run letter right after S: q;
  M ending a plain run of r: (1/k) q^r;  M ending a run of r right after S: (1/k) h(r),
  h(r) = 1/k + q (p/q)^r,  q = (k-1)/k,  p = (k-2)/k.
Sequences are given as strings over M, S, I, D (the M/S label is part of the event).
usage: weights_enum.py
"""
import itertools
from fractions import Fraction as Fr


def letters_used(w):
    nx = sum(1 for c in w if c in "MSD")
    ny = sum(1 for c in w if c in "MSI")
    return nx, ny


def obeys(w, X, Y):
    """Direct check of validity and rules 1-5 for the labelled move sequence w on (X, Y)."""
    i = j = 0
    pos = []  # (move, i_before, j_before)
    for c in w:
        pos.append((c, i, j))
        if c in "MS":
            if (X[i] == Y[j]) != (c == "M"):
                return False
            i += 1
            j += 1
        elif c == "I":
            j += 1
        else:
            i += 1
    L = len(w)
    for t in range(L - 1):
        a, b = w[t], w[t + 1]
        if {a, b} == {"I", "D"}:
            return False  # rule 1
        if a in "ID" and b == "S":
            return False  # rule 2
    for t in range(1, L - 1):
        if w[t] == "M" and {w[t - 1], w[t + 1]} == {"I", "D"}:
            return False  # rule 4 (D M I or I M D)
    # maximal runs
    t = 0
    while t < L:
        if w[t] in "ID":
            kind = w[t]
            u = t
            while u < L and w[u] == kind:
                u += 1
            run = [Y[pos[s][2]] if kind == "I" else X[pos[s][1]] for s in range(t, u)]
            if u < L and w[u] == "M":
                c, ii, jj = pos[u]
                if X[ii] in run:
                    return False  # rule 3
            if t > 0 and w[t - 1] == "S":
                c, ii, jj = pos[t - 1]
                forbidden = X[ii] if kind == "I" else Y[jj]
                if forbidden in run:
                    return False  # rule 5
            t = u
        else:
            t += 1
    return True


def prob_enum(k, w):
    nx, ny = letters_used(w)
    cnt = 0
    for X in itertools.product(range(k), repeat=nx):
        for Y in itertools.product(range(k), repeat=ny):
            cnt += obeys(w, X, Y)
    return Fr(cnt, k ** (nx + ny))


def prob_auto(k, w):
    q, p = Fr(k - 1, k), Fr(k - 2, k)
    h = lambda r: Fr(1, k) + q * (p / q) ** r
    P = Fr(1)
    L = len(w)
    t = 0
    prev = None  # last diagonal label before the current run
    while t < L:
        c = w[t]
        if c == "M":
            P *= Fr(1, k)
            t += 1
            prev = "M"
        elif c == "S":
            P *= q
            t += 1
            prev = "S"
        else:
            u = t
            while u < L and w[u] == c:
                u += 1
            r = u - t
            after_s = t > 0 and w[t - 1] == "S"
            P *= q ** r if after_s else 1
            if u < L and w[u] == "M":
                P *= Fr(1, k) * (h(r) if after_s else q ** r)
                u += 1
            t = u
    return P


def structurally_ok(w):
    L = len(w)
    for t in range(L - 1):
        if {w[t], w[t + 1]} == {"I", "D"} or (w[t] in "ID" and w[t + 1] == "S"):
            return False
    for t in range(1, L - 1):
        if w[t] == "M" and {w[t - 1], w[t + 1]} == {"I", "D"}:
            return False
    return True


if __name__ == "__main__":
    bad = 0
    # 1. single blocks, explicit
    for k in (2, 3, 4, 5):
        q, p = Fr(k - 1, k), Fr(k - 2, k)
        for r in range(1, 5 if k <= 3 else 4):
            for w, closed in (("S" + "I" * r + "M", q * Fr(1, k) * (Fr(1, k) * q ** r + q * p ** r)),
                              ("S" + "D" * r + "M", q * Fr(1, k) * (Fr(1, k) * q ** r + q * p ** r)),
                              ("S" + "I" * r, q * q ** r),
                              ("S" + "D" * r, q * q ** r),
                              ("I" * r + "M", Fr(1, k) * q ** r),
                              ("M" + "I" * r + "M", Fr(1, k) * Fr(1, k) * q ** r)):
                e, a = prob_enum(k, w), prob_auto(k, w)
                ok = e == a == closed
                bad += not ok
                print(f"k={k} {w:8s} enum={str(e):>12s} automaton={str(a):>12s} closed={str(closed):>12s} {'OK' if ok else 'MISMATCH'}")
    # 2. every structurally allowed labelled sequence of length <= 6 (k=2,3) / 5 (k=4): independence across blocks
    for k, Lmax in ((2, 6), (3, 5), (4, 4)):
        n_seq = 0
        for L in range(1, Lmax + 1):
            for w in itertools.product("MSID", repeat=L):
                w = "".join(w)
                if not structurally_ok(w):
                    continue
                n_seq += 1
                e, a = prob_enum(k, w), prob_auto(k, w)
                if e != a:
                    bad += 1
                    print(f"MISMATCH k={k} w={w}: enum {e} automaton {a}")
        print(f"k={k}: all {n_seq} structurally allowed labelled sequences of length <= {Lmax}: enum == automaton product")
    # 3. h(r) facts used by the cap
    for k in (2, 3, 4, 5, 8, 32):
        q, p = Fr(k - 1, k), Fr(k - 2, k)
        h = [Fr(1, k) + q * (p / q) ** r for r in range(1, 60)]
        mono = all(h[r] >= h[r + 1] for r in range(len(h) - 1))
        print(f"k={k}: h(1) = {h[0]} == q: {h[0] == q}; h non-increasing on 1..59: {mono}; h(r) >= 1/k: {min(h) >= Fr(1, k)}")
        bad += not (mono and h[0] == q)
    print("TOTAL MISMATCHES:", bad)
