"""Independent checks of the M_h automaton used for the gamma_2 upper bound.

1. covering: for random binary string pairs, the canonical LCS decomposition
   (computed by the full dynamic program) is a list of minimal match pairs
   followed by a null pair, and the automaton (without the complement
   reduction) never rejects that list.  This is the property Theorem 3.3 needs.
2. agreement: for random states of the table written by lcs_upper.c, this
   file's own implementation of delta (written from Figure 7 of Lueker's paper,
   separately from the C code) gives the same successor for every input T_{c,l}.

Usage: python check_cover.py <h> <pairs> <string_len> [trans_file] [samples]
"""
import random
import struct
import sys


def delta(h, state, c, l, reduce=True):
    """Figure 7 of Lueker (2009).  state = (x, y, dx, dy) as tuples of bits,
    most recent character last.  Input T_{c,l}: top = (1-c)^l c, bottom = c
    if l >= 0; top = c, bottom = (1-c)^{-l} c if l < 0."""
    x, y, dx, dy = state
    top = list(x) + ([1 - c] * l if l > 0 else []) + [c]
    bot = list(y) + ([1 - c] * (-l) if l < 0 else []) + [c]
    n, m = len(top), len(bot)
    INF = None
    L = [[INF] * (m + 1) for _ in range(n + 1)]
    L[h][h] = 0
    # column h, rows h-1 .. 0: dx holds successive differences L[i][h]-L[i-1][h], i=1..h
    for i in range(h, 0, -1):
        L[i - 1][h] = L[i][h] - dx[i - 1]
    for j in range(h, 0, -1):
        L[h][j - 1] = L[h][j] - dy[j - 1]
    for j in range(h + 1, m + 1):
        L[0][j] = L[0][j - 1]
    for i in range(h + 1, n + 1):
        L[i][0] = L[i - 1][0]
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if i <= h and j <= h:
                continue
            if top[i - 1] == bot[j - 1]:
                L[i][j] = L[i - 1][j - 1] + 1
            else:
                L[i][j] = max(L[i - 1][j], L[i][j - 1])
    i, j, matches = n, m, 0
    while i >= h and j >= h:
        if L[i][j] == L[i - 1][j]:
            i -= 1
        elif L[i][j] == L[i][j - 1]:
            j -= 1
        else:
            assert L[i][j] == L[i - 1][j - 1] + 1
            i -= 1
            j -= 1
            matches += 1
        if i == h and j == h:
            assert matches == 1
            nx = tuple(top[n - h:])
            ny = tuple(bot[m - h:])
            ndx = tuple(L[t][m] - L[t - 1][m] for t in range(n - h + 1, n + 1))
            ndy = tuple(L[n][t] - L[n][t - 1] for t in range(m - h + 1, m + 1))
            assert set(ndx) <= {0, 1} and set(ndy) <= {0, 1}
            if reduce and nx[0] == 1:
                nx = tuple(1 - b for b in nx)
                ny = tuple(1 - b for b in ny)
            return (nx, ny, ndx, ndy)
    return None


def canonical_decomposition(A, B):
    n, m = len(A), len(B)
    L = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            L[i][j] = L[i - 1][j - 1] + 1 if A[i - 1] == B[j - 1] else max(L[i - 1][j], L[i][j - 1])
    i, j, matches = n, m, []
    while i > 0 and j > 0:
        if L[i][j] == L[i - 1][j]:
            i -= 1
        elif L[i][j] == L[i][j - 1]:
            j -= 1
        else:
            matches.append((i, j))
            i -= 1
            j -= 1
    matches.reverse()
    pairs, pi, pj = [], 0, 0
    for (i, j) in matches:
        pairs.append((A[pi:i], B[pj:j]))
        pi, pj = i, j
    rest = (A[pi:], B[pj:])
    return pairs, rest, L[n][m]


def to_T(top, bot):
    """Return (c, l) with top/bottom = T_{c,l}, or raise if not minimal."""
    c = top[-1]
    assert bot[-1] == c
    if len(bot) == 1:
        assert all(t == 1 - c for t in top[:-1]), (top, bot)
        return c, len(top) - 1
    assert len(top) == 1 and all(t == 1 - c for t in bot[:-1]), (top, bot)
    return c, -(len(bot) - 1)


def key_to_state(h, key):
    mask = (1 << h) - 1
    parts = [(key >> (3 * h)) & mask, (key >> (2 * h)) & mask, (key >> h) & mask, key & mask]
    return tuple(tuple((p >> (h - 1 - t)) & 1 for t in range(h)) for p in parts)


def main():
    h = int(sys.argv[1])
    pairs_n = int(sys.argv[2])
    slen = int(sys.argv[3])
    rnd = random.Random(12345 + h)
    start = (tuple([0] * h), tuple([0] * h), tuple([1] * h), tuple([1] * h))
    for _ in range(pairs_n):
        A = [rnd.randint(0, 1) for _ in range(rnd.randint(1, slen))]
        B = [rnd.randint(0, 1) for _ in range(rnd.randint(1, slen))]
        pairs, rest, ell = canonical_decomposition(A, B)
        assert len(pairs) == ell
        assert not (set(rest[0]) & set(rest[1])), "trailing pair is not null"
        q = start
        for (top, bot) in pairs:
            c, l = to_T(top, bot)
            l = max(-(h - 1), min(h - 1, l))       # mu: keep at most h characters
            q = delta(h, q, c, l, reduce=False)
            assert q is not None, ("canonical list rejected", A, B)
    print(f"h={h}: covering held on {pairs_n} random pairs (lengths <= {slen})")
    if len(sys.argv) > 4:
        path = sys.argv[4]
        samples = int(sys.argv[5]) if len(sys.argv) > 5 else 2000
        with open(path, "rb") as f:
            hh, n = struct.unpack("<ii", f.read(8))
            assert hh == h
            nt = 2 * (2 * h - 1)
            f.seek(8 + 4 * nt * n)
            keys = struct.unpack(f"<{n}Q", f.read(8 * n))
            for s in rnd.sample(range(1, n + 1), min(samples, n)):
                f.seek(8 + 4 * nt * (s - 1))
                row = struct.unpack(f"<{nt}I", f.read(4 * nt))
                st = key_to_state(h, keys[s - 1])
                k = 0
                for c in (0, 1):
                    for l in range(-(h - 1), h):
                        mine = delta(h, st, c, l, reduce=True)
                        theirs = row[k]
                        k += 1
                        if mine is None:
                            assert theirs == 0, (s, c, l)
                        else:
                            assert theirs != 0, (s, c, l)
                            assert key_to_state(h, keys[theirs - 1]) == mine, (s, c, l)
        print(f"h={h}: C transition table agrees with the Python delta on {min(samples, n)} random states")


if __name__ == "__main__":
    main()
