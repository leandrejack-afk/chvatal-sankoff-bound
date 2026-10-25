"""Covering test for the k-letter machine M_h (Lueker 2009, Lemma 3.5 and Sec. 3.2).

For random string pairs over k letters: compute the full lcs table, the
canonical backtracking (up, else left, else diagonal), cut after each match;
check that the pieces are minimal match pairs followed by an optional null
pair; feed mu(piece) to kref.delta from the start state (no relabelling) and
check that the machine never rejects.  Some pairs use few letters and long
runs, so that pieces longer than h (the truncation mu) occur often.

Usage: python kcover.py <k> <h> <pairs> [seed]
"""
import random
import sys

import kref


def canonical_pieces(a, b):
    n, m = len(a), len(b)
    T = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            T[i][j] = T[i - 1][j - 1] + 1 if a[i - 1] == b[j - 1] else max(T[i - 1][j], T[i][j - 1])
    i, j, cuts = n, m, []
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
    for (i, j) in cuts:
        pieces.append((tuple(a[pi:i]), tuple(b[pj:j])))
        pi, pj = i, j
    tail = (tuple(a[pi:]), tuple(b[pj:]))
    return pieces, tail, T[n][m]


def rand_string(k, n, rng, mode):
    if mode == 0:
        return [rng.randrange(k) for _ in range(n)]
    letters = rng.sample(range(k), rng.randint(1, k))
    out = []
    while len(out) < n:
        out += [rng.choice(letters)] * rng.randint(1, 3 * mode + 2)
    return out[:n]


def main():
    k, h, pairs = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    rng = random.Random(int(sys.argv[4]) if len(sys.argv) > 4 else 7)
    steps = truncated = 0
    for t in range(pairs):
        mode = t % 4
        a = rand_string(k, rng.randint(1, 60), rng, mode)
        b = rand_string(k, rng.randint(1, 60), rng, mode)
        pieces, tail, L = canonical_pieces(a, b)
        assert len(pieces) == L
        assert kref.lcs_len(tail[0], tail[1]) == 0
        q = kref.start(h)
        for (pa, pb) in pieces:
            assert kref.is_mmp(pa, pb), (a, b, pa, pb)
            al, be = pa[-h:], pb[-h:]
            truncated += (len(pa) > h) or (len(pb) > h)
            q = kref.delta(h, q, al, be)
            steps += 1
            if q is None:
                print("REJECTED a canonical decomposition:", a, b)
                sys.exit(1)
    print(f"k={k} h={h}: {pairs} pairs, {steps} canonical steps ({truncated} truncated inputs), no rejection")


if __name__ == "__main__":
    main()
