"""Reference implementation, written separately from kary.c, for checking it.

- minimal match pairs by the definition in Lueker (2009), Sec. 3.1;
- the transition function of Figure 7, with a dict for the lcs table;
- relabelling to the first-occurrence form;
- reading kary.c's table file.

Strings are tuples of ints in range(k); a state is (x, y, dx, dy) with x, y
tuples of length h and dx, dy tuples of 0/1 of length h, oldest entry first.
"""
import itertools
import struct


def lcs_len(a, b):
    prev = [0] * (len(b) + 1)
    for ca in a:
        cur = [0]
        for j, cb in enumerate(b, 1):
            cur.append(prev[j - 1] + 1 if ca == cb else max(cur[j - 1], prev[j]))
        prev = cur
    return prev[-1]


def is_mmp(a, b):
    """Minimal match pair: a match pair with lcs 1 that ends on an essential match."""
    if not a or not b or a[-1] != b[-1]:
        return False
    L = lcs_len(a, b)
    return L == 1 and lcs_len(a[:-1], b) < L and lcs_len(a, b[:-1]) < L


def inputs_by_definition(k, h):
    """All mu-images: minimal match pairs whose two strings have length 1..h."""
    strs = [s for m in range(1, h + 1) for s in itertools.product(range(k), repeat=m)]
    by_last = {}
    for s in strs:
        by_last.setdefault(s[-1], []).append(s)
    out = []
    for c in range(k):
        for a in by_last[c]:
            for b in by_last[c]:
                if is_mmp(a, b):
                    out.append((a, b))
    return out


def weight_descriptor(k, h, a, b):
    """(size, ta, tb, nA, nB) normalised as in kary.c."""
    ta, tb = int(len(a) == h), int(len(b) == h)
    na, nb = len(set(a[:-1])), len(set(b[:-1]))
    if not ta and not tb:
        na = nb = 0
    elif ta and not tb:
        na = 0
    elif tb and not ta:
        nb = 0
    return (len(a) + len(b), ta, tb, na, nb)


def delta(h, q, alpha, beta):
    """Figure 7 of Lueker (2009); returns the new state or None (reject)."""
    x, y, dx, dy = q
    lcs = {(h, h): h}
    for i in range(h, 0, -1):            # lcs[i,h] - lcs[i-1,h] = dx[i-1]
        lcs[(i - 1, h)] = lcs[(i, h)] - dx[i - 1]
    for j in range(h, 0, -1):
        lcs[(h, j - 1)] = lcs[(h, j)] - dy[j - 1]
    a, b = x + tuple(alpha), y + tuple(beta)
    n, m = len(a), len(b)
    for i in range(n + 1):
        for j in range(m + 1):
            if i > h or j > h:
                if i == 0:
                    lcs[(i, j)] = lcs[(i, j - 1)]
                elif j == 0:
                    lcs[(i, j)] = lcs[(i - 1, j)]
                elif a[i - 1] == b[j - 1]:
                    lcs[(i, j)] = lcs[(i - 1, j - 1)] + 1
                else:
                    lcs[(i, j)] = max(lcs[(i, j - 1)], lcs[(i - 1, j)])
    i, j = n, m
    while i >= h and j >= h and not (i == h and j == h):
        if lcs[(i, j)] == lcs[(i - 1, j)]:
            i -= 1
        elif lcs[(i, j)] == lcs[(i, j - 1)]:
            j -= 1
        else:
            i -= 1
            j -= 1
    if (i, j) != (h, h):
        return None
    ndx = tuple(lcs[(t, m)] - lcs[(t - 1, m)] for t in range(n - h + 1, n + 1))
    ndy = tuple(lcs[(n, t)] - lcs[(n, t - 1)] for t in range(m - h + 1, m + 1))
    assert all(d in (0, 1) for d in ndx + ndy)
    return (a[-h:], b[-h:], ndx, ndy)


def canon(q):
    x, y, dx, dy = q
    lab = {}
    for c in x + y:
        if c not in lab:
            lab[c] = len(lab)
    return (tuple(lab[c] for c in x), tuple(lab[c] for c in y), dx, dy)


def start(h):
    return ((0,) * h, (0,) * h, (1,) * h, (1,) * h)


def pack(h, q):
    x, y, dx, dy = q
    vx = vy = vdx = vdy = 0
    for c in x:
        vx = vx * 4 + c
    for c in y:
        vy = vy * 4 + c
    for d in dx:
        vdx = vdx * 2 + d
    for d in dy:
        vdy = vdy * 2 + d
    return (vx << (4 * h)) | (vy << (2 * h)) | (vdx << h) | vdy


def unpack(h, key):
    dy = key & ((1 << h) - 1)
    dx = (key >> h) & ((1 << h) - 1)
    y = (key >> (2 * h)) & ((1 << (2 * h)) - 1)
    x = (key >> (4 * h)) & ((1 << (2 * h)) - 1)
    xs = tuple((x >> (2 * (h - 1 - t))) & 3 for t in range(h))
    ys = tuple((y >> (2 * (h - 1 - t))) & 3 for t in range(h))
    dxs = tuple((dx >> (h - 1 - t)) & 1 for t in range(h))
    dys = tuple((dy >> (h - 1 - t)) & 1 for t in range(h))
    return (xs, ys, dxs, dys)


def read_table(path):
    """Returns k, h, N, NI, wtypes, offsets, records (list of (tgt, wt, cnt) as arrays), keys."""
    import numpy as np
    with open(path, "rb") as f:
        assert f.read(4) == b"KLUE"
        k, h = struct.unpack("<ii", f.read(8))
        N, R = struct.unpack("<QQ", f.read(16))
        NI, NW = struct.unpack("<ii", f.read(8))
        wt = [struct.unpack("<5i", f.read(20)) for _ in range(NW)]
        off = np.frombuffer(f.read(8 * (N + 1)), dtype="<u8")
        rec = np.frombuffer(f.read(8 * R), dtype=[("tgt", "<u4"), ("wt", "<u2"), ("cnt", "<u2")])
        keys = np.frombuffer(f.read(8 * N), dtype="<u8")
    return k, h, N, NI, wt, off, rec, keys
