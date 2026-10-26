"""Naive reference for small (k, h), written from the spec with no shared code.

naive.py k h witness.bin              r, argmin class, max of T - w
naive.py k h witness.bin --class D    neighbours of class D with their class indices and values
"""
import struct
import sys

k, h, path = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
n = 2 * h
ONE = 1 << 27

rgs = []
z = [0] * n


def rec(i, used):
    if i == n:
        rgs.append(tuple(z))
        return
    for v in range(min(used + 1, k)):
        z[i] = v
        rec(i + 1, max(used, v + 1))


rec(1, 1)
assert rgs == sorted(rgs)
rank = {x: i for i, x in enumerate(rgs)}


def relabel(seq):
    seen = {}
    return tuple(seen.setdefault(c, len(seen)) for c in seq)


def class_rank(s, t):
    return min(rank[relabel(s + t)], rank[relabel(t + s)])


canon = [i for i, x in enumerate(rgs) if i <= rank[relabel(x[h:] + x[:h])]]
dense = {r: d for d, r in enumerate(canon)}
raw = open(path, "rb").read()
assert len(raw) == 4 * len(canon), f"file has {len(raw) / 4} values, classes {len(canon)}"
w = struct.unpack(f"<{len(canon)}I", raw)


def cls(s, t):
    return dense[class_rank(s, t)]


def T_of(s, t, show=False):
    if s[0] == t[0]:
        terms = [(s[1:] + (c,), t[1:] + (c2,)) for c in range(k) for c2 in range(k)]
        if show:
            for a, b in terms:
                print("  S term", "".join(map(str, a)), "".join(map(str, b)), "class", cls(a, b), "w", w[cls(a, b)])
        return ONE + sum(w[cls(a, b)] for a, b in terms) // (k * k)
    ta = [(s[1:] + (c,), t) for c in range(k)]
    tb = [(s, t[1:] + (c,)) for c in range(k)]
    A = sum(w[cls(a, b)] for a, b in ta)
    B = sum(w[cls(a, b)] for a, b in tb)
    if show:
        for name, terms in (("A", ta), ("B", tb)):
            for a, b in terms:
                print(f"  {name} term", "".join(map(str, a)), "".join(map(str, b)), "class", cls(a, b), "w", w[cls(a, b)])
        print(f"  A = {A}, B = {B}")
    return max(A, B) // k


if len(sys.argv) > 4 and sys.argv[4] == "--class":
    d = int(sys.argv[5])
    x = rgs[canon[d]]
    s, t = x[:h], x[h:]
    T = T_of(s, t, show=True)
    print(f"class {d} rank {canon[d]} s={''.join(map(str, s))} t={''.join(map(str, t))} w={w[d]} T={T} T-w={T - w[d]}")
    sys.exit(0)

best = None
mx = None
for d, r in enumerate(canon):
    x = rgs[r]
    diff = T_of(x[:h], x[h:]) - w[d]
    if best is None or diff < best[0]:
        best = (diff, d, r)
    if mx is None or diff > mx:
        mx = diff
x = rgs[best[2]]
print(f"naive k={k} h={h}: RGS {len(rgs)}, classes {len(canon)}, r = {best[0]} at class {best[1]} rank {best[2]} "
      f"s={''.join(map(str, x[:h]))} t={''.join(map(str, x[h:]))}, max T-w = {mx}")
