"""Facts about a witness, from the spec's class definition only: witness_facts.py k h witness.bin

Prints the value range, where the zero sits, and whether match-pair values ignore the shared
first letter (Lueker's third identity, which holds exactly for any vector of the form T(w') - c).
"""
import itertools
import struct
import sys

k, h, path = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
n = 2 * h
rgs, z = [], [0] * n


def rec(i, used):
    if i == n:
        rgs.append(tuple(z))
        return
    for v in range(min(used + 1, k)):
        z[i] = v
        rec(i + 1, max(used, v + 1))


rec(1, 1)
rank = {x: i for i, x in enumerate(rgs)}


def relabel(seq):
    seen = {}
    return tuple(seen.setdefault(c, len(seen)) for c in seq)


canon = [i for i, x in enumerate(rgs) if i <= rank[relabel(x[h:] + x[:h])]]
dense = {r: d for d, r in enumerate(canon)}
raw = open(path, "rb").read()
w = struct.unpack(f"<{len(raw) // 4}I", raw)


def W(s, t):
    return w[dense[min(rank[relabel(s + t)], rank[relabel(t + s)])]]


zero = tuple([0] * h + [1] * h)
zero_class = dense[min(rank[relabel(zero)], rank[relabel(zero[h:] + zero[:h])])]
bad = total = 0
if k ** (2 * h - 2) <= 200000:
    for s1 in itertools.product(range(k), repeat=h - 1):
        for t1 in itertools.product(range(k), repeat=h - 1):
            total += 1
            bad += len({W((a,) + s1, (a,) + t1) for a in range(k)}) > 1
    identity = f"match-pair values that depend on the shared first letter: {bad} of {total}"
else:
    identity = "first-letter identity: skipped (too many suffix pairs)"
lo, hi = min(w), max(w)
print(f"k={k} h={h}: min {lo} at classes {[i for i, v in enumerate(w) if v == lo]} "
      f"(class of (0^h,1^h) is {zero_class}), max {hi} = {hi / 2**27:.3f} units of 2^27; {identity}")
