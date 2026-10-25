"""Cross-check of sampled states for large k without rebuilding the whole automaton.

For h = 3 the reduced automaton has the same 316 states (as tuples) for every k >= 5.
This takes the state list from the k = 5 build, checks that it is closed under the
transitions for the given k on the sampled states, and compares, for each sampled state,
the multiset of (target orbit, input class) of the review implementation (rv_machine.py)
with that of the reference implementation (kref.py, classes by definition).

Usage: python cross_check_fast.py <k> <sample size>
"""
import random
import sys
import time

sys.path.insert(0, "../review")
sys.path.insert(0, "..")
import kref
import rv_machine as M

k, sample = int(sys.argv[1]), int(sys.argv[2])
h = 3
t0 = time.time()
states5, _, _ = M.build(5, h, sym=True)
cl = M.classes(k, h)
ref_classes = kref.inputs_by_definition(k, h)
assert sorted(cl) == sorted(ref_classes), "class lists differ"
orbits = {kref.canon(q) for q in states5}
idx = random.Random(k).sample(range(len(states5)), sample)
bad = checked = 0
for s in idx:
    q = states5[s]
    mine = []
    for p in cl:
        t = M.delta(h, q, p)
        if t is not None:
            o = kref.canon(t)
            assert o in orbits, "rv leaves the 316-state set"
            mine.append((o, p))
    ref = []
    for a, b in ref_classes:
        t = kref.delta(h, q, a, b)
        if t is not None:
            ref.append((kref.canon(t), (a, b)))
    mine.sort(); ref.sort()
    checked += len(ref)
    bad += mine != ref
print(f"k={k} h=3: {len(cl)} classes (same as by definition), sampled {sample} of the "
      f"{len(states5)} states, {checked} accepted transitions against kref: {bad} mismatching "
      f"states, all targets inside the state set ({time.time() - t0:.0f}s)", flush=True)
