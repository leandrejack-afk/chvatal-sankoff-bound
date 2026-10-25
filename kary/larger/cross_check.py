"""Cross-check two separately written implementations of M_h for k letters.

rv_machine.py (review implementation): reduced automaton by BFS, its own canonical form
and its own class list. kref.py (reference implementation): classes found from the
definition of a minimal match pair by lcs computations, its own Figure 7 and its own
canonical form. For every state of rv's reduced automaton (or a random sample), the
multiset of (target orbit, input class) must be the same in both, and the set of orbits
reached must be the same.

Usage: python cross_check.py <k> <h> [sample size]
"""
import random
import sys
import time

sys.path.insert(0, "../review")
sys.path.insert(0, "..")
import kref
import rv_machine as M

k, h = int(sys.argv[1]), int(sys.argv[2])
sample = int(sys.argv[3]) if len(sys.argv) > 3 else 0
t0 = time.time()
states, rows, cl = M.build(k, h, sym=True)
ref_classes = kref.inputs_by_definition(k, h)
assert sorted(cl) == sorted(ref_classes), "class lists differ"
orbit = lambda q: kref.canon(q)               # map any state to kref's representative
rv_orbits = {orbit(q) for q in states}
assert len(rv_orbits) == len(states), "two rv states in one orbit"
idx = list(range(len(states)))
if sample and sample < len(idx):
    idx = random.Random(1).sample(idx, sample)
bad = 0
checked = 0
for s in idx:
    q = states[s]
    mine = sorted((orbit(states[t]), cl[ci]) for t, ci in rows[s])
    ref = []
    for a, b in ref_classes:
        nq = kref.delta(h, q, a, b)
        if nq is not None:
            o = orbit(nq)
            assert o in rv_orbits, "kref reaches an orbit rv does not have"
            ref.append((o, (a, b)))
    ref.sort()
    checked += len(ref)
    if mine != ref:
        bad += 1
print(f"k={k} h={h}: {len(cl)} classes (same as by definition), {len(states)} states, "
      f"checked {len(idx)} states and {checked} accepted transitions against kref: "
      f"{bad} mismatching states ({time.time() - t0:.0f}s)", flush=True)
