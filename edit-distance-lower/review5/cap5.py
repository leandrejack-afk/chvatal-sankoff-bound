"""Review of rule 5: (1) capping direction, end to end: the optimum a over (x, t) from the closed
form (opt5.py) for R = 1, 2, 4, 8, 32, oo must be non-decreasing in R (a smaller cap over-counts
more, so it may only give a smaller bound); (2) the capped series dominate the uncapped one
coefficientwise (two-counter DP with R = 1, 2 vs R = None); (3) the author's count_check5.automaton
equals the reviewer's two-counter DP.
usage: cap5.py"""
import importlib.util
import sys

sys.dont_write_bytecode = True
import mpmath  # noqa: E402

from opt5 import optimum  # noqa: E402
from table5 import two_counter  # noqa: E402

HERE = __file__.rsplit("/", 1)[0]
spec = importlib.util.spec_from_file_location("count_check5", HERE + "/../count_check5.py")
cc5 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cc5)

mpmath.mp.dps = 30
for k in (2, 4, 32):
    vals = []
    for R in (1, 2, 4, 8, 32, None):
        t, a = optimum(k, R)
        vals.append(a)
    mono = all(vals[i] <= vals[i + 1] + mpmath.mpf(10) ** -25 for i in range(len(vals) - 1))
    print(f"k={k}: optimum a for R = 1, 2, 4, 8, 32, oo: " + ", ".join(mpmath.nstr(v, 10) for v in vals)
          + f"; non-decreasing: {mono}", flush=True)
for k, n, m in ((2, 6, 6), (2, 3, 8), (3, 4, 5), (4, 3, 5)):
    exact = two_counter(k, n, m)
    for R in (1, 2):
        cap = two_counter(k, n, m, R)
        dom = all(cap.get(c, 0) >= exact.get(c, 0) for c in set(cap) | set(exact))
        strict = any(cap.get(c, 0) > exact.get(c, 0) for c in set(cap) | set(exact))
        print(f"k={k} n={n} m={m} R={R}: capped >= uncapped at every cost: {dom} (strictly somewhere: {strict})")
for k, n in ((2, 5), (2, 6), (3, 4), (4, 4), (5, 3)):
    print(f"k={k} n={n}: author's count_check5.automaton == reviewer's two_counter: {cc5.automaton(k, n) == two_counter(k, n, n)}")
