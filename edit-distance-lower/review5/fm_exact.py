"""Review of rule 5: exact first-moment check at finite n.  For all pairs in [k]^n x [k]^m, the exact
distribution of ed (edist.c) against the exact expected number of rule 1-5 alignments by cost (the
reviewer's two-counter automaton, table5.two_counter, uncapped): P(ed <= c) <= E[N_<=c] must hold for
every c. Prints the ratio E/P per c (values near 1 are where the first moment is tight).
usage: fm_exact.py k n m
"""
import subprocess
import sys
from fractions import Fraction as Fr

sys.dont_write_bytecode = True
from table5 import two_counter  # noqa: E402

k, n, m = map(int, sys.argv[1:4])
out = subprocess.run(["./edist", str(k), str(n), str(m)], capture_output=True, text=True, check=True).stdout
cnt = {int(a): int(b) for a, b in (l.split() for l in out.split("\n") if l)}
K = k ** (n + m)
E = two_counter(k, n, m)
cp, ce, ok, line = Fr(0), Fr(0), True, []
for c in range(0, n + m + 1):
    cp += Fr(cnt.get(c, 0), K)
    ce += E.get(c, Fr(0))
    if cp > 0:
        ok &= cp <= ce
        line.append(f"c={c}: {float(ce / cp):.4g}")
print(f"k={k} n={n} m={m}: P(ed<=c) <= E[N<=c] for all c: {ok}; E/P by c: " + ", ".join(line[:14]))
