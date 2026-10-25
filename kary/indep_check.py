"""Independent end-to-end check for small (k, h), without kary.c or relabelling.

Breadth-first search over the raw states of M_h (no relabelling) with kref.delta
and the minimal match pairs found from the definition; the weighted matrix A_z
with the closed-form weights (kweights.closed_mp); its spectral radius by
power iteration in numpy; the bound 2 log(kz)/log(rho).  Compare with kary.c.

Usage: python indep_check.py <k> <h> <z>
"""
import math
import sys

import numpy as np
import scipy.sparse as sp

import kref
import kweights

k, h, z = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
inputs = kref.inputs_by_definition(k, h)
w = [float(kweights.closed_mp(k, kref.weight_descriptor(k, h, a, b), kweights.mpmath.mpf(z))) for a, b in inputs]
idx = {kref.start(h): 0}
order = [kref.start(h)]
rows, cols, vals = [], [], []
s = 0
while s < len(order):
    q = order[s]
    for (a, b), wt in zip(inputs, w):
        r = kref.delta(h, q, a, b)
        if r is None:
            continue
        if r not in idx:
            idx[r] = len(order)
            order.append(r)
        rows.append(s); cols.append(idx[r]); vals.append(wt)
    s += 1
n = len(order)
A = sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
canon_states = len({kref.canon(q) for q in order})
v = np.ones(n)
for it in range(5000):
    u = A @ v
    r = u / v
    lo, hi = r.min(), r.max()
    v = u / u.max()
    if hi / lo < 1 + 1e-13:
        break
print(f"k={k} h={h} z={z}: raw states {n} (relabelling orbits {canon_states}), "
      f"rho in [{lo:.15f}, {hi:.15f}], bound {2 * math.log(k * z) / math.log(hi):.10f}")
