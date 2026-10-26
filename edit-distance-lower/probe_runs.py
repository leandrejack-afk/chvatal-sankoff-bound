"""Edit-distance lower bounds with run-aware canonical rules (exploration).
Rules: no I-D, D-I, I-S, D-S adjacency; after a run of r insertions (deletions), a match must use a
letter different from all r inserted (deleted) letters, probability factor ((k-1)/k)^r; runs longer
than R count as R. States: diag, I_1..I_R, D_1..D_R."""
import sys
import numpy as np
from scipy.optimize import brentq, minimize_scalar

def build(k, x, t, R):
    M, S, I = x * x / k, x * x * t * (k - 1) / k, x * t
    q = (k - 1) / k
    n = 1 + 2 * R
    A = np.zeros((n, n))
    A[0, 0] = M + S; A[0, 1] = I; A[0, 1 + R] = I
    for r in range(1, R + 1):
        i, d = r, R + r
        A[i, 0] = M * q ** r; A[d, 0] = M * q ** r
        A[i, min(r + 1, R)] = I; A[d, R + min(r + 1, R)] = I
    return A

def rho(k, x, t, R):
    return max(abs(np.linalg.eigvals(build(k, x, t, R))))

def best(k, R):
    def x_of_t(t):
        return brentq(lambda x: rho(k, x, t, R) - 1, 1e-9, 64.0)
    def neg_a(t):
        x = x_of_t(t)
        return 2 * np.log(x) / np.log(t) if x > 1 else 0.0
    grid = np.linspace(0.005, 0.995, 199)
    vals = [neg_a(t) for t in grid]
    i = int(np.argmin(vals))
    r = minimize_scalar(neg_a, bounds=(grid[max(i-1,0)], grid[min(i+1,198)]), method="bounded", options={"xatol": 1e-12})
    return -r.fun

BS = {2: 0.170552, 3: 0.283660, 4: 0.359783, 5: 0.415173, 6: 0.457766, 7: 0.491836, 8: 0.519901, 16: 0.644758, 32: 0.738677}
print("k   BS2026   R=1      R=2      R=4      R=8      R=16")
for k in (2, 3, 4, 5, 6, 7, 8, 16, 32):
    row = [best(k, R) for R in (1, 2, 4, 8, 16)]
    print(f"{k:<3} {BS[k]:.6f} " + " ".join(f"{v:.6f}" for v in row), flush=True)
