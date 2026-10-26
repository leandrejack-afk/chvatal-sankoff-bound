"""Add rule 4: no D M I and no I M D (single match between opposite indels). States: diag, MD, MI,
I_1..I_R, D_1..D_R."""
import numpy as np
from scipy.optimize import brentq, minimize_scalar

def build(k, x, t, R):
    M, S, I, q = x * x / k, x * x * t * (k - 1) / k, x * t, (k - 1) / k
    n = 3 + 2 * R; A = np.zeros((n, n))
    DG, MD, MI = 0, 1, 2
    i0, d0 = 3, 3 + R  # I_r at i0 + r - 1, D_r at d0 + r - 1
    for s in (DG, MD, MI):
        A[s, DG] = M + S
        if s != MD: A[s, i0] = I        # no I right after (D-run, M)
        if s != MI: A[s, d0] = I        # no D right after (I-run, M)
    for r in range(1, R + 1):
        A[i0 + r - 1, MI] = M * q ** r; A[i0 + r - 1, i0 + min(r, R - 1)] = I
        A[d0 + r - 1, MD] = M * q ** r; A[d0 + r - 1, d0 + min(r, R - 1)] = I
    return A

def best(k, R):
    rho = lambda x, t: max(abs(np.linalg.eigvals(build(k, x, t, R))))
    def x_of_t(t): return brentq(lambda x: rho(x, t) - 1, 1e-9, 64.0)
    def neg_a(t):
        x = x_of_t(t); return 2 * np.log(x) / np.log(t) if x > 1 else 0.0
    grid = np.linspace(0.005, 0.995, 199); vals = [neg_a(t) for t in grid]; i = int(np.argmin(vals))
    r = minimize_scalar(neg_a, bounds=(grid[max(i-1,0)], grid[min(i+1,198)]), method="bounded", options={"xatol": 1e-12})
    return -r.fun

prev = {2: 0.209789, 3: 0.328355, 4: 0.408774, 5: 0.468384, 6: 0.515005, 7: 0.552816, 8: 0.584300, 16: 0.726970, 32: 0.832820}
for k in (2, 3, 4, 5, 6, 7, 8, 16, 32):
    b = best(k, 24)
    print(f"k={k}: rules 1-3 {prev[k]:.6f} -> with rule 4 {b:.6f} ({b - prev[k]:+.6f})", flush=True)
