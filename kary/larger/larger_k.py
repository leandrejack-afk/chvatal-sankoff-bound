"""Exploratory (not certified): does Lueker's M_h beat Dancik's 1994 bounds for k >= 5?
Uses the review implementation (../review/rv_machine.py) for the automaton and weights,
float power iteration for the spectral radius, golden-section search on z."""
import math, sys, time
sys.path.insert(0, '../review')
import numpy as np
import rv_machine as M

DANCIK = {3: 0.765803, 4: 0.708236, 5: 0.664428, 6: 0.629316, 7: 0.600184, 8: 0.575407,
          9: 0.553937, 10: 0.534855, 11: 0.517842, 12: 0.502591, 13: 0.488800, 14: 0.476198, 15: 0.464619}

def setup(k, h):
    states, rows, cl = M.build(k, h, sym=True)
    keys = [M.wkey(k, h, c) for c in cl]
    ukeys = sorted(set(keys))
    kidx = {kk: i for i, kk in enumerate(ukeys)}
    src, dst, wt = [], [], []
    for s, row in enumerate(rows):
        for t, ci in row:
            src.append(s); dst.append(t); wt.append(kidx[keys[ci]])
    return len(states), len(cl), np.array(src), np.array(dst), np.array(wt), ukeys

def rho(k, h, z, data, iters=3000):
    n, _, src, dst, wt, ukeys = data
    w = np.array([float(M.weight(k, kk, z)) for kk in ukeys])
    ew = w[wt]
    v = np.ones(n)
    lo = hi = 0.0
    for _ in range(iters):
        nv = np.zeros(n)
        np.add.at(nv, src, ew * v[dst])
        r = nv / v
        lo, hi = r.min(), r.max()
        v = nv / nv.max()
        if hi / lo < 1 + 1e-13:
            break
    return lo, hi

def bound(k, z, lam):
    return 2 * math.log(k * z) / math.log(lam)

def best(k, h, data):
    a, b = 0.45 / k, 0.95 / k
    g = (math.sqrt(5) - 1) / 2
    f = lambda z: bound(k, z, rho(k, h, z, data)[1])
    c, d = b - g * (b - a), a + g * (b - a)
    fc, fd = f(c), f(d)
    while b - a > 1e-6:
        if fc < fd: b, d, fd = d, c, fc; c = b - g * (b - a); fc = f(c)
        else: a, c, fc = c, d, fd; d = a + g * (b - a); fd = f(d)
    z = (a + b) / 2
    return z, f(z)

for h in (2, 3):
    for k in range(3, 16):
        t0 = time.time()
        try:
            data = setup(k, h)
        except MemoryError:
            print(k, h, 'memory'); break
        z, B = best(k, h, data)
        d = DANCIK[k]
        print(f"k={k:2d} h={h} states={data[0]:6d} classes={data[1]:6d} z={z:.5f} bound={B:.6f} dancik={d:.6f} {'BEATS' if B < d else 'no'} by {d-B:+.6f} ({time.time()-t0:.0f}s)", flush=True)
        if time.time() - t0 > 600:
            print('stopping this h: too slow'); break
