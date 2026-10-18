"""Moore partition refinement on the accepting states of M_h (rejecting state = class 0).
States in one class accept the same language of inputs, so they have the same
f_l(z, q) for every l and z; the quotient automaton has the same lambda(z).
Usage: python minimize.py <trans file>"""
import sys, struct, numpy as np, time
f = open(sys.argv[1], 'rb'); h, n = struct.unpack('<ii', f.read(8)); nt = 2 * (2 * h - 1)
T = np.fromfile(f, dtype=np.uint32, count=n * nt).reshape(n, nt).astype(np.int64)
cls = np.ones(n + 1, dtype=np.int64); cls[0] = 0
k = 1; t0 = time.time(); rounds = 0
while True:
    sig = np.concatenate([cls[1:, None], cls[T]], axis=1)
    _, inv = np.unique(sig, axis=0, return_inverse=True)
    newcls = np.zeros(n + 1, dtype=np.int64); newcls[1:] = inv.ravel() + 1
    knew = int(newcls.max()); rounds += 1
    if knew == k: break
    cls, k = newcls, knew
print(f"h={h}: {n} reachable accepting states, minimal automaton has {k} ({k/n:.3f}), {rounds} rounds, {time.time()-t0:.1f}s")
np.save(sys.argv[1].replace('.bin', '_classes.npy'), cls)
