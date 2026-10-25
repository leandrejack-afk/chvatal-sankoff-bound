"""Exact certificates for Lueker's automaton M_h with k letters (small h).

For given k, h and a decimal z: build the reduced automaton with the review
implementation (../review/rv_machine.py), find a positive vector v by power
iteration in floating point, read v as exact rationals, and evaluate
(A_z v)_q / v_q for every state in exact rational arithmetic with exact class
weights. lambda = the largest ratio, so A_z v <= lambda v holds exactly; mu = the
smallest ratio is a lower bound on the spectral radius. The bound
2 log(kz)/log(lambda) is evaluated in mpmath interval arithmetic and rounded up.

Usage: python certify_exact.py <k> <h> <z decimal>
"""
import math
import sys
import time
from fractions import Fraction

sys.path.insert(0, "../review")
import mpmath
import numpy as np
from mpmath import iv

import rv_machine as M

DANCIK = {3: "0.765803", 4: "0.708236", 5: "0.664428", 6: "0.629316", 7: "0.600184",
          8: "0.575407", 9: "0.553937", 10: "0.534855", 11: "0.517842", 12: "0.502591",
          13: "0.488800", 14: "0.476198", 15: "0.464619"}


def dec(x, d, rnd):
    """The rational x as a decimal with d places, rounded by rnd (math.ceil or math.floor)."""
    n = rnd(x * 10**d)
    sign, n = ("-", -n) if n < 0 else ("", n)
    s = str(n).rjust(d + 1, "0")
    return f"{sign}{s[:-d]}.{s[-d:]}"


def main():
    k, h, zs = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    t0 = time.time()
    states, rows, cl = M.build(k, h, sym=True)
    keys = [M.wkey(k, h, c) for c in cl]
    ukeys = sorted(set(keys))
    kidx = {kk: i for i, kk in enumerate(ukeys)}
    accepted = sum(len(r) for r in rows)
    agg = []                                   # per state: {(target, key index): count}
    for row in rows:
        d = {}
        for t, ci in row:
            kk = (t, kidx[keys[ci]])
            d[kk] = d.get(kk, 0) + 1
        agg.append(d)
    n = len(states)
    t_build = time.time() - t0

    src, dst, kid, cnt = [], [], [], []
    for s, d in enumerate(agg):
        for (t, ki), c in d.items():
            src.append(s); dst.append(t); kid.append(ki); cnt.append(c)
    src, dst, kid, cnt = map(np.array, (src, dst, kid, cnt))

    def power(zf, tol=1e-15):
        wf = np.array([float(M.weight(k, kk, Fraction(zf))) for kk in ukeys])
        coef = wf[kid] * cnt
        v = np.ones(n)
        for it in range(20000):
            nv = np.zeros(n)
            np.add.at(nv, src, coef * v[dst])
            r = nv / v
            v = nv / nv.max()
            if r.max() / r.min() < 1 + tol:
                break
        return v, r.max(), it

    if zs == "search":                         # golden-section on the float bound
        f = lambda zf: 2 * math.log(k * zf) / math.log(power(zf, 1e-12)[1])
        a, b = 0.45 / k, 0.95 / k
        g = (math.sqrt(5) - 1) / 2
        c, dd = b - g * (b - a), a + g * (b - a)
        fc, fd = f(c), f(dd)
        while b - a > 2e-5:
            if fc < fd: b, dd, fd = dd, c, fc; c = b - g * (b - a); fc = f(c)
            else: a, c, fc = c, dd, fd; dd = a + g * (b - a); fd = f(dd)
        zs = f"{(a + b) / 2:.4f}"
    z = Fraction(zs)
    assert 0 < z < Fraction(1, k)
    wexact = [M.weight(k, kk, z) for kk in ukeys]          # exact Fractions
    v, _, it = power(float(z))

    V = [Fraction(float(x)) for x in v]
    assert all(x > 0 for x in V)
    lam = mu = None
    for s, d in enumerate(agg):
        num = sum(c * wexact[ki] * V[t] for (t, ki), c in d.items())
        ratio = num / V[s]
        lam = ratio if lam is None or ratio > lam else lam
        mu = ratio if mu is None or ratio < mu else mu
    assert 0 < mu <= lam < 1

    iv.dps = 60
    mpmath.mp.dps = 60
    ivf = lambda f: iv.mpf(f.numerator) / iv.mpf(f.denominator)
    B = 2 * iv.log(k * ivf(z)) / iv.log(ivf(lam))
    Bmu = 2 * iv.log(k * ivf(z)) / iv.log(ivf(mu))
    def exact(x):                                 # an mpf endpoint as an exact rational
        m, e = mpmath.mpf(x).man_exp
        return Fraction(m) * (Fraction(2) ** e)
    hi = exact(B.b)                               # upper end of the interval for the bound
    up = Fraction(math.ceil(hi * 10**6), 10**6)
    assert up >= hi
    d = Fraction(DANCIK.get(k, "1"))
    import os
    os.makedirs("certs", exist_ok=True)
    with open(f"certs/cert_k{k}_h{h}.txt", "w") as fo:   # the certificate itself
        fo.write(f"# Exact certificate for Lueker's M_h: k={k} h={h}\n")
        fo.write(f"# states are rv_machine.canon representatives; v entries are exact doubles (hex)\n")
        fo.write(f"# lambda and mu are exact; the decimal lines are rounded up (upper) or down (lower)\n")
        fo.write(f"k {k}\nh {h}\nz {zs}\nstates {n}\n")
        fo.write(f"lambda {lam.numerator}/{lam.denominator}\nmu {mu.numerator}/{mu.denominator}\n")
        fo.write(f"lambda_upper {dec(lam, 17, math.ceil)}\nmu_lower {dec(mu, 17, math.floor)}\n"
                 f"bound_upper {dec(hi, 15, math.ceil)}\ncertified {float(up):.6f}\n")
        for q, x in zip(states, V):
            fo.write(f"v {q!r} {float(x).hex()}\n")
    print(f"k={k} h={h} z={zs} states={n} classes={len(cl)} accepted={accepted} "
          f"lambda<={dec(lam, 15, math.ceil)} (exact) mu>={dec(mu, 15, math.floor)} "
          f"bound<={dec(hi, 12, math.ceil)} best_at_z>={dec(exact(Bmu.a), 12, math.floor)} "
          f"CERTIFIED gamma_{k} <= {float(up):.6f}  dancik={DANCIK.get(k)} "
          f"{'BEATS' if up < d else 'does not beat'}  "
          f"(build {t_build:.0f}s, total {time.time() - t0:.0f}s, iters {it + 1})", flush=True)


if __name__ == "__main__":
    main()
