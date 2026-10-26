"""Brute-force exact check of a klow witness on the FULL vector (every pair of strings).

Written from Lueker J.ACM 2009 Sec. 2 and the witness format only (no code from klow.c):
  * class of (s,t): relabel s.t and t.s by first appearance (RGS), take the one that comes
    first in lexicographic order; classes are listed in lexicographic order of that RGS.
  * T(w)(as',at') = 1 + (1/k^2) sum_{c,c'} w(s'c, t'c')
    T(w)(as',bt') = max( (1/k) sum_c w(s'c, bt'), (1/k) sum_c w(as', t'c) ),  a != b
  * r_exact = min over ALL k^(2h) pairs of T(w) - w, exact rational (units 2^-27).
Usage: python3 bf_low.py k h witness.bin [claimed_r]
"""
import sys, struct
from fractions import Fraction
from math import floor

def rgs(z):
    m = {}
    return tuple(m.setdefault(x, len(m)) for x in z)

def all_rgs_lex(L, K):
    out = []
    def rec(p, mx):
        if len(p) == L:
            out.append(tuple(p)); return
        for v in range(min(mx + 2, K)):
            p.append(v); rec(p, max(mx, v)); p.pop()
    rec([0], 0)
    return out

def main():
    K, H = int(sys.argv[1]), int(sys.argv[2])
    data = open(sys.argv[3], "rb").read()
    claimed = int(sys.argv[4]) if len(sys.argv) > 4 else None
    lst = all_rgs_lex(2 * H, K)
    assert lst == sorted(lst) and len(set(lst)) == len(lst)
    rank = {z: i for i, z in enumerate(lst)}
    canon = [z for z in lst if rank[z] <= rank[rgs(z[H:] + z[:H])]]
    dense = {z: i for i, z in enumerate(canon)}
    NW = len(canon)
    P = K ** (H - 1); PH = P * K; N = PH * PH
    def digits(x):
        d = []
        for _ in range(H): d.append(x % K); x //= K
        return tuple(reversed(d))          # first letter = most significant digit
    strs = [digits(x) for x in range(PH)]
    if len(data) == 4 * N and N != NW:     # klow_plain dump: full vector, index s*k^h + t
        W = list(struct.unpack(f"<{N}I", data))
        print("plain-format witness (all pairs)")
    elif len(data) == 4 * NW:
        wv = struct.unpack(f"<{NW}I", data)
        W = [0] * N
        for s in range(PH):
            a = strs[s]
            for t in range(PH):
                b = strs[t]
                z1, z2 = rgs(a + b), rgs(b + a)
                W[s * PH + t] = wv[dense[z1] if rank[z1] <= rank[z2] else dense[z2]]
    else:
        print(f"FAIL: witness has {len(data)//4} values, class count is {NW}"); sys.exit(2)
    # exact sums
    A = [0] * (P * PH)   # A[s', y] = sum_c W[s'c, y]
    for sp in range(P):
        for y in range(PH):
            A[sp * PH + y] = sum(W[(sp * K + c) * PH + y] for c in range(K))
    B = [0] * (PH * P)   # B[x, t'] = sum_c W[x, t'c]
    for x in range(PH):
        for tp in range(P):
            base = x * PH + tp * K
            B[x * P + tp] = sum(W[base:base + K])
    U = 1 << 27; K2 = K * K
    best = None; arg = None; worst = None
    for s in range(PH):
        s0, sp = divmod(s, P)
        for t in range(PH):
            t0, tp = divmod(t, P)
            w = W[s * PH + t]
            if s0 == t0:
                S = sum(A[sp * PH + tp * K + c2] for c2 in range(K))  # sum_{c,c'} W[s'c, t'c']
                num = K2 * U + S - K2 * w                              # K^2 (T - w)
            else:
                num = K * (max(A[sp * PH + t], B[s * P + tp]) - K * w)
            if best is None or num < best: best, arg = num, (s, t)
            if worst is None or num > worst: worst = num
    r = Fraction(best, K2)
    print(f"k={K} h={H} classes={NW} pairs={N}")
    print(f"exact r = {best}/{K2} units = {float(r):.6f} units = {float(r)/U:.12f}; max T-w = {worst/K2/U:.12f}")
    s, t = arg
    print(f"argmin pair s={strs[s]} t={strs[t]}")
    if r > 0:
        print(f"bound from exact r: {floor(Fraction(2) * r / (U + r) * 10**9) / 1e9:.9f}")
    if claimed is not None:
        ok = r >= claimed
        if claimed > 0:
            print(f"bound from claimed r={claimed}: 0.{floor(Fraction(2*claimed, U + claimed) * 10**9):09d}")
        print(("PASS" if ok else "FAIL") + f": exact r - claimed r = {float(r - claimed):.4f} units")
        sys.exit(0 if ok else 1)

main()
