"""Brute-force exact check of a kedit witness on the FULL vector (every pair of strings).

From Rosenfeld arXiv 2407.18113, Sec. 2 (Facts 5, 6, Lemma 7) and the witness format only:
  T(v, v')(as', at') = (1/k^2) sum_{c,c'} v'(s'c, t'c')
  T(v, v')(as', bt') = 1 + min( (1/k) sum_c v(s'c, bt'), (1/k) sum_c v(as', t'c), (1/k^2) sum v'(s'c, t'c') )
  Lemma 7: T(v, v - r) <= v + r everywhere  =>  alpha_k <= 2r.
Computes the least such r exactly (rational, units 2^-27) over all k^(2h) pairs, and checks
T(v, v - r_claimed) <= v + r_claimed at every pair with exact arithmetic.
Usage: python3 bf_edit.py k h witness.bin [claimed_r]
"""
import sys, struct
from fractions import Fraction
from math import ceil

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

K, H = int(sys.argv[1]), int(sys.argv[2])
data = open(sys.argv[3], "rb").read()
claimed = int(sys.argv[4]) if len(sys.argv) > 4 else None
lst = all_rgs_lex(2 * H, K)
rank = {z: i for i, z in enumerate(lst)}
canon = [z for z in lst if rank[z] <= rank[rgs(z[H:] + z[:H])]]
dense = {z: i for i, z in enumerate(canon)}
NW = len(canon)
if len(data) != 4 * NW:
    print(f"FAIL: witness has {len(data)//4} values, class count is {NW}"); sys.exit(2)
wv = struct.unpack(f"<{NW}I", data)
P = K ** (H - 1); PH = P * K; N = PH * PH
def digits(x):
    d = []
    for _ in range(H): d.append(x % K); x //= K
    return tuple(reversed(d))
strs = [digits(x) for x in range(PH)]
V = [0] * N
for s in range(PH):
    for t in range(PH):
        z1, z2 = rgs(strs[s] + strs[t]), rgs(strs[t] + strs[s])
        V[s * PH + t] = wv[dense[z1] if rank[z1] <= rank[z2] else dense[z2]]
A = [sum(V[(sp * K + c) * PH + y] for c in range(K)) for sp in range(P) for y in range(PH)]  # A[s',y]
B = [sum(V[x * PH + tp * K: x * PH + tp * K + K]) for x in range(PH) for tp in range(P)]      # B[x,t']
U = 1 << 27; K2 = K * K
rstar = None; arg = None; viol = 0
for s in range(PH):
    s0, sp = divmod(s, P)
    for t in range(PH):
        t0, tp = divmod(t, P)
        v = V[s * PH + t]
        G = Fraction(sum(A[sp * PH + tp * K + c2] for c2 in range(K)), K2)
        if s0 == t0:
            need = (G - v) / 2
        else:
            mab = Fraction(min(A[sp * PH + t], B[s * P + tp]), K)
            need = min(1 * U + mab - v, (U + G - v) / 2)
        if rstar is None or need > rstar: rstar, arg = need, (s, t)
        if claimed is not None:
            rc = claimed
            if s0 == t0: Tv = G - rc
            else: Tv = U + min(Fraction(A[sp * PH + t], K), Fraction(B[s * P + tp], K), G - rc)
            if Tv > v + rc: viol += 1
print(f"k={K} h={H} classes={NW} pairs={N}")
print(f"exact least r = {rstar} units = {float(rstar):.4f} units = {float(rstar)/U:.12f}")
s, t = arg
print(f"argmax pair s={strs[s]} t={strs[t]}")
q = ceil(Fraction(2) * rstar / U * 10**9)
print(f"alpha bound from exact least r: {q // 10**9}.{q % 10**9:09d}")
if claimed is not None:
    ok = claimed >= rstar and viol == 0
    print(("PASS" if ok else "FAIL") + f": claimed r = {claimed}, claimed - exact least r = {float(claimed - rstar):.4f} units, pairs violating T(v,v-r)<=v+r: {viol}")
    sys.exit(0 if ok else 1)
