"""Check a certificate file written by certify_exact.py, from the file alone.

The file gives k, h, z, the list of states (letter-renaming representatives) and the
vector v as exact doubles. For every listed state q and every input class p of M_h this
computes the transition delta(q, p) with the review implementation
(../review/rv_machine.py), takes the representative of the target and looks it up in
the list. It checks that

- every listed state is its own representative, the start state is listed, v > 0;
- every target is listed (the set is closed under the transitions);
- the largest and the smallest of (A_z v)_q / v_q, in exact rational arithmetic with the
  exact class weights at the decimal z, are the lambda and mu written in the file;
- 2 log(kz)/log(lambda), in interval arithmetic, is at most the certified value.

It does not use the automaton build or the aggregation of certify_exact.py, and it
runs the states in parallel.

Usage: python check_cert.py <cert file> [workers]
"""
import ast
import math
import os
import sys
import time
from fractions import Fraction
from multiprocessing import Pool

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../review"))
import rv_machine as M

_K = _H = _CL = _KEYS = _IDX = None


def _init(k, h, states):
    global _K, _H, _CL, _KEYS, _IDX
    os.nice(10)
    _K, _H = k, h
    _CL = M.classes(k, h)
    _KEYS = [M.wkey(k, h, p) for p in _CL]
    _IDX = {q: i for i, q in enumerate(states)}


def _row(args):
    """Counts of (target index, weight key) over the accepted classes of one state."""
    s, q = args
    counts, missing, acc = {}, 0, 0
    for p, key in zip(_CL, _KEYS):
        t = M.delta(_H, q, p)
        if t is None:
            continue
        acc += 1
        ti = _IDX.get(M.canon(t))
        if ti is None:
            missing += 1
            continue
        counts[(ti, key)] = counts.get((ti, key), 0) + 1
    return s, counts, acc, missing, len(_CL)


def read_cert(path):
    head, states, V = {}, [], []
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        if line.startswith("v "):
            rep, hx = line[2:].rstrip("\n").rsplit(" ", 1)
            states.append(ast.literal_eval(rep))
            V.append(Fraction(float.fromhex(hx)))
        else:
            name, value = line.split(None, 1)
            head[name] = value.strip()
    return head, states, V


def dec(x, d, rnd):
    n = rnd(x * 10**d)
    s = str(n).rjust(d + 1, "0")
    return f"{s[:-d]}.{s[-d:]}"


def main():
    path = sys.argv[1]
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    t0 = time.time()
    head, states, V = read_cert(path)
    k, h = int(head["k"]), int(head["h"])
    z = Fraction(head["z"])
    lam_file, mu_file = Fraction(head["lambda"]), Fraction(head["mu"])
    certified = Fraction(head["certified"])
    n = len(states)
    assert n == int(head["states"]), "state count"
    assert len(set(states)) == n, "repeated state"
    assert all(M.canon(q) == q for q in states), "a listed state is not a representative"
    assert M.canon(M.start_state(h)) in set(states), "start state missing"
    assert all(x > 0 for x in V), "v is not positive"
    assert 0 < z < Fraction(1, k)

    rows = [None] * n
    accepted = missing = 0
    with Pool(workers, _init, (k, h, states)) as pool:
        for s, counts, acc, miss, ncl in pool.imap_unordered(_row, enumerate(states)):
            rows[s] = counts
            accepted += acc
            missing += miss
    assert missing == 0, f"{missing} transitions leave the listed states"

    w = {}
    lam = mu = None
    for s, counts in enumerate(rows):
        num = Fraction(0)
        for (t, key), c in counts.items():
            if key not in w:
                w[key] = M.weight(k, key, z)
            num += c * w[key] * V[t]
        r = num / V[s]
        lam = r if lam is None or r > lam else lam
        mu = r if mu is None or r < mu else mu
    assert lam == lam_file, "largest ratio differs from the file's lambda"
    assert mu == mu_file, "smallest ratio differs from the file's mu"
    assert 0 < mu <= lam < 1

    import mpmath
    from mpmath import iv
    iv.dps = 60
    ivf = lambda f: iv.mpf(f.numerator) / iv.mpf(f.denominator)
    B = 2 * iv.log(k * ivf(z)) / iv.log(ivf(lam))
    m, e = mpmath.mpf(B.b).man_exp
    hi = Fraction(m) * Fraction(2) ** e
    assert hi <= certified, "certified value below the bound"
    print(f"{os.path.basename(path)}: k={k} h={h} z={head['z']} states={n} classes={ncl} "
          f"accepted={accepted}, closed, v>0; lambda and mu equal the file's exact values; "
          f"lambda<={dec(lam, 15, math.ceil)} mu>={dec(mu, 15, math.floor)} "
          f"bound<={dec(hi, 12, math.ceil)} <= certified {head['certified']}: OK "
          f"({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
