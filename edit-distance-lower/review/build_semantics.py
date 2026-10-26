import os
"""Review code: does build() in certify_runs.py / certify_runs4.py encode the intended automaton?

(1) Structural: every entry of their exact matrix is compared with an independently written
    transition table keyed by state NAMES (diag, MD, MI, I1..IR, D1..DR), for several k and R.
(2) Semantic: G = e_diag^T (I - A)^{-1} 1 (sum over all paths from diag) is computed exactly from
    their matrix by Gaussian elimination in Fractions and compared with the closed form
      rules 1-3: G = (1 + 2T) / (1 - a - 2b),
      rules 1-4: solve g_d = 1 + 2T + a g_d + b g_I + b g_D, g_I = 1 + T + a g_d + b g_I,
                       g_D = 1 + T + a g_d + b g_D,
    where a = M + S, b = M W_R, T = I/(1 - I) (a run that never returns).
"""
import importlib.util
import sys
from fractions import Fraction as Fr

HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..") + "/"


def load(name):
    spec = importlib.util.spec_from_file_location(name, HERE + name + ".py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def names(R, rules):
    base = ["diag"] if rules == 3 else ["diag", "MD", "MI"]
    return base + [f"I{r}" for r in range(1, R + 1)] + [f"D{r}" for r in range(1, R + 1)]


def expected(k, x, t, R, rules):
    q = Fr(k - 1, k)
    M, S, I = x * x / k, x * x * t * q, x * t
    T = {}
    diag_like = ["diag"] if rules == 3 else ["diag", "MD", "MI"]
    for s in diag_like:
        T[(s, "diag")] = M + S
        if not (rules == 4 and s == "MD"):
            T[(s, "I1")] = I
        if not (rules == 4 and s == "MI"):
            T[(s, "D1")] = I
    for r in range(1, R + 1):
        for c in "ID":
            back = "diag" if rules == 3 else ("MI" if c == "I" else "MD")
            T[(f"{c}{r}", back)] = M * q ** r
            nxt = f"{c}{min(r + 1, R)}"
            T[(f"{c}{r}", nxt)] = T.get((f"{c}{r}", nxt), 0) + I
    return T


def solve_G(A):
    n = len(A)
    B = [[(1 if i == j else 0) - A[i][j] for j in range(n)] + [Fr(1)] for i in range(n)]
    for col in range(n):
        piv = next(r for r in range(col, n) if B[r][col] != 0)
        B[col], B[piv] = B[piv], B[col]
        pv = B[col][col]
        B[col] = [v / pv for v in B[col]]
        for r in range(n):
            if r != col and B[r][col] != 0:
                f = B[r][col]
                B[r] = [vr - f * vc for vr, vc in zip(B[r], B[col])]
    return B[0][n]


def closed_G(k, x, t, R, rules):
    q = Fr(k - 1, k)
    M, S, I = x * x / k, x * x * t * q, x * t
    W = sum((I * q) ** r for r in range(1, R + 1)) + q ** R * I ** (R + 1) / (1 - I)
    a, b, T = M + S, M * W, I / (1 - I)
    if rules == 3:
        return (1 + 2 * T) / (1 - a - 2 * b)
    # by symmetry g_I = g_D = g; g = (1 + T + a g_d)/(1 - b); g_d = 1 + 2T + a g_d + 2 b g
    # => g_d (1 - a) = 1 + 2T + 2b (1 + T + a g_d)/(1 - b)
    num = (1 + 2 * T) * (1 - b) + 2 * b * (1 + T)
    den = (1 - a) * (1 - b) - 2 * a * b
    return num / den


def main():
    mods = {3: load("certify_runs"), 4: load("certify_runs4")}
    bad = 0
    for rules in (3, 4):
        for k in (2, 3, 5, 32):
            for R in (1, 2, 3, 6):
                for x, t in ((Fr(1), Fr(1, 10)), (Fr(11, 10), Fr(1, 7)), (Fr(6, 5), Fr(1, 9))):
                    A = mods[rules].build(k, x, t, R, exact=True)
                    nm = names(R, rules)
                    assert len(A) == len(nm)
                    exp = expected(k, x, t, R, rules)
                    got = {(nm[i], nm[j]): A[i][j] for i in range(len(A)) for j in range(len(A)) if A[i][j] != 0}
                    struct_ok = got == exp
                    G1, G2 = solve_G(A), closed_G(k, x, t, R, rules)
                    ok = struct_ok and G1 == G2 and G1 > 0
                    bad += not ok
                    if not ok or (k == 2 and x == 1):
                        print(f"rules 1-{rules} k={k} R={R} x={x} t={t}: structure {'same' if struct_ok else 'DIFFERENT'}; "
                              f"G matrix == G closed form: {G1 == G2} (G = {float(G1):.12f})")
    # one R = 32 point per rule set (the certified size)
    for rules in (3, 4):
        A = mods[rules].build(2, Fr(11, 10), Fr(1, 7), 32, exact=True)
        G1, G2 = solve_G(A), closed_G(2, Fr(11, 10), Fr(1, 7), 32, rules)
        print(f"rules 1-{rules} k=2 R=32 x=11/10 t=1/7: G matrix == G closed form: {G1 == G2} ({float(G1):.15f})")
        bad += G1 != G2
    print("ALL OK" if bad == 0 else f"{bad} FAILURES")


if __name__ == "__main__":
    main()
