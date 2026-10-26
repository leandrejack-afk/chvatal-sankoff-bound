"""Review of rule 5: the optimum of a = -2 ln x / ln t over the boundary rho(A(x, t)) = 1 from the
reviewer's closed form (cert_check5.scalar: rho(A) < 1 iff f(1) > 0, a + b < 2, x t < 1), for
R = 32 and R = oo, at 40 digits, to compare with the float optima printed in certify_rule5_R32.log
and to show what the cap R = 32 costs.
usage: opt5.py
"""
import sys

sys.dont_write_bytecode = True
import mpmath  # noqa: E402

mpmath.mp.dps = 40


def f1(k, x, t, R):
    q, p = mpmath.mpf(k - 1) / k, mpmath.mpf(k - 2) / k
    h = lambda r: 1 / mpmath.mpf(k) + q * (p / q) ** r
    M, S, I = x * x / k, x * x * t * q, x * t
    a = M + S
    if R is None:
        W = I * q / (1 - I * q)
        V = (I * q / (1 - I * q)) / k + q * (I * p / (1 - I * p))
    else:
        W = sum((I * q) ** r for r in range(1, R + 1)) + q ** R * I ** (R + 1) / (1 - I)
        V = sum((I * q) ** r * h(r) for r in range(1, R + 1)) + h(R) * (I * q) ** (R + 1) / (1 - I * q)
    b, c = M * W, M * V
    return 1 - a - b + a * b - 2 * M * b - 2 * S * c, a + b, I


def ok(k, x, t, R):
    """rho(A(x, t)) < 1 (monotone in x: every entry of A grows with x)."""
    f, apb, I = f1(k, x, t, R)
    return I < 1 and apb < 2 and f > 0


def x_of_t(k, t, R):
    lo, hi = mpmath.mpf(1), (1 / t) * mpmath.mpf("0.999999999")
    # bisection on the monotone predicate rho < 1
    if not ok(k, lo, t, R):
        return lo  # no x > 1 works: a = 0 here
    assert not ok(k, hi, t, R)
    for _ in range(200):
        mid = (lo + hi) / 2
        if ok(k, mid, t, R):
            lo = mid
        else:
            hi = mid
    return lo


def a_of_t(k, t, R):
    return -2 * mpmath.log(x_of_t(k, t, R)) / mpmath.log(t)


def optimum(k, R):
    lo, hi = mpmath.mpf("0.03"), mpmath.mpf("0.5")
    g = (mpmath.sqrt(5) - 1) / 2
    c, d = hi - g * (hi - lo), lo + g * (hi - lo)
    fc, fd = a_of_t(k, c, R), a_of_t(k, d, R)
    for _ in range(110):
        if fc > fd:
            hi, d, fd = d, c, fc
            c = hi - g * (hi - lo)
            fc = a_of_t(k, c, R)
        else:
            lo, c, fc = c, d, fd
            d = lo + g * (hi - lo)
            fd = a_of_t(k, d, R)
    t = (c + d) / 2
    return t, a_of_t(k, t, R)


if __name__ == "__main__":
    for k in (2, 3, 4, 5, 6, 7, 8, 16, 32):
        t32, a32 = optimum(k, 32)
        tinf, ainf = optimum(k, None)
        print(f"k={k}: R=32 optimum a = {mpmath.nstr(a32, 13)} at t = {mpmath.nstr(t32, 10)}; "
              f"R=oo optimum a = {mpmath.nstr(ainf, 13)}; R=oo - R=32 = {mpmath.nstr(ainf - a32, 3)}", flush=True)
