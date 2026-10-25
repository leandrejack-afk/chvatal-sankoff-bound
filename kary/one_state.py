"""The one-state automaton that accepts every list of minimal match pairs.
Its matrix is the 1x1 sum of all class weights; minimizing 2 log(kz)/log(lambda) over z
gives the bound obtained by counting every list (Dancik's thesis, Fig. 4.2, 'Deken')."""
import mpmath as mp
from math import comb
mp.mp.dps = 30
def surj(m, z):   # strings that use each of m given letters, and no others
    return sum((-1) ** (m - j) * comb(m, j) / (1 - j * z) for j in range(m + 1))
def lam(k, z):    # k choices of c; disjoint letter sets S, T for a', b'
    n = k - 1
    return k * z ** 2 * sum(mp.factorial(n) / (mp.factorial(s) * mp.factorial(t) * mp.factorial(n - s - t))
                            * surj(s, z) * surj(t, z) for s in range(n + 1) for t in range(n + 1 - s))
for k in (2, 3, 4):
    f = lambda z: 2 * mp.log(k * z) / mp.log(lam(k, z))
    z = mp.findroot(lambda z: mp.diff(f, z), 0.8 / (2 * k))
    z = mp.re(z)
    print(f"k={k}: z*={mp.nstr(z, 8)} bound={mp.nstr(mp.re(f(z)), 10)}")
