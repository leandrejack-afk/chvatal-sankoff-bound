"""First-moment lower bounds on the edit-distance constant alpha_k (exploration, not certified).

P(ed(X,Y) <= a n) <= E[# canonical valid alignments of cost <= a n]. Moves: M (match, letters equal,
prob 1/k), S (substitution, letters differ, prob (k-1)/k), I, D (no condition). With weights
M = x^2/k, S = x^2 t (k-1)/k, I = D = x t (x marks a letter of either string, t a unit of cost),
the count is at most G(x,t) / (x^{2n} t^{a n}) for t <= 1 whenever the move automaton's spectral
radius rho(x,t) < 1. So alpha_k >= a whenever rho(x,t) < 1 and x^2 t^a > 1 for some x, t; the best
a is max over (x,t) of -2 log x / log t subject to rho(x,t) <= 1.
Automata: 'free' (any sequence), 'noID' (no I directly next to D, as an optimal alignment never has
an adjacent insertion and deletion), 'noID+order' (also D/I never directly after M... see code).
"""
import numpy as np
from scipy.optimize import brentq

def rho(k, x, t, auto):
    M, S, I, D = x * x / k, x * x * t * (k - 1) / k, x * t, x * t
    if auto == "free":
        return M + S + I + D
    if auto == "noID":
        # states: last move in {diag, I, D}; from I cannot go to D and vice versa
        A = np.array([[M + S, I, D],   # from diag
                      [M + S, I, 0],   # from I
                      [M + S, 0, D]])  # from D
        return max(abs(np.linalg.eigvals(A)))
    if auto == "noID_noIS":
        # also no substitution directly after an indel (I S -> S I, D S -> S D keep the cost)
        A = np.array([[M + S, I, D],   # from diag
                      [M,     I, 0],   # from I: no D, no S
                      [M,     0, D]])  # from D: no I, no S
        return max(abs(np.linalg.eigvals(A)))
    raise ValueError

def best_alpha(k, auto):
    best = 0.0
    for t in np.linspace(0.02, 0.999, 400):
        # largest x with rho <= 1 (rho increasing in x)
        f = lambda x: rho(k, x, t, auto) - 1
        hi = 1.0
        while f(hi) < 0: hi *= 2
        x = brentq(f, 1e-9, hi)
        if x > 1:
            a = -2 * np.log(x) / np.log(t)
            best = max(best, a)
    return best

for k in (2, 3, 4, 5, 6, 7, 8, 16, 32):
    print(k, {auto: round(best_alpha(k, auto), 5) for auto in ("free", "noID", "noID_noIS")})
