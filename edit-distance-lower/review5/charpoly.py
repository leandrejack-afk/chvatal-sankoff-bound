"""Symbolic check of the censored 4x4 matrix's characteristic polynomial (sympy)."""
import sympy as sp
M, S, b, c, l = sp.symbols("M S b c l")
K = sp.Matrix([[M, S, b, b], [M, S, c, c], [M, S, b, 0], [M, S, 0, b]])
cp = sp.factor((l * sp.eye(4) - K).det())
print("det(l I - K) =", cp)
a = M + S
target = sp.expand(l * (l - b) * (l**2 - (a + b) * l + (a * b - 2 * M * b - 2 * S * c)))
print("equals l (l - b) (l^2 - (a+b) l + ab - 2Mb - 2Sc):", sp.simplify(sp.expand((l * sp.eye(4) - K).det()) - target) == 0)
print("f(1) =", sp.expand(1 - (a + b) + a * b - 2 * M * b - 2 * S * c))
# rules 1-4 reduction (c = b, G0 and G1 merged) gives l^2 - (a+b) l - ab
print("c=b:", sp.factor(sp.expand(l**2 - (a + b) * l + (a * b - 2 * M * b - 2 * S * b))))
