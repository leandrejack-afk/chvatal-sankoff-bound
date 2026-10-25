"""Locate the best z from a kbig eigen run: parabola through the three best points.

Reads the final "z=... lambda in [lo, hi] bound=..." lines of kbig eigen's output,
prints each bound (from the upper end hi, as the certificate will use) and the
vertex of the parabola through the best point and its two neighbours.

Usage: python zfit.py <eigen output file>
"""
import re
import sys

pts = []
for line in open(sys.argv[1]):
    m = re.match(r"z=([0-9.]+) lambda in \[([0-9.]+), ([0-9.]+)\] bound=([0-9.]+)", line)
    if m:
        z, lo, hi, b = map(float, m.groups())
        pts.append((z, b, hi / lo - 1))
pts.sort()
for z, b, gap in pts:
    print(f"z={z:.4f} bound={b:.9f} (max/min - 1 = {gap:.1e})")
i = min(range(len(pts)), key=lambda t: pts[t][1])
if i == 0 or i == len(pts) - 1:
    print(f"minimum at the end of the range (z={pts[i][0]:.4f}): extend the range")
    sys.exit(1)
(z0, f0, _), (z1, f1, _), (z2, f2, _) = pts[i - 1], pts[i], pts[i + 1]
# vertex of the parabola through three points
num = (z1 - z0) ** 2 * (f1 - f2) - (z1 - z2) ** 2 * (f1 - f0)
den = (z1 - z0) * (f1 - f2) - (z1 - z2) * (f1 - f0)
zv = z1 - 0.5 * num / den
a = ((f2 - f1) / (z2 - z1) - (f1 - f0) / (z1 - z0)) / (z2 - z0)
fv = f1 - a * (z1 - zv) ** 2
print(f"vertex z={zv:.6f}, predicted bound {fv:.9f}, curvature {a:.3g}")
print(f"4-decimal z nearest the vertex: {round(zv, 4):.4f}")
