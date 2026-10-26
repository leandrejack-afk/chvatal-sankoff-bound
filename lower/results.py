"""Build the results tables from the CERT lines in runs/*.log (and server logs copied there).

Previous best lower bounds on gamma_k: Rosenfeld, arXiv 2407.18113 (2024), Table 1b, except
k = 4 (Heineman et al. / Soiffer et al., arXiv 2407.10925, ISIT 2025, Table 3).
Previous best upper bounds on alpha_k: Rosenfeld 2024, Table 1a.
"""
import glob
import re
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR

PREV_LOW = {3: ("0.68422", "R24"), 4: ("0.614333", "H24"), 5: ("0.56206", "R24"), 6: ("0.51850", "R24"),
            7: ("0.48712", "R24"), 8: ("0.46074", "R24"), 9: ("0.43806", "R24"), 10: ("0.41826", "R24"),
            11: ("0.40072", "R24"), 12: ("0.38504", "R24"), 13: ("0.37088", "R24"), 14: ("0.35798", "R24"),
            15: ("0.34616", "R24")}
PREV_LOW_H = {3: 10, 4: 8, 5: 7, 6: 6, 7: 6, 8: 6, 9: 6, 10: 6, 11: 6, 12: 6, 13: 6, 14: 6, 15: 6}
PREV_EDIT = {3: "0.47276", 4: "0.56578", 5: "0.6325", 6: "0.68424", 7: "0.72016", 8: "0.74896", 9: "0.77264",
             10: "0.7925", 11: "0.8095", 12: "0.82432", 13: "0.83744", 14: "0.84646", 15: "0.85608"}
PREV_EDIT_H = {3: 9, 4: 8, 5: 7, 6: 6, 7: 6, 8: 6, 9: 6, 10: 6, 11: 6, 12: 6, 13: 6, 14: 6, 15: 6}

pat = re.compile(r"CERT k=(\d+) h=(\d+) classes=(\d+) shift=27 r=(-?\d+) (gamma|alpha)_\d+ (>=|<=) ([0-9.]+)")
best = {}
for f in sorted(glob.glob("runs/*.log") + glob.glob("runs/server/*.log")):
    for line in open(f):
        m = pat.search(line)
        if not m:
            continue
        k, h, cl, r, kind, rel, val = int(m[1]), int(m[2]), int(m[3]), int(m[4]), m[5], m[6], m[7]
        key = (kind, k)
        v = Decimal(val)
        if key not in best or (v > best[key][3] if kind == "gamma" else v < best[key][3]):
            best[key] = (h, cl, r, v, f)


def low6(v):  # lower bounds: round down to 6 decimals
    return v.quantize(Decimal("0.000001"), rounding=ROUND_FLOOR)


def up6(v):  # upper bounds: round up to 6 decimals
    return v.quantize(Decimal("0.000001"), rounding=ROUND_CEILING)

print("Lower bounds on gamma_k (LCS)\n")
print("| k | previous best | h | this work | h | classes | gain |")
print("|---|---|---|---|---|---|---|")
for k in range(3, 16):
    if ("gamma", k) in best:
        h, cl, r, val, f = best[("gamma", k)]
        p, src = PREV_LOW[k]
        shown = low6(val)
        print(f"| {k} | {p} ({src}) | {PREV_LOW_H[k]} | **{shown}** | {h} | {cl:,} | {shown - Decimal(p):+.6f} |")
print("\nUpper bounds on alpha_k (edit distance)\n")
print("| k | previous best (R24) | h | this work | h | classes | gain |")
print("|---|---|---|---|---|---|---|")
for k in range(2, 16):
    if ("alpha", k) in best:
        h, cl, r, val, f = best[("alpha", k)]
        p = PREV_EDIT.get(k, "?")
        shown = up6(val)
        gain = f"{Decimal(p) - shown:+.6f}" if p != "?" else ""
        print(f"| {k} | {p} | {PREV_EDIT_H.get(k, '')} | **{shown}** | {h} | {cl:,} | {gain} |")
