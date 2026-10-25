"""Build the LaTeX rows of the main results table from the certificate logs."""
import re

DANCIK = {3: "0.765803", 4: "0.708236", 5: "0.664428", 6: "0.629316", 7: "0.600184",
          8: "0.575407", 9: "0.553937", 10: "0.534855", 11: "0.517842", 12: "0.502591",
          13: "0.488800", 14: "0.476198", 15: "0.464619"}
pat = re.compile(r"k=(\d+) h=(\d+) z=([\d.]+) states=(\d+) classes=(\d+) accepted=(\d+) .*CERTIFIED gamma_\d+ <= ([\d.]+)")


def read(path):
    out = {}
    for line in open(path):
        m = pat.search(line)
        if m:
            k, h, z, st, cl, acc, b = m.groups()
            out[(int(k), int(h))] = dict(z=z, states=int(st), classes=int(cl), acc=int(acc), bound=b)
    return out


def thin(n):
    s = f"{n:,}"
    return s.replace(",", r"\,")


rows = {}
for f in ["certify_h2.log", "certify_h3.log"] + [f"certify_h3_k{k}.log" for k in (10, 11, 12, 13, 14, 15)]:
    try:
        rows.update(read(f))
    except FileNotFoundError:
        pass
best = {3: dict(h=6, z="0.2048", states=1711682, classes=525, bound="0.739975"),
        4: dict(h=5, z="0.1496", states=290442, classes=3460, bound="0.683622")}
lines = []
for k in range(3, 16):
    d = DANCIK[k]
    h2 = rows.get((k, 2))
    h2s = h2["bound"] if h2 else "--"
    if k in best:
        b = best[k]
    else:
        r = rows.get((k, 3))
        if r is None:
            lines.append(f"% k={k}: h=3 result missing")
            continue
        b = dict(h=3, z=r["z"], states=r["states"], classes=r["classes"], bound=r["bound"])
    gain = float(d) - float(b["bound"])
    lines.append(f"{k} & {d} & {h2s} & {b['bound']} & {b['h']} & {thin(b['states'])} & {thin(b['classes'])} & {b['z']} \\\\  % gain {gain:.6f}")
print("\n".join(lines))
for k in range(5, 16):
    r = rows.get((k, 3))
    if r:
        print(f"% k={k} h=3 accepted={r['acc']}")
