"""Weights at the new (k, h) before any certificate exists (reviewer, 2026-10-25).

For (3, 7) and (4, 6):
1. kbig's input count and weight types (wtharness, make_inputs only) against kref's classes
   by definition and kref.weight_descriptor: same count, same set of types.
2. kweights.check parts 2 and 3: for every class, the number of minimal match pairs of each
   total size <= S that mu maps to it (brute force) equals the coefficients of the closed form.
3. kweights.write at a plausible z (a stub kvec header in this folder), then
   review/rv_theirweights.py (exact Fractions, the earlier referee's own weight formula).

Run from kary/: python review-big/wt/check_weights_new_h.py
"""
import os
import struct
import subprocess
import sys
from collections import defaultdict

KA = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, KA)
import kref  # noqa: E402
import kweights  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable

for k, h, S, z in [(3, 7, 18, 0.2126), (4, 6, 14, 0.1575)]:
    print(f"=== k={k} h={h}")
    types = [tuple(map(int, l.split()))[1:] for l in open(f"{HERE}/wtypes_k{k}_h{h}.txt")]
    classes = kref.inputs_by_definition(k, h)
    descs = {kref.weight_descriptor(k, h, a, b) for a, b in classes}
    print(f"classes by definition: {len(classes)}; kbig weight types: {len(types)}; "
          f"same set: {set(types) == descs}")
    # part 2 of kweights.check, without part 1 (definition vs constructive, h-independent and done before)
    counts = defaultdict(lambda: [0] * (S + 1))
    n = 0
    for a, b in kweights.mmps_constructive(k, S):
        counts[(a[-h:], b[-h:])][len(a) + len(b)] += 1
        n += 1
    assert set(counts) == set(classes), "class sets differ"
    bad = 0
    per_type = defaultdict(int)
    for cls in classes:
        desc = kref.weight_descriptor(k, h, *cls)
        if counts[cls] != kweights.closed_series(k, desc, S):
            bad += 1
        per_type[desc] += 1
    nz = min(sum(1 for c in kweights.closed_series(k, d, S) if c) for d in descs)
    print(f"{n} minimal match pairs of size <= {S}: {bad} classes whose counts differ from the closed form "
          f"(every type has >= {nz} nonzero coefficients checked)")
    # write weights at z with a stub vector header, then the exact check
    d = f"{HERE}/k{k}h{h}"
    os.makedirs(d, exist_ok=True)
    with open(f"{d}/kvec_k{k}_h{h}.bin", "wb") as f:
        f.write(struct.pack("<iiQd", k, h, 0, z))
    with open(f"{HERE}/wtypes_k{k}_h{h}.txt") as src, open(f"{d}/wtypes_k{k}_h{h}.txt", "w") as dst:
        dst.write(src.read())
    print(subprocess.run([PY, f"{KA}/kweights.py", "write", str(k), str(h), d], capture_output=True, text=True,
                         check=True).stdout.strip())
    print(subprocess.run([PY, "rv_theirweights.py", str(k), str(h), d, f"{d}/refweights.txt"], cwd=f"{KA}/review",
                         capture_output=True, text=True, check=True).stdout.strip())
