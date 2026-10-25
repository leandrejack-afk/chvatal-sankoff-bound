"""check_table.py for kbig's table files: sampled states against kref.py.

1. The inputs: kbig's count equals the number of minimal match pairs with
   both strings of length <= h, found from the definition.
2. For a random sample of states (plus the first and last 200), every
   transition is recomputed with kref.delta + kref.canon and merged by
   (target, weight type); the multiset must equal kbig's records for that state.
3. State 0 is the start state, sampled keys are in relabelled form, and the
   sampled records point at states in the table.  (kbverify checks
   tgt < N, wt < NW and cnt > 0 for every record of the file.)

The file is read through numpy.memmap, so only the sampled rows are loaded.

Usage: python check_table_big.py <k> <h> <dir> <sample size> [seed]
"""
import random
import struct
import sys
from collections import Counter
from multiprocessing import Pool

import numpy as np

import kref

k, h, d = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
sample = int(sys.argv[4])
seed = int(sys.argv[5]) if len(sys.argv) > 5 else 1

INPUTS = kref.inputs_by_definition(k, h)
DESC = [kref.weight_descriptor(k, h, a, b) for a, b in INPUTS]


def ref_row(key):
    q = kref.unpack(h, key)
    row = Counter()
    for (a, b), desc in zip(INPUTS, DESC):
        r = kref.delta(h, q, a, b)
        if r is not None:
            row[(kref.pack(h, kref.canon(r)), desc)] += 1
    return row


def read_kbig(path):
    with open(path, "rb") as f:
        assert f.read(4) == b"KBIG"
        kk, hh = struct.unpack("<ii", f.read(8))
        N, R = struct.unpack("<QQ", f.read(16))
        NI, NW = struct.unpack("<ii", f.read(8))
        wt = [struct.unpack("<5i", f.read(20)) for _ in range(NW)]
    base = 4 + 8 + 16 + 8 + 20 * NW
    rec = np.memmap(path, dtype=[("tgt", "<u4"), ("wt", "<u2"), ("cnt", "<u2")], mode="r", offset=base, shape=(R,))
    off = np.memmap(path, dtype="<u8", mode="r", offset=base + 8 * R, shape=(N + 1,))
    keys = np.memmap(path, dtype="<u8", mode="r", offset=base + 8 * R + 8 * (N + 1), shape=(N,))
    return kk, hh, N, R, NI, wt, off, rec, keys


def main():
    kk, hh, N, R, NI, wt, off, rec, keys = read_kbig(f"{d}/kbig_k{k}_h{h}.bin")
    assert (kk, hh) == (k, h)
    print(f"k={k} h={h}: states {N}, records {R}, inputs kbig {NI}, by definition {len(INPUTS)}")
    assert NI == len(INPUTS), "input count differs"
    assert int(keys[0]) == kref.pack(h, kref.start(h)), "state 0 is not the start state"
    assert int(off[0]) == 0 and int(off[N]) == R, "bad offsets"
    wdesc = [tuple(w[:5]) for w in wt]
    random.seed(seed)
    ids = set(random.sample(range(N), min(sample, N)))
    ids |= set(range(min(200, N))) | set(range(max(0, N - 200), N))
    ids = sorted(ids)
    for s in ids:
        x = int(keys[s])
        assert kref.canon(kref.unpack(h, x)) == kref.unpack(h, x), f"state {s} not in relabelled form"
    with Pool(8) as p:
        rows = p.map(ref_row, [int(keys[s]) for s in ids], chunksize=4)
    bad = 0
    ntrans = 0
    for s, row in zip(ids, rows):
        mine = Counter()
        part = np.array(rec[int(off[s]):int(off[s + 1])])
        assert len(part) == 0 or int(part["tgt"].max()) < N
        for r in part:
            mine[(int(keys[int(r["tgt"])]), wdesc[int(r["wt"])])] += int(r["cnt"])
        ntrans += sum(row.values())
        if mine != row:
            bad += 1
            if bad <= 3:
                print("MISMATCH at state", s, kref.unpack(h, int(keys[s])))
                print("  only kbig:", list((mine - row).items())[:4])
                print("  only kref:", list((row - mine).items())[:4])
    print(f"checked {len(ids)} states ({sample} random, seed {seed}, plus the first and last 200), "
          f"{ntrans} accepted transitions: {bad} mismatching states")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
