"""Check kary.c's table against kref.py.

1. The inputs: kary.c's count equals the number of minimal match pairs with
   both strings of length <= h, found from the definition.
2. For every state (or a random sample), every transition is recomputed with
   kref.delta + kref.canon and merged by (target, weight type); the multiset
   must equal kary.c's records for that state.
3. The state set is closed: every target is a state in the table, and the
   start state is state 0.

Usage: python check_table.py <k> <h> <dir> [sample size] [seed]
"""
import random
import sys
from collections import Counter
from multiprocessing import Pool

import kref

k, h, d = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
sample = int(sys.argv[4]) if len(sys.argv) > 4 else 0
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


def main():
    kk, hh, N, NI, wt, off, rec, keys = kref.read_table(f"{d}/ktab_k{k}_h{h}.bin")
    assert (kk, hh) == (k, h)
    print(f"k={k} h={h}: states {N}, inputs kary.c {NI}, by definition {len(INPUTS)}")
    assert NI == len(INPUTS), "input count differs"
    assert int(keys[0]) == kref.pack(h, kref.start(h)), "state 0 is not the start state"
    assert all(kref.canon(kref.unpack(h, int(x))) == kref.unpack(h, int(x)) for x in keys[: min(N, 20000)])
    wdesc = [tuple(w[:5]) for w in wt]
    ids = list(range(N))
    if sample and sample < N:
        random.seed(seed)
        ids = sorted(random.sample(ids, sample))
    with Pool(8) as p:
        rows = p.map(ref_row, [int(keys[s]) for s in ids], chunksize=16)
    bad = 0
    ntrans = 0
    for s, row in zip(ids, rows):
        mine = Counter()
        for r in rec[int(off[s]):int(off[s + 1])]:
            mine[(int(keys[r["tgt"]]), wdesc[r["wt"]])] += int(r["cnt"])
        ntrans += sum(row.values())
        if mine != row:
            bad += 1
            if bad <= 3:
                print("MISMATCH at state", s, kref.unpack(h, int(keys[s])))
                print("  only kary.c:", list((mine - row).items())[:4])
                print("  only kref:  ", list((row - mine).items())[:4])
    tg = rec["tgt"]
    assert int(tg.max()) < N
    print(f"checked {len(ids)} states, {ntrans} accepted transitions: {bad} mismatching states")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
