"""Compare a kbig file with kary.c's ktab file section by section, without kbig's own "klue"."""
import struct
import sys
import numpy as np

def sections(path, magic):
    with open(path, "rb") as f:
        assert f.read(4) == magic
        k, h = struct.unpack("<ii", f.read(8)); N, R = struct.unpack("<QQ", f.read(16)); NI, NW = struct.unpack("<ii", f.read(8))
        wt = f.read(20 * NW)
    base = 36 + 20 * NW
    if magic == b"KBIG":
        rec_at, off_at, key_at = base, base + 8 * R, base + 8 * R + 8 * (N + 1)
    else:
        off_at, rec_at, key_at = base, base + 8 * (N + 1), base + 8 * (N + 1) + 8 * R
    m = lambda at, n: np.memmap(path, dtype="<u8", mode="r", offset=at, shape=(n,))
    return (k, h, N, R, NI, NW, wt), m(rec_at, R), m(off_at, N + 1), m(key_at, N)

for k, h in [(3, 4), (3, 5), (3, 6), (4, 4), (4, 5), (2, 8)]:
    a = sections(f"review-big/repro/kbig_k{k}_h{h}.bin", b"KBIG")
    b = sections(f"data/ktab_k{k}_h{h}.bin", b"KLUE")
    same = a[0] == b[0]
    for i in (1, 2, 3):
        x, y = a[i], b[i]
        for s in range(0, len(x), 1 << 26):
            if not np.array_equal(x[s:s + (1 << 26)], y[s:s + (1 << 26)]):
                same = False
                break
    print(f"k={k} h={h} N={a[0][2]} R={a[0][3]}: header, records, offsets, keys identical: {same}")
