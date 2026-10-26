"""Copy a witness and add a signed delta to one uint32 entry: tamper.py src dst index delta."""
import struct
import sys

src, dst, index, delta = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
data = bytearray(open(src, "rb").read())
old = struct.unpack_from("<I", data, 4 * index)[0]
new = old + delta
if not 0 <= new < 2**32:
    sys.exit(f"entry {index} = {old}; {old} + {delta} leaves uint32 range")
struct.pack_into("<I", data, 4 * index, new)
open(dst, "wb").write(bytes(data))
print(f"{dst}: entry {index} changed {old} -> {new} (delta {delta:+d}), {len(data) // 4} values")
