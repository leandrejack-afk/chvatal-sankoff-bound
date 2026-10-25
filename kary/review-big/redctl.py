"""Red controls for kbverify, check_table_big.py and rvcmp_big (reviewer, 2026-10-25).

Makes corrupted copies of review-big/repro/kbig_k<k>_h<h>.bin (my own kbig rebuild, byte-identical
to kary.c's table after "klue") under review-big/redctl/<case>/, each with the data/ vector, weights
and wtypes copied beside it, then runs the three checkers on every case and prints what each says.

Usage (from kary/): python review-big/redctl.py 3 4
"""
import os
import shutil
import struct
import subprocess
import sys

import numpy as np

KA = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RB = os.path.join(KA, "review-big")
PY = sys.executable
k, h = int(sys.argv[1]), int(sys.argv[2])
src = f"{RB}/repro/kbig_k{k}_h{h}.bin"
names = [f"kvec_k{k}_h{h}.bin", f"weights_k{k}_h{h}.txt", f"wtypes_k{k}_h{h}.txt"]

with open(src, "rb") as f:
    assert f.read(4) == b"KBIG"
    f.read(8)
    N, R = struct.unpack("<QQ", f.read(16))
    NI, NW = struct.unpack("<ii", f.read(8))
    WT = [struct.unpack("<5i", f.read(20)) for _ in range(NW)]
BASE = 36 + 20 * NW
OFFP = BASE + 8 * R
KEYP = OFFP + 8 * (N + 1)
recdt = np.dtype([("tgt", "<u4"), ("wt", "<u2"), ("cnt", "<u2")])


def case_dir(name):
    d = f"{RB}/redctl/k{k}h{h}_{name}"
    os.makedirs(d, exist_ok=True)
    shutil.copyfile(src, f"{d}/kbig_k{k}_h{h}.bin")
    for n in names:
        shutil.copyfile(f"{KA}/data/{n}", f"{d}/{n}")
    return d


def maps(d):
    p = f"{d}/kbig_k{k}_h{h}.bin"
    rec = np.memmap(p, dtype=recdt, mode="r+", offset=BASE, shape=(R,))
    off = np.memmap(p, dtype="<u8", mode="r+", offset=OFFP, shape=(N + 1,))
    keys = np.memmap(p, dtype="<u8", mode="r+", offset=KEYP, shape=(N,))
    return rec, off, keys


# a state in the middle with at least 3 records, and its first records
off0 = np.memmap(src, dtype="<u8", mode="r", offset=OFFP, shape=(N + 1,))
rec0 = np.memmap(src, dtype=recdt, mode="r", offset=BASE, shape=(R,))
S = next(s for s in range(N // 2, N) if off0[s + 1] - off0[s] >= 3)
r0 = int(off0[S])
first = rec0[r0]
# a weight type of the same size class with a smaller weight is hard to pick blind: use any other type
other_wt = (int(first["wt"]) + 1) % NW
other_tgt = (int(first["tgt"]) + 1) % N
r_multi = next(r for r in range(R // 2, R) if rec0[r]["cnt"] >= 2)

cases = {}


def c_trunc():
    d = case_dir("trunc")
    p = f"{d}/kbig_k{k}_h{h}.bin"
    os.truncate(p, os.path.getsize(p) - 8)
    return d, "file 8 bytes short"


def c_tgt_oob():
    d = case_dir("tgt_oob")
    rec, _, _ = maps(d)
    rec[r0]["tgt"] = N
    rec.flush()
    return d, f"record {r0}: tgt = N"


def c_wt_other():
    d = case_dir("wt_other")
    rec, _, _ = maps(d)
    rec[r0]["wt"] = other_wt
    rec.flush()
    return d, f"record {r0} (state {S}): weight type {int(first['wt'])} -> {other_wt}"


def c_tgt_other():
    d = case_dir("tgt_other")
    rec, _, _ = maps(d)
    rec[r0]["tgt"] = other_tgt
    rec.flush()
    return d, f"record {r0} (state {S}): target {int(first['tgt'])} -> {other_tgt}"


def c_cnt_dec():
    d = case_dir("cnt_dec")
    rec, _, _ = maps(d)
    rec[r_multi]["cnt"] = rec[r_multi]["cnt"] - 1
    rec.flush()
    return d, f"record {r_multi}: multiplicity minus 1"


def c_off_shift():
    d = case_dir("off_shift")
    _, off, _ = maps(d)
    off[S + 1] = off[S + 1] - 1          # last record of state S now belongs to state S+1
    off.flush()
    return d, f"off[{S + 1}] - 1: one record moves from state {S} to state {S + 1}"


def c_off_nonmono():
    d = case_dir("off_nonmono")
    _, off, _ = maps(d)
    off[S + 1] = off[S] - 1 if off[S] > 0 else 0
    off[S] = off[S] + 1
    off.flush()
    return d, f"off[{S}], off[{S + 1}] made decreasing"


def c_key0():
    d = case_dir("key0")
    _, _, keys = maps(d)
    keys[0] = keys[1]
    keys.flush()
    return d, "keys[0] := keys[1]"


def c_hdr_wtype():
    d = case_dir("hdr_wtype")
    p = f"{d}/kbig_k{k}_h{h}.bin"
    with open(p, "r+b") as f:
        f.seek(36 + 20 * (NW - 1))
        w = list(WT[NW - 1])
        w[0] += 1
        f.write(struct.pack("<5i", *w))
    return d, "header weight type NW-1: size + 1"


def c_drop_record():
    """remove record r0 entirely: R - 1, later offsets - 1 (a consistent-looking file)"""
    d = f"{RB}/redctl/k{k}h{h}_drop_record"
    os.makedirs(d, exist_ok=True)
    for n in names:
        shutil.copyfile(f"{KA}/data/{n}", f"{d}/{n}")
    with open(src, "rb") as f:
        head = bytearray(f.read(BASE))
    rec = np.array(rec0)
    rec = np.concatenate([rec[:r0], rec[r0 + 1:]])
    off = np.array(off0).copy()
    off[S + 1:] -= 1
    keys = np.fromfile(src, dtype="<u8", count=N, offset=KEYP)
    struct.pack_into("<QQ", head, 12, N, R - 1)
    with open(f"{d}/kbig_k{k}_h{h}.bin", "wb") as f:
        f.write(head)
        f.write(rec.tobytes())
        f.write(off.astype("<u8").tobytes())
        f.write(keys.tobytes())
    return d, f"record {r0} (state {S}) deleted, header and offsets made consistent"


def c_vec_zero():
    d = case_dir("vec_zero")
    v = np.memmap(f"{d}/kvec_k{k}_h{h}.bin", dtype="<f8", mode="r+", offset=24, shape=(N,))
    v[S] = 0.0
    v.flush()
    return d, f"v[{S}] = 0"


def c_vec_nan():
    d = case_dir("vec_nan")
    v = np.memmap(f"{d}/kvec_k{k}_h{h}.bin", dtype="<f8", mode="r+", offset=24, shape=(N,))
    v[S] = np.nan
    v.flush()
    return d, f"v[{S}] = NaN"


def c_vec_long():
    d = case_dir("vec_long")
    with open(f"{d}/kvec_k{k}_h{h}.bin", "ab") as f:
        f.write(b"\0" * 8)
    return d, "vector file 8 bytes long"


def c_wtypes_extra():
    d = case_dir("wtypes_extra")
    with open(f"{d}/wtypes_k{k}_h{h}.txt", "a") as f:
        f.write(f"{NW} 2 0 0 0 0\n")
    return d, "wtypes file has one extra line"


def c_control():
    d = case_dir("control")
    return d, "unmodified copy"


def run(cmd, env=None):
    p = subprocess.run(cmd, capture_output=True, text=True, cwd=KA, env=env)
    out = (p.stdout + p.stderr).strip().splitlines()
    keep = [l for l in out if any(s in l for s in ("lambda_cert", "rho_lower", "MISMATCH", "mismatch", "checked",
                                                  "differing", "orbit", "bad", "short", "size", "differs",
                                                  "not", "Error", "error", "positive", "long", "more",
                                                  "monotone", "start"))]
    return p.returncode, keep


for fn in [c_control, c_trunc, c_tgt_oob, c_wt_other, c_tgt_other, c_cnt_dec, c_off_shift, c_off_nonmono,
           c_key0, c_hdr_wtype, c_drop_record, c_vec_zero, c_vec_nan, c_vec_long, c_wtypes_extra]:
    d, what = fn()
    print(f"=== {os.path.basename(d)}: {what}")
    rc, out = run([f"{KA}/kbverify", str(k), str(h), d])
    print(f"  kbverify exit {rc}: " + " | ".join(out[:3]))
    if fn in (c_vec_zero, c_vec_nan, c_vec_long, c_wtypes_extra):
        continue
    rc, out = run([PY, "check_table_big.py", str(k), str(h), d, str(N), "5"])
    print(f"  check_table_big (all {N} states) exit {rc}: " + " | ".join(l[:160] for l in out[:4]))
    rc, out = run([f"{KA}/rvcmp_big", str(k), str(h), f"{d}/kbig_k{k}_h{h}.bin"])
    print(f"  rvcmp_big (stride 1) exit {rc}: " + " | ".join(l[:200] for l in out[:3]))
