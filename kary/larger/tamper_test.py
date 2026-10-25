"""Red controls for check_cert.py: damaged copies of a certificate must be rejected.

Writes six damaged copies of certs/cert_k3_h3.txt to a temporary directory and runs
check_cert.py on each and on the undamaged file. Prints one line per case.

Usage: python tamper_test.py
"""
import ast
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "certs", "cert_k3_h3.txt")


def run(path):
    r = subprocess.run([sys.executable, os.path.join(HERE, "check_cert.py"), path, "2"],
                       capture_output=True, text=True)
    out = (r.stdout + r.stderr).strip().splitlines()
    return r.returncode, out[-1] if out else ""


def main():
    lines = open(SRC).read().splitlines(keepends=True)
    vi = [i for i, l in enumerate(lines) if l.startswith("v ")]

    def edit(fn):
        a = lines[:]
        fn(a)
        return a

    def set_field(a, name, value):
        for i, l in enumerate(a):
            if l.split(None, 1)[0] == name:
                a[i] = f"{name} {value}\n"

    def scale_v(a, j, f):
        rep, hx = a[vi[j]].rstrip("\n").rsplit(" ", 1)
        a[vi[j]] = f"{rep} {(float.fromhex(hx) * f).hex()}\n"

    def drop_state(a):
        del a[vi[100]]
        set_field(a, "states", str(len(vi) - 1))

    def rename_letters(a):                  # swap letters 0 and 1 in one state: not a representative
        rep, hx = a[vi[50]][2:].rstrip("\n").rsplit(" ", 1)
        x, y, dx, dy = ast.literal_eval(rep)
        sw = lambda t: tuple({0: 1, 1: 0}.get(c, c) for c in t)
        a[vi[50]] = f"v {(sw(x), sw(y), dx, dy)!r} {hx}\n"

    def bump_lambda(a):                     # lambda one unit larger in the numerator
        for i, l in enumerate(a):
            if l.startswith("lambda "):
                n, d = l.split()[1].split("/")
                a[i] = f"lambda {int(n) + 1}/{d}\n"

    cases = [
        ("undamaged file", lines, 0),
        ("z changed from 0.1652 to 0.1653", edit(lambda a: set_field(a, "z", "0.1653")), 1),
        ("certified value lowered to 0.754806", edit(lambda a: set_field(a, "certified", "0.754806")), 1),
        ("one state removed", edit(drop_state), 1),
        ("one state with letters 0 and 1 swapped", edit(rename_letters), 1),
        ("one entry of v multiplied by 0.999", edit(lambda a: scale_v(a, 7, 0.999)), 1),
        ("lambda in the file changed in its last digit", edit(bump_lambda), 1),
    ]
    bad = 0
    with tempfile.TemporaryDirectory() as tmp:
        for i, (name, content, expect_fail) in enumerate(cases):
            p = os.path.join(tmp, f"case{i}.txt")
            open(p, "w").writelines(content)
            code, last = run(p)
            ok = (code != 0) if expect_fail else (code == 0)
            bad += not ok
            verdict = ("rejected" if code else "accepted") + (", as it should be" if ok else ", WRONG")
            print(f"{name}: {verdict} ({last[:120]})", flush=True)
    print(f"{len(cases)} cases, {bad} wrong")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
