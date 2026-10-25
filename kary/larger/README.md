# Five to fifteen letters

Lueker's automaton improves Dančík's 1994 upper bounds (thesis, Figure 4.6) for every
alphabet size in his table, already with history length h = 3. These certificates use no
floating point in the check: the vector from the power iteration is read as exact
rationals, the class weights are exact rationals at the decimal z, and every ratio
(A_z v)_q / v_q is computed exactly (`certify_exact.py`). λ is the largest ratio and μ the
smallest, so μ ≤ ρ(A_z) ≤ λ, and the bound 2 log(kz)/log λ is evaluated in interval
arithmetic and rounded up.

| k | Dančík 1994 | h = 2 | h = 3 | classes (h = 3) | z (h = 3) |
|---|---|---|---|---|---|
| 5 | 0.664428 | 0.663835 | **0.653115** | 1,045 | 0.1063 |
| 6 | 0.629316 | 0.628671 | **0.618257** | 3,006 | 0.0912 |
| 7 | 0.600184 | 0.599487 | **0.589457** | 7,315 | 0.0803 |
| 8 | 0.575407 | 0.574670 | **0.565042** | 15,688 | 0.0719 |
| 9 | 0.553937 | 0.553170 | **0.543938** | 30,537 | 0.0652 |
| 10 | 0.534855 | 0.534271 | **0.525418** | 55,090 | 0.0597 |
| 11 | 0.517842 | 0.517460 | **0.508965** | 93,511 | 0.0552 |
| 12 | 0.502591 | 0.502359 | **0.494201** | 151,020 | 0.0514 |
| 13 | 0.488800 | 0.488683 | **0.480840** | 234,013 | 0.0480 |
| 14 | 0.476198 | 0.476210 | **0.468662** | 350,182 | 0.0452 |
| 15 | 0.464619 | 0.464765 | **0.457492** | 508,635 | 0.0426 |

The reduced automaton has 10 states for h = 2 and 316 states for h = 3, for every k from 5
to 15. With h = 2 the bound beats Dančík's for k ≤ 13 only. The same program gives exact
certificates for three and four letters at h = 2 and 3 (0.764972 and 0.754807 for k = 3,
0.707648 and 0.696851 for k = 4), the same values as the directed-rounding verifier in `../`.

## Certificates

`certs/cert_k{k}_h{h}.txt` (k = 3 to 15, h = 2 and 3) holds each certificate: k, h, z, the
list of states (letter-renaming representatives of `../review/rv_machine.py`), the vector v
with every entry an exact double written in hexadecimal, the exact rationals λ and μ, λ
rounded up and μ rounded down to 17 decimals, the bound rounded up, and the certified value.

`check_cert.py` checks a certificate file from the file alone. For every listed state and
every input class it computes the transition with `rv_machine.delta`, takes the representative
of the target and looks it up in the list, so it checks that the list contains the start
state and is closed under the transitions; it then recomputes every ratio in exact
arithmetic, checks that the largest and smallest are the λ and μ in the file, and checks the
bound in interval arithmetic. It does not use the automaton build of `certify_exact.py`. All
26 files pass (`logs/check_cert.log`). `tamper_test.py` makes six damaged copies of one
certificate (z changed, the certified value lowered, a state removed, a state not in
representative form, an entry of v changed, λ changed) and checks that each is rejected and
the undamaged file accepted (`logs/tamper_test.log`).

```
python3 check_cert.py certs/cert_k15_h3.txt 8    # 8 worker processes; about 10 minutes
```

The certificates and `logs/certify_small.log`, `logs/certify_big.log` come from the run of
25 October 2026 (`rerun_all.sh`) at the values of z found earlier; a first run the same
day gave the same exact λ, μ and v for all 26. The older logs `logs/certify_h2.log` and
`logs/certify_h3*.log` are the runs that chose z by golden-section search; they print λ and
μ rounded to the nearest 15 decimals, so the signs < and ≥ there can be off in the last
digit. The exact values in `certs/` give the same bounds.

## Checks

- The automaton and weights come from the review implementation `../review/rv_machine.py`.
  `cross_check.py` compares it with the reference implementation `../kref.py`, which finds
  the input classes from the definition of a minimal match pair: the classes agree for every
  k and h used, and every transition agrees for h = 2 (k = 5..15) and h = 3 (k = 5, 6), and
  on random samples of 60, 60, 60, 30 and 20 states for h = 3, k = 7..11
  (`logs/cross_check.log`).
- `cross_check_fast.py` (k = 12..15, h = 3) uses the fact that the 316 states are the same
  for every k from 5 to 15: it samples states from the k = 5 automaton and checks that both
  implementations give the same transitions and never leave that set
  (`logs/cross_fast_k12.log` to `logs/cross_fast_k15.log`).
- `../kcover.py`: canonical decompositions of random strings are never rejected
  (`logs/cover_largek.log`).
- `larger_k.py` is the floating-point scan that first showed the improvement
  (`logs/larger_k.log`); it is not part of the certificates.

These computations have not been repeated by anyone else. The code was written with the
help of a large language model.

## Reproducing

```
cd kary/larger
python3 certify_exact.py 5 3 0.1063     # needs numpy and mpmath; seconds
python3 certify_exact.py 15 3 0.0426    # half an hour to an hour on a laptop
python3 certify_exact.py 13 2 0.0462    # h = 2: instant
```
