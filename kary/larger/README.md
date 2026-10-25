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

The reduced automaton has 10 states for h = 2 and 316 states for h = 3, for every k ≥ 5.
With h = 2 the bound beats Dančík's for k ≤ 13 only. Full output: `logs/certify_h2.log`,
`logs/certify_h3.log`, `logs/certify_h3_k10.log` to `logs/certify_h3_k15.log`.

## Checks

- The automaton and weights come from the review implementation `../review/rv_machine.py`.
  `cross_check.py` compares it with the reference implementation `../kref.py`, which finds
  the input classes from the definition of a minimal match pair: the classes agree for every
  k and h used, and every transition agrees for h = 2 (k = 5..13) and h = 3 (k = 5, 6), and
  on random samples of 60, 60, 60, 30 and 20 states for h = 3, k = 7..11
  (`logs/cross_check.log`).
- `cross_check_fast.py` (k = 12..15, h = 3) uses the fact that the 316 states are the same
  for every k ≥ 5: it samples states from the k = 5 automaton and checks that both
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
python3 certify_exact.py 15 3 0.0426    # about half an hour
python3 certify_exact.py 13 2 0.0462    # h = 2: instant
```
