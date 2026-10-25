# Review of the disk-streamed k-letter code (kbig, kbverify, kbcertify, checks), 2026-10-25

Reviewer: fresh-context opus-reviewer, 12:46 to 17:35. All work, scratch programs and logs are in
`kary/review-big/`; nothing outside it was written (the implementer's big/, data/ and logs were
only read; my repro and red controls used my own copies).

VERDICT (code): SAFE. It was SAFE WITH EDITS at 13:32; every required edit was made by 15:30 and I
  verified each one (see "Status of the fixes"). kbig build makes exactly kary.c's table (16
  byte-identical cases, ids proven equal, re-run with the final binary), kbverify's streaming and
  directed rounding hold under every attack I made, and my own verifier (different streaming,
  emulated rounding, plus a rounding-mode-free bound) agrees with kbverify bit for bit.
VERDICT (certificate k = 3, h = 7, gamma_3 <= 0.737495): SAFE. kbcertify.py re-run, my verifier
  and the exact weight check reproduce it (lambda_cert 0x1.2f970b6938e14p-2 bit for bit), the
  SHA-256 of table, vector, weights and wtypes equal SHA256_k3_h7.txt, and finding 1 is closed for
  this table: all four quarter comparisons say "RESULT: no difference" and "exit 0", and together
  they cover 28,966,771 states (= N) and 9,520,016,086 transitions (= the build's accepted count).
VERDICT (certificate k = 4, h = 6, gamma_4 <= 0.679701): SAFE, CONDITIONAL on the three remaining
  quarter comparisons. Everything else is reproduced (lambda_cert 0x1.047f14d17bdedp-2 bit for bit,
  exact weights 28/28, SHA-256 equal to SHA256_k4_h6.txt), and quarter q0 is clean. It becomes
  SAFE when big/logs/rvcmp_k4_h6_q1..q3.log each end with "RESULT: no difference" and
  "phase P exit 0", compare 2,125,824 (q1), 2,125,729 (q2) and 2,124,800 (q3) states, and the four
  quarters' transitions add up to 27,720,529,578 (the build's accepted count); q0 had 6,920,046,139.
  Until then, do not call the k = 4 bound certified in public text.

## Work log (newest last)

All commands were run from `kary/` or `kary/review-big/` with `nice -n 10`. Logs named below are in
`kary/review-big/`.

1. Hashes of every reviewed file at the start: `sha256_at_start.txt`. My compiles of kbig.c
   (-O3 -march=native), kbverify.c and kverify.c (-O2 -ffp-model=strict -march=native) and
   rvcmp_big.c have the same `__text`, `__cstring` and `__const` sections as the binaries in kary/
   (`otool -tV` and `otool -s` compared; only UUID and signature bytes differ).
2. kbig.c changed at 12:52 (the `KBIG_QOS` variable picks the dispatch QoS class). Deleting that edit
   from the new file gives a source that compiles to exactly the `__text` I reviewed
   (`snap/kbig_v1_reconstructed.c`), so the edit changes scheduling only. The k = 3, h = 7 build
   (started 12:38) runs the earlier binary, the k = 4, h = 6 build (12:53) the new one; the new
   binary equals my compile of the current kbig.c. The team lead confirmed the edit is QoS only.
   With my compile of the current kbig.c and KBIG_QOS=user (the setting the (4,6) build and the
   eigen runs use), (3,4), (3,5), (4,3), (4,4) and (3,6) are again byte-identical to kary.c's
   tables (`repro_rv_qos.sh`, `log_repro_qos_user.txt`).
   Second edit, 13:43 (`snap/kbig.c.1343`, diff against `snap/kbig.c.1252`): only the docstring and
   the checkpoint branch of load_start change. Each z now starts from (1 - t) old[a] + t old[b], with
   a the checkpoint entry of largest z0 <= z, b of smallest z0 >= z, t = (z - z0[a]) / (z0[b] - z0[a])
   only when z0[b] > z0[a] and 0 otherwise, so t is in [0, 1], there is no division by zero when z
   equals a checkpoint z or the checkpoint repeats a z (its unused slots do), and entries stay
   positive. kary/kbig (13:43:26) equals my compile of it. Tested at (3,5) from my own checkpoint
   (z0 = 0.1948, 0.1950 x 7) with z = 0.1947, 0.1948, 0.1949, 0.1950, 0.1951 (`eig35_interp/`):
   weights 0, 0, 0.5, 0, 0 as intended, convergence to 1e-12 in 43 passes, and the z = 0.1949 vector
   gives CERTIFIED 0.743282 through kbcertify.py with my verifier equal bit for bit. The build path
   is untouched, and the vector only matters through kbverify, which checks it.
3. kbverify disassembly (`otool -tV kary/kbverify`): `fesetround(0x400000)` (FE_UPWARD) then the
   1 + 1e-30 self-test, loop `ucvtf; fmul; fmul; fadd` in the source order ((cnt * w) * v, then
   sum + term), `fdiv`; then `fesetround(0x800000)` (FE_DOWNWARD), both self-tests, the same loop;
   then `fesetround(0)`. No fused multiply-add anywhere in kbverify or kverify.
4. Byte-for-byte reproduction with the current kbig (my copy of big/repro.sh writing to
   review-big/repro, so the implementer's big/repro is untouched): `log_repro_small.txt`,
   `log_repro_big.txt`. All 16 cases identical to kary.c's tables in data/ (k = 2, h = 2..8;
   k = 3, h = 2..6; k = 4, h = 2..5), including (3,6), (4,5) and (2,8), where the hash table is
   rehashed during the build and the parallel lookup misses 15,057,583, 7,086,695 and 9,567,538
   targets respectively.
5. Three verifiers agree bit for bit (`log_small_verify.txt`, `log_headline_verify.txt`): kverify's
   recorded numbers in data/certify_*.out, kbverify on my kbig rebuilds with the data/ vectors, and
   my own verifier rvb_verify (below), at (3,2..6) and (4,2..5). (3,6) and (4,5) stream 3 chunks.
6. kbverify with the chunk cap lowered to 700, 1,000 and 4,099 records (scratch copies in mut/):
   up to 19,189 chunk boundaries at (3,5); lambda_cert and rho_lower bit-identical to the 1-chunk
   run at (3,4), (3,5), (4,3), (4,4). Cap 100 at (4,4) (a state has 635 records) is refused
   ("state 0 has too many records", exit 2). `log_kbverify_lowcap.txt`.
7. Red controls on corrupted copies of the (3,4) table, vector and wtypes (`redctl.py`,
   `log_redctl_k3h4.txt`). kbverify refuses: truncated file, tgt = N, decreasing offsets, a changed
   header weight type, v_s = 0, v_s = NaN, a vector file 8 bytes long, a wtypes file with an extra
   line. kbverify accepts, as designed, tables that are wrong but well formed (changed weight type,
   changed target, multiplicity minus 1, one record moved to the next state by an offset, one record
   deleted with header and offsets made consistent, keys[0] changed); check_table_big.py (all 5,116
   states) and rvcmp_big (stride 1) flag every one of those.
8. Weights at the new (k, h) (`wt/check_weights_new_h.py`, `log_weights_new_h.txt`): kbig's
   make_inputs (harness, no build) gives 975 inputs and 24 weight types at (3,7), 9,748 and 28 at
   (4,6), the same type sets as kref.py's classes by definition; brute-force counts of all
   787,143 (k = 3, sizes <= 18) and 6,571,684 (k = 4, sizes <= 14) minimal match pairs per class
   equal the closed-form coefficients; kweights.py write at z = 0.2126 and 0.1575 passes the earlier
   referee's exact Fraction check (24/24, 28/28).
9. The eigen phase: run on my (3,5) rebuild (`eig35/`), the table's SHA-256 is unchanged afterwards
   (it is opened read-only); with the chunk cap lowered to give 121, 5 and 3 chunks, the vectors and
   the checkpoint are bit-identical to the 1-chunk run (`log_eigen_lowcap.txt`), so the
   two-buffer reader and the buffer parity across passes are right. Its z = 0.1948 vector, put
   through kbcertify.py, gives CERTIFIED 0.743282 (same as kary.c's), and my verifier reproduces
   lambda_cert bit for bit.

10. kbig vs kary.c without "klue" (`sections_cmp.py`, `log_sections_cmp.txt`): header, Rec[],
    off[] and key[] of my kbig rebuilds read straight from the KBIG layout equal data/ktab_* at
    (3,4), (3,5), (3,6), (4,4), (4,5), (2,8). The block of kbig.c said to be copied from kary.c
    (kbig.c 75 to 235 of the 12:38 file) is textually identical to kary.c 60 to 220 (diff).
11. End-to-end rehearsal on the implementer's own kbig pipeline run at (3,6) (big/k3h6, z = 0.2048),
    with `check_cert.sh 3 6` (FULL=1), which is the script I would run on the (3,7) and (4,6)
    certificates (`log_check_cert_k3h6_dryrun.txt`): kbcertify.py re-run into review-big/c36
    reproduces "CERTIFIED: gamma_3 <= 0.739975", lambda_cert 0x1.127d82e722abcp-2; the regenerated
    weights equal big/k3h6's; rvb_verify gives the same lambda bit for bit and the a-priori route the
    same 0.739975; rv_theirweights 20/20; rvcmp_big's certificate mode (every transition recomputed
    by the earlier referee's own delta, kbig's records unused) gives 0x1.127d82e722abdp-2, one ulp
    above (401 s, 4 threads, loaded machine).
12. Covering test re-run with a new seed (`log_kcover_new.txt`): k = 3, h = 7: 4,000 pairs,
    46,478 canonical steps (2,333 truncated inputs); k = 4, h = 6: 4,000 pairs, 41,575 steps
    (3,311 truncated); no rejection.
13. The finished k = 3, h = 7 table (build.out at 13:24: N = 28,966,771, accepted transitions
    9,520,016,086, records 9,427,252,448 > 2^33, file 75,881,488,444 bytes = 36 + 20 NW + 8R +
    8(N+1) + 8N exactly; wtypes file byte-identical to my make_inputs harness), sampled with a
    different sample and block phase from big/checks_k3h7.sh (`sampled_k3h7.sh`,
    `log_sampled_k3h7.txt`):
    - rvcmp_big preamble over ALL 28,966,771 stored states: 0 share an orbit, 0 not in canonical
      form, start state index 0;
    - rvcmp_big on every 101st block of 1,024 states (phase 50): 286,720 states, 95,729,080
      transitions recomputed independently, 0 states with differing records, 0 targets missing
      (35 s on 4 threads, so a full pass is about an hour on this loaded machine);
    - check_table_big.py, 600 random states (seed 2027) plus the first and last 200: 1,000 states,
      304,413 transitions, 0 mismatching states.

## My own verifier (review-big/rvb_verify.c, rvb_emu.h, rvb_bound.py)

Written for this review, different from kbverify by design:
- parses the KBIG file with pread at fixed offsets; checks size, offsets, records, vector, wtypes
  against the header, weights (finite, up >= down > 0) and that keys[0] is the start state;
- streams records in fixed 2^21-record blocks that ignore state boundaries (rows straddle blocks),
  one contiguous state range per pthread;
- never changes the rounding mode: directed rounding is emulated in round-to-nearest with
  error-free transformations (fma for products and the division remainder, TwoSum for sums),
  in kbverify's order of operations. `rvb_emu_test.c` checks the emulation against the hardware's
  FE_UPWARD and FE_DOWNWARD on 120,000,000 random operations: 0 mismatches
  (`log_emu_test.txt`);
- also a bound that needs no directed rounding at all: with round-to-nearest, q_s <= qn_s
  (1 - 2^-53)^-(n_s + 2) for a row of n_s nonnegative terms, so lambda_apriori =
  max qn_s (1 - 2^-53)^-(nmax + 2), evaluated exactly in rvb_bound.py, which also evaluates the
  bound in 120-digit interval arithmetic and compares the printed 6-decimal value as an exact
  decimal.
At (3,6) and (4,5) its lambda equals kbverify's bit for bit, lambda_apriori is 5.7e-14 and
2.4e-13 above it, and both routes give 0.739975 and 0.683622 (`log_rvb_bound_headline.txt`).

## Code reading (against kary.c, kverify.c and Lueker 2009, Section 3.2, Figure 7)

- kbig.c lines 75 to 235 of the 12:38 file ("copied unchanged"): wtype_of, make_inputs, delta, the
  hash table (mix, rehash, find_or_add), KW and cmp_kw are textually identical to kary.c lines 60
  to 220 (`diff`, work log item 10); `pack`, the types and the double-precision weight functions
  match too. delta is Figure 7 with relabelling, as reviewed before.
- Id order (kbig.c build, lines 300 to 350 of the 12:38 file): in kary.c, ids are assigned by
  find_or_add in (state, record) order; a call for a key already present inserts nothing. kbig
  looks every target up read-only in parallel against the table as it stood at the start of the
  batch, then calls find_or_add, in the same (state, record) order, exactly for the keys that
  were absent; a key absent at the lookup but inserted earlier in the same batch is found by that
  call. So the sequence of insertions, and with it every id and every rehash, is kary.c's. Batch
  boundaries do not enter the numbering. Empirically confirmed by the 16 byte-identical tables.
- Races: during dispatch_apply nothing inserts, so `states`, `hkeys`, `hvals`, `hmask` are not
  written or reallocated while workers read them; dispatch_apply returns only after every block,
  so the sequential phase sees all of `buf`, `mult`, `ids`, `cntm`. Ids stay valid across a rehash.
  Each worker writes only its own slices.
- Overflow: N < 2^32 is enforced (find_or_add), NOT_FOUND = UINT32_MAX can never be an id;
  R, offsets, record counters and all file offsets are uint64/off_t in kbig, kbverify,
  check_table_big (numpy uint64) and rvcmp_big. At (3,7) R passes 2^32 (6.96e9 records by 13:15),
  so that path is exercised only by the new runs; no 32-bit quantity carries a record index.
  Multiplicities are uint16 (<= NI = 9,748 at (4,6), checked NI <= 65535); weight types uint8 /
  uint16 (<= 256, checked).
- File layout: header 36 + 20 NW bytes, Rec[R], off[N+1], key[N]; every reader (open_kbig,
  kbverify, check_table_big, rvcmp_big, my rvb_verify) uses the same offsets and checks the file
  size. The header is written last (fseek to byte 12 after fflush), so an interrupted build leaves
  N = R = 0 and a file every reader refuses.
- kbverify streaming: chunks are whole states [s0, s1) with off[s1] - off[s0] <= 2^27; records are
  read sequentially from Rec[0], so chunk c starts exactly at off[s0]; every record is range-checked
  as it is read; per-state sums are formed in kverify's order, so the chunking cannot change any
  number (confirmed with thousands of chunk boundaries, item 6).
- Rounding: -ffp-model=strict, FE_UPWARD and FE_DOWNWARD set and self-tested inside every block;
  all operands nonnegative, so every upward-rounded product, sum and quotient is an upper bound and
  the downward ones lower bounds; no contraction; the max/min reductions are exact. A NaN would be
  ignored by the `>` comparisons, but no NaN can arise (v > 0 finite is checked, weights > 0, no
  subtraction; an infinite weight gives lambda = inf, which kbcertify's `lam < 1` rejects).
- kbcertify.py is kcertify.py with the binary name changed (diff: 4 lines, all naming). Weights are
  regenerated from the vector's own z on every run, so a stale or edited weights file cannot leak
  into a certificate; kbverify ties the weights to the table through the wtypes file it checks
  against the header.
- The eigen phase opens the table read-only and writes only kvecs_*.bin (via .tmp and rename) and
  kvec_*_z<z>.bin; it cannot change any file the certificate reads except through the pipeline's
  explicit `cp` of the chosen vector, and any positive vector gives a valid certificate.

## Findings

1. [MAJOR for the certificates, not a code defect] A certificate is only as good as the table,
   and at the new sizes the table is checked by sampling only.
   kbverify proves A v <= lambda v for whatever well-formed table it reads; it cannot know whether
   the table is Lueker's machine. Repro (`log_redctl_k3h4.txt`): delete one record from the (3,4)
   table and fix the header and offsets (case k3h4_drop_record) -> kbverify exits 0 with the
   unchanged lambda_cert 0x1.9699a6593d1cbp-3, i.e. a "CERTIFIED" line for a matrix that is not
   the machine. Only check_table_big.py and rvcmp_big catch it, and only if they look at that
   state. The plan in big/checks_k3h7.sh checks 5,000 random states plus 400 (check_table_big.py)
   and every 20th block of 1,024 states (rvcmp_big, 5 %), so about 95 % of the rows of a 76 GB
   file written in 46 minutes on a machine under memory pressure (load 37 to 78, 7 to 10 GB
   compressed, heavy swapping during the build) are never compared with an independent
   implementation. The (4,6) plan (big/after_k4h6_build.sh: check_table_big.py on 2,000 states,
   rvcmp_big on every 40th block) covers about 2.5 %. A single wrong row invalidates the proof even when the number is still true:
   the vector is tuned to the stored matrix, and the true row would exceed lambda_cert.
   Demonstrated end to end at (3,4) (`redctl/k3h4_drop_record/`): kbig eigen on the table with one
   deleted record, then kbcertify.py, prints "CERTIFIED: gamma_3 <= 0.747887" with
   lambda_cert = 0x1.96999b7c886a8p-3 = 0.19853517030519652, which is BELOW the true machine's
   rho_lower 0.19853525123431212, i.e. below its Perron root; the same vector and weights on the
   true table give lambda_cert = 0.19999500875954684. The printed 6 decimals survived only because
   one record out of 366,370 moves lambda by 4e-7 relative.
   The evidence that the CODE is right is strong (item 4: identical tables at 16 smaller (k, h),
   the id argument above), but it does not cover I/O or memory faults in this particular file.
   The earlier certificates (3,6) and (4,5) had a full referee comparison
   (review/log_rvcmp_headline.txt: "compared all ... 0 states with differing records").
   Required before the new bounds are called certified: run rvcmp_big on the whole table
   (RVCMP_STRIDE=1; measured on the finished (3,7) table: 1/101 of it in 35 s on 4 threads on this
   loaded machine, so about an hour for all of it; at (4,6) roughly 3.5 times the transitions),
   or better its certificate mode with the final vector and
   rv_theirweights.py's weights, which recomputes every transition independently and so certifies
   the vector without using kbig's records at all (`log_rvcmp_cert_k3h5.txt` shows it running on
   a KBIG file: 0x1.e29d481fee276p-3, a few ulps above kbverify's merged sum, as expected; two
   passes, so about twice the times above). And record the SHA-256 of the table, vector and
   weights that were certified, as for the earlier runs (data/SHA256_kary.txt).

1b. [MAJOR, operational, no effect on soundness; RESOLVED at 13:40] The post-build pipeline
   deadlocked on its own wait.
   big/after_k3h7_build.sh:4 and big/checks_k3h7.sh:5 wait with
   `while pgrep -f "kbig build 3 7"`. pgrep -f matches full command lines, and the implementer's
   own background wait loops (zsh wrappers whose command line contains the same quoted string,
   pids 73677 and 81565, plus 14318) match it, and match each other; BSD pgrep excludes only its
   own ancestors (man pgrep, -a). Observed at 13:26: the build had exited at 13:24 with its final
   line in build.out, `pgrep -f "kbig build 3 7"` still listed 6852 14318 73677 81565, and
   after_k3h7_build.sh (pid 79869) and checks_k3h7.sh (79937) were still in their loops, so no
   eigen phase and no table check had started. (6852 is my own background job, whose heredoc
   contained the string; it ends by itself and does not cause the deadlock.) Reported to the team
   lead at 13:27. Fix: wait on the build's PID or on build.out's final line, never on pgrep -f with
   a string that other shells also contain. Status: the lead cleared it; the (3,7) pipeline started
   at 13:40:24 and big/after_k4h6_build.sh (13:40) waits with `kill -0 $BUILD_PID`.

2. [MINOR] rvcmp_big exits 0 when its comparison fails.
   Repro: `./rvcmp_big 3 4 review-big/redctl/k3h4_tgt_other/kbig_k3_h4.bin; echo $?` prints
   "1 states with differing records" and exit status 0; the same for the key0 case
   ("start state index -1", "9 targets not among kbig's states", exit 0). rvcmp_big.c:253-258
   (inherited from review/rvcmp.c:240-243) always returns 0 after the comparison and never checks
   dup, noncanon or the start index. big/checks_k3h7.sh:7-8 only redirects the output, so nothing
   would stop on a failure; whoever reads the logs must read the counts. Fix: return 1 when
   totm, totn, dup or noncanon is nonzero or the start index is not 0, and have the check script
   test the exit codes. (check_table_big.py does exit 1 on any mismatch.)

3. [MINOR] The byte-for-byte evidence in big/repro was produced by an earlier kbig.
   big/repro/build_k3_h6.err line 3 reads "ids assigned after the parallel lookup: 15057583
   (= states - 1)" (wrong label: (3,6) has 1,711,682 states); the current kbig.c:376 prints
   "records whose target was not yet known at the parallel lookup". The binary that built (3,7)
   was compiled at 12:38:40, after those files (12:35 to 12:37). My re-run with the current source
   (item 4) gives identical tables, so nothing is wrong, but the repro logs that get published
   should be regenerated with the final binary.

4. [NIT] kbcertify.py:40-46 (inherited from kcertify.py) checks the double `six` against the
   interval end but prints `six` rounded to 6 decimals. If the interval end lies within about
   1e-17 above a 6-decimal value d whose nearest double is above d, the assert passes and the
   printed bound d is below the proven value. Chance about 1e-11 per run; rvb_bound.py checks the
   exact decimal and confirms the printed values at (3,5), (3,6), (4,5). Fix: assert
   `Fraction(f"{six:.6f}") >= exact(B.b)`.

5. [NIT] big/pipeline.sh:41 tees kbcertify's output into certify_k<k>_h<h>.out even when kbcertify
   fails (set -o pipefail stops the script, but the partial file stays). A certify file is valid
   only if it ends with a CERTIFIED line; worth saying wherever the files are cited.

6. [NIT] kbig.c:369 (wtypes file), the start-vector writer in `project`, and the mallocs in
   read_prev_keys/project are unchecked; a failure crashes rather than writing a wrong file, so no
   soundness effect. rvcmp_big with RVCMP_STRIDE=0 divides by zero (crash). In the 13:43 load_start,
   a checkpoint whose z values were all NaN would leave a = b = -1 and read old[s * MZ - 1]
   (out of bounds); checkpoints are written by kbig itself, so this needs a corrupted file. The eigen output name
   kvec_..._z%.4f.bin collides for z values equal to 4 decimals; the file holds the exact z, so a
   certificate stays self-consistent.

No finding touches the id order, the parallel phase, the streaming chunk boundaries, the rounding
modes or the file layout: those were attacked directly (items 4 to 9 of the work log) and held.

## The k = 3, h = 7 certificate (big/k3h7/certify_k3_h7.out, 15:08)

Claimed: CERTIFIED gamma_3 <= 0.737495, z = 0.2129 (0x1.b404ea4a8c155p-3), lambda_cert =
0x1.2f970b6938e14p-2 = 0.29647462682461989, rho_lower = 0x1.2f970b68f84efp-2, N = 28,966,771,
R = 9,427,252,448. Checked so far (15:12):
- kvec_k3_h7.bin is byte-identical to kvec_k3_h7_z0.2129.bin; its header is k = 3, h = 7, N = N of
  the table, z = the double nearest 0.2129 (0 < z < 1/3).
- The bound from lambda_cert, recomputed in 120-digit interval arithmetic with exact rational z and
  lambda: 0.73749426501961209749; the printed 0.737495 is >= it as an exact decimal and 7.3e-7
  above it. At rho_lower: 0.73749426498957; lambda_cert / rho_lower - 1 = 4.95e-11, so this is
  the best this machine gives at this z, to 10 digits.
- The eigen logs (stage 1: 58 passes, 3,371 s; stage 2: 31 passes, 1,875 s, started from the
  interpolated checkpoint, weights as designed) show no error lines; big/pipeline2.sh, used for
  (4,6), differs from pipeline.sh only in the TOL1/TOL2 variables.
- Plausibility: the gains per step of h for k = 3 run 0.0102, 0.0069, 0.0046, 0.0033, 0.0025.
`review-big/check_cert.sh 3 7` (15:10 to 15:44, `log_check_cert_k3h7.txt`; no certificate mode,
as the lead asked):
- kbcertify.py re-run into review-big/c37 prints the same certificate: lambda_cert
  0x1.2f970b6938e14p-2, rho_lower 0x1.2f970b68f84efp-2, CERTIFIED gamma_3 <= 0.737495.
- The weights it regenerates at z = 0.2129 are byte-identical to big/k3h7/weights_k3_h7.txt, and
  rv_theirweights.py (the earlier referee's exact Fraction weights) accepts all 24 (up >= exact >=
  down, adjacent doubles).
- rvb_verify, my own verifier: all 28,966,771 states and all 9,427,252,448 records summed, 0 bad
  records, 0 empty rows, 0 tiny operands, state 0 is the start state, vector entries in
  [0x1.079e36ca4a0dfp-5, 1]; lambda 0x1.2f970b6938e14p-2 and rho 0x1.2f970b68f84efp-2, bit-identical
  to kbverify's, with directed rounding emulated in round-to-nearest. The rounding-mode-free a-priori
  bound (round-to-nearest maximum times (1 - 2^-53)^-(964 + 2)) is 1.06e-13 above lambda_cert and
  gives the same 0.737495; the printed 0.737495 is >= the exact upper end of the interval.
- SHA-256 (shasum -a 256) of the four certified files equal the implementer's
  big/k3h7/SHA256_k3_h7.txt: table 51d9d19bd9d4130dada1212730277a59a61b99d9b9d77649a99a11cf20845619,
  vector 766c3dda9238a664ea324aeca3051a209e4eae308671f3dad605af1c765ee35b, weights
  d8e63f1d1059206b40b10c41df857ef71f83d29ad24f98cabb804ef267229931, wtypes
  eb4373e2b0380cd6e47e29640641e0ff72879dd5d1ba3251f2b18f306129795a.
Finding 1 for this table (four quarter runs, big/logs/rvcmp_k3_h7_q0..3.log): see the comparison
status at the end of this file.

## The k = 4, h = 6 certificate (big/k4h6/certify_k4_h6.out, 16:09)

Claimed: CERTIFIED gamma_4 <= 0.679701, z = 0.157 (0x1.4189374bc6a7fp-3), lambda_cert =
0x1.047f14d17bdedp-2 = 0.25439102676952424, rho_lower = 0x1.047f14cf05242p-2, N = 8,502,177,
R = 21,142,620,884 (> 2^34; file 169,277,002,508 bytes = 36 + 20 NW + 8R + 8(N+1) + 8N exactly).
Checked without the table:
- kvec_k4_h6.bin is byte-identical to kvec_k4_h6_z0.1570.bin; header k = 4, h = 6, N of the table,
  z = the double nearest 0.157 (0 < z < 1/4); wtypes_k4_h6.txt is byte-identical to my
  make_inputs harness output (28 types, all brute-force checked in item 8).
- The bound from lambda_cert in 120-digit interval arithmetic: 0.67970046757835947; the printed
  0.679701 is >= it as an exact decimal (5.3e-7 above). At rho_lower: 0.67970046729844;
  lambda_cert / rho_lower - 1 = 5.6e-10 (stage 2 ran to 1 + 1e-9).
- Eigen logs without error lines; stage 2 started from the interpolated checkpoint as designed; the
  first attempt (big/k4h6/aborted-run1/) was stopped to switch to the aligned reader.
- The implementer's check_table_big.py: 2,400 states, 7,504,681 transitions, 0 mismatches.
- Plausibility: gains per step of h for k = 4 run 0.0108, 0.0078, 0.0054, 0.0039.
`review-big/check_cert.sh 4 6` (16:24 to 16:50, `log_check_cert_k4h6.txt`; my verifier now reads
page-aligned spans like the implementer's fixed readers, with the same outputs as before at (3,4),
(3,5), (4,4) and on the (3,6) certificate; the table is hashed with sha256nc, F_NOCACHE, which
gives the same digests as shasum on two test files, and kary/sha256nc equals my compile):
- kbcertify.py re-run into review-big/c46 prints the same certificate (lambda_cert
  0x1.047f14d17bdedp-2, rho_lower 0x1.047f14cf05242p-2, CERTIFIED gamma_4 <= 0.679701); the
  regenerated weights equal big/k4h6's; rv_theirweights.py accepts all 28.
- rvb_verify: all 8,502,177 states and 21,142,620,884 records, 0 bad records, 0 empty rows, 0 tiny
  operands, state 0 is the start state, vector entries in [0x1.744793d166f69p-4, 1]; lambda and
  rho bit-identical to kbverify's. The rounding-mode-free a-priori bound (nmax = 7,481 records in
  one row) is 8.3e-13 above lambda_cert and gives the same 0.679701.
- SHA-256 equal to big/k4h6/SHA256_k4_h6.txt: table
  549ef26d620e906f073eec9f594b8552738634faebf92777264fe916d9b9bc4f, vector
  4ee380850c2d940226bcc6547f69ec5f5b983eb11fb09cf8c1fde320af28a19e, weights
  09b8c2d2295cec4fb46ff7534570c34e18c53f39594715242c19adde4c2a5d2d, wtypes
  e39ba391b1d584ac07fc578ab8a559d8e5b988aae260d2ddb4d2a57538f05fb3.

## Status of the fixes (reported by the team lead at 13:50, checked by me)

- Finding 1, ACCEPTED, in progress: the 5 % sample run was stopped. Since 14:05 the full comparison
  is big/compare_full.sh k h: rvcmp_big four times with RVCMP_STRIDE=4 and RVCMP_PHASE=0..3, one
  after the other, each log ending with its RESULT line and "phase P exit S"; for both tables it
  starts from big/after_k4h6_build3.sh when the (4,6) build ends (the paused single run, pid 40667,
  is gone). I checked the scheme (`log_quarters_small.txt`): block b goes to phase b mod 4, so every
  state is compared exactly once; at (3,5) and (4,4) the four quarters compare 98,041 and 9,427
  states in total (= N), and on the deleted-record (3,4) table only phase 2, which holds the
  damaged state 2558, exits 1 with "RESULT: DIFFERENCES FOUND". Closing condition: all eight logs
  say "RESULT: no difference" and "phase P exit 0", and the four "compared X of N" counts add up to
  N. For (3,7) they must be 7,241,728, 7,241,728, 7,241,728 and 7,241,587 (28,288 blocks, the last
  one of 883 states in phase 3); compare_full.sh does not add them up itself.
- Finding 1b, FIXED at 13:40 (waiter shells killed; the (4,6) script waits with kill -0 on the PID;
  only one after_k4h6 script is running).
- Finding 2, FIXED (rvcmp_big.c 13:49, binary equals my compile). Compare mode now prints
  "RESULT: no difference" or "RESULT: DIFFERENCES FOUND" and exits 1 when any state differs, a target
  is missing, two stored states share an orbit, a state is not canonical or the start state is not
  state 0; RVCMP_STRIDE=0 becomes 1. On my red controls (`log_rvcmp_fix_exitcodes.txt`): control
  exit 0, tgt_other, cnt_dec, off_shift, key0, drop_record exit 1 with the RESULT line, hdr_wtype exit 1.
  Certificate mode keeps its old exits (0, or 3 on a self-test or lookup failure; key0 gives 3,
  `log_rvcmp_certmode_redctl.txt`), and it does not judge its own ratio, so I added
  `rv_certmode_bound.py` to check_cert.sh: the bound from the referee's independent transitions must
  equal the printed certificate. It does at (3,5) and on the (3,6) dry run, and it exposes the
  corrupted-table certificate of finding 1 (referee lambda 0.19999500875954704, bound 0.751291,
  against the printed 0.747887: DIFFERENT, exit 1).
- Finding 4, FIXED (kbcertify.py 13:49): the upper end of B is read exactly from B._mpi_[1], rounded
  up to 6 decimals in Fractions, and the printed decimal string is compared exactly. Repro of the old
  defect and of the fix (`finding4_demo.py`, `log_finding4_demo.txt`): 0.739975 itself has its
  nearest double 4.9e-17 above it; for an interval end 1e-30 above 0.739975 the old step printed
  0.739975 (below the proven value), the new step prints 0.739976. The fixed kbcertify.py gives the
  same output as before on my (3,5), (3,5) z = 0.1949 and (3,6) directories
  (`log_kbcertify_fixed.txt`; the only difference is a relative path in the "wrote" line).
- Finding 3, DONE: big/logs/repro_identity.log (15:29) rebuilds all 16 tables with the final kbig
  binary (SHA-256 1bedda70..., built from kbig.c 4ba8efbd..., both equal to the files in kary/ now)
  and finds them identical, with the (3,6) and (4,5) table hashes equal to data/SHA256_kary.txt;
  big/logs/kbverify_vs_kverify.log has 10 of 10 kbverify results identical to kverify's.
- I/O changes reported at 15:30 (page-aligned reads, so that F_NOCACHE really bypasses the page
  cache; snapshots `snap/kbig.c.aligned`, `snap/rvcmp_big.c.aligned`):
  - kbig.c reader_main reads [a0, a1), the 16 KiB-aligned span around the chunk's records, into a
    posix_memalign'ed buffer of maxr * 8 + 2 pages (the span is at most the chunk plus two pages,
    so it fits; a1 is clamped to the file size) and sets RBUF = RAW + (first - a0) before it marks
    the buffer ready under the mutex. RECBASE = 36 + 20 NW is a multiple of 4, so the records stay
    4-byte aligned. kary/kbig (15:22) equals my compile. With the chunk cap lowered to 121 and 3
    chunks, and at the default, the (3,5) vectors and checkpoint are bit-identical to the old
    reader's 1-chunk run (`eig35_al*/`).
  - rvcmp_big.c compare mode reads the same kind of span per block of 1,024 states (buffer grown
    when needed, rbuf = raw + skip). kary/rvcmp_big (15:24) equals my compile. All 9 red-control
    tables give the right exit (control 0; tgt_oob, wt_other, tgt_other, cnt_dec, off_shift,
    off_nonmono, key0, drop_record 1 with "RESULT: DIFFERENCES FOUND"), and the four quarters
    still add up to N at (3,5), (4,4) and on the drop_record table, where only phase 2 fails
    (`log_rvcmp_aligned.txt`).
  - kbverify.c and kbcertify.py are unchanged since my last check (12:39 and 13:49).

## Criteria check (the claims in the brief, each against my own evidence)

- "build" gives exactly kary.c's table, same numbering and records: PASS. 16 (k, h) byte-identical
  (item 4) and section-identical without klue (item 10); id argument in the code reading.
- Parallel read-only lookup, then sequential insertion in kary.c's order, no races: PASS (code
  reading; the three cases with rehashes and millions of late insertions are identical).
- Eigen phase cannot corrupt what the certificate uses: PASS (read-only open, table hash unchanged,
  multi-chunk vectors bit-identical, item 9).
- kbverify reads every record exactly once and attributes it to the right state: PASS (items 5, 6;
  my block-agnostic verifier agrees bit for bit).
- kbverify's rounding: PASS (disassembly item 3; emulated directed rounding bit-identical; the
  a-priori round-to-nearest route gives the same certified decimals).
- kbverify cross-checks the weight types: PASS (red controls hdr_wtype and wtypes_extra refused).
- kbcertify.py = kcertify.py with kbverify: PASS (diff), NIT 4 inherited.
- check_table_big.py detects wrong rows: PASS for the rows it samples (red controls, item 7).
- rvcmp_big.c = rvcmp.c plus three mechanical changes: PASS (diff: layout, class array, stride, and
  added read checks), MINOR 2 inherited.
- big/repro.sh compares byte for byte: PASS (my copy), MINOR 3 on its saved logs.
- Integer widths and file layout at the new sizes: PASS (code reading; the finished (3,7) file has
  R = 9,427,252,448 > 2^33 and its size equals the header's formula exactly).
- Certificates: (3,7) PASS in full; (4,6) PASS except the three quarter comparisons still running.

## How to check the certificates when they exist

From kary/, once big/logs/rvcmp_k3_h7_full.log ends with "RESULT: no difference":
`review-big/check_cert.sh 3 7` (and later `4 6`). It writes only into review-big/c37 (symlinks to
big/k3h7 plus its own weights) and:
- re-runs kbcertify.py (the fixed one) and checks the regenerated weights equal big/k3h7's;
- runs rvb_verify and rvb_bound.py: lambda bit-equal to lambda_cert, the rounding-mode-free a-priori
  bound, the exact decimal check of the printed value;
- runs rv_theirweights.py (exact Fractions);
- prints the SHA-256 of the certified table, vector, weights and wtypes.
`FULL=1` adds rvcmp_big's certificate mode and rv_certmode_bound.py: an independent certificate of
the same vector from recomputed transitions, whose 6-decimal bound must equal the printed one. It is
optional once the full comparison has passed (about 2 hours at (3,7) on 4 loaded threads).
Rehearsed on the implementer's (3,6) pipeline run (item 11). It reads the 76 GB table (165 GB at
(4,6)) three times (kbverify, rvb_verify, shasum); FULL=1 adds computation but no further reads of
the records. Run it when the eigen phase is not streaming the same disk.

## Attacks attempted (all run; results in the logs named)

- Rebuilt 16 tables with kbig and compared byte for byte and section by section with kary.c's
  (repro_rv.sh, sections_cmp.py): identical.
- Compiled every reviewed program myself and compared machine code with the binaries used.
- Traced the 12:52 kbig.c edit by reconstructing the pre-edit source and compiling it.
- Disassembled kbverify to confirm the rounding-mode calls, self-tests, operation order and the
  absence of FMA.
- kbverify against kverify's recorded certificates at 9 (k, h): bit-identical.
- kbverify with chunk caps of 100, 700, 1,000, 4,099 records: identical numbers, oversized state
  refused.
- Wrote an independent verifier (emulated directed rounding, verified against the hardware on
  1.2e8 operations; rounding-mode-free a-priori bound; block streaming across state boundaries):
  bit-identical lambda and rho at 9 (k, h) and on two kbig-pipeline certificates.
- 14 corrupted inputs and an unmodified control (truncation, out-of-range target, changed weight type, changed target,
  multiplicity minus 1, shifted offset, decreasing offsets, keys[0], header weight type, deleted
  record with consistent header, v = 0, v = NaN, long vector file, extra wtypes line) against
  kbverify, check_table_big.py and rvcmp_big; plus the full pipeline on the deleted-record table,
  which certifies a lambda below the true Perron root (finding 1).
- Weights at the new (k, h): input counts and type sets against the definition, brute-force series
  of 7.4 million minimal match pairs, exact Fractions at plausible z.
- Eigen phase: table hash before and after, 3 chunkings with bit-identical vectors, its vector
  through kbcertify and my verifier.
- Covering test at (3,7) and (4,6) with a new seed.
- The finished (3,7) table: header and size formula, wtypes, all-state orbit, canonical-form and
  start-state checks, 1 % independent recomputation (rvcmp_big), 1,000 states against kref.py.
- rvcmp_big's certificate mode on KBIG files at (3,5) and (3,6).
- kbcertify's final rounding: exact decimal comparison of every printed bound I produced.
- Process state after the build: found the pgrep deadlock (finding 1b).
- Both mid-review edits to kbig.c: source diff, binary against my compile, and (for the 13:43 one)
  restarts from a checkpoint with z inside, at and outside its range and with repeated z values.

## Questions (suspicions I could not demonstrate)

- kbverify reads the table through stdio without F_NOCACHE, so a certificate run pushes 76 GB
  (169 GB at (4,6)) through the page cache of a 24 GB machine that is already swapping; slowness
  for the other jobs, not wrong numbers (the implementer measured the cache effect for the other
  readers and fixed those; kbverify itself is unchanged and its numbers are reproduced by my verifier).
- (Resolved) The (4,6) table came out at 169 GB and fitted; the eigen checkpoints and certificate
  files fitted beside it.

## Commands I ran (all with nice -n 10; from kary/ unless noted)

    clang -O3 -march=native -o review-big/bin/kbig kbig.c -lm          (and kbig_v2 after the edit)
    clang -O2 -ffp-model=strict -march=native -o review-big/bin/kbverify kbverify.c -lm   (kverify too)
    clang -O3 -ffp-model=strict -march=native -o review-big/bin/rvcmp_big rvcmp_big.c -lpthread
    otool -tV / otool -s __TEXT __cstring|__const on each pair; cmp of the listings
    review-big/repro_rv.sh "2 2" ... "4 4"; review-big/repro_rv.sh "3 6" "4 5" "2 8"
    python review-big/sections_cmp.py
    review-big/small_verify.sh "3 2" "3 3" "3 4" "3 5" "4 2" "4 3" "4 4"; ... "3 6" "4 5"
    review-big/mut/kbverify_cap{100,700,1000,4099} k h review-big/t/k<k>h<h>
    clang -O2 -ffp-model=strict -march=native -o review-big/bin/rvb_emu_test review-big/rvb_emu_test.c -lm; run
    clang -O2 -ffp-contract=off -march=native -o review-big/bin/rvb_verify review-big/rvb_verify.c -lm -lpthread
    review-big/bin/rvb_verify k h <kbig> <kvec> <weights> <wtypes> 4; python review-big/rvb_bound.py k <out> <certify.out>
    python review-big/redctl.py 3 4
    python review-big/wt/check_weights_new_h.py      (wtharness: kbig.c's make_inputs only)
    review-big/bin/kbig_v2 eigen 3 5 review-big/eig35 - 1e-12 400 0.1948 0.1950; python kbcertify.py 3 5 review-big/eig35
    review-big/mut/kbig_cap{100003,2500000,4100000} eigen 3 5 ...; cmp of vectors and checkpoint
    review-big/bin/kbig_v2 eigen 3 4 review-big/redctl/k3h4_drop_record - 1e-13 2000 0.1821; kbcertify; kbverify on the true table
    (review/) python rv_theirweights.py 3 5 ../review-big/eig35 ...; ./rvcmp_big 3 5 <kbig> <kvec> <refweights>
    FULL=1 THREADS=4 review-big/check_cert.sh 3 6
    python kcover.py 3 7 4000 11; python kcover.py 4 6 4000 11
    review-big/sampled_k3h7.sh   (RVK_THREADS=4 RVCMP_STRIDE=101 RVCMP_PHASE=50 ./rvcmp_big 3 7 ...;
                                  python check_table_big.py 3 7 big/k3h7 600 2027)
    pgrep -fl "kbig build 3 7"; ps -o args= for the waiting pids; man pgrep

## Comparison status (finding 1) at 17:30, read from the logs myself

| table | quarter | states compared | expected | transitions | RESULT | exit |
|---|---|---|---|---|---|---|
| (3,7) | q0 | 7,241,728 | 7,241,728 | 2,379,806,092 | no difference | 0 |
| (3,7) | q1 | 7,241,728 | 7,241,728 | 2,377,445,667 | no difference | 0 |
| (3,7) | q2 | 7,241,728 | 7,241,728 | 2,376,055,748 | no difference | 0 |
| (3,7) | q3 | 7,241,587 | 7,241,587 | 2,386,708,579 | no difference | 0 |
| (3,7) | sum | 28,966,771 | N = 28,966,771 | 9,520,016,086 | build: 9,520,016,086 | |
| (4,6) | q0 | 2,125,824 | 2,125,824 | 6,920,046,139 | no difference | 0 |
| (4,6) | q1 | running | 2,125,824 | | | |
| (4,6) | q2 | running | 2,125,729 | | | |
| (4,6) | q3 | running | 2,124,800 | | | |
| (4,6) | sum | | N = 8,502,177 | must be 27,720,529,578 | | |

Every (3,7) log and the (4,6) q0 log also report 0 stored states sharing an orbit, 0 not in
canonical form, start state index 0, and 0 targets missing. q1 to q3 of (4,6) run side by side since
16:58 (big/compare_parallel.sh 4 6 1 2 3, separate logs, each appending its exit status); the
team lead will read them.
