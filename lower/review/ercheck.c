// ercheck.c (review): independent EXACT checker for kedit.c witnesses (edit distance, upper
// bound on alpha_k), from Rosenfeld arXiv 2407.18113 Sec. 2 (Facts 5, 6, Lemma 7) and the
// witness format only; same class enumeration as rcheck.c (no code from kedit.c).
//   matched:   least r = (G - v)/2,  G = (1/k^2) sum_{c,c'} v(s'c, t'c')
//   mismatched: least r = min( 1 + min(A, B) - v, (1 + G - v)/2 ),
//              A = (1/k) sum_c v(s'c, bt'), B = (1/k) sum_c v(as', t'c)
//   certificate r must be >= max over classes; then alpha_k <= 2r.
// Numerators below are over 2k^2 (units 2^-27).
//
// Differences from klow.c on purpose:
//  * classes are found by DFS enumeration of all RGS of length 2h in lexicographic order;
//    the rank function is self-tested against the DFS counter for EVERY RGS;
//  * the canonical orientation is chosen by comparing the two relabelled strings
//    lexicographically (memcmp), not by comparing ranks;
//  * T is evaluated exactly (no floor anywhere): numerators over k^2, every letter c of the
//    alphabet enumerated explicitly (no "fresh letter times (k-u)" shortcut), and the
//    mismatched second term uses w(as', t'c) directly (no swap-symmetric F lookup);
//  * a random sample of raw pairs (s,t) (not representatives) is evaluated too, and must
//    give exactly the same T - w as its class representative.
//
// usage: ercheck k h witness.bin [threads] [claimed_r] [samples]
// exit 0 iff the exact r >= claimed_r (when given).
#include <fcntl.h>
#include <inttypes.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#define U27 ((int64_t)1 << 27)
static int K, H, L, NT;
static uint64_t D[48][34];          // D[j][m]: completions of positions j..L-1 with m labels used
static uint64_t NR;                 // number of RGS of length L with at most K labels
static uint64_t *bm, *pc;           // canonical bitmap over ranks, prefix popcounts per word
static const uint32_t *W;
static uint64_t NWIT;

static void die(const char *m) { fprintf(stderr, "ercheck: %s\n", m); exit(9); }

static void build_D(void) {
  for (int m = 1; m <= K; m++) D[L][m] = 1;
  for (int j = L - 1; j >= 1; j--)
    for (int m = 1; m <= K; m++) {
      __uint128_t v = (__uint128_t)m * D[j + 1][m] + (m < K ? D[j + 1][m + 1] : 0);
      if (v >> 64) die("count overflow");
      D[j][m] = (uint64_t)v;
    }
  NR = D[1][1];
}

static inline void relabel(const uint8_t *in, int n, uint8_t *out) {
  uint8_t map[64];
  memset(map, 0xff, sizeof map);
  int nx = 0;
  for (int i = 0; i < n; i++) {
    uint8_t x = in[i];
    if (map[x] == 0xff) map[x] = (uint8_t)nx++;
    out[i] = map[x];
  }
}

static inline uint64_t rank_rgs(const uint8_t *z) {
  uint64_t R = 0;
  int m = 1;
  for (int i = 1; i < L; i++) {
    R += (uint64_t)z[i] * D[i + 1][m];
    if (z[i] == m) m++;
  }
  return R;
}

// canonical relabelled string of the pair (x, y), |x| = |y| = H
static inline void canon(const uint8_t *x, const uint8_t *y, uint8_t *z) {
  uint8_t a[48], b[48], z2[48];
  memcpy(a, x, H); memcpy(a + H, y, H);
  memcpy(b, y, H); memcpy(b + H, x, H);
  relabel(a, L, z);
  relabel(b, L, z2);
  if (memcmp(z2, z, L) < 0) memcpy(z, z2, L);
}

static inline uint64_t dense(uint64_t R) {
  uint64_t wd = R >> 6, b = R & 63;
  return pc[wd] + (uint64_t)__builtin_popcountll(bm[wd] & ((b ? (~0ULL >> (64 - b)) : 0)));
}

static inline int64_t wl(const uint8_t *x, const uint8_t *y) {
  uint8_t z[48];
  canon(x, y, z);
  uint64_t R = rank_rgs(z);
  if (!((bm[R >> 6] >> (R & 63)) & 1)) die("lookup hit a non-canonical rank");
  uint64_t d = dense(R);
  if (d >= NWIT) die("dense index out of range");
  return (int64_t)W[d];
}

// exact 2k^2 * (least r at the pair (s,t)) in units of 2^-27
static int64_t num_at(const uint8_t *s, const uint8_t *t) {
  uint8_t a[48], b[48];
  int64_t v = wl(s, t), KK = (int64_t)K * K, S = 0;
  memcpy(a, s + 1, H - 1);
  memcpy(b, t + 1, H - 1);
  for (int c = 0; c < K; c++) {
    a[H - 1] = (uint8_t)c;
    for (int c2 = 0; c2 < K; c2++) { b[H - 1] = (uint8_t)c2; S += wl(a, b); }
  }
  if (s[0] == t[0]) return S - KK * v;                       // 2k^2 (G - v)/2
  int64_t SA = 0, SB = 0;
  for (int c = 0; c < K; c++) {
    a[H - 1] = (uint8_t)c; SA += wl(a, t);                   // v(s'c, bt')
    b[H - 1] = (uint8_t)c; SB += wl(s, b);                   // v(as', t'c), direct
  }
  int64_t mab = SA < SB ? SA : SB;
  int64_t n1 = 2 * KK * U27 + 2 * (int64_t)K * mab - 2 * KK * v;  // 2k^2 (1 + min(A,B) - v)
  int64_t n2 = KK * U27 + S - KK * v;                             // 2k^2 (1 + G - v)/2
  return n1 < n2 ? n1 : n2;
}

// ---- pass 1: DFS over all RGS in lex order: self-test the rank, mark canonical ----
static uint64_t dfs_counter, rank_fail;
static uint8_t zz[48];
static void dfs1(int i, int m) {  // m = labels used in zz[0..i-1]
  if (i == L) {
    uint64_t R = dfs_counter++;
    if (rank_rgs(zz) != R) rank_fail++;
    uint8_t sw[48], zs[48];
    memcpy(sw, zz + H, H); memcpy(sw + H, zz, H);
    relabel(sw, L, zs);
    if (memcmp(zz, zs, L) <= 0) bm[R >> 6] |= 1ULL << (R & 63);
    return;
  }
  int top = m < K ? m : m - 1;  // labels 0..m-1 or the new label m (if m < K)
  for (int v = 0; v <= top; v++) { zz[i] = (uint8_t)v; dfs1(i + 1, v == m ? m + 1 : m); }
}

// ---- pass 2: every canonical RGS (= every class), exact T - w ----
#define PFX 5
typedef struct { int id; uint64_t ncls; int64_t mn, mx; uint64_t mn_d; uint8_t mn_z[48]; uint64_t nmatched; } job_t;
static uint8_t prefixes[100000][PFX];
static int nprefix, plen;
static void gen_prefix(uint8_t *p, int i, int m) {
  if (i == plen) { memcpy(prefixes[nprefix++], p, plen); return; }
  int top = m < K ? m : m - 1;
  for (int v = 0; v <= top; v++) { p[i] = (uint8_t)v; gen_prefix(p, i + 1, v == m ? m + 1 : m); }
}
static void dfs2(job_t *j, uint8_t *z, int i, int m) {
  if (i == L) {
    uint8_t sw[48], zs[48];
    memcpy(sw, z + H, H); memcpy(sw + H, z, H);
    relabel(sw, L, zs);
    if (memcmp(z, zs, L) > 0) return;  // not canonical
    j->ncls++;
    if (z[H] == z[0]) j->nmatched++;
    int64_t v = num_at(z, z + H);
    if (v < j->mn) { j->mn = v; j->mn_d = dense(rank_rgs(z)); memcpy(j->mn_z, z, L); }
    if (v > j->mx) j->mx = v;
    return;
  }
  int top = m < K ? m : m - 1;
  for (int v = 0; v <= top; v++) { z[i] = (uint8_t)v; dfs2(j, z, i + 1, v == m ? m + 1 : m); }
}
static void *work2(void *arg) {
  job_t *j = (job_t *)arg;
  uint8_t z[48];
  for (int p = j->id; p < nprefix; p += NT) {
    memcpy(z, prefixes[p], plen);
    int m = 0;
    for (int i = 0; i < plen; i++) if (z[i] == m) m++;
    dfs2(j, z, plen, m);
  }
  return NULL;
}

static uint64_t rng = 0x9E3779B97F4A7C15ULL;
static inline uint64_t next64(void) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng; }

int main(int argc, char **argv) {
  if (argc < 4) { fprintf(stderr, "usage: ercheck k h witness.bin [threads] [claimed_r] [samples]\n"); return 1; }
  K = atoi(argv[1]); H = atoi(argv[2]); L = 2 * H;
  NT = argc > 4 ? atoi(argv[4]) : 4;
  int have_claim = argc > 5;
  int64_t claimed = have_claim ? strtoll(argv[5], 0, 10) : 0;
  long samples = argc > 6 ? atol(argv[6]) : 200000;
  if (K < 2 || K > 32 || H < 2 || L > 46) die("bad k or h");
  build_D();
  int fd = open(argv[3], O_RDONLY);
  if (fd < 0) die("cannot open witness");
  struct stat st; fstat(fd, &st);
  if (st.st_size % 4) die("witness size not a multiple of 4");
  NWIT = (uint64_t)st.st_size / 4;
  W = mmap(0, st.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
  if (W == MAP_FAILED) die("mmap");
  uint64_t nw = (NR + 63) / 64;
  bm = calloc(nw + 1, 8); pc = calloc(nw + 1, 8);
  if (!bm || !pc) die("oom");
  dfs1(1, 1);  // zz[0] = 0
  if (dfs_counter != NR) die("DFS count differs from the completion-count total");
  uint64_t ncan = 0;
  for (uint64_t i = 0; i < nw; i++) { pc[i] = ncan; ncan += (uint64_t)__builtin_popcountll(bm[i]); }
  pc[nw] = ncan;
  printf("k=%d h=%d rgs=%" PRIu64 " (DFS %" PRIu64 ", rank self-test failures %" PRIu64 ") classes=%" PRIu64 " witness values=%" PRIu64 "\n",
         K, H, NR, dfs_counter, rank_fail, ncan, NWIT);
  if (rank_fail) die("rank self-test failed");
  if (ncan != NWIT) die("class count differs from the witness length");
  plen = L < PFX ? L : PFX;
  uint8_t p[PFX]; p[0] = 0;
  gen_prefix(p, 1, 1);
  job_t *jobs = calloc(NT, sizeof(job_t));
  pthread_t th[64];
  for (int i = 0; i < NT; i++) { jobs[i].id = i; jobs[i].mn = INT64_MAX; jobs[i].mx = INT64_MIN; pthread_create(&th[i], 0, work2, &jobs[i]); }
  uint64_t tot = 0, matched = 0; int64_t mn = INT64_MAX, mx = INT64_MIN; uint64_t mn_d = 0; uint8_t mn_z[48] = {0};
  for (int i = 0; i < NT; i++) {
    pthread_join(th[i], 0);
    tot += jobs[i].ncls; matched += jobs[i].nmatched;
    if (jobs[i].mn < mn || (jobs[i].mn == mn && jobs[i].mn_d < mn_d)) { mn = jobs[i].mn; mn_d = jobs[i].mn_d; memcpy(mn_z, jobs[i].mn_z, L); }
    if (jobs[i].mx > mx) mx = jobs[i].mx;
  }
  if (tot != NWIT) die("pass 2 visited a different number of classes");
  int64_t KK = (int64_t)K * K;
  printf("classes visited %" PRIu64 " (matched first letters %" PRIu64 ")\n", tot, matched);
  printf("exact least certificate r = max over classes = %" PRId64 "/%" PRId64 " units = %.6f units = %.12f (min over classes %.12f)\n", mx, 2 * KK,
         (double)mx / (2 * KK), (double)mx / (2 * KK) / U27, (double)mn / (2 * KK) / U27);
  printf("(class of the minimum, not the maximum) %" PRIu64 ": s=", mn_d);
  for (int i = 0; i < H; i++) printf("%d", mn_z[i]);
  printf(" t=");
  for (int i = 0; i < H; i++) printf("%d", mn_z[H + i]);
  printf("\n");
  // random raw pairs: T - w must equal the representative's exactly
  long bad = 0;
  for (long n = 0; n < samples; n++) {
    uint8_t s[24], t[24], z[48];
    for (int i = 0; i < H; i++) { s[i] = (uint8_t)(next64() % K); t[i] = (uint8_t)(next64() % K); }
    if (next64() & 1) t[0] = s[0];  // half of them with matched first letters
    canon(s, t, z);
    int64_t a = num_at(s, t), b = num_at(z, z + H);
    if (a != b) { if (bad < 5) fprintf(stderr, "sample mismatch %" PRId64 " vs %" PRId64 "\n", a, b); bad++; }
    if (a > mx) { fprintf(stderr, "sample above class maximum\n"); bad++; }
  }
  printf("random raw pairs: %ld, mismatches with their class representative: %ld\n", samples, bad);
  {
    // alpha_k <= 2r, r = mx/(2KK*2^27): ceil(1e9 * mx / (KK*2^27))
    __int128 den = (__int128)KK * U27, q = ((__int128)mx * 1000000000 + den - 1) / den;
    printf("bound from exact least r: alpha_%d <= %" PRId64 ".%09" PRId64 "\n", K, (int64_t)(q / 1000000000), (int64_t)(q % 1000000000));
  }
  int ok = 1;
  if (have_claim) {
    ok = (__int128)claimed * 2 * KK >= (__int128)mx && bad == 0;
    printf("%s: claimed r = %" PRId64 ", claimed - exact least r = %.4f units\n", ok ? "PASS" : "FAIL", claimed,
           (double)claimed - (double)mx / (2 * KK));
  }
  return ok ? 0 : 1;
}
