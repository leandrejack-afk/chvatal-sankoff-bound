// kedit.c: upper bounds on the average edit distance constant alpha_k (Rosenfeld,
// arXiv 2407.18113, Section 2: Facts 5 and 6, Lemma 7), with the same class storage as
// klow.c (alphabet permutations and the swap of the two strings).
//
// Edit distance (unit-cost substitutions, insertions, deletions) is unchanged by
// permuting the alphabet and by swapping the strings, so V_n(s, t) is constant on the
// classes of klow.c. With F(x, y) = ceil( sum_c w(xc, y) / k ) for |x| = h - 1, |y| = h,
// and for a state (as, bt) with tails s, t of length h - 1:
//   G = min( ceil( sum_c F(t, sc) / k ), ceil( sum_c F(s, tc) / k ) )
//       >= (1/k^2) sum_{c,c'} w(sc, tc')
//   A = F(s, bt) >= (1/k) sum_c w(sc, bt),   B = F(t, as) >= (1/k) sum_c w(as, tc)
// Rosenfeld's operator with second argument v - r is
//   a == b:  T = G - r          a != b:  T = 1 + min(A, B, G - r)
// and Lemma 7 says: if T(v, v - r) <= v + r everywhere then alpha_k <= 2r. For a fixed
// v the least such r per class is
//   a == b:  ceil((G - v)/2)    a != b:  min( 1 + min(A, B) - v, ceil((1 + G - v)/2) )
// and the certificate is r = the maximum over classes (all roundings go up). The
// iteration that finds v is v <- T(v, v - rho), rho the current estimate of the growth
// per step, shifted to minimum 0 after every step.
//
// usage: kedit k h iters [threads] [dump_prefix]
#include <inttypes.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define SHIFT 27
#define MAXL 40
#define MAXK 32
typedef uint32_t val_t;

static int K, H, LW, LF, NT;
static uint64_t cntW[MAXL + 1][MAXK], cntF[MAXL + 1][MAXK];
static uint64_t NWF, NF, NW;  // all RGS of length 2h, all RGS of length 2h-1, kept classes
static uint64_t *bits, *cum;  // canonical bitmap over w-space ranks, prefix counts per word
static val_t *w, *F;
static int64_t RHO;  // growth estimate used by the iteration
// warm start: the witness of a run with h - 1, read through its own class space
static int HP;
static uint64_t cntP[MAXL + 1][MAXK];
static uint64_t *bitsP, *cumP;
static val_t *wprev;

static double now(void) {
  struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

static void build_counts(uint64_t cnt[][MAXK], int L) {
  for (int m = 0; m < K; m++) cnt[L][m] = 1;
  for (int i = L - 1; i >= 1; i--)
    for (int m = 0; m < K; m++) {
      uint64_t a, b = m + 1 < K ? cnt[i + 1][m + 1] : 0;
      if (__builtin_mul_overflow((uint64_t)(m + 1), cnt[i + 1][m], &a) || __builtin_add_overflow(a, b, &cnt[i][m])) {
        fprintf(stderr, "too many classes for k=%d h=%d (counts overflow 64 bits)\n", K, H);
        exit(1);
      }
    }
}

// Rank of the RGS of p[0..n1-1] . q[0..n2-1] (raw letters < MAXK + 2).
static inline uint64_t rank2(const uint8_t *p, int n1, const uint8_t *q, int n2,
                             uint64_t cnt[][MAXK]) {
  uint8_t map[MAXK + 2];
  memset(map, 0xff, sizeof map);
  int next = 0, m = -1, i = 0;
  uint64_t R = 0;
  for (int j = 0; j < n1; j++, i++) {
    int x = p[j], r = map[x];
    if (r == 0xff) r = map[x] = (uint8_t)next++;
    if (i) R += (uint64_t)r * cnt[i + 1][m];
    if (r > m) m = r;
  }
  for (int j = 0; j < n2; j++, i++) {
    int x = q[j], r = map[x];
    if (r == 0xff) r = map[x] = (uint8_t)next++;
    if (i) R += (uint64_t)r * cnt[i + 1][m];
    if (r > m) m = r;
  }
  return R;
}

static inline uint64_t dense_of_rank(uint64_t R) {
  uint64_t wd = R >> 6, b = R & 63;
  uint64_t mask = b ? (~0ULL >> (64 - b)) : 0;
  return cum[wd] + (uint64_t)__builtin_popcountll(bits[wd] & mask);
}

// Dense index of the class of the pair (a, b), |a| = |b| = h.
static inline uint64_t wdense(const uint8_t *a, const uint8_t *b) {
  uint64_t r1 = rank2(a, H, b, H, cntW), r2 = rank2(b, H, a, H, cntW);
  return dense_of_rank(r1 < r2 ? r1 : r2);
}

// Index of the class of the pair (x, y), |x| = h - 1, |y| = h.
static inline uint64_t fidx(const uint8_t *x, const uint8_t *y) {
  return rank2(x, H - 1, y, H, cntF);
}

static void unrank(uint64_t R, int L, uint64_t cnt[][MAXK], uint8_t *r, int8_t *mx) {
  r[0] = 0; mx[0] = 0;
  int m = 0;
  for (int i = 1; i < L; i++) {
    int top = m + 1 < K ? m + 1 : K - 1, v;
    for (v = 0; v <= top; v++) {
      uint64_t c = cnt[i + 1][v > m ? v : m];
      if (R < c) break;
      R -= c;
    }
    r[i] = (uint8_t)v;
    if (v > m) m = v;
    mx[i] = (int8_t)m;
  }
}

static inline int next_rgs(uint8_t *r, int8_t *mx, int L) {
  for (int i = L - 1; i >= 1; i--) {
    if (r[i] <= mx[i - 1] && r[i] + 1 < K) {
      r[i]++;
      int m = mx[i - 1] > r[i] ? mx[i - 1] : r[i];
      mx[i] = (int8_t)m;
      for (int j = i + 1; j < L; j++) { r[j] = 0; mx[j] = (int8_t)m; }
      return 1;
    }
  }
  return 0;
}

typedef struct {
  int id, phase;
  uint64_t lo, hi;
  uint64_t cnt;
  uint64_t vmin, vmax;
  int64_t dmin, dmax;
  uint64_t dmin_idx;
  int64_t need;  // max over the job's classes of the least r (certificate)
  uint64_t need_idx;
} job_t;

static inline int64_t ceil2(int64_t x) { return x >= 0 ? (x + 1) / 2 : -((-x) / 2); }

static void *work(void *arg) {
  job_t *j = (job_t *)arg;
  uint8_t r[MAXL], a[MAXL], b[MAXL];
  int8_t mx[MAXL];
  const uint64_t one = (uint64_t)1 << SHIFT;
  if (j->hi <= j->lo) { j->cnt = 0; j->vmin = UINT64_MAX; j->vmax = 0; j->dmin = INT64_MAX; j->dmax = INT64_MIN; j->need = INT64_MIN; return NULL; }
  if (j->phase == 0) {  // canonical bitmap over w-space RGS ranks [lo, hi), lo and hi multiples of 64
    unrank(j->lo, LW, cntW, r, mx);
    uint64_t word = 0, c = 0;
    for (uint64_t R = j->lo; R < j->hi; R++) {
      uint64_t R2 = rank2(r + H, H, r, H, cntW);
      if (R <= R2) { word |= 1ULL << (R & 63); c++; }
      if ((R & 63) == 63 || R + 1 == j->hi) { bits[R >> 6] = word; word = 0; }
      if (R + 1 < j->hi) next_rgs(r, mx, LW);
    }
    j->cnt = c;
  } else if (j->phase == 1) {  // F from w, F-space ranks [lo, hi)
    unrank(j->lo, LF, cntF, r, mx);
    for (uint64_t R = j->lo; R < j->hi; R++) {
      int u = mx[LF - 1] + 1;  // letters present in x.y
      memcpy(a, r, H - 1);
      const uint8_t *y = r + H - 1;
      uint64_t sum = 0;
      for (int c = 0; c < u; c++) { a[H - 1] = (uint8_t)c; sum += w[wdense(a, y)]; }
      if (u < K) { a[H - 1] = (uint8_t)u; sum += (uint64_t)(K - u) * w[wdense(a, y)]; }
      F[R] = (val_t)((sum + K - 1) / K);
      if (R + 1 < j->hi) next_rgs(r, mx, LF);
    }
  } else if (j->phase == 3) {  // warm start: w(s, t) = wprev(s[0..h-2], t[0..h-2])
    unrank(j->lo, LW, cntW, r, mx);
    uint64_t d = dense_of_rank(j->lo), cnt = 0;
    for (uint64_t R = j->lo; R < j->hi; R++) {
      if ((bits[R >> 6] >> (R & 63)) & 1) {
        const uint8_t *s = r, *t = r + H;
        uint64_t q1 = rank2(s, HP, t, HP, cntP), q2 = rank2(t, HP, s, HP, cntP);
        uint64_t q = q1 < q2 ? q1 : q2, wd = q >> 6, b = q & 63;
        uint64_t mask = b ? (~0ULL >> (64 - b)) : 0;
        w[d] = wprev[cumP[wd] + (uint64_t)__builtin_popcountll(bitsP[wd] & mask)];
        d++; cnt++;
      }
      if (R + 1 < j->hi) next_rgs(r, mx, LW);
    }
    j->cnt = cnt;
  } else {  // phase 2: v <- T(v, v - rho) through F, and the certificate of the input v
    unrank(j->lo, LW, cntW, r, mx);
    uint64_t d = dense_of_rank(j->lo), cnt = 0;
    uint64_t vmin = UINT64_MAX, vmax = 0, dmi = 0, ni = 0;
    int64_t dmin = INT64_MAX, dmax = INT64_MIN, need_max = INT64_MIN;
    const int64_t ONE = (int64_t)one;
    for (uint64_t R = j->lo; R < j->hi; R++) {
      if ((bits[R >> 6] >> (R & 63)) & 1) {
        const uint8_t *s = r, *t = r + H;
        uint32_t seen = 0;  // letters present in the tails s'.t'
        for (int i = 1; i < H; i++) seen |= 1u << s[i];
        for (int i = 1; i < H; i++) seen |= 1u << t[i];
        int u = __builtin_popcount(seen);
        uint64_t sum1 = 0, sum2 = 0;
        memcpy(a, s + 1, H - 1);  // a = s' c
        memcpy(b, t + 1, H - 1);  // b = t' c
        for (int c = 0; c < K + 1 && seen; c++) {
          if (seen & (1u << c)) {
            seen &= ~(1u << c);
            a[H - 1] = (uint8_t)c; sum1 += F[fidx(t + 1, a)];
            b[H - 1] = (uint8_t)c; sum2 += F[fidx(s + 1, b)];
          }
        }
        if (u < K) {
          a[H - 1] = (uint8_t)(MAXK); sum1 += (uint64_t)(K - u) * F[fidx(t + 1, a)];
          b[H - 1] = (uint8_t)(MAXK); sum2 += (uint64_t)(K - u) * F[fidx(s + 1, b)];
        }
        int64_t g1 = (int64_t)((sum1 + K - 1) / K), g2 = (int64_t)((sum2 + K - 1) / K);
        int64_t G = g1 < g2 ? g1 : g2, wv = w[d], tv, nd;
        if (t[0] == s[0]) {
          tv = G - RHO;
          nd = ceil2(G - wv);
        } else {
          int64_t A = F[fidx(s + 1, t)], B = F[fidx(t + 1, s)], mab = A < B ? A : B;
          tv = ONE + (mab < G - RHO ? mab : G - RHO);
          int64_t n1 = ONE + mab - wv, n2 = ceil2(ONE + G - wv);
          nd = n1 < n2 ? n1 : n2;
        }
        int64_t st = tv + RHO;  // >= 0; the shift by rho is undone by the normalization
        if (st < 0 || st > 0xffffffffLL) { fprintf(stderr, "range error %lld\n", (long long)st); exit(2); }
        int64_t dd = tv - wv;
        if (dd < dmin) { dmin = dd; dmi = d; }
        if (dd > dmax) dmax = dd;
        if (nd > need_max) { need_max = nd; ni = d; }
        if ((uint64_t)st < vmin) vmin = (uint64_t)st;
        if ((uint64_t)st > vmax) vmax = (uint64_t)st;
        w[d] = (val_t)st;
        d++; cnt++;
      }
      if (R + 1 < j->hi) next_rgs(r, mx, LW);
    }
    j->cnt = cnt; j->vmin = vmin; j->vmax = vmax; j->dmin = dmin; j->dmax = dmax; j->dmin_idx = dmi;
    j->need = need_max; j->need_idx = ni;
  }
  return NULL;
}

static job_t jobs[512];
static pthread_t th[512];
static int last_nj;

static void run_phase(int phase, uint64_t total, int align64) {
  int nj = NT * 8;  // more jobs than threads for balance, run in waves of NT
  if (nj > 512) nj = 512;
  if ((uint64_t)nj > total) nj = 1;
  last_nj = nj;
  for (int i = 0; i < nj; i++) {
    uint64_t lo = total * i / nj, hi = total * (i + 1) / nj;
    if (align64) { lo &= ~63ULL; hi = (i + 1 == nj) ? total : (hi & ~63ULL); }
    jobs[i] = (job_t){.id = i, .phase = phase, .lo = lo, .hi = hi};
  }
  for (int base = 0; base < nj; base += NT) {
    int n = nj - base < NT ? nj - base : NT;
    for (int i = 0; i < n; i++) pthread_create(&th[i], NULL, work, &jobs[base + i]);
    for (int i = 0; i < n; i++) pthread_join(th[i], NULL);
  }
}

static void bound_str(int64_t r, char *out) {  // ceil(10^9 * 2r / 2^SHIFT)
  if (r < 0) { strcpy(out, "<0 (impossible: check the operator)"); return; }
  __int128 num = (__int128)2 * r * 1000000000, den = (__int128)1 << SHIFT;
  int64_t q = (int64_t)((num + den - 1) / den);
  sprintf(out, "%" PRId64 ".%09" PRId64, q / 1000000000, q % 1000000000);
}

int main(int argc, char **argv) {
  if (argc < 4) { fprintf(stderr, "usage: kedit k h iters [threads] [dump_prefix]\n"); return 1; }
  K = atoi(argv[1]); H = atoi(argv[2]);
  int iters = atoi(argv[3]);
  NT = argc > 4 ? atoi(argv[4]) : 8;
  const char *dump = argc > 5 ? argv[5] : NULL;
  if (K < 2 || K > 30 || H < 2 || 2 * H > MAXL) { fprintf(stderr, "bad k or h\n"); return 1; }
  LW = 2 * H; LF = 2 * H - 1;
  build_counts(cntW, LW); build_counts(cntF, LF);
  NWF = cntW[1][0]; NF = cntF[1][0];
  uint64_t nwords = (NWF + 63) / 64;
  bits = calloc(nwords + 1, 8); cum = calloc(nwords + 1, 8);
  F = malloc(NF * sizeof(val_t));
  double t0 = now();
  run_phase(0, NWF, 1);
  NW = 0;
  for (uint64_t i = 0; i < nwords; i++) { cum[i] = NW; NW += (uint64_t)__builtin_popcountll(bits[i]); }
  cum[nwords] = NW;
  w = calloc(NW, sizeof(val_t));
  if (!bits || !cum || !F || !w) { fprintf(stderr, "out of memory\n"); return 1; }
  const char *init = getenv("KLOW_INIT");
  if (init) {  // the witness of the run with h - 1
    HP = H - 1;
    int LP = 2 * HP;
    build_counts(cntP, LP);
    uint64_t nfp = cntP[1][0], nwp_words = (nfp + 63) / 64, nwp = 0;
    bitsP = calloc(nwp_words + 1, 8); cumP = calloc(nwp_words + 1, 8);
    uint8_t r[MAXL]; int8_t mx[MAXL];
    unrank(0, LP, cntP, r, mx);
    for (uint64_t R = 0; R < nfp; R++) {
      if (R <= rank2(r + HP, HP, r, HP, cntP)) bitsP[R >> 6] |= 1ULL << (R & 63);
      if (R + 1 < nfp) next_rgs(r, mx, LP);
    }
    for (uint64_t i = 0; i < nwp_words; i++) { cumP[i] = nwp; nwp += (uint64_t)__builtin_popcountll(bitsP[i]); }
    cumP[nwp_words] = nwp;
    wprev = malloc(nwp * sizeof(val_t));
    FILE *f = fopen(init, "rb");
    if (!f || fread(wprev, sizeof(val_t), nwp, f) != nwp || fgetc(f) != EOF) { fprintf(stderr, "bad KLOW_INIT file (want %" PRIu64 " values)\n", nwp); return 5; }
    fclose(f);
    run_phase(3, NWF, 1);
    uint64_t cnt = 0;
    for (int i = 0; i < last_nj; i++) if (jobs[i].hi > jobs[i].lo) cnt += jobs[i].cnt;
    if (cnt != NW) { fprintf(stderr, "warm start count mismatch\n"); return 5; }
    free(wprev); free(bitsP); free(cumP);
    printf("warm start from %s (h=%d, %" PRIu64 " classes)\n", init, HP, nwp);
  }
  printf("k=%d h=%d rgs=%" PRIu64 " classes=%" PRIu64 " F=%" PRIu64 " shift=%d threads=%d (setup %.1fs)\n",
         K, H, NWF, NW, NF, SHIFT, NT, now() - t0);
  fflush(stdout);
  int64_t best = INT64_MAX; int best_it = 0;
  int64_t dmin = 0, dmax = 0, need = 0;
  int cert = 0;
  RHO = 0;
  for (int it = 1; ; it++) {
    double ta = now();
    if (cert && dump) {  // the witness: the input vector of this step
      char fn[1024]; snprintf(fn, sizeof fn, "%s_w.bin", dump);
      FILE *f = fopen(fn, "wb");
      if (!f || fwrite(w, sizeof(val_t), NW, f) != NW || fclose(f)) { fprintf(stderr, "dump failed\n"); return 4; }
    }
    run_phase(1, NF, 0);
    run_phase(2, NWF, 1);
    uint64_t vmin = UINT64_MAX, vmax = 0, cnt = 0, ni = 0;
    dmin = INT64_MAX; dmax = INT64_MIN; need = INT64_MIN;
    for (int i = 0; i < last_nj; i++) {
      if (jobs[i].hi <= jobs[i].lo) continue;
      cnt += jobs[i].cnt;
      if (jobs[i].vmin < vmin) vmin = jobs[i].vmin;
      if (jobs[i].vmax > vmax) vmax = jobs[i].vmax;
      if (jobs[i].dmin < dmin) dmin = jobs[i].dmin;
      if (jobs[i].dmax > dmax) dmax = jobs[i].dmax;
      if (jobs[i].need > need) { need = jobs[i].need; ni = jobs[i].need_idx; }
    }
    if (cnt != NW) { fprintf(stderr, "count mismatch %" PRIu64 " vs %" PRIu64 "\n", cnt, NW); return 3; }
    for (uint64_t i = 0; i < NW; i++) w[i] -= (val_t)vmin;
    if (need < best) { best = need; best_it = it; }
    char bs[64]; bound_str(need, bs);
    if (it <= 3 || it % 5 == 0 || cert)
      printf("it %4d  growth %.10f..%.10f  rho %.10f  cert r %.10f  bound %s  (argmax %" PRIu64 ", %.1fs)\n", it,
             (double)dmin / (1 << SHIFT), (double)dmax / (1 << SHIFT), (double)RHO / (1 << SHIFT),
             (double)need / (1 << SHIFT), bs, ni, now() - ta);
    fflush(stdout);
    if (cert) break;
    RHO = (dmin + dmax) / 2;
    if (RHO < 0) RHO = 0;
    if (it >= iters || (it >= 20 && dmax - dmin <= 2 && need - dmax <= 4)) cert = 1;
  }
  char bs[64]; bound_str(best, bs);
  printf("best certificate over all steps = %" PRId64 " / 2^%d at step %d (%s)\n", best, SHIFT, best_it, bs);
  bound_str(need, bs);
  printf("CERT k=%d h=%d classes=%" PRIu64 " shift=%d r=%" PRId64 " alpha_%d <= %s%s\n", K, H, NW, SHIFT,
         need, K, bs, dump ? " (witness dumped)" : "");
  return 0;
}
