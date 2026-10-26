// klow_plain.c: Lueker's lower bound for the Chvatal-Sankoff constant gamma_k
// (J. ACM 56(3) 2009, Section 2, Lemma 2.3), any alphabet size k, plain layout.
//
// State: a pair (s, t) of strings of length h over {0..k-1}, stored at index
// s * k^h + t (s[0] is the most significant base-k digit). One step of
//   T(w)(as', at') = 1 + (1/k^2) sum_{c,c'} w(s'c, t'c')
//   T(w)(as', bt') = max( (1/k) sum_c w(s'c, bt'), (1/k) sum_c w(as', t'c) ),  a != b
// is iterated from w = 0. Values are fixed point (unit 2^-SHIFT) and every
// division rounds down, so the computed T is never above the true T. For the
// input vector w of any step, r = min_i (T(w)_i - w_i) satisfies T(w,w) >= w + r,
// hence gamma_k >= 2r/(1+r) (Lemma 2.3). No symmetry is used: this program is the
// simple reference for small (k, h).
//
// usage: klow_plain k h iters [threads] [dumpfile]
#include <inttypes.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SHIFT 27
typedef uint32_t val_t;

static int K, H, NT;
static uint64_t P, PH, N; // k^(h-1), k^h, k^(2h)
static val_t *w, *tw, *Fp;
static int FROUTE;

typedef struct {
  int id;
  uint64_t s_lo, s_hi;
  uint64_t tmin, tmax;  // min, max of T(w)
  int64_t dmin, dmax;   // min, max of T(w) - w
} job_t;

// F-route (same rounding as klow.c): F(x, y) = floor(sum_c w(xc, y) / k), |x| = h-1, |y| = h.
static void *fworker(void *arg) {
  job_t *j = (job_t *)arg;  // s_lo..s_hi range over x in [0, k^(h-1))
  for (uint64_t x = j->s_lo; x < j->s_hi; x++)
    for (uint64_t y = 0; y < PH; y++) {
      uint64_t sum = 0;
      for (int c = 0; c < K; c++) sum += w[(x * K + c) * PH + y];
      Fp[x * PH + y] = (val_t)(sum / K);
    }
  return NULL;
}

static inline uint64_t tval_froute(uint64_t s, uint64_t t) {
  const uint64_t one = (uint64_t)1 << SHIFT;
  uint64_t s0 = s / P, sr = s % P, t0 = t / P, tr = t % P;
  if (s0 == t0) {
    uint64_t sum1 = 0, sum2 = 0;  // the two rounding orders; take the larger
    for (int c = 0; c < K; c++) { sum1 += Fp[tr * PH + sr * K + c]; sum2 += Fp[sr * PH + tr * K + c]; }
    return one + (sum1 > sum2 ? sum1 : sum2) / K;
  }
  uint64_t A = Fp[sr * PH + t], B = Fp[tr * PH + s];
  return A > B ? A : B;
}

static void *worker(void *arg) {
  job_t *j = (job_t *)arg;
  const uint64_t one = (uint64_t)1 << SHIFT, kk = (uint64_t)K * K;
  uint64_t tmin = UINT64_MAX, tmax = 0;
  int64_t dmin = INT64_MAX, dmax = INT64_MIN;
  for (uint64_t s = j->s_lo; s < j->s_hi; s++) {
    uint64_t s0 = s / P, sr = s % P;
    for (uint64_t t = 0; t < PH; t++) {
      uint64_t t0 = t / P, tr = t % P, v;
      if (FROUTE) v = tval_froute(s, t);
      else if (s0 == t0) {
        uint64_t sum = 0;
        for (int c = 0; c < K; c++) {
          const val_t *row = w + (sr * K + c) * PH + tr * K;
          for (int c2 = 0; c2 < K; c2++) sum += row[c2];
        }
        v = one + sum / kk;
      } else {
        uint64_t a = 0, b = 0;
        for (int c = 0; c < K; c++) a += w[(sr * K + c) * PH + t];
        const val_t *row = w + s * PH + tr * K;
        for (int c = 0; c < K; c++) b += row[c];
        v = (a > b ? a : b) / K;
      }
      uint64_t i = s * PH + t;
      tw[i] = (val_t)v;
      if (v > 0xffffffffu) { fprintf(stderr, "overflow at %" PRIu64 "\n", i); exit(2); }
      int64_t d = (int64_t)v - (int64_t)w[i];
      if (v < tmin) tmin = v;
      if (v > tmax) tmax = v;
      if (d < dmin) dmin = d;
      if (d > dmax) dmax = d;
    }
  }
  j->tmin = tmin; j->tmax = tmax; j->dmin = dmin; j->dmax = dmax;
  return NULL;
}

static void bound_str(int64_t r, char *out) {
  // floor(10^9 * 2r/(2^SHIFT + r)) as 0.xxxxxxxxx
  if (r <= 0) { strcpy(out, "<=0"); return; }
  __int128 num = (__int128)2 * r * 1000000000;
  __int128 den = ((__int128)1 << SHIFT) + r;
  int64_t q = (int64_t)(num / den);
  sprintf(out, "0.%09" PRId64, q);
}

int main(int argc, char **argv) {
  if (argc < 4) { fprintf(stderr, "usage: klow_plain k h iters [threads] [dumpfile]\n"); return 1; }
  K = atoi(argv[1]); H = atoi(argv[2]);
  int iters = atoi(argv[3]);
  NT = argc > 4 ? atoi(argv[4]) : 8;
  const char *dump = argc > 5 ? argv[5] : NULL;
  P = 1; for (int i = 0; i < H - 1; i++) P *= K;
  PH = P * K; N = PH * PH;
  FROUTE = getenv("KLOW_FROUTE") && atoi(getenv("KLOW_FROUTE"));
  w = calloc(N, sizeof(val_t)); tw = malloc(N * sizeof(val_t));
  if (FROUTE) Fp = malloc(P * PH * sizeof(val_t));
  if (!w || !tw) { fprintf(stderr, "out of memory for %" PRIu64 " entries\n", N); return 1; }
  printf("k=%d h=%d states=%" PRIu64 " shift=%d threads=%d froute=%d\n", K, H, N, SHIFT, NT, FROUTE);
  pthread_t th[256]; job_t jb[256];
  int64_t best_r = INT64_MIN; int best_it = -1;
  for (int it = 1; it <= iters; it++) {
    if (FROUTE) {
      for (int i = 0; i < NT; i++) {
        jb[i].s_lo = P * i / NT; jb[i].s_hi = P * (i + 1) / NT;
        pthread_create(&th[i], NULL, fworker, &jb[i]);
      }
      for (int i = 0; i < NT; i++) pthread_join(th[i], NULL);
    }
    for (int i = 0; i < NT; i++) {
      jb[i].id = i; jb[i].s_lo = PH * i / NT; jb[i].s_hi = PH * (i + 1) / NT;
      pthread_create(&th[i], NULL, worker, &jb[i]);
    }
    uint64_t tmin = UINT64_MAX, tmax = 0; int64_t dmin = INT64_MAX, dmax = INT64_MIN;
    for (int i = 0; i < NT; i++) {
      pthread_join(th[i], NULL);
      if (jb[i].tmin < tmin) tmin = jb[i].tmin;
      if (jb[i].tmax > tmax) tmax = jb[i].tmax;
      if (jb[i].dmin < dmin) dmin = jb[i].dmin;
      if (jb[i].dmax > dmax) dmax = jb[i].dmax;
    }
    char bs[64]; bound_str(dmin, bs);
    if (dmin > best_r) { best_r = dmin; best_it = it; }
    if (it <= 3 || it % 5 == 0 || it == iters)
      printf("it %4d  rmin %.10f  rmax %.10f  spread %.6f  bound %s\n", it,
             (double)dmin / (1 << SHIFT), (double)dmax / (1 << SHIFT),
             (double)(tmax - tmin) / (1 << SHIFT), bs);
    fflush(stdout);
    if (it == iters && dump) { // dump the input vector of the last step (the witness)
      FILE *f = fopen(dump, "wb");
      if (!f || fwrite(w, sizeof(val_t), N, f) != N || fclose(f)) { fprintf(stderr, "dump failed\n"); return 4; }
    }
    for (uint64_t i = 0; i < N; i++) w[i] = tw[i] - (val_t)tmin;
  }
  char bs[64]; bound_str(best_r, bs);
  printf("best r = %" PRId64 " / 2^%d at step %d: gamma_%d >= %s\n", best_r, SHIFT, best_it, K, bs);
  return 0;
}
