// kedit_plain.c: the iteration of kedit.c on the plain layout (every pair of strings
// stored, no symmetry), with the same arithmetic, for checking kedit.c on small (k, h).
// usage: kedit_plain k h iters
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SHIFT 27
typedef uint32_t val_t;

static inline int64_t ceil2(int64_t x) { return x >= 0 ? (x + 1) / 2 : -((-x) / 2); }

static void bound_str(int64_t r, char *out) {
  __int128 num = (__int128)2 * r * 1000000000, den = (__int128)1 << SHIFT;
  int64_t q = (int64_t)((num + den - 1) / den);
  sprintf(out, "%" PRId64 ".%09" PRId64, q / 1000000000, q % 1000000000);
}

int main(int argc, char **argv) {
  if (argc < 4) { fprintf(stderr, "usage: kedit_plain k h iters\n"); return 1; }
  int K = atoi(argv[1]), H = atoi(argv[2]), iters = atoi(argv[3]);
  uint64_t P = 1; for (int i = 0; i < H - 1; i++) P *= K;
  uint64_t PH = P * K, N = PH * PH;
  val_t *w = calloc(N, sizeof(val_t)), *tw = malloc(N * sizeof(val_t)), *F = malloc(P * PH * sizeof(val_t));
  const int64_t ONE = (int64_t)1 << SHIFT;
  int64_t RHO = 0, best = INT64_MAX, need = 0, dmin = 0, dmax = 0;
  int best_it = 0, cert = 0;
  for (int it = 1; ; it++) {
    for (uint64_t x = 0; x < P; x++)  // F(x, y) = ceil(sum_c w(xc, y) / k)
      for (uint64_t y = 0; y < PH; y++) {
        uint64_t sum = 0;
        for (int c = 0; c < K; c++) sum += w[(x * K + c) * PH + y];
        F[x * PH + y] = (val_t)((sum + K - 1) / K);
      }
    int64_t vmin = INT64_MAX;
    dmin = INT64_MAX; dmax = INT64_MIN; need = INT64_MIN;
    for (uint64_t s = 0; s < PH; s++)
      for (uint64_t t = 0; t < PH; t++) {
        uint64_t s0 = s / P, sr = s % P, t0 = t / P, tr = t % P;
        uint64_t sum1 = 0, sum2 = 0;
        for (int c = 0; c < K; c++) { sum1 += F[tr * PH + sr * K + c]; sum2 += F[sr * PH + tr * K + c]; }
        int64_t g1 = (int64_t)((sum1 + K - 1) / K), g2 = (int64_t)((sum2 + K - 1) / K);
        int64_t G = g1 < g2 ? g1 : g2, wv = w[s * PH + t], tv, nd;
        if (s0 == t0) {
          tv = G - RHO;
          nd = ceil2(G - wv);
        } else {
          int64_t A = F[sr * PH + t], B = F[tr * PH + s], mab = A < B ? A : B;
          tv = ONE + (mab < G - RHO ? mab : G - RHO);
          int64_t n1 = ONE + mab - wv, n2 = ceil2(ONE + G - wv);
          nd = n1 < n2 ? n1 : n2;
        }
        int64_t st = tv + RHO, dd = tv - wv;
        if (dd < dmin) dmin = dd;
        if (dd > dmax) dmax = dd;
        if (nd > need) need = nd;
        if (st < vmin) vmin = st;
        tw[s * PH + t] = (val_t)st;
      }
    for (uint64_t i = 0; i < N; i++) w[i] = tw[i] - (val_t)vmin;
    if (need < best) { best = need; best_it = it; }
    char bs[64]; bound_str(need, bs);
    if (it <= 3 || it % 5 == 0 || cert)
      printf("it %4d  growth %.10f..%.10f  rho %.10f  cert r %.10f  bound %s\n", it,
             (double)dmin / (1 << SHIFT), (double)dmax / (1 << SHIFT), (double)RHO / (1 << SHIFT),
             (double)need / (1 << SHIFT), bs);
    if (cert) break;
    RHO = (dmin + dmax) / 2;
    if (RHO < 0) RHO = 0;
    if (it >= iters || (it >= 20 && dmax - dmin <= 2 && need - dmax <= 4)) cert = 1;
  }
  char bs[64]; bound_str(best, bs);
  printf("best certificate over all steps = %" PRId64 " / 2^%d at step %d (%s)\n", best, SHIFT, best_it, bs);
  return 0;
}
