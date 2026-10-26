// mc_ed.c (review code): Monte Carlo of ed(X, Y)/n for independent uniform X, Y in [k]^n.
// Exact Levenshtein distance by the O(n^2) two-row DP; xoshiro256** RNG.
// Prints mean, sample sd, standard error, min and max of ed/n, and how many samples fall at or
// below each given bound (a first-moment lower bound a must have P(ed <= a n) -> 0, so for large
// n essentially none should).
// usage: mc_ed k n N seed [bound ...]
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

static uint64_t s[4];
static inline uint64_t rotl(const uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static uint64_t nxt(void) {
  const uint64_t r = rotl(s[1] * 5, 7) * 9, t = s[1] << 17;
  s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t; s[3] = rotl(s[3], 45);
  return r;
}
static int rnd(int k) {  // unbiased
  uint64_t lim = UINT64_MAX - UINT64_MAX % (uint64_t)k, v;
  do v = nxt(); while (v >= lim);
  return (int)(v % (uint64_t)k);
}

int main(int argc, char **argv) {
  int k = atoi(argv[1]), n = atoi(argv[2]), N = atoi(argv[3]);
  uint64_t seed = strtoull(argv[4], 0, 10);
  int nb = argc - 5; double *bnd = malloc(sizeof(double) * (nb > 0 ? nb : 1));
  for (int b = 0; b < nb; b++) bnd[b] = atof(argv[5 + b]);
  long *below = calloc(nb > 0 ? nb : 1, sizeof(long));
  s[0] = seed ^ 0x9E3779B97F4A7C15ULL; s[1] = seed * 0xBF58476D1CE4E5B9ULL + 7; s[2] = 0x94D049BB133111EBULL; s[3] = seed + 99;
  for (int z = 0; z < 32; z++) nxt();
  unsigned char *X = malloc(n), *Y = malloc(n);
  int *prev = malloc(sizeof(int) * (n + 1)), *cur = malloc(sizeof(int) * (n + 1));
  double sum = 0, sum2 = 0, mn = 1e9, mx = -1;
  for (int it = 0; it < N; it++) {
    for (int i = 0; i < n; i++) X[i] = (unsigned char)rnd(k);
    for (int j = 0; j < n; j++) Y[j] = (unsigned char)rnd(k);
    for (int j = 0; j <= n; j++) prev[j] = j;
    for (int i = 1; i <= n; i++) {
      cur[0] = i;
      unsigned char xi = X[i - 1];
      for (int j = 1; j <= n; j++) {
        int b = prev[j - 1] + (xi != Y[j - 1]);
        int d = prev[j] + 1; if (d < b) b = d;
        int in = cur[j - 1] + 1; if (in < b) b = in;
        cur[j] = b;
      }
      int *tmp = prev; prev = cur; cur = tmp;
    }
    double r = (double)prev[n] / n;
    sum += r; sum2 += r * r; if (r < mn) mn = r; if (r > mx) mx = r;
    for (int b = 0; b < nb; b++) if (prev[n] <= bnd[b] * n) below[b]++;
  }
  double mean = sum / N, sd = sqrt((sum2 - N * mean * mean) / (N - 1));
  printf("k=%d n=%d N=%d: mean ed/n = %.5f  sd = %.5f  se = %.5f  min = %.5f  max = %.5f",
         k, n, N, mean, sd, sd / sqrt(N), mn, mx);
  for (int b = 0; b < nb; b++) printf("  #(ed<=%.6f n)=%ld", bnd[b], below[b]);
  printf("\n");
  return 0;
}
