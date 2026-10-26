// mc5.c (review of rule 5): Monte Carlo of ed(X, Y)/n for independent uniform X, Y in [k]^n.
// Full O(n^2) Levenshtein DP (two rows), splitmix64 letters. Prints mean, sd, min, max of ed/n.
// usage: mc5 k n samples seed
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
static unsigned long long rs;
static unsigned long long rnd64(void) { rs += 0x9E3779B97F4A7C15ULL; unsigned long long z = rs; z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL; z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL; return z ^ (z >> 31); }
int main(int argc, char **argv) {
  int k = atoi(argv[1]), n = atoi(argv[2]), S = atoi(argv[3]); rs = strtoull(argv[4], 0, 10);
  unsigned char *X = malloc(n), *Y = malloc(n);
  int *prev = malloc(sizeof(int) * (n + 1)), *cur = malloc(sizeof(int) * (n + 1));
  double sum = 0, sum2 = 0, mn = 1e9, mx = -1;
  for (int s = 0; s < S; s++) {
    for (int i = 0; i < n; i++) { X[i] = rnd64() % k; Y[i] = rnd64() % k; }
    for (int j = 0; j <= n; j++) prev[j] = j;
    for (int i = 1; i <= n; i++) {
      cur[0] = i;
      unsigned char xi = X[i - 1];
      for (int j = 1; j <= n; j++) {
        int c = prev[j - 1] + (xi != Y[j - 1]);
        int d = prev[j] + 1; if (d < c) c = d;
        d = cur[j - 1] + 1; if (d < c) c = d;
        cur[j] = c;
      }
      int *tmp = prev; prev = cur; cur = tmp;
    }
    double v = (double)prev[n] / n;
    sum += v; sum2 += v * v; if (v < mn) mn = v; if (v > mx) mx = v;
  }
  double mean = sum / S, sd = sqrt(sum2 / S - mean * mean);
  printf("k=%d n=%d samples=%d: mean ed/n = %.5f (sd %.5f, se %.5f), min %.5f, max %.5f\n", k, n, S, mean, sd, sd / sqrt(S), mn, mx);
  return 0;
}
