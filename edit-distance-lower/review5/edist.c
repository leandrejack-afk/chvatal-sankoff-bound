// edist.c (review of rule 5): exact distribution of ed(X, Y) over all pairs in [k]^n x [k]^m.
// Prints "c count" lines. usage: edist k n m
#include <stdio.h>
#include <stdlib.h>
int main(int argc, char **argv) {
  int k = atoi(argv[1]), n = atoi(argv[2]), m = atoi(argv[3]);
  long long NX = 1, NY = 1; for (int i = 0; i < n; i++) NX *= k; for (int j = 0; j < m; j++) NY *= k;
  int X[32], Y[32], prev[33], cur[33]; long long cnt[70] = {0};
  for (long long xi = 0; xi < NX; xi++) {
    long long t = xi; for (int i = n - 1; i >= 0; i--) { X[i] = t % k; t /= k; }
    for (long long yi = 0; yi < NY; yi++) {
      long long u = yi; for (int j = m - 1; j >= 0; j--) { Y[j] = u % k; u /= k; }
      for (int j = 0; j <= m; j++) prev[j] = j;
      for (int i = 1; i <= n; i++) {
        cur[0] = i;
        for (int j = 1; j <= m; j++) {
          int c = prev[j-1] + (X[i-1] != Y[j-1]);
          if (prev[j] + 1 < c) c = prev[j] + 1;
          if (cur[j-1] + 1 < c) c = cur[j-1] + 1;
          cur[j] = c;
        }
        for (int j = 0; j <= m; j++) prev[j] = cur[j];
      }
      cnt[prev[m]]++;
    }
  }
  for (int c = 0; c <= n + m; c++) if (cnt[c]) printf("%d %lld\n", c, cnt[c]);
  return 0;
}
