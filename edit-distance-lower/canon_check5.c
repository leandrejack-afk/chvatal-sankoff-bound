// canon_check5.c: as canon_check3.c plus rule 5 (after a substitution on (a, b), the following
// insertion run avoids a and the following deletion run avoids b; matching that letter instead
// would be cheaper by 1, so this holds in every optimal alignment). Exhaustive check that every
// pair (X, Y) has an optimal alignment obeying rules 1 to 5.
// During an insertion run i does not change, and during a deletion run j does not change, so the
// substitution letter a = X[i-1] (b = Y[j-1]) is still at hand when the run grows.
// usage: canon_check5 k n m
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define INF 1000000
#define NS 8
int main(int argc, char **argv) {
  int k = atoi(argv[1]), n = atoi(argv[2]), m = atoi(argv[3]);
  long long NX = 1, NY = 1; for (int i = 0; i < n; i++) NX *= k; for (int i = 0; i < m; i++) NY *= k;
  int X[32], Y[32]; long long bad = 0, total = 0;
  int K2 = 1 << k;  // run letter sets
  // states: 0 diagonal after M (or start), 1 I-run after a diagonal-type state other than S,
  // 2 D-run likewise, 3 M right after a D-run, 4 M right after an I-run, 5 diagonal after S,
  // 6 I-run right after S, 7 D-run right after S
  static int E[33][33];
  int *C = malloc(sizeof(int) * (n + 1) * (m + 1) * NS * K2);
#define AT(i, j, s, mask) C[((((i) * (m + 1) + (j)) * NS + (s)) * K2) + (mask)]
  for (long long xi = 0; xi < NX; xi++) {
    long long t = xi; for (int i = n - 1; i >= 0; i--) { X[i] = t % k; t /= k; }
    for (long long yi = 0; yi < NY; yi++) {
      long long u = yi; for (int j = m - 1; j >= 0; j--) { Y[j] = u % k; u /= k; }
      for (int i = 0; i <= n; i++) for (int j = 0; j <= m; j++) {
        if (i == 0 && j == 0) { E[i][j] = 0; continue; }
        int best = INF;
        if (i > 0 && j > 0) { int c = E[i-1][j-1] + (X[i-1] != Y[j-1]); if (c < best) best = c; }
        if (j > 0 && E[i][j-1] + 1 < best) best = E[i][j-1] + 1;
        if (i > 0 && E[i-1][j] + 1 < best) best = E[i-1][j] + 1;
        E[i][j] = best;
      }
      for (int z = 0; z < (n + 1) * (m + 1) * NS * K2; z++) C[z] = INF;
      AT(0, 0, 0, 0) = 0;
      for (int i = 0; i <= n; i++) for (int j = 0; j <= m; j++) for (int s = 0; s < NS; s++) for (int mask = 0; mask < K2; mask++) {
        int c = AT(i, j, s, mask); if (c >= INF) continue;
        int irun = (s == 1 || s == 6), drun = (s == 2 || s == 7), inrun = irun || drun;
        if (i < n && j < m) {
          if (X[i] == Y[j]) {  // M: after a run, the letter must avoid the run's letters
            if (!(inrun && (mask >> X[i] & 1))) { int ns = drun ? 3 : irun ? 4 : 0; if (c < AT(i+1, j+1, ns, 0)) AT(i+1, j+1, ns, 0) = c; }
          } else if (!inrun) {  // S only after a diagonal-type move or at the start
            if (c + 1 < AT(i+1, j+1, 5, 0)) AT(i+1, j+1, 5, 0) = c + 1;
          }
        }
        if (j < m && !drun && s != 3) {  // insertion of Y[j]
          int afterS = (s == 5 || s == 6);
          if (!(afterS && Y[j] == X[i-1])) {
            int nm = (irun ? mask : 0) | (1 << Y[j]), ns = afterS ? 6 : 1;
            if (c + 1 < AT(i, j+1, ns, nm)) AT(i, j+1, ns, nm) = c + 1;
          }
        }
        if (i < n && !irun && s != 4) {  // deletion of X[i]
          int afterS = (s == 5 || s == 7);
          if (!(afterS && X[i] == Y[j-1])) {
            int nm = (drun ? mask : 0) | (1 << X[i]), ns = afterS ? 7 : 2;
            if (c + 1 < AT(i+1, j, ns, nm)) AT(i+1, j, ns, nm) = c + 1;
          }
        }
      }
      int cc = INF; for (int s = 0; s < NS; s++) for (int mask = 0; mask < K2; mask++) if (AT(n, m, s, mask) < cc) cc = AT(n, m, s, mask);
      total++;
      if (cc != E[n][m]) { bad++; if (bad <= 3) { printf("mismatch ed=%d canon=%d X=", E[n][m], cc); for (int i = 0; i < n; i++) printf("%d", X[i]); printf(" Y="); for (int j = 0; j < m; j++) printf("%d", Y[j]); printf("\n"); } }
    }
  }
  printf("k=%d n=%d m=%d: %lld pairs, %lld mismatches\n", k, n, m, total, bad);
  return bad != 0;
}
