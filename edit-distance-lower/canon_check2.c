// canon_check2.c: exhaustive check that every pair (X, Y) has an optimal alignment obeying:
// no I-D, D-I, I-S, D-S adjacency; S only on unequal letters, M only on equal letters; and after a
// maximal run of insertions (deletions), a match uses a letter that differs from every letter of
// that run. DP state: (i, j, last move type, set of letters in the current run).
// usage: canon_check2 k n m
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define INF 1000000
int main(int argc, char **argv) {
  int k = atoi(argv[1]), n = atoi(argv[2]), m = atoi(argv[3]);
  long long NX = 1, NY = 1; for (int i = 0; i < n; i++) NX *= k; for (int i = 0; i < m; i++) NY *= k;
  int X[32], Y[32]; long long bad = 0, total = 0;
  int K2 = 1 << k;  // run letter sets
  static int E[33][33];
  int *C = malloc(sizeof(int) * (n + 1) * (m + 1) * 3 * K2);
#define AT(i, j, s, mask) C[((((i) * (m + 1) + (j)) * 3 + (s)) * K2) + (mask)]
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
      for (int z = 0; z < (n + 1) * (m + 1) * 3 * K2; z++) C[z] = INF;
      AT(0, 0, 0, 0) = 0;
      for (int i = 0; i <= n; i++) for (int j = 0; j <= m; j++) for (int s = 0; s < 3; s++) for (int mask = 0; mask < K2; mask++) {
        int c = AT(i, j, s, mask); if (c >= INF) continue;
        if (i < n && j < m) {
          if (X[i] == Y[j]) {  // M: after a run, the letter must avoid the run's letters
            if (!(s != 0 && (mask >> X[i] & 1)) && c < AT(i+1, j+1, 0, 0)) AT(i+1, j+1, 0, 0) = c;
          } else if (s == 0) {  // S only after a diagonal move or at the start
            if (c + 1 < AT(i+1, j+1, 0, 0)) AT(i+1, j+1, 0, 0) = c + 1;
          }
        }
        if (j < m && s != 2) { int nm = (s == 1 ? mask : 0) | (1 << Y[j]); if (c + 1 < AT(i, j+1, 1, nm)) AT(i, j+1, 1, nm) = c + 1; }
        if (i < n && s != 1) { int nm = (s == 2 ? mask : 0) | (1 << X[i]); if (c + 1 < AT(i+1, j, 2, nm)) AT(i+1, j, 2, nm) = c + 1; }
      }
      int cc = INF; for (int s = 0; s < 3; s++) for (int mask = 0; mask < K2; mask++) if (AT(n, m, s, mask) < cc) cc = AT(n, m, s, mask);
      total++;
      if (cc != E[n][m]) { bad++; if (bad <= 3) { printf("mismatch ed=%d canon=%d X=", E[n][m], cc); for (int i = 0; i < n; i++) printf("%d", X[i]); printf(" Y="); for (int j = 0; j < m; j++) printf("%d", Y[j]); printf("\n"); } }
    }
  }
  printf("k=%d n=%d m=%d: %lld pairs, %lld mismatches\n", k, n, m, total, bad);
  return bad != 0;
}
