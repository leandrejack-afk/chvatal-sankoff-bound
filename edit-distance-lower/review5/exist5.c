// exist5.c (review of rule 5, independent of canon_check5.c).
// Two checks per pair (X, Y):
//  (A) rule 5 holds in EVERY optimal alignment: with prefix distances E and suffix distances F, for
//      every substitution cell (i, j) (x_i != y_j, 1-indexed) and every s >= 1 with y_{j+s} == x_i,
//      E[i-1][j-1] + 1 + s + F[i][j+s] > ed (an alignment containing S I^s with the forbidden letter
//      is never optimal); likewise x_{i+s} == y_j for deletions: E[i-1][j-1] + 1 + s + F[i+s][j] > ed.
//      No rules 1-4 involved.
//  (B) min cost over alignments obeying rules 1-5 == ed, by a DP whose run states store the run's
//      START index and whether the run follows an S (no letter masks); rules 3 and 5 are checked by
//      scanning the run when it ends (at the closing match, or at the end of the strings).
// Mutants (argv[1] = 'm'): (A) uses the wrong letter (y_{j+s} == y_j, x_{i+s} == x_i), (B) adds the
// bogus rule "an insertion run after S also avoids b = y_j"; both must report violations.
// usage: exist5 [m] exh k n m            exhaustive over all pairs
//        exist5 [m] rnd k n m samples seed   random pairs
//        exist5 [m] mut k n M samples seed e  X random, Y = X after e random edits (M must be MAXL-safe; ignored)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define INF 1000000
#define MAXL 160
static int k, n, m, SD, mutant = 0;
static int X[MAXL], Y[MAXL];
static int E[MAXL + 1][MAXL + 1], F[MAXL + 2][MAXL + 2];
// (B) DP arrays: diag states 0 G0, 1 G1, 2 MD, 3 MI; runs: I-run after non-S (start s), I-run after S,
// D-run after non-S, D-run after S
static int Cd[MAXL + 1][MAXL + 1][4];
static int *Cr;  // [i][j][kind 0..3][start 0..MAXL]
#define CR(i, j, kd, s) Cr[((((i) * (m + 1) + (j)) * 4 + (kd)) * SD) + (s)]

static unsigned long long rs;
static unsigned long long rnd64(void) { rs += 0x9E3779B97F4A7C15ULL; unsigned long long z = rs; z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL; z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL; return z ^ (z >> 31); }

static void dists(void) {
  for (int i = 0; i <= n; i++) for (int j = 0; j <= m; j++) {
    if (i == 0) { E[i][j] = j; continue; }
    if (j == 0) { E[i][j] = i; continue; }
    int c = E[i-1][j-1] + (X[i-1] != Y[j-1]);
    if (E[i-1][j] + 1 < c) c = E[i-1][j] + 1;
    if (E[i][j-1] + 1 < c) c = E[i][j-1] + 1;
    E[i][j] = c;
  }
  // F[i][j] = ed(X[i..n-1], Y[j..m-1]) (0-indexed suffixes)
  for (int i = n; i >= 0; i--) for (int j = m; j >= 0; j--) {
    if (i == n) { F[i][j] = m - j; continue; }
    if (j == m) { F[i][j] = n - i; continue; }
    int c = F[i+1][j+1] + (X[i] != Y[j]);
    if (F[i+1][j] + 1 < c) c = F[i+1][j] + 1;
    if (F[i][j+1] + 1 < c) c = F[i][j+1] + 1;
    F[i][j] = c;
  }
}

// (A): returns number of violating (S, s) witnesses that lie on an optimal alignment
static long checkA(void) {
  int ed = E[n][m]; long viol = 0;
  for (int i = 1; i <= n; i++) for (int j = 1; j <= m; j++) {
    int a = X[i-1], b = Y[j-1];
    if (a == b) continue;
    int pre = E[i-1][j-1] + 1;
    for (int s = 1; j + s <= m; s++) {  // insertions of y_{j+1..j+s}; 1-indexed y_{j+s} = Y[j+s-1]
      int letter = Y[j + s - 1];
      int forb = mutant ? b : a;
      if (letter == forb && pre + s + F[i][j + s] == ed) viol++;
    }
    for (int s = 1; i + s <= n; s++) {
      int letter = X[i + s - 1];
      int forb = mutant ? a : b;
      if (letter == forb && pre + s + F[i + s][j] == ed) viol++;
    }
  }
  return viol;
}

static inline void relax(int *slot, int c) { if (c < *slot) *slot = c; }

// does the run (kind: 0 I plain, 1 I after S, 2 D plain, 3 D after S) from start s0 to the current
// (i, j) satisfy rule 5 (only after S)? Scans the run letters.
static int rule5_ok(int kd, int s0, int i, int j) {
  if (kd == 1) {  // I-run Y[s0..j-1] after S on (X[i-1], Y[s0-1])
    int a = X[i-1], b = Y[s0-1];
    for (int u = s0; u < j; u++) { if (Y[u] == a) return 0; if (mutant && Y[u] == b) return 0; }
  } else if (kd == 3) {  // D-run X[s0..i-1] after S on (X[s0-1], Y[j-1])
    int b = Y[j-1];
    for (int u = s0; u < i; u++) if (X[u] == b) return 0;
  }
  return 1;
}

static int checkB(void) {
  for (int i = 0; i <= n; i++) for (int j = 0; j <= m; j++) {
    for (int d = 0; d < 4; d++) Cd[i][j][d] = INF;
    for (int kd = 0; kd < 4; kd++) for (int s = 0; s < SD; s++) CR(i, j, kd, s) = INF;
  }
  Cd[0][0][0] = 0;
  for (int i = 0; i <= n; i++) for (int j = 0; j <= m; j++) {
    // diagonal-type states
    for (int d = 0; d < 4; d++) {
      int c = Cd[i][j][d]; if (c >= INF) continue;
      if (i < n && j < m) {
        if (X[i] == Y[j]) relax(&Cd[i+1][j+1][0], c);
        else relax(&Cd[i+1][j+1][1], c + 1);
      }
      if (j < m && d != 2) relax(&CR(i, j+1, d == 1 ? 1 : 0, j), c + 1);  // I-run starting at Y[j]
      if (i < n && d != 3) relax(&CR(i+1, j, d == 1 ? 3 : 2, i), c + 1);  // D-run starting at X[i]
    }
    // run states
    for (int kd = 0; kd < 4; kd++) {
      int ins = kd < 2;
      int smax = ins ? m : n;
      for (int s0 = 0; s0 <= smax; s0++) {
        int c = CR(i, j, kd, s0); if (c >= INF) continue;
        if (ins && j < m) relax(&CR(i, j+1, kd, s0), c + 1);
        if (!ins && i < n) relax(&CR(i+1, j, kd, s0), c + 1);
        if (i < n && j < m && X[i] == Y[j] && rule5_ok(kd, s0, i, j)) {  // closing match, rule 3 scan
          int ok = 1, letter = X[i];
          if (ins) { for (int u = s0; u < j; u++) if (Y[u] == letter) { ok = 0; break; } }
          else { for (int u = s0; u < i; u++) if (X[u] == letter) { ok = 0; break; } }
          if (ok) relax(&Cd[i+1][j+1][ins ? 3 : 2], c);
        }
      }
    }
  }
  int best = INF;
  for (int d = 0; d < 4; d++) if (Cd[n][m][d] < best) best = Cd[n][m][d];
  for (int kd = 0; kd < 4; kd++) for (int s0 = 0; s0 < SD; s0++) {
    int c = CR(n, m, kd, s0);
    if (c < best && rule5_ok(kd, s0, n, m)) best = c;  // run ending at the end of the strings
  }
  return best;
}

static long badA = 0, badB = 0, total = 0;
static void one(void) {
  dists();
  long va = checkA();
  int cb = checkB();
  total++;
  if (va) { badA++; if (badA <= 2) { printf("  (A) violation: X="); for (int i = 0; i < n; i++) printf("%d", X[i]); printf(" Y="); for (int j = 0; j < m; j++) printf("%d", Y[j]); printf(" ed=%d\n", E[n][m]); } }
  if (cb != E[n][m]) { badB++; if (badB <= 2) { printf("  (B) mismatch: X="); for (int i = 0; i < n; i++) printf("%d", X[i]); printf(" Y="); for (int j = 0; j < m; j++) printf("%d", Y[j]); printf(" ed=%d canon=%d\n", E[n][m], cb); } }
}

int main(int argc, char **argv) {
  int a = 1;
  if (argc > 1 && strcmp(argv[1], "m") == 0) { mutant = 1; a++; }
  const char *mode = argv[a];
  k = atoi(argv[a+1]); n = atoi(argv[a+2]); m = atoi(argv[a+3]);
  if (n > MAXL || m > MAXL) return 2;
  if (strcmp(mode, "mut") == 0) m = MAXL;  /* allocate for the largest possible m */
  SD = (n > m ? n : m) + 1;
  Cr = malloc(sizeof(int) * (size_t)(n + 1) * (m + 1) * 4 * SD);
  if (strcmp(mode, "exh") == 0) {
    long long NX = 1, NY = 1; for (int i = 0; i < n; i++) NX *= k; for (int j = 0; j < m; j++) NY *= k;
    for (long long xi = 0; xi < NX; xi++) {
      long long t = xi; for (int i = n - 1; i >= 0; i--) { X[i] = t % k; t /= k; }
      for (long long yi = 0; yi < NY; yi++) {
        long long u = yi; for (int j = m - 1; j >= 0; j--) { Y[j] = u % k; u /= k; }
        one();
      }
    }
  } else if (strcmp(mode, "mut") == 0) {  // X random of length n, Y = X after e random edits (m ignored)
    long samples = atol(argv[a+4]); rs = strtoull(argv[a+5], 0, 10); int e = atoi(argv[a+6]);
    int nn = n;
    for (long s = 0; s < samples; s++) {
      n = nn;
      int T[2 * MAXL]; int L = n;
      for (int i = 0; i < n; i++) { X[i] = rnd64() % k; T[i] = X[i]; }
      for (int z = 0; z < e; z++) {
        int op = rnd64() % 3;
        if (op == 0 && L > 0) { int p = rnd64() % L; T[p] = rnd64() % k; }
        else if (op == 1 && L < MAXL) { int p = rnd64() % (L + 1); for (int u = L; u > p; u--) T[u] = T[u-1]; T[p] = rnd64() % k; L++; }
        else if (op == 2 && L > 1) { int p = rnd64() % L; for (int u = p; u < L - 1; u++) T[u] = T[u+1]; L--; }
      }
      m = L; for (int j = 0; j < m; j++) Y[j] = T[j];
      one();
    }
  } else {
    long samples = atol(argv[a+4]); rs = strtoull(argv[a+5], 0, 10);
    for (long s = 0; s < samples; s++) {
      for (int i = 0; i < n; i++) X[i] = rnd64() % k;
      for (int j = 0; j < m; j++) Y[j] = rnd64() % k;
      one();
    }
  }
  printf("%s%s k=%d n=%d m=%d: %ld pairs; (A) pairs with an optimal alignment violating rule 5: %ld; "
         "(B) pairs with min rule 1-5 cost != ed: %ld\n", mutant ? "MUTANT " : "", mode, k, n, m, total, badA, badB);
  return 0;
}
