// count_canon.c (review code): exhaustive over all pairs (X, Y) in [k]^n x [k]^m.
// For each pair: (1) Levenshtein distance ed; (2) the number of canonical alignments of each cost,
// by a DP whose state is (i, j, t, p): t = 0 last move diagonal (or start), 1 inside an I-run,
// 2 inside a D-run, 3 M right after a D-run, 4 M right after an I-run; p = index where the
// current run started (in Y for I-runs, in X for D-runs), so rule 3 is tested by scanning the run
// letters directly (no letter masks, unlike canon_check2/3).  Rules: 1 no I-D/D-I adjacency;
// 2 no S right after I or D; 3 an M right after a run avoids every run letter; 4 (if rules=4) no
// D M I and no I M D.
// Output: per cost c, the total count sum_pairs N_c(pair) (compare with k^(n+m) * automaton),
// the number of pairs with ed = c (first-moment check P(ed <= c) <= E[N_{<=c}]), and the
// number of pairs whose minimum canonical cost differs from ed (existence check).
// usage: count_canon rules k n m
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MX 12
#define MC (2 * MX + 1)
static int K, RULES, n, m, X[MX], Y[MX];
// cnt[i][j][t][p][c], flat, sized to (n, m)
static uint64_t *cntv; static int DP_, DC_;
#define CNT(i, j, t, p, c) cntv[((((size_t)(i) * (m + 1) + (j)) * 5 + (t)) * DP_ + (p)) * DC_ + (c)]

static int avoid_I(int p, int j, int letter) {  // letter not in Y[p..j-1]
  for (int u = p; u < j; u++) if (Y[u] == letter) return 0;
  return 1;
}
static int avoid_D(int p, int i, int letter) {  // letter not in X[p..i-1]
  for (int u = p; u < i; u++) if (X[u] == letter) return 0;
  return 1;
}

int main(int argc, char **argv) {
  RULES = atoi(argv[1]); K = atoi(argv[2]); n = atoi(argv[3]); m = atoi(argv[4]);
  if (n > MX || m > MX) return 2;
  long long NX = 1, NY = 1;
  for (int i = 0; i < n; i++) NX *= K;
  for (int j = 0; j < m; j++) NY *= K;
  int maxc = n + m;
  DP_ = (n > m ? n : m) + 1; DC_ = maxc + 2;
  size_t NC = (size_t)(n + 1) * (m + 1) * 5 * DP_ * DC_;
  cntv = calloc(NC, sizeof(uint64_t));
  static unsigned long long total[MC], ed_hist[MC];
  long long mism = 0, pairs = 0;
  static int E[MX + 1][MX + 1];
  for (long long xi = 0; xi < NX; xi++) {
    long long tt = xi; for (int i = 0; i < n; i++) { X[i] = tt % K; tt /= K; }
    for (long long yi = 0; yi < NY; yi++) {
      long long uu = yi; for (int j = 0; j < m; j++) { Y[j] = uu % K; uu /= K; }
      for (int i = 0; i <= n; i++) for (int j = 0; j <= m; j++) {
        if (!i || !j) { E[i][j] = i + j; continue; }
        int b = E[i - 1][j - 1] + (X[i - 1] != Y[j - 1]);
        if (E[i - 1][j] + 1 < b) b = E[i - 1][j] + 1;
        if (E[i][j - 1] + 1 < b) b = E[i][j - 1] + 1;
        E[i][j] = b;
      }
      memset(cntv, 0, NC * sizeof(uint64_t));
      CNT(0, 0, 0, 0, 0) = 1;
      for (int i = 0; i <= n; i++) for (int j = 0; j <= m; j++) for (int t = 0; t < 5; t++) for (int p = 0; p < DP_; p++) {
        int need_p = (t == 1 || t == 2);
        if (!need_p && p) continue;
        for (int c = 0; c <= maxc; c++) {
          uint64_t w = CNT(i, j, t, p, c); if (!w) continue;
          if (i < n && j < m) {
            if (X[i] == Y[j]) {
              int ok = 1, nt = 0;
              if (t == 1) { ok = avoid_I(p, j, X[i]); nt = RULES == 4 ? 4 : 0; }
              if (t == 2) { ok = avoid_D(p, i, Y[j]); nt = RULES == 4 ? 3 : 0; }
              if (ok) CNT(i + 1, j + 1, nt, 0, c) += w;
            } else if (!need_p) {
              CNT(i + 1, j + 1, 0, 0, c + 1) += w;
            }
          }
          if (j < m && t != 2 && !(RULES == 4 && t == 3)) CNT(i, j + 1, 1, t == 1 ? p : j, c + 1) += w;
          if (i < n && t != 1 && !(RULES == 4 && t == 4)) CNT(i + 1, j, 2, t == 2 ? p : i, c + 1) += w;
        }
      }
      int minc = -1;
      for (int c = 0; c <= maxc; c++) {
        uint64_t s = 0;
        for (int t = 0; t < 5; t++) for (int p = 0; p < DP_; p++) s += CNT(n, m, t, p, c);
        total[c] += s;
        if (s && minc < 0) minc = c;
      }
      ed_hist[E[n][m]]++;
      if (minc != E[n][m]) mism++;
      pairs++;
    }
  }
  printf("rules=%d k=%d n=%d m=%d pairs=%lld existence-mismatches=%lld\n", RULES, K, n, m, pairs, mism);
  unsigned long long cum_ed = 0, cum_cnt = 0; int fm_ok = 1;
  printf("totals");
  for (int c = 0; c <= maxc; c++) printf(" %llu", total[c]);
  printf("\n");
  for (int c = 0; c <= maxc; c++) {
    cum_ed += ed_hist[c]; cum_cnt += total[c];
    if (cum_ed > cum_cnt) fm_ok = 0;
  }
  printf("first-moment P(ed<=c) <= E[N_<=c] for all c: %s\n", fm_ok ? "yes" : "NO");
  return mism != 0 || !fm_ok;
}
