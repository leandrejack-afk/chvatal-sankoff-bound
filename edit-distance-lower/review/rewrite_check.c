// rewrite_check.c (review code): run the README's REWRITING PROCEDURE on random optimal alignments
// and check every claim of the proof step by step.
//   - an optimal alignment is drawn by backtracking the Levenshtein table with a tie-break policy
//     (0 uniform random, 1 prefer I, 2 prefer D, 3 prefer indels, 4 prefer diagonals, 5 prefer
//     late diagonals = indels first when walking backwards is policy 3, so 5 = random but biased
//     toward runs), giving adversarial non-canonical starting points;
//   - while some rule is violated, ONE violation is picked at random and rewritten:
//       rule 2: I S -> S I, D S -> S D (the new diagonal pair must be unequal, else ERROR)
//       rule 3: I^r M -> I^s M I^(r-s) for a RANDOM matching offset s (D likewise)
//       rule 4: D M I -> S S, I M D -> S S (both new pairs must be unequal, else ERROR)
//     rule 1 (I D / D I adjacency) must never occur (ERROR if seen);
//   - after every rewrite: alignment valid (letters consumed in order, M equal, S unequal), cost
//     unchanged, potential (#indels, sum over diagonal moves of #indels before it) strictly smaller
//     lexicographically (ERROR otherwise);
//   - at the fixpoint all rules hold and cost == ed.
// usage: rewrite_check rules(3|4) k nmin nmax trials seed [skew]
//   skew 0: uniform letters; 1: letter 0 with prob 1/2; 2: random period-p strings with noise.
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static uint64_t s[4];
static inline uint64_t rotl(const uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static uint64_t nxt(void) {
  const uint64_t r = rotl(s[1] * 5, 7) * 9, t = s[1] << 17;
  s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t; s[3] = rotl(s[3], 45);
  return r;
}
static int rnd(int n) { return (int)(nxt() % (uint64_t)n); }

#ifndef MUT
#define MUT 0  // mutation tests: 1 rule-3 rewrite to the run start (wrong offset), 3 bogus D M M I -> S S S
#endif
#define MAXN 64
#define MAXL (2 * MAXN + 2)
static int K, RULES;
static int X[MAXN], Y[MAXN], n, m;
static int E[MAXN + 1][MAXN + 1];
static char A[MAXL];  // move types
static int L;

static void fill_string(int *Z, int len, int skew) {
  if (skew == 0) { for (int i = 0; i < len; i++) Z[i] = rnd(K); return; }
  if (skew == 1) { for (int i = 0; i < len; i++) Z[i] = rnd(2) ? 0 : rnd(K); return; }
  int p = 1 + rnd(3), base[3];
  for (int i = 0; i < p; i++) base[i] = rnd(K);
  for (int i = 0; i < len; i++) Z[i] = rnd(5) ? base[i % p] : rnd(K);
}

static void edtable(void) {
  for (int i = 0; i <= n; i++) for (int j = 0; j <= m; j++) {
    if (!i) { E[i][j] = j; continue; }
    if (!j) { E[i][j] = i; continue; }
    int b = E[i - 1][j - 1] + (X[i - 1] != Y[j - 1]);
    if (E[i - 1][j] + 1 < b) b = E[i - 1][j] + 1;
    if (E[i][j - 1] + 1 < b) b = E[i][j - 1] + 1;
    E[i][j] = b;
  }
}

// backtrack an optimal alignment with a tie-break policy; builds A in forward order
static void backtrack(int policy) {
  char rev[MAXL]; int r = 0, i = n, j = m;
  while (i > 0 || j > 0) {
    int opts[3], no = 0;  // 0 diag, 1 I (consume Y), 2 D (consume X)
    if (i > 0 && j > 0 && E[i][j] == E[i - 1][j - 1] + (X[i - 1] != Y[j - 1])) opts[no++] = 0;
    if (j > 0 && E[i][j] == E[i][j - 1] + 1) opts[no++] = 1;
    if (i > 0 && E[i][j] == E[i - 1][j] + 1) opts[no++] = 2;
    int c = opts[rnd(no)];
    int pol = policy == 5 ? rnd(5) : policy;
    for (int q = 0; q < no; q++) {
      if (pol == 1 && opts[q] == 1) c = 1;
      if (pol == 2 && opts[q] == 2) c = 2;
      if (pol == 3 && opts[q] != 0) { c = opts[q]; }
      if (pol == 4 && opts[q] == 0) c = 0;
    }
    if (c == 0) { rev[r++] = X[i - 1] == Y[j - 1] ? 'M' : 'S'; i--; j--; }
    else if (c == 1) { rev[r++] = 'I'; j--; }
    else { rev[r++] = 'D'; i--; }
  }
  L = r;
  for (int q = 0; q < r; q++) A[q] = rev[r - 1 - q];
}

// positions before each move; returns cost, or -1 if invalid
static int PI[MAXL + 1], PJ[MAXL + 1];
static int scan(void) {
  int i = 0, j = 0, cost = 0;
  for (int q = 0; q < L; q++) {
    PI[q] = i; PJ[q] = j;
    char c = A[q];
    if (c == 'M' || c == 'S') {
      if (i >= n || j >= m) return -1;
      if ((c == 'M') != (X[i] == Y[j])) return -1;
      cost += c == 'S'; i++; j++;
    } else if (c == 'I') { if (j >= m) return -1; cost++; j++; }
    else if (c == 'D') { if (i >= n) return -1; cost++; i++; }
    else return -1;
  }
  PI[L] = i; PJ[L] = j;
  if (i != n || j != m) return -1;
  return cost;
}
static void potential(long *ind, long *phi) {
  long c = 0, p = 0;
  for (int q = 0; q < L; q++) {
    if (A[q] == 'I' || A[q] == 'D') c++;
    else p += c;
  }
  *ind = c; *phi = p;
}

// violation list: (kind, index). kind 1: I-D adjacency, 2: indel then S, 3: run then M with a
// matching run letter (index = index of the M), 4: D M I or I M D (index = index of the M)
static int VK[4 * MAXL], VI[4 * MAXL], NV;
static int run_start(int q) {  // q: index of an M preceded by an indel run; returns run start
  char c = A[q - 1]; int p = q - 1;
  while (p - 1 >= 0 && A[p - 1] == c) p--;
  return p;
}
static void violations(void) {
  NV = 0;
  for (int q = 0; q + 1 < L; q++) {
    char a = A[q], b = A[q + 1];
    if ((a == 'I' && b == 'D') || (a == 'D' && b == 'I')) { VK[NV] = 1; VI[NV++] = q; }
    if ((a == 'I' || a == 'D') && b == 'S') { VK[NV] = 2; VI[NV++] = q; }
  }
  for (int q = 1; q < L; q++) if (A[q] == 'M' && (A[q - 1] == 'I' || A[q - 1] == 'D')) {
    int p = run_start(q), bad = 0;
    for (int u = p; u < q; u++) {
      int letter = A[u] == 'I' ? Y[PJ[u]] : X[PI[u]];
      int ml = X[PI[q]];
      if (letter == ml) bad = 1;
    }
    if (bad) { VK[NV] = 3; VI[NV++] = q; }
    if (RULES == 4 && q + 1 < L && ((A[q - 1] == 'D' && A[q + 1] == 'I') || (A[q - 1] == 'I' && A[q + 1] == 'D'))) {
      VK[NV] = 4; VI[NV++] = q;
    }
    if (MUT == 3 && q + 2 < L && A[q + 1] == 'M' && ((A[q - 1] == 'D' && A[q + 2] == 'I') || (A[q - 1] == 'I' && A[q + 2] == 'D'))) {
      VK[NV] = 5; VI[NV++] = q;
    }
  }
}

static long errors = 0, total_rewrites = 0, max_rewrites = 0, trials_done = 0;
static long cnt_kind[6];
static void fail(const char *why) {
  errors++;
  if (errors <= 5) {
    printf("ERROR %s: X=", why); for (int i = 0; i < n; i++) printf("%d", X[i]);
    printf(" Y="); for (int j = 0; j < m; j++) printf("%d", Y[j]);
    printf(" A=%.*s\n", L, A);
  }
}

static int one_trial(int policy) {
  edtable();
  backtrack(policy);
  int cost0 = scan();
  if (cost0 != E[n][m]) { fail("backtrack not optimal"); return 0; }
  long rw = 0;
  for (;;) {
    if (scan() != cost0) { fail("cost/validity changed"); return 0; }
    violations();
    if (NV == 0) break;
    int v = rnd(NV), q = VI[v];
    long i0, p0, i1, p1; potential(&i0, &p0);
    char save[MAXL]; int Ls = L; memcpy(save, A, L);
    cnt_kind[VK[v]]++;
    if (VK[v] == 1) { fail("rule 1 violated in an optimal alignment"); return 0; }
    if (VK[v] == 2) {  // indel at q, S at q+1 -> S at q, indel at q+1
      int i = PI[q], j = PJ[q];
      if (X[i] == Y[j]) { fail("rule-2 swap creates a match (alignment not optimal)"); return 0; }
      char ind = A[q]; A[q] = 'S'; A[q + 1] = ind;
    } else if (VK[v] == 3) {  // run [p, q) then M at q
      int p = run_start(q), ml = X[PI[q]], offs[MAXL], no = 0;
      for (int u = p; u < q; u++) {
        int letter = A[u] == 'I' ? Y[PJ[u]] : X[PI[u]];
        if (letter == ml) offs[no++] = u;
      }
      int u = offs[rnd(no)];  // M moves to index u, run letters u..q-1 shift right
      if (MUT == 1) u = p;
      char ind = A[p];
      for (int z = q; z > u; z--) A[z] = ind;
      A[u] = 'M';
    } else if (VK[v] == 5) {  // MUT 3: bogus D M M I -> S S S
      A[q - 1] = 'S'; A[q] = 'S'; A[q + 1] = 'S';
      memmove(A + q + 2, A + q + 3, L - q - 3); L--;
    } else {  // rule 4 at M index q: A[q-1] A[q] A[q+1] = D M I or I M D -> S S
      int i = PI[q - 1], j = PJ[q - 1];
      // the two new diagonal pairs are (X[i], Y[j]) and (X[i+1], Y[j+1])
      if (X[i] == Y[j] || X[i + 1] == Y[j + 1]) { fail("rule-4 rewrite creates a match (not optimal)"); return 0; }
      A[q - 1] = 'S'; A[q] = 'S';
      memmove(A + q + 1, A + q + 2, L - q - 2); L--;
    }
    rw++;
    int c1 = scan();
    if (c1 != cost0) { memcpy(A, save, Ls); L = Ls; fail("rewrite changed cost or broke validity"); return 0; }
    potential(&i1, &p1);
    if (!(i1 < i0 || (i1 == i0 && p1 < p0))) { fail("potential did not drop"); return 0; }
    if (rw > 100000) { fail("no termination"); return 0; }
  }
  // fixpoint: all rules hold (violations() found none, including rule 1)
  total_rewrites += rw; if (rw > max_rewrites) max_rewrites = rw;
  return 1;
}

int main(int argc, char **argv) {
  if (argc < 7) { fprintf(stderr, "usage: rewrite_check rules k nmin nmax trials seed [skew]\n"); return 2; }
  RULES = atoi(argv[1]); K = atoi(argv[2]);
  int nmin = atoi(argv[3]), nmax = atoi(argv[4]); long trials = atol(argv[5]);
  uint64_t seed = strtoull(argv[6], 0, 10); int skew = argc > 7 ? atoi(argv[7]) : 0;
  s[0] = seed ^ 0x9E3779B97F4A7C15ULL; s[1] = seed * 0xBF58476D1CE4E5B9ULL + 1; s[2] = 0x94D049BB133111EBULL; s[3] = seed + 12345;
  for (int z = 0; z < 20; z++) nxt();
  for (long tr = 0; tr < trials; tr++) {
    n = nmin + rnd(nmax - nmin + 1); m = nmin + rnd(nmax - nmin + 1);
    fill_string(X, n, skew == 3 ? rnd(3) : skew); fill_string(Y, m, skew == 3 ? rnd(3) : skew);
    one_trial(rnd(6));
    trials_done++;
  }
  printf("rules=%d k=%d n,m in [%d,%d] skew=%d trials=%ld errors=%ld rewrites: total %ld max %ld "
         "(rule2 %ld, rule3 %ld, rule4 %ld, rule1-seen %ld)\n",
         RULES, K, nmin, nmax, skew, trials_done, errors, total_rewrites, max_rewrites,
         cnt_kind[2], cnt_kind[3], cnt_kind[4], cnt_kind[1]);
  return errors != 0;
}
