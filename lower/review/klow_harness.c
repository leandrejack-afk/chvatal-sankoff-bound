// klow_harness.c (review): tests klow.c's OWN class machinery by including the source
// unchanged (main renamed). Build: clang -O2 -Dmain=klow_main -o bin/klow_harness klow_harness.c -lpthread
//
//   klow_harness exh k h      exhaustive tests (small k, h)
//   klow_harness rnd k h n    random tests on big shapes (no bitmap, no big allocation)
#include "../klow.c"
#undef main

static uint64_t hs = 0x243F6A8885A308D3ULL;
static inline uint64_t rnd64(void) { hs ^= hs << 13; hs ^= hs >> 7; hs ^= hs << 17; return hs; }
static long fails;
#define CHECK(c, ...) do { if (!(c)) { if (fails < 10) { fprintf(stderr, __VA_ARGS__); fputc('\n', stderr); } fails++; } } while (0)

static void setup(int k, int h) {
  K = k; H = h; LW = 2 * H; LF = 2 * H - 1;
  memset(cntW, 0, sizeof cntW); memset(cntF, 0, sizeof cntF);
  build_counts(cntW, LW); build_counts(cntF, LF);
  NWF = cntW[1][0]; NF = cntF[1][0];
}

// independent RGS check and lex compare
static int is_rgs(const uint8_t *r, int L) {
  int m = -1;
  for (int i = 0; i < L; i++) { if (r[i] > m + 1 || r[i] >= K) return 0; if (r[i] > m) m = r[i]; }
  return r[0] == 0;
}

// walk every RGS of length L with unrank/next_rgs, compare with unrank(R) and rank2
static void walk(int L, uint64_t N, uint64_t cnt[][MAXK], const char *name) {
  uint8_t r[MAXL], q[MAXL], prev[MAXL]; int8_t mx[MAXL], mq[MAXL];
  unrank(0, L, cnt, r, mx);
  uint64_t R = 0;
  for (;;) {
    CHECK(is_rgs(r, L), "%s: invalid RGS at %llu", name, (unsigned long long)R);
    unrank(R, L, cnt, q, mq);
    CHECK(!memcmp(q, r, L) && !memcmp(mq, mx, L), "%s: unrank(%llu) differs from the walk", name, (unsigned long long)R);
    CHECK(rank2(r, L, r, 0, cnt) == R, "%s: rank2 (one piece) differs at %llu", name, (unsigned long long)R);
    for (int cut = 0; cut <= L; cut += (L > 6 ? 3 : 1))
      CHECK(rank2(r, cut, r + cut, L - cut, cnt) == R, "%s: rank2 split %d differs at %llu", name, cut, (unsigned long long)R);
    int m = -1; for (int i = 0; i < L; i++) { if (r[i] > m) m = r[i]; CHECK(mx[i] == m, "%s: mx wrong", name); }
    if (R) CHECK(memcmp(prev, r, L) < 0, "%s: not increasing at %llu", name, (unsigned long long)R);
    memcpy(prev, r, L);
    R++;
    if (!next_rgs(r, mx, L)) break;
  }
  CHECK(R == N, "%s: walk visited %llu of %llu", name, (unsigned long long)R, (unsigned long long)N);
}

static void decode(uint64_t x, int n, uint8_t *s) { for (int i = n - 1; i >= 0; i--) { s[i] = x % K; x /= K; } }
static uint64_t encode(const uint8_t *s, int n) { uint64_t x = 0; for (int i = 0; i < n; i++) x = x * K + s[i]; return x; }

// brute-force orbit key of a pair under S_K x swap: min over all permutations of both orientations
static int perms[40320][8], nperm;
static void gen_perms(int *p, int i) {
  if (i == K) { memcpy(perms[nperm++], p, sizeof(int) * K); return; }
  for (int j = i; j < K; j++) { int t = p[i]; p[i] = p[j]; p[j] = t; gen_perms(p, i + 1); t = p[i]; p[i] = p[j]; p[j] = t; }
}
static uint64_t orbit_key(const uint8_t *s, const uint8_t *t, int n1, int n2, int swap_ok) {
  uint64_t best = UINT64_MAX; uint8_t a[MAXL];
  for (int p = 0; p < nperm; p++) {
    for (int i = 0; i < n1; i++) a[i] = perms[p][s[i]];
    for (int i = 0; i < n2; i++) a[n1 + i] = perms[p][t[i]];
    uint64_t c = encode(a, n1 + n2); if (c < best) best = c;
    if (swap_ok) {
      for (int i = 0; i < n2; i++) a[i] = perms[p][t[i]];
      for (int i = 0; i < n1; i++) a[n2 + i] = perms[p][s[i]];
      c = encode(a, n1 + n2); if (c < best) best = c;
    }
  }
  return best;
}

static void exhaustive(int k, int h) {
  setup(k, h);
  walk(LW, NWF, cntW, "W");
  walk(LF, NF, cntF, "F");
  // canonical bitmap with the program's own phase 0
  NT = 2;
  uint64_t nwords = (NWF + 63) / 64;
  bits = calloc(nwords + 1, 8); cum = calloc(nwords + 1, 8);
  run_phase(0, NWF, 1);
  NW = 0;
  for (uint64_t i = 0; i < nwords; i++) { cum[i] = NW; NW += (uint64_t)__builtin_popcountll(bits[i]); }
  cum[nwords] = NW;
  int p[8]; for (int i = 0; i < K; i++) p[i] = i;
  nperm = 0;
  uint64_t PH = 1; for (int i = 0; i < H; i++) PH *= K;
  uint64_t N = PH * PH;
  double fact = 1; for (int i = 2; i <= K; i++) fact *= i;
  int doperm = K <= 8 && fact * (double)N <= 3e8;
  if (doperm) gen_perms(p, 0);
  // pair -> dense index (program) and pair -> orbit key (brute force); partitions must agree
  uint64_t *dk = malloc(N * 8), *ok = malloc(N * 8);
  uint8_t s[MAXL], t[MAXL], ps[MAXL], pt[MAXL];
  for (uint64_t x = 0; x < N; x++) {
    decode(x / PH, H, s); decode(x % PH, H, t);
    uint64_t d = wdense(s, t);
    CHECK(d < NW, "dense out of range");
    CHECK(d == wdense(t, s), "swap changes the class");
    // a random permutation, letters remapped (also to labels other than 0..K-1 order)
    int q[32]; for (int i = 0; i < K; i++) q[i] = i;
    for (int i = K - 1; i > 0; i--) { int j = rnd64() % (i + 1), tt = q[i]; q[i] = q[j]; q[j] = tt; }
    for (int i = 0; i < H; i++) { ps[i] = q[s[i]]; pt[i] = q[t[i]]; }
    CHECK(d == wdense(ps, pt), "permutation changes the class");
    dk[x] = d;
    ok[x] = doperm ? orbit_key(s, t, H, H, 1) : 0;
  }
  if (doperm) {
    // same orbit <=> same dense index: map orbit key -> dense must be a bijection
    // O(N) check via arrays: orbit key -> first dense seen, dense -> first key seen
    uint64_t maxkey = N;  // keys are pair codes < N
    uint64_t *k2d = malloc(maxkey * 8), *d2k = malloc(NW * 8);
    for (uint64_t i = 0; i < maxkey; i++) k2d[i] = UINT64_MAX;
    for (uint64_t i = 0; i < NW; i++) d2k[i] = UINT64_MAX;
    uint64_t norb = 0;
    for (uint64_t x = 0; x < N; x++) {
      if (k2d[ok[x]] == UINT64_MAX) { k2d[ok[x]] = dk[x]; norb++; }
      else CHECK(k2d[ok[x]] == dk[x], "one orbit split across classes");
      if (d2k[dk[x]] == UINT64_MAX) d2k[dk[x]] = ok[x];
      else CHECK(d2k[dk[x]] == ok[x], "two orbits merged into one class");
    }
    CHECK(norb == NW, "orbit count %llu != classes %llu", (unsigned long long)norb, (unsigned long long)NW);
    free(k2d); free(d2k);
  }
  // F space: fidx invariant under permutations, distinct on alphabet orbits, in range
  uint64_t P = PH / K, NFp = P * PH;
  uint64_t *fk = malloc(NFp * 8), *fo = malloc(NFp * 8);
  for (uint64_t x = 0; x < NFp; x++) {
    decode(x / PH, H - 1, s); decode(x % PH, H, t);
    uint64_t f = fidx(s, t);
    CHECK(f < NF, "fidx out of range");
    int q[32]; for (int i = 0; i < K; i++) q[i] = i;
    for (int i = K - 1; i > 0; i--) { int j = rnd64() % (i + 1), tt = q[i]; q[i] = q[j]; q[j] = tt; }
    for (int i = 0; i < H; i++) { ps[i] = q[s[i]]; pt[i] = q[t[i]]; }
    CHECK(f == fidx(ps, pt), "fidx not permutation invariant");
    fk[x] = f; fo[x] = doperm ? orbit_key(s, t, H - 1, H, 0) : 0;
  }
  if (doperm) {
    uint64_t *k2f = malloc(NFp * 8), *f2k = malloc(NF * 8), norb = 0;
    for (uint64_t i = 0; i < NFp; i++) k2f[i] = UINT64_MAX;
    for (uint64_t i = 0; i < NF; i++) f2k[i] = UINT64_MAX;
    for (uint64_t x = 0; x < NFp; x++) {
      if (k2f[fo[x]] == UINT64_MAX) { k2f[fo[x]] = fk[x]; norb++; } else CHECK(k2f[fo[x]] == fk[x], "F orbit split");
      if (f2k[fk[x]] == UINT64_MAX) f2k[fk[x]] = fo[x]; else CHECK(f2k[fk[x]] == fo[x], "F orbits merged");
    }
    CHECK(norb == NF, "F orbit count %llu != NF %llu", (unsigned long long)norb, (unsigned long long)NF);
    free(k2f); free(f2k);
  }
  // fresh-letter label MAXK must act like any absent letter, in both lookups
  for (int trial = 0; trial < 20000; trial++) {
    uint8_t x[MAXL], y[MAXL], a[MAXL];
    for (int i = 0; i < H - 1; i++) x[i] = rnd64() % K;
    for (int i = 0; i < H; i++) y[i] = rnd64() % K;
    uint32_t seen = 0; for (int i = 0; i < H - 1; i++) seen |= 1u << x[i]; for (int i = 0; i < H; i++) seen |= 1u << y[i];
    (void)seen;
    // y' = y with last letter replaced by MAXK vs by an absent letter, as used in phase 2: fidx(t', s'c)
    memcpy(a, y, H);
    uint32_t seen2 = 0; for (int i = 0; i < H - 1; i++) seen2 |= 1u << x[i]; for (int i = 0; i < H - 1; i++) seen2 |= 1u << y[i];
    if (__builtin_popcount(seen2) == K) continue;
    a[H - 1] = MAXK; uint64_t f1 = fidx(x, a);
    for (int c = 0; c < K; c++) if (!(seen2 & (1u << c))) { a[H - 1] = (uint8_t)c; CHECK(fidx(x, a) == f1, "fresh label differs"); }
  }
  printf("exh k=%d h=%d: W rgs %llu, F rgs %llu, classes %llu, pairs %llu, perms %d -> %s (%ld failures)\n", k, h,
         (unsigned long long)NWF, (unsigned long long)NF, (unsigned long long)NW, (unsigned long long)N, nperm,
         fails ? "FAIL" : "ok", fails);
  free(dk); free(ok); free(fk); free(fo); free(bits); free(cum);
}

static void random_big(int k, int h, long n) {
  setup(k, h);
  uint8_t r[MAXL], q[MAXL], s[MAXL], t[MAXL], ps[MAXL], pt[MAXL]; int8_t mx[MAXL], mq[MAXL];
  for (long i = 0; i < n; i++) {
    int L = (i & 1) ? LF : LW; uint64_t N = (i & 1) ? NF : NWF;
    uint64_t (*cnt)[MAXK] = (i & 1) ? cntF : cntW;
    uint64_t R = rnd64() % N;
    if (i % 1000 == 0) R = (i % 2000 == 0) ? 0 : N - 1;  // the ends
    unrank(R, L, cnt, r, mx);
    CHECK(is_rgs(r, L), "invalid RGS from unrank");
    CHECK(rank2(r, L / 2, r + L / 2, L - L / 2, cnt) == R, "rank2(unrank(R)) != R");
    if (R + 1 < N) {
      memcpy(q, r, L); memcpy(mq, mx, L);
      CHECK(next_rgs(q, mq, L) == 1, "next_rgs ended early");
      uint8_t u[MAXL]; int8_t mu[MAXL];
      unrank(R + 1, L, cnt, u, mu);
      CHECK(!memcmp(u, q, L) && !memcmp(mu, mq, L), "next_rgs != unrank(R+1)");
      CHECK(memcmp(r, q, L) < 0, "successor not larger");
    } else {
      memcpy(q, r, L); memcpy(mq, mx, L);
      CHECK(next_rgs(q, mq, L) == 0, "next_rgs past the end");
    }
    // class key min(rank s.t, rank t.s) invariant under permutation and swap, on raw pairs
    for (int j = 0; j < H; j++) { s[j] = rnd64() % K; t[j] = rnd64() % K; }
    if (rnd64() & 1) t[0] = s[0];
    int pm[32]; for (int j = 0; j < K; j++) pm[j] = j;
    for (int j = K - 1; j > 0; j--) { int z = rnd64() % (j + 1), tt = pm[j]; pm[j] = pm[z]; pm[z] = tt; }
    for (int j = 0; j < H; j++) { ps[j] = pm[s[j]]; pt[j] = pm[t[j]]; }
    uint64_t a1 = rank2(s, H, t, H, cntW), a2 = rank2(t, H, s, H, cntW);
    uint64_t b1 = rank2(pt, H, ps, H, cntW), b2 = rank2(ps, H, pt, H, cntW);
    CHECK((a1 < a2 ? a1 : a2) == (b1 < b2 ? b1 : b2), "class key not invariant");
    CHECK(a1 < NWF && a2 < NWF, "rank out of range");
  }
  printf("rnd k=%d h=%d: rgs %llu F %llu, %ld samples -> %s (%ld failures)\n", k, h, (unsigned long long)NWF,
         (unsigned long long)NF, n, fails ? "FAIL" : "ok", fails);
}

int main(int argc, char **argv) {
  if (argc < 4) return 1;
  int k = atoi(argv[2]), h = atoi(argv[3]);
  if (!strcmp(argv[1], "exh")) exhaustive(k, h);
  else random_big(k, h, argc > 4 ? atol(argv[4]) : 1000000);
  return fails ? 1 : 0;
}
