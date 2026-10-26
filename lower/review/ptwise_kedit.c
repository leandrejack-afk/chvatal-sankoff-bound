// ptwise_kedit.c (review): POINTWISE test of kedit.c's operator and per-class certificate on
// arbitrary (random) vectors and a random rho. Includes kedit.c unchanged (main renamed), runs
// the program's own phase 1 + phase 2, and compares with exact values computed here:
//   T_exact(v, v - rho): a == b: G - rho;  a != b: 1 + min(A, B, G - rho)
//   (G, A, B exact rationals, every letter enumerated, B via v(as', t'c) directly)
//   require k^2 * (stored - rho) >= k^2 * T_exact (never below) and gap < 3 units;
//   the program's max per-class certificate (job need) must be >= the exact max least r.
// Build: clang -O2 -Dmain=kedit_main -o bin/ptwise_kedit ptwise_kedit.c -lpthread
// usage: ptwise_kedit k h [seed] [maxval_bits]
#include "../kedit.c"
#undef main

static uint64_t st64 = 88172645463325252ULL;
static inline uint64_t rnd(void) { st64 ^= st64 << 13; st64 ^= st64 >> 7; st64 ^= st64 << 17; return st64; }

int main(int argc, char **argv) {
  if (argc < 3) return 1;
  K = atoi(argv[1]); H = atoi(argv[2]);
  if (argc > 3) st64 ^= strtoull(argv[3], 0, 10) * 0x9E3779B97F4A7C15ULL;
  int vb = argc > 4 ? atoi(argv[4]) : 30;
  NT = 2; LW = 2 * H; LF = 2 * H - 1;
  build_counts(cntW, LW); build_counts(cntF, LF);
  NWF = cntW[1][0]; NF = cntF[1][0];
  uint64_t nwords = (NWF + 63) / 64;
  bits = calloc(nwords + 1, 8); cum = calloc(nwords + 1, 8);
  F = malloc(NF * sizeof(val_t));
  run_phase(0, NWF, 1);
  NW = 0;
  for (uint64_t i = 0; i < nwords; i++) { cum[i] = NW; NW += (uint64_t)__builtin_popcountll(bits[i]); }
  cum[nwords] = NW;
  w = malloc(NW * sizeof(val_t));
  val_t *w0 = malloc(NW * sizeof(val_t));
  for (uint64_t i = 0; i < NW; i++) w0[i] = w[i] = (val_t)(rnd() & ((1ULL << vb) - 1));
  RHO = (int64_t)(rnd() % ((uint64_t)1 << 27));   // a random growth estimate in [0, 1)
  run_phase(1, NF, 0);
  run_phase(2, NWF, 1);   // w <- T(w0, w0 - RHO) + RHO, in place
  int64_t need = INT64_MIN;
  for (int i = 0; i < last_nj; i++) if (jobs[i].hi > jobs[i].lo && jobs[i].need > need) need = jobs[i].need;
  uint8_t r[MAXL], a[MAXL], b[MAXL]; int8_t mx[MAXL];
  unrank(0, LW, cntW, r, mx);
  const int64_t KK = (int64_t)K * K, U = (int64_t)1 << SHIFT;
  uint64_t d = 0, below = 0, toohigh = 0; int64_t maxgap = 0, rstar = INT64_MIN;  // rstar: 2KK * exact least r
  for (uint64_t R = 0; R < NWF; R++) {
    if ((bits[R >> 6] >> (R & 63)) & 1) {
      const uint8_t *s = r, *t = r + H;
      int64_t v = w0[d], S = 0, SA = 0, SB = 0;
      memcpy(a, s + 1, H - 1); memcpy(b, t + 1, H - 1);
      for (int c = 0; c < K; c++) { a[H - 1] = (uint8_t)c; for (int c2 = 0; c2 < K; c2++) { b[H - 1] = (uint8_t)c2; S += w0[wdense(a, b)]; } }
      int64_t texact, lr;   // KK * T_exact(v, v - RHO); 2KK * least r
      if (s[0] == t[0]) {
        texact = S - KK * RHO;
        lr = S - KK * v;
      } else {
        for (int c = 0; c < K; c++) { a[H - 1] = (uint8_t)c; SA += w0[wdense(a, t)]; b[H - 1] = (uint8_t)c; SB += w0[wdense(s, b)]; }
        int64_t mab = SA < SB ? SA : SB;   // K * min(A, B)
        int64_t o1 = K * mab, o2 = S - KK * RHO;
        texact = KK * U + (o1 < o2 ? o1 : o2);
        int64_t n1 = 2 * KK * U + 2 * K * mab - 2 * KK * v, n2 = KK * U + S - KK * v;
        lr = n1 < n2 ? n1 : n2;
      }
      if (lr > rstar) rstar = lr;
      int64_t tk = KK * ((int64_t)w[d] - RHO);
      if (tk < texact) below++;
      int64_t gap = tk - texact;
      if (gap > maxgap) maxgap = gap;
      if (gap >= 3 * KK) toohigh++;
      d++;
    }
    if (R + 1 < NWF) next_rgs(r, mx, LW);
  }
  int certok = (__int128)need * 2 * KK >= (__int128)rstar;
  printf("ptwise kedit k=%d h=%d classes=%llu rho=%.6f: T below exact at %llu classes; max excess %.3f units; >=3 units at %llu; "
         "cert need %lld vs exact least %.3f units (excess %.3f) -> %s\n",
         K, H, (unsigned long long)NW, (double)RHO / U, (unsigned long long)below, (double)maxgap / KK,
         (unsigned long long)toohigh, (long long)need, (double)rstar / (2 * KK), (double)need - (double)rstar / (2 * KK),
         (below == 0 && certok && d == NW) ? "ok" : "FAIL");
  return !(below == 0 && certok);
}
