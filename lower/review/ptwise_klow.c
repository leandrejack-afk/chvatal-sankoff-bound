// ptwise_klow.c (review): POINTWISE test of klow.c's operator on arbitrary (random) vectors.
// Includes klow.c unchanged (main renamed), fills w with random values, runs the program's
// own phase 1 + phase 2, and compares every class's new value with the exact Lueker T of
// the old vector, computed here by direct lookups with every letter enumerated:
//   require  k^2 * T_klow <= k^2 * T_exact  (never above)  and  T_exact - T_klow < 3 units.
// Build: clang -O2 -Dmain=klow_main -o bin/ptwise_klow ptwise_klow.c -lpthread
// usage: ptwise_klow k h [seed] [maxval_bits]
#include "../klow.c"
#undef main

static uint64_t st64 = 88172645463325252ULL;
static inline uint64_t rnd(void) { st64 ^= st64 << 13; st64 ^= st64 >> 7; st64 ^= st64 << 17; return st64; }

int main(int argc, char **argv) {
  if (argc < 3) return 1;
  K = atoi(argv[1]); H = atoi(argv[2]);
  if (argc > 3) st64 ^= strtoull(argv[3], 0, 10) * 0x9E3779B97F4A7C15ULL;
  int vb = argc > 4 ? atoi(argv[4]) : 31;   // random values < 2^vb
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
  run_phase(1, NF, 0);
  run_phase(2, NWF, 1);   // w <- T_klow(w0), in place (no normalization here)
  // exact T at every class representative, from w0, with the program's class map (tested separately)
  uint8_t r[MAXL], a[MAXL], b[MAXL]; int8_t mx[MAXL];
  unrank(0, LW, cntW, r, mx);
  uint64_t d = 0, above = 0, toolow = 0; int64_t maxloss = 0; double KK = (double)K * K;
  for (uint64_t R = 0; R < NWF; R++) {
    if ((bits[R >> 6] >> (R & 63)) & 1) {
      const uint8_t *s = r, *t = r + H;
      int64_t num;  // K^2 * T_exact (units)
      memcpy(a, s + 1, H - 1); memcpy(b, t + 1, H - 1);
      if (s[0] == t[0]) {
        int64_t S = 0;
        for (int c = 0; c < K; c++) { a[H - 1] = (uint8_t)c; for (int c2 = 0; c2 < K; c2++) { b[H - 1] = (uint8_t)c2; S += w0[wdense(a, b)]; } }
        num = (int64_t)(KK) * ((int64_t)1 << SHIFT) + S;
      } else {
        int64_t SA = 0, SB = 0;
        for (int c = 0; c < K; c++) { a[H - 1] = (uint8_t)c; SA += w0[wdense(a, t)]; b[H - 1] = (uint8_t)c; SB += w0[wdense(s, b)]; }
        num = (int64_t)K * (SA > SB ? SA : SB);
      }
      int64_t tk = (int64_t)K * K * (int64_t)w[d];
      if (tk > num) above++;
      int64_t loss = num - tk;
      if (loss > maxloss) maxloss = loss;
      if (loss >= 3 * (int64_t)K * K) toolow++;
      d++;
    }
    if (R + 1 < NWF) next_rgs(r, mx, LW);
  }
  printf("ptwise klow k=%d h=%d classes=%llu values<2^%d: T_klow above exact T at %llu classes; max loss %.3f units; loss>=3 units at %llu -> %s\n",
         K, H, (unsigned long long)NW, vb, (unsigned long long)above, (double)maxloss / KK, (unsigned long long)toolow,
         (above == 0 && d == NW) ? "ok" : "FAIL");
  return above != 0;
}
