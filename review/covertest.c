/* Reviewer tool: adversarial covering test for M_h (independent code).
 * For many string pairs P: full LCS DP, canonical backtrack (up, then left, then diagonal),
 * cut after each match; check each piece is T_{c,l} (minimal match pair) and the tail is null;
 * feed the list (mu-truncated) to (a) the unreduced automaton and (b) the complement-reduced
 * automaton with complement tracking; neither may reject, and (b)'s state must equal
 * reduce((a)'s state).  Also checks the complement symmetry on every visited state.
 * usage: covertest <h> <pairs> <maxlen> <seed> */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

static int h;
static uint64_t rng = 88172645463325252ULL;
static uint64_t xr(void) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng; }

static int my_delta(uint64_t key, int c, int l, uint64_t *out, int reduce) {
    uint32_t m = (1u << h) - 1;
    uint32_t x = (uint32_t)(key >> (3 * h)) & m, y = (uint32_t)(key >> (2 * h)) & m;
    uint32_t dx = (uint32_t)(key >> h) & m, dy = (uint32_t)key & m;
    int A[64], B[64], n = 0, mm = 0;
    for (int t = h - 1; t >= 0; t--) A[n++] = (x >> t) & 1;
    for (int t = h - 1; t >= 0; t--) B[mm++] = (y >> t) & 1;
    for (int t = 0; t < l; t++) A[n++] = 1 - c;
    for (int t = 0; t < -l; t++) B[mm++] = 1 - c;
    A[n++] = c; B[mm++] = c;
    int V[64][64];
    V[h][h] = 77;
    for (int i = h; i >= 1; i--) V[i - 1][h] = V[i][h] - (int)((dx >> (h - i)) & 1);
    for (int j = h; j >= 1; j--) V[h][j - 1] = V[h][j] - (int)((dy >> (h - j)) & 1);
    for (int i = 0; i <= n; i++) for (int j = 0; j <= mm; j++) {
        if (i <= h && j <= h) continue;
        if (i == 0) V[i][j] = V[0][j - 1];
        else if (j == 0) V[i][j] = V[i - 1][0];
        else if (A[i - 1] == B[j - 1]) V[i][j] = V[i - 1][j - 1] + 1;
        else V[i][j] = V[i - 1][j] > V[i][j - 1] ? V[i - 1][j] : V[i][j - 1];
    }
    int i = n, j = mm, hit = 0;
    while (i >= h && j >= h) {
        if (i == h && j == h) { hit = 1; break; }
        if (V[i][j] == V[i - 1][j]) i--; else if (V[i][j] == V[i][j - 1]) j--; else { i--; j--; }
    }
    if (!hit) return 0;
    uint32_t nx = 0, ny = 0, ndx = 0, ndy = 0;
    for (int t = n - h; t < n; t++) nx = (nx << 1) | (uint32_t)A[t];
    for (int t = mm - h; t < mm; t++) ny = (ny << 1) | (uint32_t)B[t];
    for (int r = n - h + 1; r <= n; r++) ndx = (ndx << 1) | (uint32_t)(V[r][mm] - V[r - 1][mm]);
    for (int q = mm - h + 1; q <= mm; q++) ndy = (ndy << 1) | (uint32_t)(V[n][q] - V[n][q - 1]);
    if (reduce && (nx & (1u << (h - 1)))) { nx ^= m; ny ^= m; }
    *out = ((uint64_t)nx << (3 * h)) | ((uint64_t)ny << (2 * h)) | ((uint64_t)ndx << h) | ndy;
    return 1;
}
static uint64_t compl_xy(uint64_t k) { uint32_t m = (1u << h) - 1; return k ^ (((uint64_t)m << (3 * h)) | ((uint64_t)m << (2 * h))); }
static uint64_t reduce_key(uint64_t k) { return ((k >> (4 * h - 1)) & 1) ? compl_xy(k) : k; }

#define MAXL 400
static int Ls[MAXL + 1][MAXL + 1];
static long npairs, nsteps, nfail, nsymfail, ntrunc;

static void gen(int *S, int len, int mode) {
    double p = 0.5;
    switch (mode) {
    case 0: for (int i = 0; i < len; i++) S[i] = xr() & 1; break;                       /* uniform */
    case 1: p = 0.08; /* fallthrough */
    case 2: if (mode == 2) p = 0.2;                                                      /* long runs */
        { int b = xr() & 1; for (int i = 0; i < len; i++) { if ((xr() % 10000) < p * 10000) b ^= 1; S[i] = b; } } break;
    case 3: { int per = 1 + xr() % 6; for (int i = 0; i < len; i++) S[i] = ((i / per) & 1); } break; /* periodic blocks */
    case 4: { int b = 0; int i = 0; while (i < len) { int r = 1 + xr() % (3 * h); for (int t = 0; t < r && i < len; t++) S[i++] = b; b ^= 1; } } break;
    default: for (int i = 0; i < len; i++) S[i] = (xr() % 100) < 85; break;           /* biased */
    }
}

int main(int argc, char **argv) {
    h = atoi(argv[1]); long pairs = atol(argv[2]); int maxlen = atoi(argv[3]); rng ^= (uint64_t)atoll(argv[4]) * 0x9E3779B97F4A7C15ULL;
    uint32_t m = (1u << h) - 1;
    uint64_t start = ((uint64_t)m << h) | m;
    int A[MAXL], B[MAXL];
    for (long t = 0; t < pairs; t++) {
        int la = 1 + xr() % maxlen, lb = 1 + xr() % maxlen;
        if (t % 7 == 0) lb = 1 + xr() % 4;                 /* very unequal lengths */
        int mode = xr() % 6, mode2 = (xr() % 3) ? mode : (int)(xr() % 6);
        gen(A, la, mode); gen(B, lb, mode2);
        for (int i = 0; i <= la; i++) Ls[i][0] = 0;
        for (int j = 0; j <= lb; j++) Ls[0][j] = 0;
        for (int i = 1; i <= la; i++) for (int j = 1; j <= lb; j++)
            Ls[i][j] = (A[i - 1] == B[j - 1]) ? Ls[i - 1][j - 1] + 1 : (Ls[i - 1][j] > Ls[i][j - 1] ? Ls[i - 1][j] : Ls[i][j - 1]);
        int mi[MAXL], mj[MAXL], k = 0, i = la, j = lb;
        while (i > 0 && j > 0) {
            if (Ls[i][j] == Ls[i - 1][j]) i--; else if (Ls[i][j] == Ls[i][j - 1]) j--; else { mi[k] = i; mj[k] = j; k++; i--; j--; }
        }
        if (k != Ls[la][lb]) { printf("lcs count mismatch\n"); return 1; }
        /* tail must be null */
        { int pa = k ? mi[0] : 0, pb = k ? mj[0] : 0, a0 = 0, a1 = 0, b0 = 0, b1 = 0;
          for (int q = pa; q < la; q++) { if (A[q]) a1 = 1; else a0 = 1; }
          for (int q = pb; q < lb; q++) { if (B[q]) b1 = 1; else b0 = 1; }
          if ((a0 && b0) || (a1 && b1)) { printf("tail not null\n"); nfail++; } }
        uint64_t qa = start, qb = start; int flip = 0;  /* qb is the representative; flip=1 means actual = complement(qb) */
        int pa = 0, pb = 0;
        for (int r = k - 1; r >= 0; r--) {
            int ea = mi[r], eb = mj[r];                    /* piece = A[pa..ea-1], B[pb..eb-1] */
            int ta = ea - pa, tb = eb - pb, c = A[ea - 1];
            if (B[eb - 1] != c || (ta > 1 && tb > 1)) { printf("not a match pair / both padded\n"); nfail++; break; }
            for (int q = pa; q < ea - 1; q++) if (A[q] == c) { printf("not minimal (top)\n"); nfail++; }
            for (int q = pb; q < eb - 1; q++) if (B[q] == c) { printf("not minimal (bottom)\n"); nfail++; }
            int l = (ta > 1) ? ta - 1 : -(tb - 1);
            if (l > h - 1) { l = h - 1; ntrunc++; } if (l < -(h - 1)) { l = -(h - 1); ntrunc++; }
            uint64_t na, nb;
            if (!my_delta(qa, c, l, &na, 0)) { nfail++; printf("UNREDUCED REJECTED canonical list (h=%d)\n", h); break; }
            int cc = flip ? 1 - c : c;
            if (!my_delta(qb, cc, l, &nb, 0)) { nfail++; printf("REDUCED REJECTED\n"); break; }
            /* symmetry: delta(compl q, compl p) == compl delta(q, p) */
            { uint64_t s1; int ok1 = my_delta(compl_xy(qa), 1 - c, l, &s1, 0); if (!ok1 || s1 != compl_xy(na)) nsymfail++; }
            int nflip = (nb >> (4 * h - 1)) & 1; uint64_t nbr = nflip ? compl_xy(nb) : nb;
            flip ^= nflip; qb = nbr;
            if (reduce_key(na) != qb) { nsymfail++; }
            qa = na; pa = ea; pb = eb; nsteps++;
        }
        npairs++;
    }
    printf("h=%d pairs=%ld steps=%ld truncated-inputs=%ld failures=%ld symmetry/reduction-mismatches=%ld\n", h, npairs, nsteps, ntrunc, nfail, nsymfail);
    return 0;
}
