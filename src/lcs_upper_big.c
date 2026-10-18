/*
 * Large-h version of lcs_upper.c: same automaton M_h (Lueker 2009, Fig. 7),
 * same state numbering (breadth-first discovery order), same file format,
 * but
 *   - transitions of a batch of states are computed in parallel, then ids are
 *     assigned sequentially in (state, c, l) order, so the table is identical
 *     to the one lcs_upper.c writes;
 *   - the hash table stores only 32-bit ids and compares keys through the
 *     state array (about 12 bytes per state in RAM);
 *   - the eigen phase streams the transition table from disk in chunks.
 *
 * Usage:
 *   lcs_upper_big build <h> <dir> [log2 hash capacity]
 *   lcs_upper_big eigen <h> <dir> <start-vector or -> <z> [<z> ...]
 *   lcs_upper_big project <h> <dir>    (write a start vector for h from eigen_h<h-1>.bin)
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <dispatch/dispatch.h>

static int H;
static uint32_t MASK;
static int NT, LL;

static double now(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

static inline uint64_t pack(uint32_t x, uint32_t y, uint32_t dx, uint32_t dy) {
    return ((uint64_t)x << (3 * H)) | ((uint64_t)y << (2 * H)) |
           ((uint64_t)dx << H) | (uint64_t)dy;
}

#define REJECT UINT64_MAX

static uint64_t delta(uint64_t key, int c, int l) {
    uint32_t x = (uint32_t)(key >> (3 * H)) & MASK;
    uint32_t y = (uint32_t)(key >> (2 * H)) & MASK;
    uint32_t dx = (uint32_t)(key >> H) & MASK;
    uint32_t dy = (uint32_t)key & MASK;
    int a[40], b[40];
    int T[40][40];
    int aLen = H, bLen = H;
    for (int i = H; i >= 1; i--) { a[i] = x & 1; x >>= 1; }
    for (int i = H; i >= 1; i--) { b[i] = y & 1; y >>= 1; }
    T[H][H] = H;
    for (int ii = H; ii >= 1; ii--) { T[ii - 1][H] = T[ii][H] - (int)(dx & 1); dx >>= 1; }
    for (int jj = H; jj >= 1; jj--) { T[H][jj - 1] = T[H][jj] - (int)(dy & 1); dy >>= 1; }
    if (l > 0) { for (int t = 0; t < l; t++) a[++aLen] = 1 - c; }
    else       { for (int t = 0; t < -l; t++) b[++bLen] = 1 - c; }
    a[++aLen] = c;
    b[++bLen] = c;
    for (int j = H + 1; j <= bLen; j++) T[0][j] = T[0][j - 1];
    for (int i = H + 1; i <= aLen; i++) T[i][0] = T[i - 1][0];
    for (int i = 1; i <= aLen; i++)
        for (int j = (i > H) ? 1 : H + 1; j <= bLen; j++)
            T[i][j] = (a[i] == b[j]) ? T[i - 1][j - 1] + 1
                    : (T[i][j - 1] > T[i - 1][j] ? T[i][j - 1] : T[i - 1][j]);
    int i = aLen, j = bLen, matches = 0, ok = 0;
    while (i >= H && j >= H) {
        if (T[i][j] == T[i - 1][j]) i--;
        else if (T[i][j] == T[i][j - 1]) j--;
        else { i--; j--; matches++; }
        if (i == H && j == H) { ok = 1; break; }
    }
    if (!ok) return REJECT;
    if (matches != 1) { fprintf(stderr, "unexpected match count\n"); exit(3); }
    uint32_t nx = 0, ny = 0, ndx = 0, ndy = 0;
    for (int t = aLen - H + 1; t <= aLen; t++) nx = (nx << 1) | (uint32_t)a[t];
    for (int t = bLen - H + 1; t <= bLen; t++) ny = (ny << 1) | (uint32_t)b[t];
    for (int t = aLen - H + 1; t <= aLen; t++) {
        int d = T[t][bLen] - T[t - 1][bLen];
        if (d != 0 && d != 1) { fprintf(stderr, "bad difference\n"); exit(3); }
        ndx = (ndx << 1) | (uint32_t)d;
    }
    for (int t = bLen - H + 1; t <= bLen; t++) {
        int d = T[aLen][t] - T[aLen][t - 1];
        if (d != 0 && d != 1) { fprintf(stderr, "bad difference\n"); exit(3); }
        ndy = (ndy << 1) | (uint32_t)d;
    }
    if (nx >> (H - 1)) { nx = MASK - nx; ny = MASK - ny; }
    return pack(nx, ny, ndx, ndy);
}

static inline uint64_t mix(uint64_t k) {
    k ^= k >> 33; k *= 0xff51afd7ed558ccdULL;
    k ^= k >> 33; k *= 0xc4ceb9fe1a85ec53ULL;
    k ^= k >> 33; return k;
}

static uint32_t *htab;          /* 0 = empty, else state id */
static uint64_t hmask;
static uint64_t *states;         /* states[id] */
static uint64_t nstates, scap;

static uint32_t find_or_add(uint64_t key) {
    uint64_t p = mix(key) & hmask;
    for (;;) {
        uint32_t id = htab[p];
        if (id == 0) {
            nstates++;
            if (nstates >= scap) { fprintf(stderr, "state array full\n"); exit(2); }
            if (nstates > 0xfffffff0ULL) { fprintf(stderr, "too many states\n"); exit(2); }
            states[nstates] = key;
            htab[p] = (uint32_t)nstates;
            return (uint32_t)nstates;
        }
        if (states[id] == key) return id;
        p = (p + 1) & hmask;
    }
}

static void build(const char *dir, int lg) {
    uint64_t hcap = 1ULL << lg;
    hmask = hcap - 1;
    htab = calloc(hcap, sizeof(uint32_t));
    scap = hcap * 7 / 10;          /* keep load factor <= 0.7 */
    states = malloc(scap * sizeof(uint64_t));
    if (!htab || !states) { fprintf(stderr, "out of memory\n"); exit(2); }
    char path[1024];
    snprintf(path, sizeof path, "%s/trans_h%d.bin", dir, H);
    FILE *f = fopen(path, "wb");
    if (!f) { perror(path); exit(2); }
    int32_t hdr[2] = {H, 0};
    fwrite(hdr, sizeof hdr, 1, f);
    nstates = 0;
    find_or_add(pack(0, 0, MASK, MASK));
    const uint64_t B = 1 << 21;
    uint64_t *succ = malloc(B * NT * sizeof(uint64_t));
    uint32_t *rows = malloc(B * NT * sizeof(uint32_t));
    double t0 = now();
    uint64_t done = 0;
    while (done < nstates) {
        uint64_t lo = done + 1, hi = nstates;          /* states known so far */
        if (hi - lo + 1 > B) hi = lo + B - 1;
        uint64_t cnt = hi - lo + 1;
        dispatch_apply(cnt, dispatch_get_global_queue(QOS_CLASS_USER_INITIATED, 0), ^(size_t k) {
            uint64_t key = states[lo + k];
            uint64_t *o = succ + k * NT;
            int m = 0;
            for (int c = 0; c < 2; c++)
                for (int l = -(H - 1); l <= H - 1; l++)
                    o[m++] = delta(key, c, l);
        });
        for (uint64_t k = 0; k < cnt * NT; k++)
            rows[k] = succ[k] == REJECT ? 0 : find_or_add(succ[k]);
        fwrite(rows, sizeof(uint32_t), cnt * NT, f);
        done = hi;
        fprintf(stderr, "  processed %llu, discovered %llu (%.0fs)\n",
                (unsigned long long)done, (unsigned long long)nstates, now() - t0);
    }
    fwrite(&states[1], sizeof(uint64_t), nstates, f);
    fseek(f, 0, SEEK_SET);
    hdr[1] = (int32_t)nstates;
    fwrite(hdr, sizeof hdr, 1, f);
    fclose(f);
    printf("h=%d accepting states=%llu build time %.1fs -> %s\n", H,
           (unsigned long long)nstates, now() - t0, path);
}

/* ---------------- streaming eigen phase -------------------------------- */
static uint64_t N;
static double *V, *VN;
static char tpath[1024];

static void weights(double z, double *w) {
    for (int l = -(H - 1); l <= H - 1; l++) {
        int al = l < 0 ? -l : l;
        w[l + H - 1] = (al == H - 1) ? pow(z, H + 1) / (1.0 - z) : pow(z, al + 2);
    }
}

static uint32_t *chunkbuf;
static const uint64_t CH = 1 << 24;          /* states per streamed chunk */

static void iterate(const double *w, double *rmin, double *rmax) {
    FILE *f = fopen(tpath, "rb");
    fseek(f, 8, SEEK_SET);
    double mn = INFINITY, mx = 0;
    for (uint64_t lo = 1; lo <= N; lo += CH) {
        uint64_t cnt = (N - lo + 1 < CH) ? N - lo + 1 : CH;
        if (fread(chunkbuf, sizeof(uint32_t), cnt * NT, f) != cnt * NT) { fprintf(stderr, "short read\n"); exit(2); }
        const int nch = 256;
        double *mins = malloc(nch * sizeof(double)), *maxs = malloc(nch * sizeof(double));
        uint64_t per = (cnt + nch - 1) / nch;
        dispatch_apply(nch, dispatch_get_global_queue(QOS_CLASS_USER_INITIATED, 0), ^(size_t ci) {
            double a = INFINITY, b = 0;
            uint64_t s0 = ci * per, s1 = s0 + per;
            if (s1 > cnt) s1 = cnt;
            for (uint64_t k = s0; k < s1; k++) {
                const uint32_t *r = chunkbuf + k * NT;
                double sum = 0;
                for (int q = 0; q < LL; q++) sum += w[q] * V[r[q]];
                for (int q = 0; q < LL; q++) sum += w[q] * V[r[LL + q]];
                uint64_t s = lo + k;
                VN[s] = sum;
                double ratio = sum / V[s];
                if (ratio < a) a = ratio;
                if (ratio > b) b = ratio;
            }
            mins[ci] = a; maxs[ci] = b;
        });
        for (int i = 0; i < nch; i++) { if (mins[i] < mn) mn = mins[i]; if (maxs[i] > mx) mx = maxs[i]; }
        free(mins); free(maxs);
    }
    fclose(f);
    *rmin = mn; *rmax = mx;
}

static void dump_vector(const char *dir, double z) {
    char path[1024];
    snprintf(path, sizeof path, "%s/eigen_h%d.bin", dir, H);
    FILE *f = fopen(path, "wb");
    int32_t hh = H; uint64_t n = N;
    fwrite(&hh, sizeof hh, 1, f);
    fwrite(&n, sizeof n, 1, f);
    fwrite(&z, sizeof z, 1, f);
    fwrite(V + 1, sizeof(double), N, f);
    fclose(f);
    fprintf(stderr, "  wrote %s\n", path);
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: see source\n"); return 1; }
    H = atoi(argv[2]);
    MASK = (1u << H) - 1;
    NT = 2 * (2 * H - 1);
    LL = 2 * H - 1;
    const char *dir = argv[3];
    if (!strcmp(argv[1], "build")) {
        int lg = argc > 4 ? atoi(argv[4]) : 30;
        build(dir, lg);
        return 0;
    }
    snprintf(tpath, sizeof tpath, "%s/trans_h%d.bin", dir, H);
    FILE *f = fopen(tpath, "rb");
    if (!f) { perror(tpath); return 2; }
    int32_t hdr[2];
    if (fread(hdr, sizeof hdr, 1, f) != 1 || hdr[0] != H) return 2;
    N = (uint64_t)hdr[1];
    V = malloc((N + 1) * sizeof(double));
    VN = malloc((N + 1) * sizeof(double));
    V[0] = VN[0] = 0;
    if (!strcmp(argv[1], "project")) {
        /* start vector for h from the (h-1) vector: each h-state is mapped to the
         * (h-1)-state obtained by dropping the oldest character of x and y and the
         * oldest difference of dx and dy (then complement-reduced); if that state
         * is not in the (h-1) table the entry is 1. */
        uint64_t *keys = malloc(N * sizeof(uint64_t));
        fseek(f, 8 + (uint64_t)4 * NT * N, SEEK_SET);
        if (fread(keys, sizeof(uint64_t), N, f) != N) return 2;
        fclose(f);
        int hp = H - 1, NTp = 2 * (2 * hp - 1);
        char pp[1024];
        snprintf(pp, sizeof pp, "%s/trans_h%d.bin", dir, hp);
        FILE *g = fopen(pp, "rb");
        int32_t hd2[2];
        if (fread(hd2, sizeof hd2, 1, g) != 1) return 2;
        uint64_t Np = (uint64_t)hd2[1];
        uint64_t *kp = malloc(Np * sizeof(uint64_t));
        fseek(g, 8 + (uint64_t)4 * NTp * Np, SEEK_SET);
        if (fread(kp, sizeof(uint64_t), Np, g) != Np) return 2;
        fclose(g);
        snprintf(pp, sizeof pp, "%s/eigen_h%d.bin", dir, hp);
        g = fopen(pp, "rb");
        int32_t eh; uint64_t en; double ez;
        if (fread(&eh, 4, 1, g) != 1 || fread(&en, 8, 1, g) != 1 || fread(&ez, 8, 1, g) != 1) return 2;
        double *vp = malloc(Np * sizeof(double));
        if (fread(vp, sizeof(double), Np, g) != Np) return 2;
        fclose(g);
        /* hash the (h-1) keys */
        uint64_t cap = 1; while (cap < 2 * Np) cap <<= 1;
        uint64_t *hk = malloc(cap * 8); uint32_t *hv = malloc(cap * 4);
        memset(hk, 0xff, cap * 8);
        for (uint64_t i = 0; i < Np; i++) {
            uint64_t p = mix(kp[i]) & (cap - 1);
            while (hk[p] != UINT64_MAX) p = (p + 1) & (cap - 1);
            hk[p] = kp[i]; hv[p] = (uint32_t)i;
        }
        uint32_t mp = (1u << hp) - 1;
        uint64_t hit = 0;
        for (uint64_t s = 0; s < N; s++) {
            uint64_t k = keys[s];
            uint32_t x = (uint32_t)(k >> (3 * H)) & mp, y = (uint32_t)(k >> (2 * H)) & mp;
            uint32_t dx = (uint32_t)(k >> H) & mp, dy = (uint32_t)k & mp;
            if (x >> (hp - 1)) { x = mp - x; y = mp - y; }
            uint64_t kk = ((uint64_t)x << (3 * hp)) | ((uint64_t)y << (2 * hp)) | ((uint64_t)dx << hp) | dy;
            uint64_t p = mix(kk) & (cap - 1);
            double val = 1.0;
            while (hk[p] != UINT64_MAX) { if (hk[p] == kk) { val = vp[hv[p]]; hit++; break; } p = (p + 1) & (cap - 1); }
            V[s + 1] = val;
        }
        fprintf(stderr, "projected %llu of %llu states onto h-1 states\n", (unsigned long long)hit, (unsigned long long)N);
        char op[1024];
        snprintf(op, sizeof op, "%s/start_h%d.bin", dir, H);
        FILE *o = fopen(op, "wb");
        fwrite(V + 1, sizeof(double), N, o);
        fclose(o);
        return 0;
    }
    fclose(f);
    chunkbuf = malloc(CH * NT * sizeof(uint32_t));
    if (strcmp(argv[4], "-")) {
        FILE *g = fopen(argv[4], "rb");
        if (fread(V + 1, sizeof(double), N, g) != N) { fprintf(stderr, "bad start vector\n"); return 2; }
        fclose(g);
    } else {
        for (uint64_t s = 1; s <= N; s++) V[s] = 1.0;
    }
    printf("h=%d states=%llu\n", H, (unsigned long long)N);
    for (int ai = 5; ai < argc; ai++) {
        double z = atof(argv[ai]);
        double w[64];
        weights(z, w);
        double rmin, rmax, t0 = now();
        int it = 0;
        for (;;) {
            iterate(w, &rmin, &rmax);
            it++;
            double mx = 0;
            for (uint64_t s = 1; s <= N; s++) if (VN[s] > mx) mx = VN[s];
            for (uint64_t s = 1; s <= N; s++) V[s] = VN[s] / mx;
            fprintf(stderr, "    z=%.6f it=%d ratios [%.12f, %.12f] bound<=%.7f (%.0fs)\n", z, it, rmin, rmax,
                    2.0 * log(2.0 * z) / log(rmax), now() - t0);
            if (rmax / rmin <= 1 + 1e-10 || it >= 5000) break;
        }
        printf("z=%.6f lambda=%.15f bound=%.9f iterations=%d\n", z, rmax, 2.0 * log(2.0 * z) / log(rmax), it);
        fflush(stdout);
        dump_vector(dir, z);
    }
    return 0;
}
