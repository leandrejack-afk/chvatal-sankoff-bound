/*
 * Upper bounds on the binary Chvatal-Sankoff constant gamma_2 by Lueker's
 * method (J. ACM 56(3), 2009, Section 3.2): the history-h automaton M_h that
 * reads minimal match pairs T_{c,l} and rejects lists that are provably not
 * canonical.  For z in (0, 1/2) let A_z be its weighted transition matrix and
 * lambda(z) its Perron root; then gamma_2 <= 2 log(2z) / log lambda(z).
 *
 * This program
 *   1. builds the reachable accepting states of M_h by breadth-first search,
 *      with the complement symmetry <x,y,dx,dy> ~ <~x,~y,dx,dy>;
 *   2. writes the transition table to <dir>/trans_h<h>.bin;
 *   3. for each requested z, runs power iteration and reports
 *      lambda(z) = max_s (A_z v)_s / v_s and the bound; writes the vector.
 *
 * Usage:
 *   lcs_upper build <h> <dir>
 *   lcs_upper eigen <h> <dir> <z> [<z> ...]
 *   lcs_upper scan  <h> <dir> <zlo> <zhi>        (golden-section on z)
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
static int NT;          /* transitions per state: 2 * (2H - 1) */

static double now(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

static inline uint64_t pack(uint32_t x, uint32_t y, uint32_t dx, uint32_t dy) {
    return ((uint64_t)x << (3 * H)) | ((uint64_t)y << (2 * H)) |
           ((uint64_t)dx << H) | (uint64_t)dy;
}

/* delta(q, T_{c,l}); returns 0 on rejection, else 1 with *out = reduced state.
 * Transcribes Figure 7 of the paper (and Lueker's DPTable/State classes). */
static int delta(uint64_t key, int c, int l, uint64_t *out) {
    uint32_t x = (uint32_t)(key >> (3 * H)) & MASK;
    uint32_t y = (uint32_t)(key >> (2 * H)) & MASK;
    uint32_t dx = (uint32_t)(key >> H) & MASK;
    uint32_t dy = (uint32_t)key & MASK;
    int a[40], b[40];
    int T[40][40];
    int aLen = H, bLen = H;
    for (int i = H; i >= 1; i--) { a[i] = x & 1; x >>= 1; }
    for (int i = H; i >= 1; i--) { b[i] = y & 1; y >>= 1; }
#ifdef DEBUG
    for (int i = 0; i < 40; i++) for (int j = 0; j < 40; j++) T[i][j] = -1000000;
#endif
    T[H][H] = H;
    for (int ii = H; ii >= 1; ii--) { T[ii - 1][H] = T[ii][H] - (int)(dx & 1); dx >>= 1; }
    for (int jj = H; jj >= 1; jj--) { T[H][jj - 1] = T[H][jj] - (int)(dy & 1); dy >>= 1; }
    if (l > 0) { for (int t = 0; t < l; t++) a[++aLen] = 1 - c; }
    else       { for (int t = 0; t < -l; t++) b[++bLen] = 1 - c; }
    a[++aLen] = c;
    b[++bLen] = c;
    for (int j = H + 1; j <= bLen; j++) T[0][j] = T[0][j - 1];
    for (int i = H + 1; i <= aLen; i++) T[i][0] = T[i - 1][0];
    for (int i = 1; i <= aLen; i++) {
        for (int j = (i > H) ? 1 : H + 1; j <= bLen; j++) {
#ifdef DEBUG
            if (T[i - 1][j - 1] < -1000 || T[i - 1][j] < -1000 || T[i][j - 1] < -1000) {
                fprintf(stderr, "read of unknown cell at %d %d\n", i, j); exit(3);
            }
#endif
            if (a[i] == b[j]) T[i][j] = T[i - 1][j - 1] + 1;
            else T[i][j] = T[i][j - 1] > T[i - 1][j] ? T[i][j - 1] : T[i - 1][j];
        }
    }
    int i = aLen, j = bLen, matches = 0;
    int ok = 0;
    while (i >= H && j >= H) {
        if (T[i][j] == T[i - 1][j]) i--;
        else if (T[i][j] == T[i][j - 1]) j--;
        else {
#ifdef DEBUG
            if (T[i][j] != T[i - 1][j - 1] + 1) { fprintf(stderr, "inconsistent table\n"); exit(3); }
#endif
            i--; j--; matches++;
        }
        if (i == H && j == H) { ok = 1; break; }
    }
    if (!ok) return 0;
    if (matches != 1) { fprintf(stderr, "unexpected match count %d\n", matches); exit(3); }
    uint32_t nx = 0, ny = 0, ndx = 0, ndy = 0;
    for (int t = aLen - H + 1; t <= aLen; t++) nx = (nx << 1) | (uint32_t)a[t];
    for (int t = bLen - H + 1; t <= bLen; t++) ny = (ny << 1) | (uint32_t)b[t];
    for (int t = aLen - H + 1; t <= aLen; t++) {
        int d = T[t][bLen] - T[t - 1][bLen];
        if (d != 0 && d != 1) { fprintf(stderr, "vertical difference %d\n", d); exit(3); }
        ndx = (ndx << 1) | (uint32_t)d;
    }
    for (int t = bLen - H + 1; t <= bLen; t++) {
        int d = T[aLen][t] - T[aLen][t - 1];
        if (d != 0 && d != 1) { fprintf(stderr, "horizontal difference %d\n", d); exit(3); }
        ndy = (ndy << 1) | (uint32_t)d;
    }
    if (nx >> (H - 1)) { nx = MASK - nx; ny = MASK - ny; }
    *out = pack(nx, ny, ndx, ndy);
    return 1;
}

/* ---------------- open-addressing hash: key -> state id (1-based) ---------- */
static uint64_t *hkeys;
static uint32_t *hvals;
static uint64_t hcap, hmask;

static inline uint64_t mix(uint64_t k) {
    k ^= k >> 33; k *= 0xff51afd7ed558ccdULL;
    k ^= k >> 33; k *= 0xc4ceb9fe1a85ec53ULL;
    k ^= k >> 33; return k;
}
#define EMPTY UINT64_MAX

static uint64_t *states;     /* states[id] = key, id from 1 */
static uint64_t nstates, scap;

static uint32_t find_or_add(uint64_t key) {
    uint64_t p = mix(key) & hmask;
    for (;;) {
        if (hkeys[p] == key) return hvals[p];
        if (hkeys[p] == EMPTY) {
            if (nstates + 1 >= scap) {
                scap *= 2;
                states = realloc(states, scap * sizeof(uint64_t));
                if (!states) { fprintf(stderr, "out of memory (states)\n"); exit(2); }
            }
            nstates++;
            if (nstates >= hcap / 2) { fprintf(stderr, "hash table too small\n"); exit(2); }
            states[nstates] = key;
            hkeys[p] = key;
            hvals[p] = (uint32_t)nstates;
            return (uint32_t)nstates;
        }
        p = (p + 1) & hmask;
    }
}

static void build(const char *dir, int log2cap) {
    hcap = 1ULL << log2cap; hmask = hcap - 1;
    hkeys = malloc(hcap * sizeof(uint64_t));
    hvals = malloc(hcap * sizeof(uint32_t));
    if (!hkeys || !hvals) { fprintf(stderr, "out of memory (hash)\n"); exit(2); }
    memset(hkeys, 0xff, hcap * sizeof(uint64_t));
    scap = 1 << 20; nstates = 0;
    states = malloc(scap * sizeof(uint64_t));
    char path[1024];
    snprintf(path, sizeof path, "%s/trans_h%d.bin", dir, H);
    FILE *f = fopen(path, "wb");
    if (!f) { perror(path); exit(2); }
    int32_t hdr[2] = {H, 0};
    fwrite(hdr, sizeof hdr, 1, f);
    find_or_add(pack(0, 0, MASK, MASK));   /* start state q0 = <0^h,0^h,1^h,1^h> */
    uint32_t *row = malloc(NT * sizeof(uint32_t));
    double t0 = now();
    for (uint64_t s = 1; s <= nstates; s++) {
        uint64_t key = states[s];
        int k = 0;
        for (int c = 0; c < 2; c++) {
            for (int l = -(H - 1); l <= H - 1; l++) {
                uint64_t nk;
                row[k++] = delta(key, c, l, &nk) ? find_or_add(nk) : 0;
            }
        }
        fwrite(row, sizeof(uint32_t), NT, f);
        if ((s & ((1 << 22) - 1)) == 0)
            fprintf(stderr, "  processed %llu, discovered %llu (%.0fs)\n",
                    (unsigned long long)s, (unsigned long long)nstates, now() - t0);
    }
    /* append the state keys so the table is self-describing */
    fwrite(&states[1], sizeof(uint64_t), nstates, f);
    fseek(f, 0, SEEK_SET);
    hdr[1] = (int32_t)nstates;
    fwrite(hdr, sizeof hdr, 1, f);
    fclose(f);
    printf("h=%d accepting states=%llu build time %.1fs -> %s\n", H,
           (unsigned long long)nstates, now() - t0, path);
}

/* ---------------- eigenvalue part -------------------------------------- */
static uint32_t *trans;   /* trans[(s-1)*NT + k] */
static uint64_t N;

static void load(const char *dir) {
    char path[1024];
    snprintf(path, sizeof path, "%s/trans_h%d.bin", dir, H);
    FILE *f = fopen(path, "rb");
    if (!f) { perror(path); exit(2); }
    int32_t hdr[2];
    if (fread(hdr, sizeof hdr, 1, f) != 1 || hdr[0] != H) { fprintf(stderr, "bad header\n"); exit(2); }
    N = (uint64_t)hdr[1];
    trans = malloc(N * NT * sizeof(uint32_t));
    if (!trans) { fprintf(stderr, "out of memory (trans)\n"); exit(2); }
    if (fread(trans, sizeof(uint32_t), N * NT, f) != N * NT) { fprintf(stderr, "short read\n"); exit(2); }
    fclose(f);
}

static void weights(double z, double *w) {
    /* w[k] for k = l + (H-1), l in [-(H-1), H-1]:
     * z^{|l|+2} if |l| < H-1, z^{H+1}/(1-z) if |l| = H-1 */
    for (int l = -(H - 1); l <= H - 1; l++) {
        int al = l < 0 ? -l : l;
        double v = pow(z, al + 2);
        if (al == H - 1) v = pow(z, H + 1) / (1.0 - z);
        w[l + H - 1] = v;
    }
}

/* one multiplication vnew = A_z v; returns min and max ratios */
static void iterate(const double *w, const double *v, double *vn, double *rmin, double *rmax) {
    const int chunks = 256;
    double *mins = malloc(chunks * sizeof(double)), *maxs = malloc(chunks * sizeof(double));
    uint64_t per = (N + chunks - 1) / chunks;
    int L = 2 * H - 1;
    dispatch_apply(chunks, dispatch_get_global_queue(QOS_CLASS_USER_INITIATED, 0), ^(size_t ci) {
        double mn = INFINITY, mx = 0;
        uint64_t lo = 1 + ci * per, hi = lo + per;
        if (hi > N + 1) hi = N + 1;
        for (uint64_t s = lo; s < hi; s++) {
            const uint32_t *r = trans + (s - 1) * NT;
            double sum = 0;
            for (int k = 0; k < L; k++) sum += w[k] * v[r[k]];
            for (int k = 0; k < L; k++) sum += w[k] * v[r[L + k]];
            vn[s] = sum;
            double q = sum / v[s];
            if (q < mn) mn = q;
            if (q > mx) mx = q;
        }
        mins[ci] = mn; maxs[ci] = mx;
    });
    double mn = INFINITY, mx = 0;
    for (int i = 0; i < chunks; i++) { if (mins[i] < mn) mn = mins[i]; if (maxs[i] > mx) mx = maxs[i]; }
    free(mins); free(maxs);
    *rmin = mn; *rmax = mx;
}

static double *V, *VN;

static double lambda_of(double z, double tol, int verbose) {
    double w[64];
    weights(z, w);
    double rmin, rmax;
    int it = 0;
    double t0 = now();
    for (;;) {
        iterate(w, V, VN, &rmin, &rmax);
        it++;
        /* normalise by the maximum entry */
        double mx = 0;
        for (uint64_t s = 1; s <= N; s++) if (VN[s] > mx) mx = VN[s];
        for (uint64_t s = 1; s <= N; s++) V[s] = VN[s] / mx;
        if (verbose && (it % 25 == 0))
            fprintf(stderr, "    z=%.6f it=%d ratios [%.12f, %.12f] (%.0fs)\n", z, it, rmin, rmax, now() - t0);
        if (rmax / rmin <= tol) break;
        if (it > 20000) { fprintf(stderr, "no convergence\n"); break; }
    }
    if (verbose)
        fprintf(stderr, "  z=%.6f converged in %d iterations: lambda in [%.15f, %.15f]\n", z, it, rmin, rmax);
    return rmax;
}

static double bound_of(double z, double lam) { return 2.0 * log(2.0 * z) / log(lam); }

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
    if (H < 2 || H > 15) { fprintf(stderr, "h out of range\n"); return 1; }
    MASK = (1u << H) - 1;
    NT = 2 * (2 * H - 1);
    const char *dir = argv[3];
    if (!strcmp(argv[1], "build")) {
        int lg = 20;
        while ((1ULL << lg) < (uint64_t)(4.0 * pow(8.3, H - 1))) lg++;
        if (argc > 4) lg = atoi(argv[4]);
        build(dir, lg);
        return 0;
    }
    load(dir);
    V = malloc((N + 1) * sizeof(double));
    VN = malloc((N + 1) * sizeof(double));
    V[0] = VN[0] = 0;
    for (uint64_t s = 1; s <= N; s++) V[s] = 1.0;
    printf("h=%d states=%llu\n", H, (unsigned long long)N);
    if (!strcmp(argv[1], "eigen")) {
        for (int i = 4; i < argc; i++) {
            double z = atof(argv[i]);
            double lam = lambda_of(z, 1 + 1e-10, 1);
            printf("z=%.6f lambda=%.15f bound=%.9f\n", z, lam, bound_of(z, lam));
            fflush(stdout);
        }
        dump_vector(dir, atof(argv[argc - 1]));
    } else if (!strcmp(argv[1], "scan")) {
        double a = atof(argv[4]), b = atof(argv[5]);
        const double g = (sqrt(5.0) - 1) / 2;
        double c = b - g * (b - a), d = a + g * (b - a);
        double fc = bound_of(c, lambda_of(c, 1 + 1e-10, 1));
        double fd = bound_of(d, lambda_of(d, 1 + 1e-10, 1));
        while (b - a > 2e-5) {
            if (fc < fd) { b = d; d = c; fd = fc; c = b - g * (b - a); fc = bound_of(c, lambda_of(c, 1 + 1e-10, 1)); }
            else { a = c; c = d; fc = fd; d = a + g * (b - a); fd = bound_of(d, lambda_of(d, 1 + 1e-10, 1)); }
            printf("  [%.6f, %.6f] best so far %.9f\n", a, b, fc < fd ? fc : fd);
            fflush(stdout);
        }
        double zb = fc < fd ? c : d;
        double lam = lambda_of(zb, 1 + 1e-12, 1);
        printf("best z=%.6f lambda=%.15f bound=%.9f\n", zb, lam, bound_of(zb, lam));
        dump_vector(dir, zb);
    }
    return 0;
}
