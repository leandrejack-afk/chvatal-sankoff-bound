/*
 * kary.c: Lueker's canonicity machine M_h (J. ACM 56(3), 2009, Section 3.2 and
 * Figure 7) for an alphabet of size k, 2 <= k <= 4.
 *
 * States <x, y, dx, dy>: x, y in Sigma^h are the last h characters of the two
 * strings, dx, dy in {0,1}^h the last h differences in the last column and the
 * last row of the lcs table.  Start state <0^h, 0^h, 1^h, 1^h>.
 *
 * Inputs: the minimal match pairs (a'c, b'c) with c not in a'b' and
 * letters(a'), letters(b') disjoint, truncated by mu to at most h characters
 * per side.  The weight of a truncated class is z^(|alpha|+|beta|) times the
 * generating function of the prefixes mu removed (see weight_double and
 * kweights.py); it depends only on the "weight type" (size, ta, tb, nA, nB).
 *
 * States are taken modulo relabelling of the alphabet (letters renamed in
 * order of first occurrence in x, then y).  Relabelling commutes with delta
 * and preserves weights, so the matrix on orbits is well defined and a
 * certificate for it lifts to one for A_z.  Transitions from one state with
 * the same target and the same weight type are merged into one record with a
 * multiplicity.
 *
 * Usage:
 *   kary build k h dir
 *   kary eigen k h dir z [tol]        power iteration at z, writes the vector
 *   kary scan  k h dir zlo zhi        golden-section search on z, writes the vector
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <dispatch/dispatch.h>

#define HMAX 10
#define NMAX (2 * HMAX + 2)

static int K, H;
static uint64_t M1, M2;            /* masks of h and 2h bits */

typedef struct { uint8_t c, alen, blen, wt; uint8_t a[HMAX], b[HMAX]; } Input;
typedef struct { int32_t size, ta, tb, na, nb; } WType;
typedef struct { uint32_t tgt; uint16_t wt; uint16_t cnt; } Rec;

static Input *inputs;
static int NI;
static WType wtypes[256];
static int NW;

static double now(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

static inline uint64_t pack(uint64_t x, uint64_t y, uint64_t dx, uint64_t dy) {
    return (x << (4 * H)) | (y << (2 * H)) | (dx << H) | dy;
}

static int wtype_of(int size, int ta, int tb, int na, int nb) {
    if (!ta && !tb) { na = 0; nb = 0; }
    else if (ta && !tb) na = 0;          /* weight depends on nb only */
    else if (!ta && tb) nb = 0;          /* weight depends on na only */
    for (int i = 0; i < NW; i++)
        if (wtypes[i].size == size && wtypes[i].ta == ta && wtypes[i].tb == tb &&
            wtypes[i].na == na && wtypes[i].nb == nb) return i;
    if (NW == 256) { fprintf(stderr, "too many weight types\n"); exit(2); }
    wtypes[NW] = (WType){size, ta, tb, na, nb};
    return NW++;
}

/* all inputs mu(p): (a'c, b'c), |a'|, |b'| <= h-1, c not in a'b', disjoint letters */
static void make_inputs(void) {
    int cap = 1 << 12;
    inputs = malloc(cap * sizeof(Input));
    NI = 0;
    for (int c = 0; c < K; c++) {
        int O[4], no = 0;
        for (int t = 0; t < K; t++) if (t != c) O[no++] = t;
        int ns = 0;
        for (int m = 0, p = 1; m <= H - 1; m++, p *= no) ns += p;
        uint8_t (*S)[HMAX] = malloc((size_t)ns * HMAX);
        uint8_t *SL = malloc(ns), *SM = malloc(ns);
        int idx = 0;
        for (int m = 0; m <= H - 1; m++) {
            int cnt = 1;
            for (int t = 0; t < m; t++) cnt *= no;
            for (int r = 0; r < cnt; r++) {
                int q = r;
                uint8_t mask = 0;
                for (int t = m - 1; t >= 0; t--) { S[idx][t] = (uint8_t)O[q % no]; mask |= (uint8_t)(1 << S[idx][t]); q /= no; }
                SL[idx] = (uint8_t)m; SM[idx] = mask; idx++;
            }
        }
        for (int i = 0; i < ns; i++)
            for (int j = 0; j < ns; j++) {
                if (SM[i] & SM[j]) continue;
                if (NI == cap) { cap *= 2; inputs = realloc(inputs, cap * sizeof(Input)); }
                Input *in = &inputs[NI++];
                memset(in, 0, sizeof *in);
                in->c = (uint8_t)c;
                in->alen = (uint8_t)(SL[i] + 1);
                in->blen = (uint8_t)(SL[j] + 1);
                memcpy(in->a, S[i], SL[i]); in->a[SL[i]] = (uint8_t)c;
                memcpy(in->b, S[j], SL[j]); in->b[SL[j]] = (uint8_t)c;
                in->wt = (uint8_t)wtype_of(in->alen + in->blen, in->alen == H, in->blen == H,
                                           __builtin_popcount(SM[i]), __builtin_popcount(SM[j]));
            }
        free(S); free(SL); free(SM);
    }
}

/* delta(q, mu(p)) as in Figure 7, followed by relabelling; 0 = reject */
static int delta(uint64_t key, const Input *in, uint64_t *out) {
    uint64_t x = (key >> (4 * H)) & M2, y = (key >> (2 * H)) & M2;
    uint64_t dx = (key >> H) & M1, dy = key & M1;
    int a[NMAX], b[NMAX];
    int T[NMAX][NMAX];
    for (int i = H; i >= 1; i--) { a[i] = (int)(x & 3); x >>= 2; }
    for (int j = H; j >= 1; j--) { b[j] = (int)(y & 3); y >>= 2; }
    T[H][H] = H;
    for (int i = H; i >= 1; i--) { T[i - 1][H] = T[i][H] - (int)(dx & 1); dx >>= 1; }
    for (int j = H; j >= 1; j--) { T[H][j - 1] = T[H][j] - (int)(dy & 1); dy >>= 1; }
    int aLen = H + in->alen, bLen = H + in->blen;
    for (int t = 0; t < in->alen; t++) a[H + 1 + t] = in->a[t];
    for (int t = 0; t < in->blen; t++) b[H + 1 + t] = in->b[t];
    for (int j = H + 1; j <= bLen; j++) T[0][j] = T[0][j - 1];
    for (int i = H + 1; i <= aLen; i++) T[i][0] = T[i - 1][0];
    for (int i = 1; i <= aLen; i++)
        for (int j = (i > H) ? 1 : H + 1; j <= bLen; j++) {
            if (a[i] == b[j]) T[i][j] = T[i - 1][j - 1] + 1;
            else T[i][j] = T[i][j - 1] > T[i - 1][j] ? T[i][j - 1] : T[i - 1][j];
        }
    int i = aLen, j = bLen, matches = 0, ok = 0;
    while (i >= H && j >= H) {
        if (T[i][j] == T[i - 1][j]) i--;
        else if (T[i][j] == T[i][j - 1]) j--;
        else { i--; j--; matches++; }
        if (i == H && j == H) { ok = 1; break; }
    }
    if (!ok) return 0;
    if (matches != 1) { fprintf(stderr, "unexpected match count %d\n", matches); exit(3); }
    int map[4] = {-1, -1, -1, -1}, nl = 0;
    for (int t = aLen - H + 1; t <= aLen; t++) if (map[a[t]] < 0) map[a[t]] = nl++;
    for (int t = bLen - H + 1; t <= bLen; t++) if (map[b[t]] < 0) map[b[t]] = nl++;
    uint64_t nx = 0, ny = 0, ndx = 0, ndy = 0;
    for (int t = aLen - H + 1; t <= aLen; t++) nx = (nx << 2) | (uint64_t)map[a[t]];
    for (int t = bLen - H + 1; t <= bLen; t++) ny = (ny << 2) | (uint64_t)map[b[t]];
    for (int t = aLen - H + 1; t <= aLen; t++) {
        int d = T[t][bLen] - T[t - 1][bLen];
        if (d != 0 && d != 1) { fprintf(stderr, "vertical difference %d\n", d); exit(3); }
        ndx = (ndx << 1) | (uint64_t)d;
    }
    for (int t = bLen - H + 1; t <= bLen; t++) {
        int d = T[aLen][t] - T[aLen][t - 1];
        if (d != 0 && d != 1) { fprintf(stderr, "horizontal difference %d\n", d); exit(3); }
        ndy = (ndy << 1) | (uint64_t)d;
    }
    *out = pack(nx, ny, ndx, ndy);
    return 1;
}

/* ---------------- open-addressing hash: key -> state id (0-based) -------- */
static uint64_t *hkeys;
static uint32_t *hvals;
static uint64_t hcap, hmask;
static uint64_t *states;
static uint64_t nstates, scap;
#define EMPTY UINT64_MAX

static inline uint64_t mix(uint64_t k) {
    k ^= k >> 33; k *= 0xff51afd7ed558ccdULL;
    k ^= k >> 33; k *= 0xc4ceb9fe1a85ec53ULL;
    k ^= k >> 33; return k;
}

static void rehash(void) {
    uint64_t ncap = hcap * 2;
    uint64_t *nk = malloc(ncap * sizeof(uint64_t));
    uint32_t *nv = malloc(ncap * sizeof(uint32_t));
    if (!nk || !nv) { fprintf(stderr, "out of memory (rehash)\n"); exit(2); }
    memset(nk, 0xff, ncap * sizeof(uint64_t));
    uint64_t nm = ncap - 1;
    for (uint64_t p = 0; p < hcap; p++) {
        if (hkeys[p] == EMPTY) continue;
        uint64_t q = mix(hkeys[p]) & nm;
        while (nk[q] != EMPTY) q = (q + 1) & nm;
        nk[q] = hkeys[p]; nv[q] = hvals[p];
    }
    free(hkeys); free(hvals);
    hkeys = nk; hvals = nv; hcap = ncap; hmask = nm;
}

static uint32_t find_or_add(uint64_t key) {
    uint64_t p = mix(key) & hmask;
    for (;;) {
        if (hkeys[p] == key) return hvals[p];
        if (hkeys[p] == EMPTY) {
            if (nstates == scap) {
                scap *= 2;
                states = realloc(states, scap * sizeof(uint64_t));
                if (!states) { fprintf(stderr, "out of memory (states)\n"); exit(2); }
            }
            if (nstates >= UINT32_MAX - 1) { fprintf(stderr, "too many states\n"); exit(2); }
            uint32_t id = (uint32_t)nstates++;
            states[id] = key;
            hkeys[p] = key; hvals[p] = id;
            if (nstates * 2 > hcap) rehash();
            return id;
        }
        p = (p + 1) & hmask;
    }
}

typedef struct { uint64_t key; uint32_t wt; } KW;
static int cmp_kw(const void *pa, const void *pb) {
    const KW *a = pa, *b = pb;
    if (a->key != b->key) return a->key < b->key ? -1 : 1;
    return (int)a->wt - (int)b->wt;
}

static void path_of(char *buf, size_t n, const char *dir, const char *stem, const char *ext) {
    snprintf(buf, n, "%s/%s_k%d_h%d.%s", dir, stem, K, H, ext);
}

static void build(const char *dir) {
    make_inputs();
    fprintf(stderr, "k=%d h=%d inputs per state=%d weight types=%d\n", K, H, NI, NW);
    hcap = 1 << 20; hmask = hcap - 1;
    hkeys = malloc(hcap * sizeof(uint64_t));
    hvals = malloc(hcap * sizeof(uint32_t));
    memset(hkeys, 0xff, hcap * sizeof(uint64_t));
    scap = 1 << 20; nstates = 0;
    states = malloc(scap * sizeof(uint64_t));
    find_or_add(pack(0, 0, M1, M1));
    uint64_t ocap = 1 << 20, rcap = 1 << 24, nrec = 0;
    uint64_t *offs = malloc(ocap * sizeof(uint64_t));
    Rec *recs = malloc(rcap * sizeof(Rec));
    offs[0] = 0;
    if (NI > 65535) { fprintf(stderr, "too many inputs for 16-bit multiplicities\n"); exit(2); }
    uint64_t B = (1ULL << 28) / ((uint64_t)NI * sizeof(KW));
    if (B > 8192) B = 8192;
    if (B < 64) B = 64;
    KW *buf = malloc(B * NI * sizeof(KW));
    uint32_t *cntm = malloc(B * sizeof(uint32_t));
    uint16_t *mult = malloc(B * NI * sizeof(uint16_t));
    if (!buf || !cntm || !mult || !offs || !recs) { fprintf(stderr, "out of memory (buffers)\n"); exit(2); }
    uint64_t accepted = 0;
    double t0 = now(), tlast = t0;
    for (uint64_t lo = 0; lo < nstates;) {
        uint64_t hi = lo + B < nstates ? lo + B : nstates;
        const uint64_t *st = states;
        dispatch_apply(hi - lo, dispatch_get_global_queue(QOS_CLASS_UTILITY, 0), ^(size_t r) {
            uint64_t key = st[lo + r];
            KW *kb = buf + r * NI;
            uint16_t *mb = mult + r * NI;
            int m = 0;
            for (int t = 0; t < NI; t++) {
                uint64_t nk;
                if (delta(key, &inputs[t], &nk)) { kb[m].key = nk; kb[m].wt = inputs[t].wt; m++; }
            }
            qsort(kb, m, sizeof(KW), cmp_kw);
            int u = 0;
            for (int t = 0; t < m; t++) {
                if (u > 0 && kb[u - 1].key == kb[t].key && kb[u - 1].wt == kb[t].wt) { mb[u - 1]++; continue; }
                kb[u] = kb[t]; mb[u] = 1; u++;
            }
            cntm[r] = (uint32_t)u;
        });
        for (uint64_t s = lo; s < hi; s++) {
            uint64_t r = s - lo;
            const KW *kb = buf + r * NI;
            const uint16_t *mb = mult + r * NI;
            if (nrec + cntm[r] > rcap) {
                while (nrec + cntm[r] > rcap) rcap *= 2;
                recs = realloc(recs, rcap * sizeof(Rec));
                if (!recs) { fprintf(stderr, "out of memory (records)\n"); exit(2); }
            }
            for (uint32_t t = 0; t < cntm[r]; t++) {
                uint32_t id = find_or_add(kb[t].key);
                recs[nrec].tgt = id; recs[nrec].wt = (uint16_t)kb[t].wt; recs[nrec].cnt = mb[t];
                accepted += mb[t];
                nrec++;
            }
            if (s + 2 > ocap) { ocap *= 2; offs = realloc(offs, ocap * sizeof(uint64_t)); }
            offs[s + 1] = nrec;
        }
        lo = hi;
        if (now() - tlast > 30) {
            tlast = now();
            fprintf(stderr, "  processed %llu, discovered %llu, records %llu (%.0fs)\n",
                    (unsigned long long)lo, (unsigned long long)nstates, (unsigned long long)nrec, tlast - t0);
        }
    }
    char path[1024];
    path_of(path, sizeof path, dir, "ktab", "bin");
    FILE *f = fopen(path, "wb");
    if (!f) { perror(path); exit(2); }
    int32_t hk[2] = {K, H};
    uint64_t hn[2] = {nstates, nrec};
    int32_t hi2[2] = {NI, NW};
    fwrite("KLUE", 4, 1, f);
    fwrite(hk, sizeof hk, 1, f);
    fwrite(hn, sizeof hn, 1, f);
    fwrite(hi2, sizeof hi2, 1, f);
    fwrite(wtypes, sizeof(WType), NW, f);
    fwrite(offs, sizeof(uint64_t), nstates + 1, f);
    fwrite(recs, sizeof(Rec), nrec, f);
    fwrite(states, sizeof(uint64_t), nstates, f);
    fclose(f);
    char wpath[1024];
    path_of(wpath, sizeof wpath, dir, "wtypes", "txt");
    FILE *fw = fopen(wpath, "w");
    for (int i = 0; i < NW; i++)
        fprintf(fw, "%d %d %d %d %d %d\n", i, wtypes[i].size, wtypes[i].ta, wtypes[i].tb, wtypes[i].na, wtypes[i].nb);
    fclose(fw);
    printf("k=%d h=%d inputs=%d accepting states=%llu accepted transitions=%llu records=%llu build %.1fs -> %s\n",
           K, H, NI, (unsigned long long)nstates, (unsigned long long)accepted,
           (unsigned long long)nrec, now() - t0, path);
}

/* ---------------- eigenvalue part ---------------------------------------- */
static uint64_t N, R;
static uint64_t *OFF;
static Rec *RECS;

static void load(const char *dir) {
    char path[1024];
    path_of(path, sizeof path, dir, "ktab", "bin");
    FILE *f = fopen(path, "rb");
    if (!f) { perror(path); exit(2); }
    char magic[4];
    int32_t hk[2], hi2[2];
    uint64_t hn[2];
    if (fread(magic, 4, 1, f) != 1 || memcmp(magic, "KLUE", 4) || fread(hk, sizeof hk, 1, f) != 1 ||
        fread(hn, sizeof hn, 1, f) != 1 || fread(hi2, sizeof hi2, 1, f) != 1 || hk[0] != K || hk[1] != H) {
        fprintf(stderr, "bad header\n"); exit(2);
    }
    N = hn[0]; R = hn[1]; NI = hi2[0]; NW = hi2[1];
    if (fread(wtypes, sizeof(WType), NW, f) != (size_t)NW) { fprintf(stderr, "short wtypes\n"); exit(2); }
    OFF = malloc((N + 1) * sizeof(uint64_t));
    RECS = malloc(R * sizeof(Rec));
    if (!OFF || !RECS) { fprintf(stderr, "out of memory (table)\n"); exit(2); }
    if (fread(OFF, sizeof(uint64_t), N + 1, f) != N + 1 || fread(RECS, sizeof(Rec), R, f) != R) {
        fprintf(stderr, "short table\n"); exit(2);
    }
    fclose(f);
}

/* exact weight of a type, in double (the certificate uses kweights.py instead) */
static double g_strings(int a, int s, double z) {
    /* generating function of strings over a + s letters that use all s of them */
    double sum = 0, binom = 1;
    for (int j = 0; j <= s; j++) {
        if (j > 0) binom = binom * (s - j + 1) / j;
        double term = binom / (1.0 - (a + j) * z);
        sum += ((s - j) & 1) ? -term : term;
    }
    return sum;
}

static double weight_double(const WType *w, double z) {
    double v = pow(z, w->size);
    if (w->ta && !w->tb) v /= 1.0 - (K - 1 - w->nb) * z;
    else if (!w->ta && w->tb) v /= 1.0 - (K - 1 - w->na) * z;
    else if (w->ta && w->tb) {
        int f = K - 1 - w->na - w->nb;
        double tot = 0;
        for (int s = 0; s <= f; s++)
            for (int t = 0; s + t <= f; t++) {
                double mult = tgamma(f + 1) / (tgamma(s + 1) * tgamma(t + 1) * tgamma(f - s - t + 1));
                tot += mult * g_strings(w->na, s, z) * g_strings(w->nb, t, z);
            }
        v *= tot;
    }
    return v;
}

static void iterate(const double *w, const double *v, double *vn, double *rmin, double *rmax) {
    const int chunks = 512;
    double mins[512], maxs[512];
    double *pmin = mins, *pmax = maxs;
    uint64_t per = (N + chunks - 1) / chunks;
    dispatch_apply(chunks, dispatch_get_global_queue(QOS_CLASS_UTILITY, 0), ^(size_t ci) {
        double mn = INFINITY, mx = 0;
        uint64_t lo = ci * per, hi = lo + per;
        if (hi > N) hi = N;
        for (uint64_t s = lo; s < hi; s++) {
            double sum = 0;
            for (uint64_t r = OFF[s]; r < OFF[s + 1]; r++) sum += RECS[r].cnt * w[RECS[r].wt] * v[RECS[r].tgt];
            vn[s] = sum;
            double q = sum / v[s];
            if (q < mn) mn = q;
            if (q > mx) mx = q;
        }
        pmin[ci] = mn; pmax[ci] = mx;
    });
    double mn = INFINITY, mx = 0;
    for (int i = 0; i < chunks; i++) { if (mins[i] < mn) mn = mins[i]; if (maxs[i] > mx) mx = maxs[i]; }
    *rmin = mn; *rmax = mx;
}

static double *V, *VN;

static double lambda_of(double z, double tol, int verbose, int *iters) {
    double w[256];
    for (int i = 0; i < NW; i++) w[i] = weight_double(&wtypes[i], z);
    double rmin, rmax;
    int it = 0;
    double t0 = now();
    for (;;) {
        iterate(w, V, VN, &rmin, &rmax);
        it++;
        double mx = 0;
        for (uint64_t s = 0; s < N; s++) if (VN[s] > mx) mx = VN[s];
        for (uint64_t s = 0; s < N; s++) V[s] = VN[s] / mx;
        if (verbose && it % 50 == 0)
            fprintf(stderr, "    z=%.6f it=%d ratios [%.15f, %.15f] (%.0fs)\n", z, it, rmin, rmax, now() - t0);
        if (rmax / rmin <= tol) break;
        if (it > 50000) { fprintf(stderr, "no convergence\n"); break; }
    }
    if (verbose)
        fprintf(stderr, "  z=%.8f %d iterations: lambda in [%.15f, %.15f] (%.0fs)\n", z, it, rmin, rmax, now() - t0);
    if (iters) *iters = it;
    return rmax;
}

static double bound_of(double z, double lam) { return 2.0 * log(K * z) / log(lam); }

static void dump_vector(const char *dir, double z) {
    char path[1024];
    path_of(path, sizeof path, dir, "kvec", "bin");
    FILE *f = fopen(path, "wb");
    int32_t hk[2] = {K, H};
    fwrite(hk, sizeof hk, 1, f);
    fwrite(&N, sizeof N, 1, f);
    fwrite(&z, sizeof z, 1, f);
    fwrite(V, sizeof(double), N, f);
    fclose(f);
    fprintf(stderr, "  wrote %s\n", path);
}

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: kary build|eigen|scan k h dir ...\n"); return 1; }
    K = atoi(argv[2]);
    H = atoi(argv[3]);
    if (K < 2 || K > 4 || H < 2 || H > HMAX) { fprintf(stderr, "k or h out of range\n"); return 1; }
    M1 = (1ULL << H) - 1;
    M2 = (1ULL << (2 * H)) - 1;
    const char *dir = argv[4];
    if (!strcmp(argv[1], "build")) { build(dir); return 0; }
    load(dir);
    V = malloc(N * sizeof(double));
    VN = malloc(N * sizeof(double));
    for (uint64_t s = 0; s < N; s++) V[s] = 1.0;
    printf("k=%d h=%d states=%llu records=%llu weight types=%d\n", K, H, (unsigned long long)N,
           (unsigned long long)R, NW);
    if (!strcmp(argv[1], "eigen")) {
        double z = atof(argv[5]);
        double tol = argc > 6 ? 1 + atof(argv[6]) : 1 + 1e-12;
        int it;
        double lam = lambda_of(z, tol, 1, &it);
        printf("z=%.8f lambda=%.15f bound=%.9f iterations=%d\n", z, lam, bound_of(z, lam), it);
        dump_vector(dir, z);
    } else if (!strcmp(argv[1], "scan")) {
        double a = atof(argv[5]), b = atof(argv[6]);
        const double g = (sqrt(5.0) - 1) / 2;
        double c = b - g * (b - a), d = a + g * (b - a);
        double fc = bound_of(c, lambda_of(c, 1 + 1e-10, 0, NULL));
        double fd = bound_of(d, lambda_of(d, 1 + 1e-10, 0, NULL));
        while (b - a > 2e-5) {
            if (fc < fd) { b = d; d = c; fd = fc; c = b - g * (b - a); fc = bound_of(c, lambda_of(c, 1 + 1e-10, 0, NULL)); }
            else { a = c; c = d; fc = fd; d = a + g * (b - a); fd = bound_of(d, lambda_of(d, 1 + 1e-10, 0, NULL)); }
            printf("  [%.6f, %.6f] best so far %.9f\n", a, b, fc < fd ? fc : fd);
            fflush(stdout);
        }
        double zb = fc < fd ? c : d;
        zb = round(zb * 1e4) / 1e4;
        int it;
        double lam = lambda_of(zb, 1 + 1e-12, 1, &it);
        printf("best z=%.4f lambda=%.15f bound=%.9f iterations=%d\n", zb, lam, bound_of(zb, lam), it);
        dump_vector(dir, zb);
    }
    return 0;
}
