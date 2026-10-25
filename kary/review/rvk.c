/* rvk.c: referee's own C implementation of Lueker's machine M_h for k <= 4 letters.
 * Written from Lueker, J. ACM 56(3) 2009, Fig. 7 (same logic as rv_machine.py, own encoding).
 *
 *   rvk classes K H                       print "la lb nA nB" per input class
 *   rvk build K H SYM WFILE OUTPREFIX     build by BFS, power-iterate, certify with directed rounding
 *   rvk delta K H SYM < states            for each packed state (hex), print accepted (class, canon target)
 *
 * State packing (uint64): x letters at bits 2t (t=0..h-1), y letters at 2h+2t,
 * dx bit t (= lcs[t+1,h]-lcs[t,h]) at 4h+t, dy bit t at 5h+t.
 * Orbit form: letters renamed by first occurrence in x1 y1 x2 y2 ... (interleaved).
 * clang -O3 -ffp-model=strict -o rvk rvk.c -lpthread
 */
#include <fenv.h>
#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#pragma STDC FENV_ACCESS ON

#define MAXH 7
typedef struct { int la, lb, a[MAXH], b[MAXH], key; } Cls;

static int K, H, NC, SYM, NTHR = 8;
static Cls *cls;
static int keyid[MAXH + 1][MAXH + 1][5][5], nkeys;

static double now(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

static void enum_classes(void) {
    int cap = 1 << 16;
    cls = malloc(sizeof(Cls) * cap);
    NC = 0;
    memset(keyid, -1, sizeof keyid);
    nkeys = 0;
    for (int c = 0; c < K; c++) {
        int oth[4], no = 0;
        for (int l = 0; l < K; l++)
            if (l != c) oth[no++] = l;
        /* all strings over oth of length 0..H-1, by length then lexicographic */
        static int strs[200000][MAXH];
        static int slen[200000];
        int ns = 0;
        for (int L = 0; L < H; L++) {
            long tot = 1;
            for (int t = 0; t < L; t++) tot *= no;
            for (long code = 0; code < tot; code++) {
                long cc = code;
                for (int t = L - 1; t >= 0; t--) { strs[ns][t] = oth[cc % no]; cc /= no; }
                slen[ns++] = L;
            }
        }
        for (int i = 0; i < ns; i++) {
            int ma = 0;
            for (int t = 0; t < slen[i]; t++) ma |= 1 << strs[i][t];
            for (int j = 0; j < ns; j++) {
                int mb = 0;
                for (int t = 0; t < slen[j]; t++) mb |= 1 << strs[j][t];
                if (ma & mb) continue;
                Cls *p = &cls[NC++];
                p->la = slen[i] + 1;
                p->lb = slen[j] + 1;
                memcpy(p->a, strs[i], sizeof(int) * slen[i]);
                memcpy(p->b, strs[j], sizeof(int) * slen[j]);
                p->a[slen[i]] = c;
                p->b[slen[j]] = c;
                int nA = __builtin_popcount(ma), nB = __builtin_popcount(mb);
                if (keyid[p->la][p->lb][nA][nB] < 0) keyid[p->la][p->lb][nA][nB] = nkeys++;
                p->key = keyid[p->la][p->lb][nA][nB];
            }
        }
    }
}

static inline uint64_t canon(uint64_t q) {
    if (!SYM) return q;
    int mp[4] = {-1, -1, -1, -1}, nxt = 0;
    uint64_t r = q & ~((1ULL << (4 * H)) - 1);
    for (int t = 0; t < H; t++) {
        int cx = (q >> (2 * t)) & 3, cy = (q >> (2 * H + 2 * t)) & 3;
        if (mp[cx] < 0) mp[cx] = nxt++;
        if (mp[cy] < 0) mp[cy] = nxt++;
        r |= (uint64_t)mp[cx] << (2 * t);
        r |= (uint64_t)mp[cy] << (2 * H + 2 * t);
    }
    return r;
}

/* Fig. 7. Returns 1 and the new (uncanonicalised) state, or 0 for the rejecting state. */
static inline int delta(uint64_t q, const Cls *p, uint64_t *out) {
    const int h = H;
    int a[2 * MAXH + 1], b[2 * MAXH + 1], L[2 * MAXH + 1][2 * MAXH + 1];
    for (int t = 0; t < h; t++) {
        a[t + 1] = (q >> (2 * t)) & 3;
        b[t + 1] = (q >> (2 * h + 2 * t)) & 3;
    }
    for (int t = 0; t < p->la; t++) a[h + 1 + t] = p->a[t];
    for (int t = 0; t < p->lb; t++) b[h + 1 + t] = p->b[t];
    int n = h + p->la, m = h + p->lb;
    L[h][h] = h;
    for (int i = h; i >= 1; i--) L[i - 1][h] = L[i][h] - (int)((q >> (4 * h + i - 1)) & 1);
    for (int j = h; j >= 1; j--) L[h][j - 1] = L[h][j] - (int)((q >> (5 * h + j - 1)) & 1);
    for (int i = 0; i <= n; i++)
        for (int j = 0; j <= m; j++) {
            if (i > h || j > h) {
                if (i == 0) L[i][j] = L[i][j - 1];
                else if (j == 0) L[i][j] = L[i - 1][j];
                else if (a[i] == b[j]) L[i][j] = L[i - 1][j - 1] + 1;
                else L[i][j] = L[i][j - 1] > L[i - 1][j] ? L[i][j - 1] : L[i - 1][j];
            }
        }
    int i = n, j = m;
    while (i >= h && j >= h) {
        if (L[i][j] == L[i - 1][j]) i--;
        else if (L[i][j] == L[i][j - 1]) j--;
        else { i--; j--; }
        if (i == h && j == h) {
            uint64_t r = 0;
            for (int t = 0; t < h; t++) {
                r |= (uint64_t)a[n - h + 1 + t] << (2 * t);
                r |= (uint64_t)b[m - h + 1 + t] << (2 * h + 2 * t);
                r |= (uint64_t)(L[n - h + t + 1][m] - L[n - h + t][m]) << (4 * h + t);
                r |= (uint64_t)(L[n][m - h + t + 1] - L[n][m - h + t]) << (5 * h + t);
            }
            *out = r;
            return 1;
        }
    }
    return 0;
}

/* ---------- hash table: packed state -> index ---------- */
static uint64_t *hkey;
static uint32_t *hval;
static uint64_t hmask;
static inline uint64_t mix(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}
static uint64_t *states;
static uint64_t nst, stcap;
static inline uint32_t lookup_insert(uint64_t s) {
    uint64_t pos = mix(s) & hmask;
    for (;;) {
        if (hkey[pos] == 0) {
            hkey[pos] = s + 1;
            if (nst == stcap) { stcap *= 2; states = realloc(states, stcap * 8); }
            states[nst] = s;
            hval[pos] = (uint32_t)nst;
            return (uint32_t)nst++;
        }
        if (hkey[pos] == s + 1) return hval[pos];
        pos = (pos + 1) & hmask;
    }
}
static inline int64_t lookup(uint64_t s) {
    uint64_t pos = mix(s) & hmask;
    for (;;) {
        if (hkey[pos] == 0) return -1;
        if (hkey[pos] == s + 1) return hval[pos];
        pos = (pos + 1) & hmask;
    }
}

/* ---------- BFS ---------- */
#define CHUNK 2048
static uint64_t *buf;
static uint64_t c0, c1;
static void *work_delta(void *arg) {
    long tid = (long)arg;
    for (uint64_t s = c0 + tid; s < c1; s += NTHR) {
        uint64_t q = states[s], t;
        uint64_t *row = buf + (s - c0) * NC;
        for (int ci = 0; ci < NC; ci++) row[ci] = delta(q, &cls[ci], &t) ? canon(t) : ~0ULL;
    }
    return NULL;
}

static uint32_t *tr;     /* (target << CB) | class */
static uint64_t ntr, trcap, *rowst;
static int CB;

static void build(void) {
    hmask = (1ULL << 25) - 1; /* 33.5M slots: enough for every run here */
    hkey = calloc(hmask + 1, 8);
    hval = malloc((hmask + 1) * 4);
    stcap = 1 << 20;
    states = malloc(stcap * 8);
    nst = 0;
    uint64_t s0 = 0;
    for (int t = 0; t < H; t++) s0 |= 1ULL << (4 * H + t) | 1ULL << (5 * H + t); /* x=y=0^h, dx=dy=1^h */
    lookup_insert(canon(s0));
    buf = malloc((size_t)CHUNK * NC * 8);
    trcap = 1 << 24;
    tr = malloc(trcap * 4);
    uint64_t rcap = 1 << 20;
    rowst = malloc(rcap * 8);
    ntr = 0;
    CB = 0;
    while ((1 << CB) < NC) CB++;
    uint64_t qi = 0;
    pthread_t th[64];
    while (qi < nst) {
        c0 = qi;
        c1 = nst < qi + CHUNK ? nst : qi + CHUNK;
        for (long t = 0; t < NTHR; t++) pthread_create(&th[t], NULL, work_delta, (void *)t);
        for (long t = 0; t < NTHR; t++) pthread_join(th[t], NULL);
        for (uint64_t s = c0; s < c1; s++) {
            if (s + 1 >= rcap) { rcap *= 2; rowst = realloc(rowst, rcap * 8); }
            rowst[s] = ntr;
            uint64_t *row = buf + (s - c0) * NC;
            for (int ci = 0; ci < NC; ci++) {
                if (row[ci] == ~0ULL) continue;
                uint32_t ti = lookup_insert(row[ci]);
                if (ntr == trcap) { trcap *= 2; tr = realloc(tr, trcap * 4); }
                if (((uint64_t)ti << CB) >> 32) { fprintf(stderr, "target index overflow\n"); exit(1); }
                tr[ntr++] = (ti << CB) | (uint32_t)ci;
            }
        }
        qi = c1;
    }
    rowst[nst] = ntr;
}

/* ---------- power iteration and directed-rounding certificate ---------- */
static double *wnear, *wup, *wdn, *v, *w;
static double part_hi[64], part_lo[64];
static int mode; /* 0 nearest spmv, 1 upward max ratio, 2 downward min ratio */
static int selftest_fail;
static void *work_spmv(void *arg) {
    long tid = (long)arg;
    uint64_t lo = nst * tid / NTHR, hi = nst * (tid + 1) / NTHR;
    const double *wt = mode == 1 ? wup : mode == 2 ? wdn : wnear;
    if (mode == 1) fesetround(FE_UPWARD);
    if (mode == 2) fesetround(FE_DOWNWARD);
    volatile double one = 1.0, tiny = 0x1p-60, three = 3.0, ten = 10.0;
    if (mode == 1 && !(one + tiny > one && one - tiny == one && one / three > 0x1.5555555555555p-2)) selftest_fail = 1;
    if (mode == 2 && !(one - tiny < one && one + tiny == one && one / ten < 0x1.999999999999ap-4)) selftest_fail = 1;
    if (mode == 0 && !(fegetround() == FE_TONEAREST && one + tiny == one && one - tiny == one)) selftest_fail = 1;
    double mx = 0, mn = INFINITY;
    for (uint64_t s = lo; s < hi; s++) {
        double sum = 0;
        for (uint64_t e = rowst[s]; e < rowst[s + 1]; e++) sum += wt[tr[e] & ((1u << CB) - 1)] * v[tr[e] >> CB];
        if (mode == 0) w[s] = sum;
        double r = sum / v[s];
        if (r > mx) mx = r;
        if (r < mn) mn = r;
    }
    part_hi[tid] = mx;
    part_lo[tid] = mn;
    fesetround(FE_TONEAREST);
    return NULL;
}
static void spmv(int md, double *hi, double *lo) {
    pthread_t th[64];
    mode = md;
    for (long t = 0; t < NTHR; t++) pthread_create(&th[t], NULL, work_spmv, (void *)t);
    for (long t = 0; t < NTHR; t++) pthread_join(th[t], NULL);
    *hi = 0;
    *lo = INFINITY;
    for (int t = 0; t < NTHR; t++) {
        if (part_hi[t] > *hi) *hi = part_hi[t];
        if (part_lo[t] < *lo) *lo = part_lo[t];
    }
}

static void read_weights(const char *fn) {
    FILE *f = fopen(fn, "r");
    if (!f) { perror(fn); exit(1); }
    wnear = malloc(NC * 8); wup = malloc(NC * 8); wdn = malloc(NC * 8);
    for (int ci = 0; ci < NC; ci++)
        if (fscanf(f, "%la %la %la", &wdn[ci], &wnear[ci], &wup[ci]) != 3) { fprintf(stderr, "weights short\n"); exit(1); }
    fclose(f);
    for (int ci = 0; ci < NC; ci++)
        if (!(wdn[ci] > 0 && wdn[ci] <= wnear[ci] && wnear[ci] <= wup[ci])) { fprintf(stderr, "bad weight %d\n", ci); exit(1); }
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage\n"); return 1; }
    K = atoi(argv[2]);
    H = atoi(argv[3]);
    if (getenv("RVK_THREADS")) NTHR = atoi(getenv("RVK_THREADS"));
    enum_classes();
    if (!strcmp(argv[1], "classes")) {
        for (int ci = 0; ci < NC; ci++) {
            int ma = 0, mb = 0;
            for (int t = 0; t < cls[ci].la - 1; t++) ma |= 1 << cls[ci].a[t];
            for (int t = 0; t < cls[ci].lb - 1; t++) mb |= 1 << cls[ci].b[t];
            printf("%d %d %d %d\n", cls[ci].la, cls[ci].lb, __builtin_popcount(ma), __builtin_popcount(mb));
        }
        return 0;
    }
    SYM = atoi(argv[4]);
    if (!strcmp(argv[1], "delta")) {
        char line[256];
        while (fgets(line, sizeof line, stdin)) {
            uint64_t q = strtoull(line, NULL, 16), t;
            printf("%llx:", (unsigned long long)q);
            for (int ci = 0; ci < NC; ci++)
                if (delta(q, &cls[ci], &t)) printf(" %d:%llx", ci, (unsigned long long)canon(t));
            printf("\n");
        }
        return 0;
    }
    /* build */
    double t0 = now();
    build();
    double t1 = now();
    /* records: distinct (target, weight key) per state */
    uint64_t rec = 0;
    {
        uint64_t *tmp = malloc(sizeof(uint64_t) * NC);
        for (uint64_t s = 0; s < nst; s++) {
            int nn = 0;
            for (uint64_t e = rowst[s]; e < rowst[s + 1]; e++)
                tmp[nn++] = ((uint64_t)(tr[e] >> CB) << 16) | (uint64_t)cls[tr[e] & ((1u << CB) - 1)].key;
            for (int x = 1; x < nn; x++) { /* insertion sort is fine for <= NC entries */
                uint64_t key = tmp[x];
                int y = x - 1;
                while (y >= 0 && tmp[y] > key) { tmp[y + 1] = tmp[y]; y--; }
                tmp[y + 1] = key;
            }
            for (int x = 0; x < nn; x++)
                if (x == 0 || tmp[x] != tmp[x - 1]) rec++;
        }
        free(tmp);
    }
    uint64_t empty_rows = 0;
    for (uint64_t s = 0; s < nst; s++) empty_rows += rowst[s + 1] == rowst[s];
    printf("k=%d h=%d sym=%d: states %llu, classes %d, accepted %llu, records %llu, empty rows %llu, build %.1fs\n",
           K, H, SYM, (unsigned long long)nst, NC, (unsigned long long)ntr, (unsigned long long)rec,
           (unsigned long long)empty_rows, t1 - t0);
    fflush(stdout);
    if (argc < 6) return 0;
    read_weights(argv[5]);
    v = malloc(nst * 8);
    w = malloc(nst * 8);
    for (uint64_t s = 0; s < nst; s++) v[s] = 1.0;
    double hi = 0, lo = 0;
    int it;
    for (it = 0; it < 5000; it++) {
        spmv(0, &hi, &lo);
        double mx = 0;
        for (uint64_t s = 0; s < nst; s++) if (w[s] > mx) mx = w[s];
        for (uint64_t s = 0; s < nst; s++) v[s] = w[s] / mx;
        if (hi / lo - 1 < 2e-13) break;
    }
    printf("power iteration: %d its, last nearest ratios [%.17g, %.17g]\n", it, lo, hi);
    for (uint64_t s = 0; s < nst; s++)
        if (!(v[s] > 0) || !isfinite(v[s])) { printf("NONPOSITIVE v at %llu\n", (unsigned long long)s); return 2; }
    double lamup, dummy, rholo;
    selftest_fail = 0;
    spmv(1, &lamup, &dummy);
    spmv(2, &dummy, &rholo);
    double ln, ll;
    spmv(0, &ln, &ll);
    if (selftest_fail) { printf("ROUNDING SELF-TEST FAILED\n"); return 3; }
    printf("lambda_up %a = %.17g\nrho_low   %a = %.17g\nnearest   [%.17g, %.17g]\n", lamup, lamup, rholo, rholo, ll, ln);
    printf("time total %.1fs\n", now() - t0);
    char fn[512];
    snprintf(fn, sizeof fn, "%s.states", argv[6 < argc ? 6 : 5]);
    if (argc > 6) {
        FILE *f = fopen(fn, "wb"); fwrite(states, 8, nst, f); fclose(f);
        snprintf(fn, sizeof fn, "%s.vec", argv[6]);
        f = fopen(fn, "wb"); fwrite(v, 8, nst, f); fclose(f);
    }
    return 0;
}
