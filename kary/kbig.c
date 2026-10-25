/*
 * kbig.c: kary.c for tables larger than memory (k = 3, h = 7 and k = 4, h = 6).
 *
 * The machine, the inputs, the weight types, the state numbering and the
 * records are those of kary.c (make_inputs, delta, wtype_of, cmp_kw and the
 * hash table are copied from it unchanged), so for every (k, h) the table holds
 * the same data.  Only the way it is produced and stored differs:
 *
 *   - records go to disk as they are produced, in state order;
 *   - each worker sorts and merges the transitions of its state as kary.c does
 *     and looks the targets up in the hash table (read only, in parallel); the
 *     keys not yet known are then inserted one at a time, in the same
 *     (state, record) order in which kary.c calls find_or_add, so every id is
 *     the one kary.c gives;
 *   - the eigen phase streams the records from disk (reader thread, two
 *     buffers) and iterates up to 8 values of z in each pass over the table.
 *
 * File kbig_k<k>_h<h>.bin:
 *   "KBIG", int32 k, h, uint64 N, R, int32 NI, NW, WType[NW],
 *   Rec[R] (state order), uint64 off[N+1], uint64 key[N].
 * kary.c's ktab file holds the same sections with off[] before the records;
 * "kbig klue" writes that file, so that small cases can be compared with
 * kary.c's output byte for byte.
 *
 * Usage:
 *   kbig build   k h dir
 *   kbig klue    k h dir                 kbig_k<k>_h<h>.bin -> ktab_k<k>_h<h>.bin
 *   kbig project k h dir prevdir         start vector from prevdir's (h-1) table and vector
 *   kbig eigen   k h dir start tol maxit z1 [z2 ... z8]
 *       start: "-" (all ones), a vector file (raw doubles or kvec format) or a
 *       checkpoint kvecs_k<k>_h<h>.bin (each z starts from the interpolation of the
 *       two checkpoint vectors whose z bracket it).
 *       Stops when max/min ratio <= 1 + tol for every z, or after maxit passes.
 *       Writes the checkpoint every 10 passes and at the end, and
 *       kvec_k<k>_h<h>_z<z>.bin (kary.c's vector format) for every z.
 *
 * KBIG_QOS=user runs the worker threads at user-initiated QoS instead of utility.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <fcntl.h>
#include <unistd.h>
#include <pthread.h>
#include <sys/stat.h>
#include <dispatch/dispatch.h>

#define HMAX 10
#define NMAX (2 * HMAX + 2)
#define MZ 8

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

/* ---------------- copied unchanged from kary.c ---------------------------- */
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
/* ---------------- end of the part copied from kary.c ---------------------- */

#define NOT_FOUND UINT32_MAX

/* worker threads run at utility QoS unless KBIG_QOS=user (user-initiated, as the binary
 * h = 10 run did); the results do not depend on it */
static dispatch_queue_t workq(void) {
    const char *q = getenv("KBIG_QOS");
    return dispatch_get_global_queue(q && !strcmp(q, "user") ? QOS_CLASS_USER_INITIATED : QOS_CLASS_UTILITY, 0);
}

/* read-only lookup, safe while no thread inserts */
static inline uint32_t lookup(uint64_t key) {
    uint64_t p = mix(key) & hmask;
    for (;;) {
        uint64_t k = hkeys[p];
        if (k == key) return hvals[p];
        if (k == EMPTY) return NOT_FOUND;
        p = (p + 1) & hmask;
    }
}

static void path_of(char *buf, size_t n, const char *dir, const char *stem, const char *ext) {
    snprintf(buf, n, "%s/%s_k%d_h%d.%s", dir, stem, K, H, ext);
}

static uint64_t header_bytes(int nw) { return 4 + 8 + 16 + 8 + (uint64_t)nw * sizeof(WType); }

static void xwrite(FILE *f, const void *p, size_t size, size_t n, const char *what) {
    if (fwrite(p, size, n, f) != n) { fprintf(stderr, "write failed (%s)\n", what); exit(2); }
}

static void build(const char *dir) {
    make_inputs();
    fprintf(stderr, "k=%d h=%d inputs per state=%d weight types=%d\n", K, H, NI, NW);
    if (NI > 65535) { fprintf(stderr, "too many inputs for 16-bit multiplicities\n"); exit(2); }
    hcap = 1 << 20; hmask = hcap - 1;
    hkeys = malloc(hcap * sizeof(uint64_t));
    hvals = malloc(hcap * sizeof(uint32_t));
    memset(hkeys, 0xff, hcap * sizeof(uint64_t));
    scap = 1 << 20; nstates = 0;
    states = malloc(scap * sizeof(uint64_t));
    find_or_add(pack(0, 0, M1, M1));
    char path[1024];
    path_of(path, sizeof path, dir, "kbig", "bin");
    FILE *f = fopen(path, "wb");
    if (!f) { perror(path); exit(2); }
    setvbuf(f, NULL, _IOFBF, 1 << 24);
    fcntl(fileno(f), F_NOCACHE, 1);
    int32_t hk[2] = {K, H};
    uint64_t hn[2] = {0, 0};
    int32_t hi2[2] = {NI, NW};
    xwrite(f, "KBIG", 4, 1, "magic");
    xwrite(f, hk, sizeof hk, 1, "header");
    xwrite(f, hn, sizeof hn, 1, "header");
    xwrite(f, hi2, sizeof hi2, 1, "header");
    xwrite(f, wtypes, sizeof(WType), NW, "weight types");
    uint64_t ocap = 1 << 20, nrec = 0;
    uint64_t *offs = malloc(ocap * sizeof(uint64_t));
    offs[0] = 0;
    uint64_t B = (1ULL << 28) / ((uint64_t)NI * sizeof(KW));
    if (B > 8192) B = 8192;
    if (B < 64) B = 64;
    KW *buf = malloc(B * NI * sizeof(KW));
    uint32_t *ids = malloc(B * NI * sizeof(uint32_t));
    uint16_t *mult = malloc(B * NI * sizeof(uint16_t));
    uint32_t *cntm = malloc(B * sizeof(uint32_t));
    Rec *rbuf = malloc(B * NI * sizeof(Rec));
    if (!buf || !ids || !mult || !cntm || !rbuf || !offs) { fprintf(stderr, "out of memory (buffers)\n"); exit(2); }
    uint64_t accepted = 0, inserted_late = 0;
    double t0 = now(), tlast = t0;
    for (uint64_t lo = 0; lo < nstates;) {
        uint64_t hi = lo + B < nstates ? lo + B : nstates;
        const uint64_t *st = states;          /* no insertion (hence no realloc) until the workers finish */
        dispatch_apply(hi - lo, workq(), ^(size_t r) {
            uint64_t key = st[lo + r];
            KW *kb = buf + r * NI;
            uint16_t *mb = mult + r * NI;
            uint32_t *ib = ids + r * NI;
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
            for (int t = 0; t < u; t++) ib[t] = lookup(kb[t].key);
            cntm[r] = (uint32_t)u;
        });
        uint64_t n = 0;
        for (uint64_t s = lo; s < hi; s++) {
            uint64_t r = s - lo;
            const KW *kb = buf + r * NI;
            const uint16_t *mb = mult + r * NI;
            const uint32_t *ib = ids + r * NI;
            for (uint32_t t = 0; t < cntm[r]; t++) {
                uint32_t id = ib[t];
                if (id == NOT_FOUND) { id = find_or_add(kb[t].key); inserted_late++; }
                rbuf[n].tgt = id; rbuf[n].wt = (uint16_t)kb[t].wt; rbuf[n].cnt = mb[t];
                accepted += mb[t];
                n++;
            }
            if (s + 2 > ocap) {
                ocap *= 2;
                offs = realloc(offs, ocap * sizeof(uint64_t));
                if (!offs) { fprintf(stderr, "out of memory (offsets)\n"); exit(2); }
            }
            offs[s + 1] = nrec + n;
        }
        xwrite(f, rbuf, sizeof(Rec), n, "records");
        nrec += n;
        lo = hi;
        if (now() - tlast > 30) {
            tlast = now();
            fprintf(stderr, "  processed %llu, discovered %llu, records %llu (%.0fs)\n",
                    (unsigned long long)lo, (unsigned long long)nstates, (unsigned long long)nrec, tlast - t0);
        }
    }
    xwrite(f, offs, sizeof(uint64_t), nstates + 1, "offsets");
    xwrite(f, states, sizeof(uint64_t), nstates, "keys");
    if (fflush(f)) { perror("fflush"); exit(2); }
    if (fseek(f, 4 + sizeof hk, SEEK_SET)) { perror("fseek"); exit(2); }
    hn[0] = nstates; hn[1] = nrec;
    xwrite(f, hn, sizeof hn, 1, "header");
    if (fclose(f)) { perror("fclose"); exit(2); }
    char wpath[1024];
    path_of(wpath, sizeof wpath, dir, "wtypes", "txt");
    FILE *fw = fopen(wpath, "w");
    for (int i = 0; i < NW; i++)
        fprintf(fw, "%d %d %d %d %d %d\n", i, wtypes[i].size, wtypes[i].ta, wtypes[i].tb, wtypes[i].na, wtypes[i].nb);
    fclose(fw);
    printf("k=%d h=%d inputs=%d accepting states=%llu accepted transitions=%llu records=%llu build %.1fs -> %s\n",
           K, H, NI, (unsigned long long)nstates, (unsigned long long)accepted,
           (unsigned long long)nrec, now() - t0, path);
    fprintf(stderr, "  records whose target was not yet known at the parallel lookup: %llu\n", (unsigned long long)inserted_late);
}

/* ---------------- reading a KBIG file ------------------------------------- */
static uint64_t N, R;
static uint64_t *OFF;
static uint64_t RECBASE;           /* byte offset of Rec[0] in the file */

static FILE *open_kbig(const char *dir, char *path, size_t n) {
    path_of(path, n, dir, "kbig", "bin");
    FILE *f = fopen(path, "rb");
    if (!f) { perror(path); exit(2); }
    char magic[4];
    int32_t hk[2], hi2[2];
    uint64_t hn[2];
    if (fread(magic, 4, 1, f) != 1 || memcmp(magic, "KBIG", 4) || fread(hk, sizeof hk, 1, f) != 1 ||
        fread(hn, sizeof hn, 1, f) != 1 || fread(hi2, sizeof hi2, 1, f) != 1 || hk[0] != K || hk[1] != H) {
        fprintf(stderr, "bad header in %s\n", path); exit(2);
    }
    N = hn[0]; R = hn[1]; NI = hi2[0]; NW = hi2[1];
    if (NW < 1 || NW > 256 || fread(wtypes, sizeof(WType), NW, f) != (size_t)NW) { fprintf(stderr, "bad weight types\n"); exit(2); }
    RECBASE = header_bytes(NW);
    struct stat sb;
    if (fstat(fileno(f), &sb) || (uint64_t)sb.st_size != RECBASE + 8 * R + 8 * (N + 1) + 8 * N) {
        fprintf(stderr, "file size does not match the header\n"); exit(2);
    }
    return f;
}

static void load_offsets(FILE *f) {
    OFF = malloc((N + 1) * sizeof(uint64_t));
    if (!OFF) { fprintf(stderr, "out of memory (offsets)\n"); exit(2); }
    if (fseeko(f, (off_t)(RECBASE + 8 * R), SEEK_SET) || fread(OFF, sizeof(uint64_t), N + 1, f) != N + 1) {
        fprintf(stderr, "short offsets\n"); exit(2);
    }
    if (OFF[0] != 0 || OFF[N] != R) { fprintf(stderr, "bad offsets\n"); exit(2); }
    for (uint64_t s = 0; s < N; s++) if (OFF[s + 1] < OFF[s]) { fprintf(stderr, "offsets not monotone\n"); exit(2); }
}

static uint64_t *load_keys(FILE *f) {
    uint64_t *keys = malloc(N * sizeof(uint64_t));
    if (!keys) { fprintf(stderr, "out of memory (keys)\n"); exit(2); }
    if (fseeko(f, (off_t)(RECBASE + 8 * R + 8 * (N + 1)), SEEK_SET) || fread(keys, sizeof(uint64_t), N, f) != N) {
        fprintf(stderr, "short keys\n"); exit(2);
    }
    return keys;
}

/* kary.c's ktab format, for comparison with kary.c */
static void klue(const char *dir) {
    char path[1024];
    FILE *f = open_kbig(dir, path, sizeof path);
    load_offsets(f);
    uint64_t *keys = load_keys(f);
    char op[1024];
    path_of(op, sizeof op, dir, "ktab", "bin");
    FILE *o = fopen(op, "wb");
    if (!o) { perror(op); exit(2); }
    int32_t hk[2] = {K, H}, hi2[2] = {NI, NW};
    uint64_t hn[2] = {N, R};
    xwrite(o, "KLUE", 4, 1, "magic");
    xwrite(o, hk, sizeof hk, 1, "header");
    xwrite(o, hn, sizeof hn, 1, "header");
    xwrite(o, hi2, sizeof hi2, 1, "header");
    xwrite(o, wtypes, sizeof(WType), NW, "weight types");
    xwrite(o, OFF, sizeof(uint64_t), N + 1, "offsets");
    if (fseeko(f, (off_t)RECBASE, SEEK_SET)) { perror("fseeko"); exit(2); }
    const uint64_t CH = 1 << 24;
    Rec *rb = malloc(CH * sizeof(Rec));
    for (uint64_t done = 0; done < R;) {
        uint64_t c = R - done < CH ? R - done : CH;
        if (fread(rb, sizeof(Rec), c, f) != c) { fprintf(stderr, "short records\n"); exit(2); }
        xwrite(o, rb, sizeof(Rec), c, "records");
        done += c;
    }
    xwrite(o, keys, sizeof(uint64_t), N, "keys");
    if (fclose(o)) { perror("fclose"); exit(2); }
    fclose(f);
    printf("wrote %s (N=%llu, R=%llu)\n", op, (unsigned long long)N, (unsigned long long)R);
}

/* ---------------- start vector from the (h-1) vector ---------------------- */
static uint64_t *read_prev_keys(const char *prevdir, int hp, uint64_t *np) {
    /* kary.c's ktab or a kbig file for (K, hp) */
    char p[1024];
    snprintf(p, sizeof p, "%s/ktab_k%d_h%d.bin", prevdir, K, hp);
    FILE *f = fopen(p, "rb");
    char magic[4];
    int32_t hk[2], hi2[2];
    uint64_t hn[2];
    if (f) {
        if (fread(magic, 4, 1, f) != 1 || memcmp(magic, "KLUE", 4) || fread(hk, sizeof hk, 1, f) != 1 ||
            fread(hn, sizeof hn, 1, f) != 1 || fread(hi2, sizeof hi2, 1, f) != 1 || hk[0] != K || hk[1] != hp) {
            fprintf(stderr, "bad header in %s\n", p); exit(2);
        }
        uint64_t n = hn[0], r = hn[1];
        uint64_t *keys = malloc(n * sizeof(uint64_t));
        off_t at = (off_t)(header_bytes(hi2[1]) + 8 * (n + 1) + 8 * r);
        if (fseeko(f, at, SEEK_SET) || fread(keys, sizeof(uint64_t), n, f) != n) { fprintf(stderr, "short keys in %s\n", p); exit(2); }
        fclose(f);
        *np = n;
        return keys;
    }
    snprintf(p, sizeof p, "%s/kbig_k%d_h%d.bin", prevdir, K, hp);
    f = fopen(p, "rb");
    if (!f) { fprintf(stderr, "no table for h=%d in %s\n", hp, prevdir); exit(2); }
    if (fread(magic, 4, 1, f) != 1 || memcmp(magic, "KBIG", 4) || fread(hk, sizeof hk, 1, f) != 1 ||
        fread(hn, sizeof hn, 1, f) != 1 || fread(hi2, sizeof hi2, 1, f) != 1 || hk[0] != K || hk[1] != hp) {
        fprintf(stderr, "bad header in %s\n", p); exit(2);
    }
    uint64_t n = hn[0], r = hn[1];
    uint64_t *keys = malloc(n * sizeof(uint64_t));
    off_t at = (off_t)(header_bytes(hi2[1]) + 8 * r + 8 * (n + 1));
    if (fseeko(f, at, SEEK_SET) || fread(keys, sizeof(uint64_t), n, f) != n) { fprintf(stderr, "short keys in %s\n", p); exit(2); }
    fclose(f);
    *np = n;
    return keys;
}

static void project(const char *dir, const char *prevdir) {
    char path[1024];
    FILE *f = open_kbig(dir, path, sizeof path);
    uint64_t *keys = load_keys(f);
    fclose(f);
    int hp = H - 1;
    uint64_t np;
    uint64_t *pk = read_prev_keys(prevdir, hp, &np);
    char p[1024];
    snprintf(p, sizeof p, "%s/kvec_k%d_h%d.bin", prevdir, K, hp);
    FILE *g = fopen(p, "rb");
    int32_t ek[2]; uint64_t en; double ez;
    if (!g || fread(ek, sizeof ek, 1, g) != 1 || fread(&en, 8, 1, g) != 1 || fread(&ez, 8, 1, g) != 1 ||
        ek[0] != K || ek[1] != hp || en != np) { fprintf(stderr, "bad vector %s\n", p); exit(2); }
    double *pv = malloc(np * sizeof(double));
    if (fread(pv, sizeof(double), np, g) != np) { fprintf(stderr, "short vector %s\n", p); exit(2); }
    fclose(g);
    uint64_t cap = 1;
    while (cap < 2 * np) cap <<= 1;
    uint64_t *hk = malloc(cap * 8);
    uint32_t *hv = malloc(cap * 4);
    memset(hk, 0xff, cap * 8);
    for (uint64_t i = 0; i < np; i++) {
        uint64_t q = mix(pk[i]) & (cap - 1);
        while (hk[q] != EMPTY) q = (q + 1) & (cap - 1);
        hk[q] = pk[i]; hv[q] = (uint32_t)i;
    }
    double *v = malloc(N * sizeof(double));
    uint64_t hit = 0;
    double sum = 0;
    uint64_t m1 = (1ULL << hp) - 1;
    for (uint64_t s = 0; s < N; s++) {
        uint64_t key = keys[s];
        uint64_t x = (key >> (4 * H)) & M2, y = (key >> (2 * H)) & M2;
        uint64_t dx = (key >> H) & M1, dy = key & M1;
        /* drop the oldest character of x and y and the oldest difference of dx and dy */
        int cx[HMAX], cy[HMAX];
        for (int t = hp - 1; t >= 0; t--) { cx[t] = (int)(x & 3); x >>= 2; }
        for (int t = hp - 1; t >= 0; t--) { cy[t] = (int)(y & 3); y >>= 2; }
        int map[4] = {-1, -1, -1, -1}, nl = 0;
        for (int t = 0; t < hp; t++) if (map[cx[t]] < 0) map[cx[t]] = nl++;
        for (int t = 0; t < hp; t++) if (map[cy[t]] < 0) map[cy[t]] = nl++;
        uint64_t nx = 0, ny = 0;
        for (int t = 0; t < hp; t++) nx = (nx << 2) | (uint64_t)map[cx[t]];
        for (int t = 0; t < hp; t++) ny = (ny << 2) | (uint64_t)map[cy[t]];
        uint64_t kk = (nx << (4 * hp)) | (ny << (2 * hp)) | ((dx & m1) << hp) | (dy & m1);
        uint64_t q = mix(kk) & (cap - 1);
        v[s] = -1;
        while (hk[q] != EMPTY) {
            if (hk[q] == kk) { v[s] = pv[hv[q]]; hit++; sum += v[s]; break; }
            q = (q + 1) & (cap - 1);
        }
    }
    double mean = hit ? sum / hit : 1.0;
    for (uint64_t s = 0; s < N; s++) if (v[s] < 0) v[s] = mean;
    path_of(p, sizeof p, dir, "start", "bin");
    FILE *o = fopen(p, "wb");
    xwrite(o, v, sizeof(double), N, "start vector");
    fclose(o);
    printf("projected %llu of %llu states onto the h=%d vector (others %.6g) -> %s\n",
           (unsigned long long)hit, (unsigned long long)N, hp, mean, p);
}

/* ---------------- streaming eigen phase ----------------------------------- */
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

/* exact weight of a type, in double (the certificate uses kweights.py instead) */
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

typedef struct { uint64_t s0, s1; } Chunk;
static Chunk *CHK;
static uint64_t NCH;
static int TFD;
static Rec *RBUF[2];               /* records of the chunk held by each buffer */
static char *RAW[2];               /* page-aligned read buffers; RBUF[i] points into RAW[i] */
static uint64_t FILESIZE;
#define PAGE 16384                 /* F_NOCACHE reads bypass the page cache only when page-aligned */
static int READY[2];
static uint64_t HELD[2];           /* chunk index held by each buffer */
static uint64_t NEXT_CHUNK;
static int STOP;
static pthread_mutex_t MU = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t CV = PTHREAD_COND_INITIALIZER;

static void *reader_main(void *arg) {
    (void)arg;
    int which = 0;
    for (;;) {
        pthread_mutex_lock(&MU);
        while (READY[which] && !STOP) pthread_cond_wait(&CV, &MU);
        if (STOP) { pthread_mutex_unlock(&MU); return NULL; }
        uint64_t c = NEXT_CHUNK;
        NEXT_CHUNK = (c + 1 == NCH) ? 0 : c + 1;      /* wrap around for the next pass */
        pthread_mutex_unlock(&MU);
        uint64_t r0 = OFF[CHK[c].s0], r1 = OFF[CHK[c].s1];
        /* read the page-aligned span around the chunk's records, so that the read is uncached */
        uint64_t first = RECBASE + r0 * sizeof(Rec), last = RECBASE + r1 * sizeof(Rec);
        uint64_t a0 = first & ~(uint64_t)(PAGE - 1), a1 = (last + PAGE - 1) & ~(uint64_t)(PAGE - 1);
        if (a1 > FILESIZE) a1 = FILESIZE;
        uint64_t bytes = a1 - a0, done = 0;
        while (done < bytes) {
            uint64_t want = bytes - done;
            if (want > (1ULL << 30)) want = 1ULL << 30;   /* macOS caps one read at INT_MAX */
            ssize_t got = pread(TFD, RAW[which] + done, want, (off_t)(a0 + done));
            if (got <= 0) { perror("pread"); exit(2); }
            done += (uint64_t)got;
        }
        RBUF[which] = (Rec *)(RAW[which] + (first - a0));
        pthread_mutex_lock(&MU);
        HELD[which] = c; READY[which] = 1;
        pthread_cond_broadcast(&CV);
        pthread_mutex_unlock(&MU);
        which ^= 1;
    }
}

static double *V, *VN;             /* N x MZ, interleaved: V[s * MZ + j] */
static double W[256][MZ];
static int NZ;
static double Z[MZ];

static void iterate(double *rmin, double *rmax) {
    static int which = 0;
    for (int j = 0; j < MZ; j++) { rmin[j] = INFINITY; rmax[j] = 0; }
    for (uint64_t c = 0; c < NCH; c++) {
        pthread_mutex_lock(&MU);
        while (!READY[which]) pthread_cond_wait(&CV, &MU);
        pthread_mutex_unlock(&MU);
        if (HELD[which] != c) { fprintf(stderr, "reader out of step\n"); exit(2); }
        const Rec *rb = RBUF[which];
        uint64_t s0 = CHK[c].s0, s1 = CHK[c].s1, rbase = OFF[s0];
        const int nb = 256;
        uint64_t per = (s1 - s0 + nb - 1) / nb;
        double bmin[256][MZ], bmax[256][MZ];
        double (*pmin)[MZ] = bmin, (*pmax)[MZ] = bmax;
        dispatch_apply(nb, workq(), ^(size_t bi) {
            double mn[MZ], mx[MZ];
            for (int j = 0; j < MZ; j++) { mn[j] = INFINITY; mx[j] = 0; }
            uint64_t a = s0 + bi * per, b = a + per;
            if (b > s1) b = s1;
            for (uint64_t s = a; s < b; s++) {
                double acc[MZ] = {0};
                for (uint64_t r = OFF[s]; r < OFF[s + 1]; r++) {
                    Rec x = rb[r - rbase];
                    const double *vt = V + (uint64_t)x.tgt * MZ;
                    const double *w = W[x.wt];
                    double cn = x.cnt;
                    for (int j = 0; j < MZ; j++) acc[j] += cn * w[j] * vt[j];
                }
                for (int j = 0; j < MZ; j++) {
                    VN[s * MZ + j] = acc[j];
                    double q = acc[j] / V[s * MZ + j];
                    if (q < mn[j]) mn[j] = q;
                    if (q > mx[j]) mx[j] = q;
                }
            }
            for (int j = 0; j < MZ; j++) { pmin[bi][j] = mn[j]; pmax[bi][j] = mx[j]; }
        });
        for (int bi = 0; bi < nb; bi++)
            for (int j = 0; j < MZ; j++) {
                if (bmin[bi][j] < rmin[j]) rmin[j] = bmin[bi][j];
                if (bmax[bi][j] > rmax[j]) rmax[j] = bmax[bi][j];
            }
        pthread_mutex_lock(&MU);
        READY[which] = 0;
        pthread_cond_broadcast(&CV);
        pthread_mutex_unlock(&MU);
        which ^= 1;
    }
}

static void normalize(void) {
    double mx[MZ];
    for (int j = 0; j < MZ; j++) mx[j] = 0;
    for (uint64_t s = 0; s < N; s++)
        for (int j = 0; j < MZ; j++) if (VN[s * MZ + j] > mx[j]) mx[j] = VN[s * MZ + j];
    const double *pm = mx;
    const int nb = 256;
    uint64_t per = (N + nb - 1) / nb;
    dispatch_apply(nb, workq(), ^(size_t bi) {
        uint64_t a = bi * per, b = a + per;
        if (b > N) b = N;
        for (uint64_t s = a; s < b; s++)
            for (int j = 0; j < MZ; j++) V[s * MZ + j] = VN[s * MZ + j] / pm[j];
    });
}

static void save_checkpoint(const char *dir) {
    char p[1024], tmp[1100];
    path_of(p, sizeof p, dir, "kvecs", "bin");
    snprintf(tmp, sizeof tmp, "%s.tmp", p);
    FILE *f = fopen(tmp, "wb");
    if (!f) { perror(tmp); exit(2); }
    int32_t hk[3] = {K, H, NZ};
    xwrite(f, "KVM8", 4, 1, "checkpoint");
    xwrite(f, hk, sizeof hk, 1, "checkpoint");
    xwrite(f, &N, sizeof N, 1, "checkpoint");
    xwrite(f, Z, sizeof(double), MZ, "checkpoint");
    xwrite(f, V, sizeof(double), N * MZ, "checkpoint");
    if (fclose(f)) { perror("fclose"); exit(2); }
    if (rename(tmp, p)) { perror("rename"); exit(2); }
}

static void save_vectors(const char *dir) {
    double *v = malloc(N * sizeof(double));
    for (int j = 0; j < NZ; j++) {
        char p[1024];
        snprintf(p, sizeof p, "%s/kvec_k%d_h%d_z%.4f.bin", dir, K, H, Z[j]);
        for (uint64_t s = 0; s < N; s++) v[s] = V[s * MZ + j];
        FILE *f = fopen(p, "wb");
        if (!f) { perror(p); exit(2); }
        int32_t hk[2] = {K, H};
        xwrite(f, hk, sizeof hk, 1, "vector");
        xwrite(f, &N, sizeof N, 1, "vector");
        xwrite(f, &Z[j], sizeof(double), 1, "vector");
        xwrite(f, v, sizeof(double), N, "vector");
        if (fclose(f)) { perror("fclose"); exit(2); }
        fprintf(stderr, "  wrote %s\n", p);
    }
    free(v);
}

static void load_start(const char *start) {
    if (!strcmp(start, "-")) {
        for (uint64_t i = 0; i < N * MZ; i++) V[i] = 1.0;
        return;
    }
    FILE *g = fopen(start, "rb");
    if (!g) { perror(start); exit(2); }
    struct stat sb;
    fstat(fileno(g), &sb);
    char magic[4] = {0};
    if (fread(magic, 4, 1, g) == 1 && !memcmp(magic, "KVM8", 4)) {
        int32_t hk[3]; uint64_t n; double z0[MZ];
        if (fread(hk, sizeof hk, 1, g) != 1 || fread(&n, 8, 1, g) != 1 || fread(z0, sizeof(double), MZ, g) != MZ ||
            hk[0] != K || hk[1] != H || n != N || hk[2] < 1 || hk[2] > MZ) { fprintf(stderr, "bad checkpoint\n"); exit(2); }
        double *old = malloc(N * MZ * sizeof(double));
        if (fread(old, sizeof(double), N * MZ, g) != N * MZ) { fprintf(stderr, "short checkpoint\n"); exit(2); }
        /* each z starts from the linear interpolation, in z, of the two checkpoint vectors
         * that bracket it (the nearest one if none is on one side); entries stay positive */
        int lo[MZ], hi[MZ];
        double tw[MZ];
        for (int j = 0; j < MZ; j++) {
            int a = -1, b = -1;
            for (int i = 0; i < hk[2]; i++) {
                if (z0[i] <= Z[j] && (a < 0 || z0[i] > z0[a])) a = i;
                if (z0[i] >= Z[j] && (b < 0 || z0[i] < z0[b])) b = i;
            }
            if (a < 0) a = b;
            if (b < 0) b = a;
            lo[j] = a; hi[j] = b;
            tw[j] = z0[b] > z0[a] ? (Z[j] - z0[a]) / (z0[b] - z0[a]) : 0.0;
            if (j < NZ) fprintf(stderr, "  z=%.4f starts from the checkpoint vectors at z=%.4f and z=%.4f (weight %.3f on the second)\n",
                                Z[j], z0[a], z0[b], tw[j]);
        }
        for (uint64_t s = 0; s < N; s++)
            for (int j = 0; j < MZ; j++)
                V[s * MZ + j] = (1 - tw[j]) * old[s * MZ + lo[j]] + tw[j] * old[s * MZ + hi[j]];
        free(old);
        fclose(g);
        return;
    }
    double *v = malloc(N * sizeof(double));
    if ((uint64_t)sb.st_size == 8 * N) {
        fseek(g, 0, SEEK_SET);
    } else if ((uint64_t)sb.st_size == 24 + 8 * N) {
        fseek(g, 24, SEEK_SET);
    } else { fprintf(stderr, "start vector has the wrong size\n"); exit(2); }
    if (fread(v, sizeof(double), N, g) != N) { fprintf(stderr, "short start vector\n"); exit(2); }
    fclose(g);
    for (uint64_t s = 0; s < N; s++) {
        if (!(v[s] > 0) || !isfinite(v[s])) { fprintf(stderr, "start vector entry %llu not positive\n", (unsigned long long)s); exit(2); }
        for (int j = 0; j < MZ; j++) V[s * MZ + j] = v[s];
    }
    free(v);
}

static double bound_of(double z, double lam) { return 2.0 * log(K * z) / log(lam); }

static void eigen(const char *dir, const char *start, double tol, int maxit) {
    char path[1024];
    FILE *f = open_kbig(dir, path, sizeof path);
    load_offsets(f);
    fclose(f);
    for (int j = NZ; j < MZ; j++) Z[j] = Z[NZ - 1];     /* unused slots repeat the last z */
    for (int j = 0; j < MZ; j++)
        if (!(Z[j] > 0 && Z[j] < 1.0 / K)) { fprintf(stderr, "z out of range\n"); exit(1); }
    for (int i = 0; i < NW; i++)
        for (int j = 0; j < MZ; j++) W[i][j] = weight_double(&wtypes[i], Z[j]);
    const uint64_t CAP = 1ULL << 27;                    /* records per buffer (1 GiB) */
    uint64_t ccap = 1024;
    CHK = malloc(ccap * sizeof(Chunk));
    NCH = 0;
    for (uint64_t s = 0; s < N;) {
        uint64_t e = s + 1;
        while (e < N && OFF[e + 1] - OFF[s] <= CAP) e++;
        if (OFF[e] - OFF[s] > CAP) { fprintf(stderr, "a state has more than %llu records\n", (unsigned long long)CAP); exit(2); }
        if (NCH == ccap) { ccap *= 2; CHK = realloc(CHK, ccap * sizeof(Chunk)); }
        CHK[NCH++] = (Chunk){s, e};
        s = e;
    }
    uint64_t maxr = 0;
    for (uint64_t c = 0; c < NCH; c++) if (OFF[CHK[c].s1] - OFF[CHK[c].s0] > maxr) maxr = OFF[CHK[c].s1] - OFF[CHK[c].s0];
    if (posix_memalign((void **)&RAW[0], PAGE, maxr * sizeof(Rec) + 2 * PAGE) ||
        posix_memalign((void **)&RAW[1], PAGE, maxr * sizeof(Rec) + 2 * PAGE)) { fprintf(stderr, "out of memory (buffers)\n"); exit(2); }
    RBUF[0] = (Rec *)RAW[0]; RBUF[1] = (Rec *)RAW[1];
    V = malloc(N * MZ * sizeof(double));
    VN = malloc(N * MZ * sizeof(double));
    if (!RBUF[0] || !RBUF[1] || !V || !VN) { fprintf(stderr, "out of memory (eigen)\n"); exit(2); }
    load_start(start);
    {
        struct stat sb;
        if (stat(path, &sb)) { perror(path); exit(2); }
        FILESIZE = (uint64_t)sb.st_size;
    }
    TFD = open(path, O_RDONLY);
    if (TFD < 0) { perror(path); exit(2); }
    fcntl(TFD, F_NOCACHE, 1);
    pthread_t th;
    pthread_create(&th, NULL, reader_main, NULL);
    printf("k=%d h=%d states=%llu records=%llu weight types=%d chunks=%llu, %d values of z, tol %g\n", K, H,
           (unsigned long long)N, (unsigned long long)R, NW, (unsigned long long)NCH, NZ, tol);
    fflush(stdout);
    double rmin[MZ], rmax[MZ], t0 = now();
    int it = 0;
    for (;;) {
        iterate(rmin, rmax);
        it++;
        normalize();
        int conv = 1;
        for (int j = 0; j < NZ; j++) if (rmax[j] / rmin[j] > 1 + tol) conv = 0;
        fprintf(stderr, "  pass %d (%.0fs):", it, now() - t0);
        for (int j = 0; j < NZ; j++)
            fprintf(stderr, " [z=%.4f %.12f %.9f %.1e]", Z[j], rmax[j], bound_of(Z[j], rmax[j]), rmax[j] / rmin[j] - 1);
        fprintf(stderr, "\n");
        if (it % 10 == 0) save_checkpoint(dir);
        if (conv || it >= maxit) break;
    }
    save_checkpoint(dir);
    save_vectors(dir);
    for (int j = 0; j < NZ; j++)
        printf("z=%.4f lambda in [%.15f, %.15f] bound=%.9f passes=%d\n", Z[j], rmin[j], rmax[j],
               bound_of(Z[j], rmax[j]), it);
    printf("eigen time %.0fs\n", now() - t0);
    pthread_mutex_lock(&MU); STOP = 1; pthread_cond_broadcast(&CV); pthread_mutex_unlock(&MU);
}

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: kbig build|klue|project|eigen k h dir ...\n"); return 1; }
    K = atoi(argv[2]);
    H = atoi(argv[3]);
    if (K < 2 || K > 4 || H < 2 || H > HMAX) { fprintf(stderr, "k or h out of range\n"); return 1; }
    M1 = (1ULL << H) - 1;
    M2 = (1ULL << (2 * H)) - 1;
    const char *dir = argv[4];
    if (!strcmp(argv[1], "build")) { build(dir); return 0; }
    if (!strcmp(argv[1], "klue")) { klue(dir); return 0; }
    if (!strcmp(argv[1], "project")) {
        if (argc != 6) { fprintf(stderr, "usage: kbig project k h dir prevdir\n"); return 1; }
        project(dir, argv[5]);
        return 0;
    }
    if (!strcmp(argv[1], "eigen")) {
        if (argc < 9) { fprintf(stderr, "usage: kbig eigen k h dir start tol maxit z1 [... z8]\n"); return 1; }
        NZ = argc - 8;
        if (NZ > MZ) { fprintf(stderr, "at most %d values of z\n", MZ); return 1; }
        for (int j = 0; j < NZ; j++) Z[j] = atof(argv[8 + j]);
        eigen(dir, argv[5], atof(argv[6]), atoi(argv[7]));
        return 0;
    }
    fprintf(stderr, "unknown command %s\n", argv[1]);
    return 1;
}
