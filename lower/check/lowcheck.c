/*
 * lowcheck: independent checker for k-ary Chvatal-Sankoff lower-bound certificates.
 *
 * Method: Lueker (J. ACM 2009), Section 2, Lemma 2.3, and its k-letter form in
 * arXiv 2407.18113, Section 3, Facts 9 and 10. If T(w) >= w + r at every pair of
 * strings of length h, then gamma_k >= 2r/(1+r). Witness values are fixed point
 * with unit 2^-27, one little-endian uint32 per class of pairs.
 *
 * Usage:
 *   lowcheck k h witness.bin [threads] [claimed_r]   check a certificate
 *   lowcheck --brute k h witness.bin [threads]       check, then exhaustive cross-checks
 *   lowcheck --count k h [threads]                   class count only, no witness
 *
 * Build: clang -O3 -march=native -lpthread lowcheck.c -o lowcheck
 *
 * Class map. A pair (s, t) is read as z = s.t of length n = 2h. RGS(z) relabels
 * letters by first appearance. CNT[L][m] counts the ways to fill L more entries of
 * an RGS when m labels are already in use (at most k labels in total), so the
 * lexicographic rank of an RGS z is sum over i >= 1 of z[i] * CNT[n-1-i][m_i],
 * with m_i the number of labels in z[0..i-1]: every smaller value at position i is
 * an old label and leaves CNT[n-1-i][m_i] completions. The class rank of (s, t) is
 * min(rank(RGS(s.t)), rank(RGS(t.s))). Bit r of a bitmap over all ranks is set
 * when the RGS of rank r is a class representative, that is
 * r <= rank(RGS(swap(z))). The dense witness position of a class is the number
 * of set bits below its rank: a stored count per 512-bit block plus popcounts
 * inside the block. A lookup costs O(h) arithmetic and three memory reads.
 */

#include <errno.h>
#include <fcntl.h>
#include <inttypes.h>
#include <pthread.h>
#include <stdarg.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

#if !defined(__BYTE_ORDER__) || __BYTE_ORDER__ != __ORDER_LITTLE_ENDIAN__
#error "witness files are little-endian uint32; this build reads them in place"
#endif

#define KMAX 16
#define NMAX 32
#define SHIFT 27
#define ONE ((uint64_t)1 << SHIFT)
#define WORDS_PER_BLOCK 8
#define CHUNK ((uint64_t)1 << 16) /* ranks per work item, a multiple of 512 */
#define PAIR_CHUNK ((uint64_t)1 << 14)
#define BRUTE_MAX_PAIRS ((uint64_t)200000000)
#define MAX_THREADS 256

typedef unsigned __int128 u128;

static int K, H, N;
static uint64_t CNT[NMAX + 1][KMAX + 2];
static uint64_t TOTAL; /* RGS sequences of length N with at most K labels */
static uint64_t *BITS;
static uint64_t *BLOCKCNT;
static uint64_t NWORDS, NBLOCKS, CLASSES;
static const uint32_t *W;
static int64_t *TCLS; /* T per class, kept only in --brute mode */
static int THREADS = 8;
static double T0;
static int HAVE_CLAIM;
static long long CLAIMED;

static const char LETTERS[] = "0123456789abcdefghijklmnopqrstuv";

static void die(const char *fmt, ...)
{
    va_list ap;
    if (HAVE_CLAIM)
        printf("claimed r = %lld: FAIL (stopped on an error before r was computed)\n", CLAIMED);
    fflush(stdout);
    va_start(ap, fmt);
    fputs("lowcheck: ERROR: ", stderr);
    vfprintf(stderr, fmt, ap);
    fputc('\n', stderr);
    va_end(ap);
    exit(2);
}

static double now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static void *xcalloc(size_t count, size_t size, const char *what)
{
    void *p = calloc(count, size);
    if (!p)
        die("cannot allocate %s (%zu x %zu bytes)", what, count, size);
    return p;
}

/* ---------- ranks of restricted growth strings ---------- */

static void build_counts(void)
{
    memset(CNT, 0, sizeof CNT);
    for (int m = 1; m <= K; m++)
        CNT[0][m] = 1;
    for (int L = 1; L < N; L++) {
        int mmax = K < N - L ? K : N - L; /* at most N-L labels are in use with L entries left */
        for (int m = 1; m <= mmax; m++) {
            uint64_t a, b = (m < K) ? CNT[L - 1][m + 1] : 0;
            if (__builtin_mul_overflow((uint64_t)m, CNT[L - 1][m], &a) ||
                __builtin_add_overflow(a, b, &CNT[L][m]))
                die("rank table overflows 64 bits for k=%d h=%d", K, H);
        }
    }
    TOTAL = CNT[N - 1][1];
}

static inline void rank_step(int8_t *lab, int *m, uint64_t *r, int c, int pos)
{
    int l = lab[c];
    if (l < 0) {
        l = *m;
        lab[c] = (int8_t)l;
    }
    *r += (uint64_t)l * CNT[N - 1 - pos][*m];
    if (l == *m)
        (*m)++;
}

/* Rank of RGS(a.b) among all RGS of length N with at most K labels; a, b have H letters. */
static inline uint64_t rank_of(const uint8_t *a, const uint8_t *b)
{
    int8_t lab[KMAX];
    memset(lab, -1, sizeof lab);
    uint64_t r = 0;
    int m = 0;
    for (int i = 0; i < H; i++)
        rank_step(lab, &m, &r, a[i], i);
    for (int i = 0; i < H; i++)
        rank_step(lab, &m, &r, b[i], H + i);
    return r;
}

static void unrank(uint64_t r, uint8_t *z)
{
    z[0] = 0;
    int m = 1;
    for (int i = 1; i < N; i++) {
        uint64_t c = CNT[N - 1 - i][m];
        uint64_t q = r / c;
        if (q < (uint64_t)m) {
            z[i] = (uint8_t)q;
            r -= q * c;
        } else {
            z[i] = (uint8_t)m;
            r -= (uint64_t)m * c;
            m++;
        }
    }
}

static void prefix_max(const uint8_t *z, uint8_t *pm)
{
    pm[0] = z[0];
    for (int i = 1; i < N; i++)
        pm[i] = z[i] > pm[i - 1] ? z[i] : pm[i - 1];
}

/* Advance z to the next RGS in lexicographic order; pm[i] = max(z[0..i]). */
static int next_rgs(uint8_t *z, uint8_t *pm)
{
    for (int i = N - 1; i >= 1; i--) {
        if (z[i] <= pm[i - 1] && z[i] + 1 < K) {
            z[i]++;
            pm[i] = z[i] > pm[i - 1] ? z[i] : pm[i - 1];
            for (int j = i + 1; j < N; j++) {
                z[j] = 0;
                pm[j] = pm[i];
            }
            return 1;
        }
    }
    return 0;
}

static inline int is_canon(uint64_t rank)
{
    return (int)((BITS[rank >> 6] >> (rank & 63)) & 1);
}

static inline uint64_t dense_of(uint64_t rank)
{
    uint64_t wi = rank >> 6, b = rank >> 9;
    uint64_t d = BLOCKCNT[b];
    for (uint64_t j = b * WORDS_PER_BLOCK; j < wi; j++)
        d += (uint64_t)__builtin_popcountll(BITS[j]);
    return d + (uint64_t)__builtin_popcountll(BITS[wi] & ((UINT64_C(1) << (rank & 63)) - 1));
}

static inline uint64_t class_rank(const uint8_t *s, const uint8_t *t)
{
    uint64_t r1 = rank_of(s, t), r2 = rank_of(t, s);
    return r1 < r2 ? r1 : r2;
}

static inline uint64_t wclass(const uint8_t *s, const uint8_t *t)
{
    return W[dense_of(class_rank(s, t))];
}

/* ---------- thread pool over numbered chunks ---------- */

static atomic_uint_fast64_t NEXT;
static uint64_t NCHUNKS;

static uint64_t take_chunk(void)
{
    return (uint64_t)atomic_fetch_add(&NEXT, 1);
}

static void run_parallel(void *(*worker)(void *), void *args, size_t argsize, uint64_t nchunks)
{
    pthread_t tid[MAX_THREADS];
    atomic_store(&NEXT, 0);
    NCHUNKS = nchunks;
    for (int i = 0; i < THREADS; i++) {
        int e = pthread_create(&tid[i], NULL, worker, (char *)args + (size_t)i * argsize);
        if (e)
            die("pthread_create: %s", strerror(e));
    }
    for (int i = 0; i < THREADS; i++)
        pthread_join(tid[i], NULL);
}

/* ---------- canonical bitmap ---------- */

static void *build_worker(void *arg)
{
    (void)arg;
    uint8_t z[NMAX], pm[NMAX];
    for (;;) {
        uint64_t c = take_chunk();
        if (c >= NCHUNKS)
            break;
        uint64_t lo = c * CHUNK, hi = lo + CHUNK < TOTAL ? lo + CHUNK : TOTAL;
        unrank(lo, z);
        prefix_max(z, pm);
        uint64_t word = 0;
        for (uint64_t r = lo; r < hi; r++) {
            if (r <= rank_of(z + H, z))
                word |= UINT64_C(1) << (r & 63);
            if ((r & 63) == 63 || r == hi - 1) {
                BITS[r >> 6] = word;
                word = 0;
            }
            if (r + 1 < hi && !next_rgs(z, pm))
                die("internal: RGS enumeration ended early at rank %" PRIu64, r);
        }
    }
    return NULL;
}

static void build_index(void)
{
    NWORDS = (TOTAL + 63) / 64;
    NBLOCKS = (NWORDS + WORDS_PER_BLOCK - 1) / WORDS_PER_BLOCK;
    BITS = xcalloc(NBLOCKS * WORDS_PER_BLOCK, sizeof(uint64_t), "canonical bitmap");
    BLOCKCNT = xcalloc(NBLOCKS, sizeof(uint64_t), "block counts");
    char dummy[MAX_THREADS];
    run_parallel(build_worker, dummy, 1, (TOTAL + CHUNK - 1) / CHUNK);
    uint64_t acc = 0;
    for (uint64_t b = 0; b < NBLOCKS; b++) {
        BLOCKCNT[b] = acc;
        for (uint64_t j = b * WORDS_PER_BLOCK; j < (b + 1) * WORDS_PER_BLOCK; j++)
            acc += (uint64_t)__builtin_popcountll(BITS[j]);
    }
    CLASSES = acc;
}

/* ---------- independent class count (Burnside) ---------- */

/*
 * Orbits of S_K x {identity, swap} on pairs of strings of length H. A letter
 * permutation with m1 fixed letters fixes m1^(2H) pairs; composed with the swap,
 * (s, t) is fixed iff t = sigma(s) and sigma^2 fixes every letter of s, which gives
 * (m1 + 2*m2)^H pairs, m2 being the number of 2-cycles. Summed over cycle types.
 */
static u128 FACT[KMAX + 1];
static u128 BURN_SUM;
static int BURN_OVERFLOW;

static u128 pow_checked(u128 b, int e)
{
    u128 r = 1;
    while (e-- > 0)
        if (__builtin_mul_overflow(r, b, &r))
            BURN_OVERFLOW = 1;
    return r;
}

static void burn_rec(int rem, int p, int *mult)
{
    if (rem == 0) {
        u128 denom = 1;
        for (int i = 1; i <= K; i++)
            denom *= pow_checked((u128)i, mult[i]) * FACT[mult[i]];
        u128 size = FACT[K] / denom, term = 0, add = 0;
        u128 f = pow_checked((u128)mult[1], 2 * H) + pow_checked((u128)(mult[1] + 2 * mult[2]), H);
        if (__builtin_mul_overflow(size, f, &term) || __builtin_add_overflow(BURN_SUM, term, &add))
            BURN_OVERFLOW = 1;
        BURN_SUM = add;
        return;
    }
    if (p == 0)
        return;
    for (int m = 0; m * p <= rem; m++) {
        mult[p] = m;
        burn_rec(rem - m * p, p - 1, mult);
    }
    mult[p] = 0;
}

/* Returns 1 and the orbit count, or 0 if 128-bit arithmetic would overflow. */
static int burnside(uint64_t *orbits)
{
    int mult[KMAX + 1] = {0};
    FACT[0] = 1;
    for (int i = 1; i <= KMAX; i++)
        FACT[i] = FACT[i - 1] * (u128)i;
    BURN_SUM = 0;
    BURN_OVERFLOW = 0;
    burn_rec(K, K, mult);
    u128 den = 2 * FACT[K];
    if (BURN_OVERFLOW || BURN_SUM % den != 0 || BURN_SUM / den > UINT64_MAX)
        return 0;
    *orbits = (uint64_t)(BURN_SUM / den);
    return 1;
}

/* ---------- the certificate check ---------- */

static int64_t compute_T(const uint8_t *z)
{
    const uint8_t *s = z, *t = z + H;
    uint8_t a[NMAX], b[NMAX];
    if (s[0] == t[0]) {
        memcpy(a, s + 1, (size_t)(H - 1));
        memcpy(b, t + 1, (size_t)(H - 1));
        uint64_t S = 0;
        for (int c = 0; c < K; c++) {
            a[H - 1] = (uint8_t)c;
            for (int c2 = 0; c2 < K; c2++) {
                b[H - 1] = (uint8_t)c2;
                S += wclass(a, b);
            }
        }
        return (int64_t)(ONE + S / (uint64_t)(K * K));
    }
    uint64_t A = 0, B = 0;
    memcpy(a, s + 1, (size_t)(H - 1));
    for (int c = 0; c < K; c++) {
        a[H - 1] = (uint8_t)c;
        A += wclass(a, t);
    }
    memcpy(b, t + 1, (size_t)(H - 1));
    for (int c = 0; c < K; c++) {
        b[H - 1] = (uint8_t)c;
        B += wclass(s, b);
    }
    return (int64_t)((A > B ? A : B) / (uint64_t)K);
}

typedef struct {
    int64_t min_diff, max_diff;
    uint64_t min_rank, min_dense, count;
    char pad[64];
} CheckResult;

static uint64_t PROGRESS_EVERY;

static void *check_worker(void *arg)
{
    CheckResult *res = arg;
    res->min_diff = INT64_MAX;
    res->max_diff = INT64_MIN;
    res->min_rank = UINT64_MAX;
    res->min_dense = 0;
    res->count = 0;
    uint8_t z[NMAX], pm[NMAX];
    for (;;) {
        uint64_t c = take_chunk();
        if (c >= NCHUNKS)
            break;
        if (PROGRESS_EVERY && c && c % PROGRESS_EVERY == 0)
            fprintf(stderr, "lowcheck: %" PRIu64 "%% of ranks dispatched, %.0f s\n",
                    (c * 100 + NCHUNKS / 2) / NCHUNKS, now() - T0);
        uint64_t lo = c * CHUNK, hi = lo + CHUNK < TOTAL ? lo + CHUNK : TOTAL;
        unrank(lo, z);
        prefix_max(z, pm);
        uint64_t d = BLOCKCNT[lo >> 9];
        for (uint64_t r = lo; r < hi; r++) {
            if (is_canon(r)) {
                int64_t T = compute_T(z);
                int64_t diff = T - (int64_t)W[d];
                if (diff < res->min_diff || (diff == res->min_diff && r < res->min_rank)) {
                    res->min_diff = diff;
                    res->min_rank = r;
                    res->min_dense = d;
                }
                if (diff > res->max_diff)
                    res->max_diff = diff;
                if (TCLS)
                    TCLS[d] = T;
                d++;
                res->count++;
            }
            if (r + 1 < hi && !next_rgs(z, pm))
                die("internal: RGS enumeration ended early at rank %" PRIu64, r);
        }
    }
    return NULL;
}

static void run_check(CheckResult *out)
{
    CheckResult *res = xcalloc((size_t)THREADS, sizeof *res, "thread results");
    uint64_t nchunks = (TOTAL + CHUNK - 1) / CHUNK;
    PROGRESS_EVERY = nchunks >= 2000 ? nchunks / 20 : 0;
    run_parallel(check_worker, res, sizeof *res, nchunks);
    *out = res[0];
    for (int i = 1; i < THREADS; i++) {
        CheckResult *x = &res[i];
        if (x->min_diff < out->min_diff || (x->min_diff == out->min_diff && x->min_rank < out->min_rank)) {
            out->min_diff = x->min_diff;
            out->min_rank = x->min_rank;
            out->min_dense = x->min_dense;
        }
        if (x->max_diff > out->max_diff)
            out->max_diff = x->max_diff;
        out->count += x->count;
    }
    free(res);
}

/* ---------- --brute: exhaustive cross-checks for small k^(2h) ---------- */

static uint64_t POWK[NMAX + 1];
static uint64_t NPAIRS;
static uint8_t PERMS[200][KMAX];
static int NPERMS, FULL_GROUP;
static uint32_t *WF;
static _Atomic uint8_t *HIT;
static _Atomic uint64_t *ORBIT_OF;

static void decode_pair(uint64_t p, uint8_t *s, uint8_t *t)
{
    uint64_t sv = p / POWK[H], tv = p % POWK[H];
    for (int i = H - 1; i >= 0; i--) {
        s[i] = (uint8_t)(sv % (uint64_t)K);
        sv /= (uint64_t)K;
        t[i] = (uint8_t)(tv % (uint64_t)K);
        tv /= (uint64_t)K;
    }
}

static uint64_t encode_pair(const uint8_t *s, const uint8_t *t)
{
    uint64_t sv = 0, tv = 0;
    for (int i = 0; i < H; i++) {
        sv = sv * (uint64_t)K + s[i];
        tv = tv * (uint64_t)K + t[i];
    }
    return sv * POWK[H] + tv;
}

/* T at pair index p straight from the definition, on the full expanded vector WF. */
static int64_t brute_T(uint64_t p)
{
    uint64_t kh = POWK[H], kh1 = POWK[H - 1], k = (uint64_t)K;
    uint64_t sv = p / kh, tv = p % kh;
    uint64_t s1 = (sv % kh1) * k, t1 = (tv % kh1) * k;
    if (sv / kh1 == tv / kh1) {
        uint64_t S = 0;
        for (uint64_t c = 0; c < k; c++)
            for (uint64_t c2 = 0; c2 < k; c2++)
                S += WF[(s1 + c) * kh + t1 + c2];
        return (int64_t)(ONE + S / (k * k));
    }
    uint64_t A = 0, B = 0;
    for (uint64_t c = 0; c < k; c++) {
        A += WF[(s1 + c) * kh + tv];
        B += WF[sv * kh + t1 + c];
    }
    return (int64_t)((A > B ? A : B) / k);
}

typedef struct {
    uint64_t bad_index, bad_swap, bad_perm, bad_orbit, orbit_reps;
    char pad[64];
} MapStats;

static void *map_worker(void *arg)
{
    MapStats *st = arg;
    uint8_t s[NMAX], t[NMAX], ps[NMAX], pt[NMAX];
    for (;;) {
        uint64_t c = take_chunk();
        if (c >= NCHUNKS)
            break;
        uint64_t lo = c * PAIR_CHUNK, hi = lo + PAIR_CHUNK < NPAIRS ? lo + PAIR_CHUNK : NPAIRS;
        for (uint64_t p = lo; p < hi; p++) {
            decode_pair(p, s, t);
            uint64_t cr = class_rank(s, t);
            if (cr >= TOTAL || !is_canon(cr)) {
                st->bad_index++;
                continue;
            }
            uint64_t ci = dense_of(cr);
            if (ci >= CLASSES) {
                st->bad_index++;
                continue;
            }
            atomic_store_explicit(&HIT[ci], 1, memory_order_relaxed);
            WF[p] = W[ci];
            if (dense_of(class_rank(t, s)) != ci)
                st->bad_swap++;
            uint64_t canon = UINT64_MAX;
            for (int j = 0; j < NPERMS; j++) {
                for (int i = 0; i < H; i++) {
                    ps[i] = PERMS[j][s[i]];
                    pt[i] = PERMS[j][t[i]];
                }
                if (dense_of(class_rank(ps, pt)) != ci)
                    st->bad_perm++;
                if (FULL_GROUP) {
                    uint64_t e1 = encode_pair(ps, pt), e2 = encode_pair(pt, ps);
                    if (e1 < canon)
                        canon = e1;
                    if (e2 < canon)
                        canon = e2;
                }
            }
            if (FULL_GROUP) {
                if (canon == p)
                    st->orbit_reps++;
                uint64_t expect = UINT64_MAX;
                if (!atomic_compare_exchange_strong(&ORBIT_OF[ci], &expect, canon) && expect != canon)
                    st->bad_orbit++;
            }
        }
    }
    return NULL;
}

typedef struct {
    int64_t min_d, max_d;
    uint64_t bad_T, first_bad;
    char pad[64];
} PairStats;

static void *pair_worker(void *arg)
{
    PairStats *st = arg;
    st->min_d = INT64_MAX;
    st->max_d = INT64_MIN;
    st->bad_T = 0;
    st->first_bad = UINT64_MAX;
    uint8_t s[NMAX], t[NMAX];
    for (;;) {
        uint64_t c = take_chunk();
        if (c >= NCHUNKS)
            break;
        uint64_t lo = c * PAIR_CHUNK, hi = lo + PAIR_CHUNK < NPAIRS ? lo + PAIR_CHUNK : NPAIRS;
        for (uint64_t p = lo; p < hi; p++) {
            int64_t T = brute_T(p);
            int64_t d = T - (int64_t)WF[p];
            if (d < st->min_d)
                st->min_d = d;
            if (d > st->max_d)
                st->max_d = d;
            decode_pair(p, s, t);
            if (T != TCLS[dense_of(class_rank(s, t))]) {
                if (p < st->first_bad)
                    st->first_bad = p;
                st->bad_T++;
            }
        }
    }
    return NULL;
}

static int next_perm(uint8_t *a, int n)
{
    int i = n - 2;
    while (i >= 0 && a[i] >= a[i + 1])
        i--;
    if (i < 0)
        return 0;
    int j = n - 1;
    while (a[j] <= a[i])
        j--;
    uint8_t x = a[i];
    a[i] = a[j];
    a[j] = x;
    for (int l = i + 1, r = n - 1; l < r; l++, r--) {
        x = a[l];
        a[l] = a[r];
        a[r] = x;
    }
    return 1;
}

static uint64_t splitmix(uint64_t *x)
{
    uint64_t z = (*x += UINT64_C(0x9E3779B97F4A7C15));
    z = (z ^ (z >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94D049BB133111EB);
    return z ^ (z >> 31);
}

static const char *verdict(int ok)
{
    return ok ? "PASS" : "FAIL";
}

static int brute(const CheckResult *cr)
{
    POWK[0] = 1;
    for (int i = 1; i <= N; i++) {
        if (POWK[i - 1] > BRUTE_MAX_PAIRS)
            die("--brute needs k^(2h) <= %" PRIu64, BRUTE_MAX_PAIRS);
        POWK[i] = POWK[i - 1] * (uint64_t)K;
    }
    NPAIRS = POWK[N];
    if (NPAIRS > BRUTE_MAX_PAIRS)
        die("--brute needs k^(2h) <= %" PRIu64 " (here %" PRIu64 ")", BRUTE_MAX_PAIRS, NPAIRS);
    int all_ok = 1;
    printf("\n--brute: exhaustive cross-checks over all %" PRIu64 " pairs (k^(2h))\n", NPAIRS);

    /* 1. rank = position in lexicographic order of all RGS, found by scanning every
          string of length 2h in base-k order (which is lexicographic order). */
    {
        uint8_t d[NMAX] = {0}, it[NMAX], itpm[NMAX], u[NMAX];
        unrank(0, it);
        prefix_max(it, itpm);
        uint64_t idx = 0, bad = 0;
        int more = 1;
        for (uint64_t x = 0; x < NPAIRS; x++) {
            int rgs = d[0] == 0, mx = 0;
            for (int i = 1; i < N && rgs; i++) {
                if (d[i] > mx + 1)
                    rgs = 0;
                else if (d[i] > mx)
                    mx = d[i];
            }
            if (rgs) {
                unrank(idx, u);
                if (rank_of(d, d + H) != idx || memcmp(u, d, (size_t)N) != 0 || !more ||
                    memcmp(it, d, (size_t)N) != 0)
                    bad++;
                more = next_rgs(it, itpm);
                idx++;
            }
            for (int i = N - 1; i >= 0; i--) {
                if (++d[i] < K)
                    break;
                d[i] = 0;
            }
        }
        int ok = bad == 0 && idx == TOTAL && !more;
        all_ok &= ok;
        printf("rank test: %" PRIu64 " RGS found by scanning all strings (table says %" PRIu64
               "); rank, unrank and successor disagreements: %" PRIu64 "  %s\n",
               idx, TOTAL, bad, verdict(ok));
    }

    /* 2. class map: every pair lands on a set bit, permutations and the swap keep the
          class, every class is hit, and (k <= 5) class == orbit of S_k x swap. */
    {
        uint8_t id[KMAX];
        for (int i = 0; i < K; i++)
            id[i] = (uint8_t)i;
        if (K <= 5) {
            FULL_GROUP = 1;
            NPERMS = 0;
            do
                memcpy(PERMS[NPERMS++], id, KMAX);
            while (next_perm(id, K));
        } else {
            FULL_GROUP = 0;
            NPERMS = 0;
            for (int i = 0; i + 1 < K; i++) {
                for (int j = 0; j < K; j++)
                    PERMS[NPERMS][j] = (uint8_t)j;
                PERMS[NPERMS][i] = (uint8_t)(i + 1);
                PERMS[NPERMS][i + 1] = (uint8_t)i;
                NPERMS++;
            }
            uint64_t seed = 12345;
            for (int q = 0; q < 8; q++, NPERMS++) {
                for (int j = 0; j < K; j++)
                    PERMS[NPERMS][j] = (uint8_t)j;
                for (int j = K - 1; j > 0; j--) {
                    int r = (int)(splitmix(&seed) % (uint64_t)(j + 1));
                    uint8_t x = PERMS[NPERMS][j];
                    PERMS[NPERMS][j] = PERMS[NPERMS][r];
                    PERMS[NPERMS][r] = x;
                }
            }
        }
        WF = xcalloc(NPAIRS, sizeof *WF, "expanded witness");
        HIT = xcalloc(CLASSES, sizeof *HIT, "class hit flags");
        if (FULL_GROUP) {
            ORBIT_OF = xcalloc(CLASSES, sizeof *ORBIT_OF, "orbit table");
            for (uint64_t i = 0; i < CLASSES; i++)
                atomic_store(&ORBIT_OF[i], UINT64_MAX);
        }
        MapStats *st = xcalloc((size_t)THREADS, sizeof *st, "thread stats");
        run_parallel(map_worker, st, sizeof *st, (NPAIRS + PAIR_CHUNK - 1) / PAIR_CHUNK);
        MapStats sum = {0};
        for (int i = 0; i < THREADS; i++) {
            sum.bad_index += st[i].bad_index;
            sum.bad_swap += st[i].bad_swap;
            sum.bad_perm += st[i].bad_perm;
            sum.bad_orbit += st[i].bad_orbit;
            sum.orbit_reps += st[i].orbit_reps;
        }
        free(st);
        uint64_t unhit = 0;
        for (uint64_t i = 0; i < CLASSES; i++)
            unhit += atomic_load(&HIT[i]) == 0;
        int ok = sum.bad_index == 0;
        all_ok &= ok;
        printf("class index: pairs whose class rank is not a set bit or out of range: %" PRIu64 "  %s\n",
               sum.bad_index, verdict(ok));
        ok = sum.bad_swap == 0;
        all_ok &= ok;
        printf("swap: pairs whose swap lands in another class: %" PRIu64 "  %s\n", sum.bad_swap, verdict(ok));
        ok = sum.bad_perm == 0;
        all_ok &= ok;
        if (FULL_GROUP)
            printf("permutations: all %d of S_%d on every pair, class changes: %" PRIu64 "  %s\n",
                   NPERMS, K, sum.bad_perm, verdict(ok));
        else
            printf("permutations: %d adjacent transpositions (they generate S_%d) and 8 random ones on every"
                   " pair, class changes: %" PRIu64 "  %s\n",
                   K - 1, K, sum.bad_perm, verdict(ok));
        ok = unhit == 0;
        all_ok &= ok;
        printf("surjective: classes hit by no pair: %" PRIu64 " of %" PRIu64 "  %s\n", unhit, CLASSES, verdict(ok));
        if (FULL_GROUP) {
            ok = sum.bad_orbit == 0 && sum.orbit_reps == CLASSES;
            all_ok &= ok;
            printf("orbits: brute-force orbits of S_%d x swap = %" PRIu64 ", classes = %" PRIu64
                   ", classes spanning two orbits = %" PRIu64 "  %s\n",
                   K, sum.orbit_reps, CLASSES, sum.bad_orbit, verdict(ok));
            free((void *)ORBIT_OF);
        }
        free((void *)HIT);
    }

    /* 3. T at every pair from the definition, on the expanded vector. */
    {
        PairStats *st = xcalloc((size_t)THREADS, sizeof *st, "thread stats");
        run_parallel(pair_worker, st, sizeof *st, (NPAIRS + PAIR_CHUNK - 1) / PAIR_CHUNK);
        int64_t mn = INT64_MAX, mx = INT64_MIN;
        uint64_t bad = 0, first = UINT64_MAX;
        for (int i = 0; i < THREADS; i++) {
            if (st[i].min_d < mn)
                mn = st[i].min_d;
            if (st[i].max_d > mx)
                mx = st[i].max_d;
            bad += st[i].bad_T;
            if (st[i].first_bad < first)
                first = st[i].first_bad;
        }
        free(st);
        int ok = bad == 0;
        all_ok &= ok;
        printf("T per pair: pairs whose T differs from their class representative's T: %" PRIu64 "  %s\n",
               bad, verdict(ok));
        if (bad) {
            uint8_t s[NMAX], t[NMAX];
            decode_pair(first, s, t);
            printf("  first such pair: s=");
            for (int i = 0; i < H; i++)
                putchar(LETTERS[s[i]]);
            printf(" t=");
            for (int i = 0; i < H; i++)
                putchar(LETTERS[t[i]]);
            putchar('\n');
        }
        ok = mn == cr->min_diff;
        all_ok &= ok;
        printf("min over all pairs of T - w = %" PRId64 ", class-based r = %" PRId64 "  %s\n", mn, cr->min_diff,
               verdict(ok));
        ok = mx == cr->max_diff;
        all_ok &= ok;
        printf("max over all pairs of T - w = %" PRId64 ", class-based max = %" PRId64 "  %s\n", mx, cr->max_diff,
               verdict(ok));
    }
    free(WF);
    printf("--brute overall: %s\n", verdict(all_ok));
    return all_ok;
}

/* ---------- driver ---------- */

static long parse_long(const char *s, const char *what, long lo, long hi)
{
    char *end;
    errno = 0;
    long v = strtol(s, &end, 10);
    if (errno || *end || end == s || v < lo || v > hi)
        die("bad %s '%s' (expected an integer in [%ld, %ld])", what, s, lo, hi);
    return v;
}

static void usage(void)
{
    fputs("usage: lowcheck k h witness.bin [threads] [claimed_r]\n"
          "       lowcheck --brute k h witness.bin [threads]\n"
          "       lowcheck --count k h [threads]\n",
          stderr);
    exit(2);
}

static void print_pair_of_rank(uint64_t rank)
{
    uint8_t z[NMAX];
    unrank(rank, z);
    printf("s=");
    for (int i = 0; i < H; i++)
        putchar(LETTERS[z[i]]);
    printf(" t=");
    for (int i = H; i < N; i++)
        putchar(LETTERS[z[i]]);
}

static void setup(int k, int h)
{
    K = k;
    H = h;
    N = 2 * h;
    build_counts();
    double t = now();
    build_index();
    double mb = (double)(NBLOCKS * WORDS_PER_BLOCK * 8 + NBLOCKS * 8) / 1048576.0;
    printf("ranks: %" PRIu64 " RGS of length %d with at most %d labels\n", TOTAL, N, K);
    printf("index memory: %.1f MB (bitmap and block counts), built in %.2f s\n", mb, now() - t);
    uint64_t orbits;
    if (burnside(&orbits)) {
        printf("classes computed: %" PRIu64 " (Burnside orbit count of S_%d x swap: %" PRIu64 ", %s)\n", CLASSES, K,
               orbits, orbits == CLASSES ? "agrees" : "DISAGREES");
        if (orbits != CLASSES)
            die("class count %" PRIu64 " disagrees with the Burnside count %" PRIu64, CLASSES, orbits);
    } else {
        printf("classes computed: %" PRIu64 " (Burnside count skipped: exceeds 128 bits)\n", CLASSES);
    }
}

int main(int argc, char **argv)
{
    int mode_brute = 0, mode_count = 0, ai = 1;
    if (argc > 1 && strcmp(argv[1], "--brute") == 0) {
        mode_brute = 1;
        ai = 2;
    } else if (argc > 1 && strcmp(argv[1], "--count") == 0) {
        mode_count = 1;
        ai = 2;
    }
    int nargs = argc - ai;
    if (mode_count ? (nargs < 2 || nargs > 3) : mode_brute ? (nargs < 3 || nargs > 4) : (nargs < 3 || nargs > 5))
        usage();
    int k = (int)parse_long(argv[ai], "k", 2, KMAX);
    int h = (int)parse_long(argv[ai + 1], "h", 1, NMAX / 2);
    const char *path = mode_count ? NULL : argv[ai + 2];
    int ti = mode_count ? ai + 2 : ai + 3;
    if (argc > ti)
        THREADS = (int)parse_long(argv[ti], "threads", 1, MAX_THREADS);
    if (!mode_count && !mode_brute && argc > ti + 1) {
        char *end;
        errno = 0;
        long long v = strtoll(argv[ti + 1], &end, 10);
        if (errno || *end || end == argv[ti + 1])
            die("bad claimed_r '%s'", argv[ti + 1]);
        CLAIMED = v;
        HAVE_CLAIM = 1;
    }

    T0 = now();
    printf("lowcheck k=%d h=%d%s%s threads=%d\n", k, h, path ? " witness=" : "", path ? path : "", THREADS);
    setup(k, h);
    if (mode_count)
        return 0;

    int fd = open(path, O_RDONLY);
    if (fd < 0)
        die("cannot open %s: %s", path, strerror(errno));
    struct stat sb;
    if (fstat(fd, &sb) != 0)
        die("cannot stat %s: %s", path, strerror(errno));
    uint64_t bytes = (uint64_t)sb.st_size;
    if (bytes % 4 != 0 || bytes / 4 != CLASSES)
        die("witness %s holds %" PRIu64 " bytes (%" PRIu64 " values%s), but there are %" PRIu64
            " classes for k=%d h=%d",
            path, bytes, bytes / 4, bytes % 4 ? " plus a partial value" : "", CLASSES, k, h);
    printf("witness values: %" PRIu64 " (equals the class count)\n", bytes / 4);
    void *map = mmap(NULL, (size_t)bytes, PROT_READ, MAP_SHARED, fd, 0);
    if (map == MAP_FAILED)
        die("cannot mmap %s: %s", path, strerror(errno));
    W = map;
    if (mode_brute)
        TCLS = xcalloc(CLASSES, sizeof *TCLS, "per-class T");

    double t = now();
    CheckResult res;
    run_check(&res);
    double tc = now() - t;
    if (res.count != CLASSES)
        die("internal: checked %" PRIu64 " classes, expected %" PRIu64, res.count, CLASSES);

    int64_t r = res.min_diff;
    printf("r = %" PRId64 " (units of 2^-27, min over classes of T - w)\n", r);
    printf("min of T - w at class %" PRIu64 " (rank %" PRIu64 "): ", res.min_dense, res.min_rank);
    print_pair_of_rank(res.min_rank);
    printf(" w=%" PRIu32 "\n", W[res.min_dense]);
    printf("max of T - w = %" PRId64 "\n", res.max_diff);
    if (r > 0) {
        u128 q = (u128)2 * (u128)r * 1000000000u / ((u128)ONE + (u128)r);
        printf("bound: gamma_%d >= %" PRIu64 ".%09" PRIu64 " (2r/(2^27 + r) rounded down to 9 decimals)\n", k,
               (uint64_t)(q / 1000000000u), (uint64_t)(q % 1000000000u));
    } else {
        printf("bound: none (r <= 0, no bound follows)\n");
    }
    printf("time: check %.2f s, total %.2f s\n", tc, now() - T0);

    int rc = 0;
    if (mode_brute)
        rc = brute(&res) ? 0 : 1;
    if (HAVE_CLAIM) {
        int pass = r >= CLAIMED;
        printf("claimed r = %lld: %s (computed r - claimed r = %lld)\n", CLAIMED, pass ? "PASS" : "FAIL",
               (long long)r - CLAIMED);
        rc = pass ? 0 : 1;
    }
    munmap(map, (size_t)bytes);
    close(fd);
    return rc;
}
