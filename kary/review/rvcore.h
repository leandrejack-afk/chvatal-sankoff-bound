/* rvcore.h: classes, canon and delta, copied verbatim from rvk.c (referee code). */
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

