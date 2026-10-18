/* Reviewer tool (independent of the reviewed code): re-derive every transition of a
 * trans_h<h>.bin table with my own implementation of Lueker's delta (Fig. 7, with the
 * Java semantics of DPTable.isCanonical: accept iff the canonical backtrack visits (h,h)),
 * and check   row[s][c,l] == 0  <=>  delta rejects,   keys[row-1] == delta(key_s,c,l).
 * Also checks key[1] is the start state and keys are distinct-format (top bits zero).
 * Streams rows with F_NOCACHE; keeps only the key array in memory.
 * usage: fullcheck <trans file> [stride]   (stride>1 checks states s = 1, 1+stride, ...) */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <dispatch/dispatch.h>

static int h;
#define UNK (-(1 << 20))

/* returns 1 and *out on accept, 0 on reject.  Strings are 0-based here: A[0..n-1]. */
static int my_delta(uint64_t key, int c, int l, uint64_t *out) {
    uint32_t m = (1u << h) - 1;
    uint32_t x = (uint32_t)(key >> (3 * h)) & m, y = (uint32_t)(key >> (2 * h)) & m;
    uint32_t dx = (uint32_t)(key >> h) & m, dy = (uint32_t)key & m;
    int A[64], B[64], n = 0, mm = 0;
    for (int t = h - 1; t >= 0; t--) A[n++] = (x >> t) & 1;      /* oldest first */
    for (int t = h - 1; t >= 0; t--) B[mm++] = (y >> t) & 1;
    for (int t = 0; t < l; t++) A[n++] = 1 - c;
    for (int t = 0; t < -l; t++) B[mm++] = 1 - c;
    A[n++] = c; B[mm++] = c;
    /* table V[i][j], i in 0..n, j in 0..mm; corner (h,h) has value 1000 (arbitrary, Lemma 3.21) */
    static __thread int V[64][64];
    for (int i = 0; i <= n; i++) for (int j = 0; j <= mm; j++) V[i][j] = UNK;
    V[h][h] = 1000;
    /* dx: MSB is difference V[1][h]-V[0][h], LSB is V[h][h]-V[h-1][h] */
    for (int i = h; i >= 1; i--) V[i - 1][h] = V[i][h] - (int)((dx >> (h - i)) & 1);
    for (int j = h; j >= 1; j--) V[h][j - 1] = V[h][j] - (int)((dy >> (h - j)) & 1);
    for (int i = 0; i <= n; i++) {
        for (int j = 0; j <= mm; j++) {
            if (i <= h && j <= h) continue;
            int val;
            if (i == 0) val = V[0][j - 1];
            else if (j == 0) val = V[i - 1][0];
            else if (A[i - 1] == B[j - 1]) val = V[i - 1][j - 1] + 1;
            else { int u = V[i - 1][j], w = V[i][j - 1]; val = u > w ? u : w; }
            if (val < -1000) { fprintf(stderr, "read unknown cell\n"); exit(7); }
            V[i][j] = val;
        }
    }
    int i = n, j = mm, diag = 0, hit = 0;
    while (i >= h && j >= h) {
        if (i == h && j == h) { hit = 1; break; }
        if (V[i][j] == V[i - 1][j]) i--;
        else if (V[i][j] == V[i][j - 1]) j--;
        else { if (V[i][j] != V[i - 1][j - 1] + 1) { fprintf(stderr, "inconsistent\n"); exit(7); } i--; j--; diag++; }
    }
    if (!hit) return 0;
    if (diag != 1) { fprintf(stderr, "diag=%d\n", diag); exit(7); }
    uint32_t nx = 0, ny = 0, ndx = 0, ndy = 0;
    for (int t = n - h; t < n; t++) nx = (nx << 1) | (uint32_t)A[t];
    for (int t = mm - h; t < mm; t++) ny = (ny << 1) | (uint32_t)B[t];
    for (int r = n - h + 1; r <= n; r++) { int d = V[r][mm] - V[r - 1][mm]; if (d < 0 || d > 1) exit(8); ndx = (ndx << 1) | (uint32_t)d; }
    for (int q = mm - h + 1; q <= mm; q++) { int d = V[n][q] - V[n][q - 1]; if (d < 0 || d > 1) exit(8); ndy = (ndy << 1) | (uint32_t)d; }
    if (nx & (1u << (h - 1))) { nx ^= m; ny ^= m; }               /* complement symmetry only */
    *out = ((uint64_t)nx << (3 * h)) | ((uint64_t)ny << (2 * h)) | ((uint64_t)ndx << h) | ndy;
    return 1;
}

int main(int argc, char **argv) {
    int fd = open(argv[1], O_RDONLY);
    if (fd < 0) { perror("open"); return 1; }
    fcntl(fd, F_NOCACHE, 1);
    long stride = argc > 2 ? atol(argv[2]) : 1;
    int32_t hdr[2];
    if (pread(fd, hdr, 8, 0) != 8) return 1;
    h = hdr[0]; uint64_t N = (uint64_t)hdr[1];
    int L = 2 * h - 1, NT = 2 * L;
    uint64_t *keys = malloc(N * 8);
    uint64_t koff = 8 + (uint64_t)4 * NT * N, got = 0;
    while (got < N * 8) { ssize_t r = pread(fd, (char *)keys + got, (N * 8 - got) > (1u << 30) ? (1u << 30) : (N * 8 - got), koff + got); if (r <= 0) { perror("pread keys"); return 1; } got += r; }
    uint32_t m = (1u << h) - 1;
    uint64_t start = ((uint64_t)m << h) | m;
    printf("h=%d N=%llu key[1]=%llx start-state=%llx %s\n", h, (unsigned long long)N, (unsigned long long)keys[0], (unsigned long long)start, keys[0] == start ? "OK" : "BAD");
    long badtop = 0;
    for (uint64_t s = 0; s < N; s++) if (keys[s] >> (4 * h)) badtop++;
    const uint64_t CH = 1 << 20;
    uint32_t *buf = malloc(CH * NT * 4);
    __block long mism = 0, checked = 0, acc = 0;
    for (uint64_t base = 1; base <= N; base += CH) {
        uint64_t cnt = (N - base + 1 < CH) ? N - base + 1 : CH;
        uint64_t bytes = cnt * NT * 4, off = 8 + (base - 1) * NT * 4, d = 0;
        while (d < bytes) { ssize_t r = pread(fd, (char *)buf + d, bytes - d, off + d); if (r <= 0) { perror("pread"); return 1; } d += r; }
        const int P = 64;
        long *mm_ = calloc(P, sizeof(long)), *ck = calloc(P, sizeof(long)), *ac = calloc(P, sizeof(long));
        dispatch_apply(P, dispatch_get_global_queue(QOS_CLASS_UTILITY, 0), ^(size_t p) {
            for (uint64_t k = p; k < cnt; k += P) {
                uint64_t s = base + k;
                if ((s - 1) % stride) continue;
                const uint32_t *row = buf + k * NT;
                for (int c = 0; c < 2; c++) for (int l = -(h - 1); l <= h - 1; l++) {
                    uint64_t o; int a = my_delta(keys[s - 1], c, l, &o);
                    uint32_t t = row[c * L + l + h - 1];
                    int okk = a ? (t >= 1 && t <= N && keys[t - 1] == o) : (t == 0);
                    if (!okk) mm_[p]++;
                    if (a) ac[p]++;
                    ck[p]++;
                }
            }
        });
        for (int p = 0; p < P; p++) { mism += mm_[p]; checked += ck[p]; acc += ac[p]; }
        free(mm_); free(ck); free(ac);
    }
    printf("checked %ld transitions (stride %ld): mismatches=%ld accepted=%ld keys-with-high-bits=%ld\n", checked, stride, mism, acc, badtop);
    return 0;
}
