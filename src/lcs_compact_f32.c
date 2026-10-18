/*
 * Compact transition table for M_h and a double-buffered streaming power
 * iteration.  For every accepting state and each character c the accepted
 * paddings l form one interval [lo_c, hi_c] containing 0 (checked during the
 * conversion), so a state is stored as a 32-bit header
 *   (lo_0 + h - 1) | (hi_0 + h - 1) << 5 | (lo_1 + h - 1) << 10 | (hi_1 + h - 1) << 15
 * followed by its accepted targets in (c, l) order.
 *
 * Usage:
 *   lcs_compact convert <h> <dir>       trans_h<h>.bin -> ctrans_h<h>.{hdr,tgt,idx}
 *   lcs_compact eigen <h> <dir> <start vector or -> <z> [<z> ...]
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
#include <dispatch/dispatch.h>

static int H, NT, LL;
static uint64_t N;
static const uint64_t BLK = 1 << 20;      /* states per index block */

static double now(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

static void convert(const char *dir) {
    char p[1024];
    snprintf(p, sizeof p, "%s/trans_h%d.bin", dir, H);
    FILE *in = fopen(p, "rb");
    int32_t hdr[2];
    if (!in || fread(hdr, sizeof hdr, 1, in) != 1 || hdr[0] != H) { fprintf(stderr, "bad input\n"); exit(2); }
    N = (uint64_t)hdr[1];
    snprintf(p, sizeof p, "%s/ctrans_h%d.hdr", dir, H); FILE *fh = fopen(p, "wb");
    snprintf(p, sizeof p, "%s/ctrans_h%d.tgt", dir, H); FILE *ft = fopen(p, "wb");
    snprintf(p, sizeof p, "%s/ctrans_h%d.idx", dir, H); FILE *fi = fopen(p, "wb");
    fwrite(&N, sizeof N, 1, fh);
    const uint64_t CH = 1 << 22;
    uint32_t *buf = malloc(CH * NT * sizeof(uint32_t));
    uint32_t *hb = malloc(CH * sizeof(uint32_t));
    uint32_t *tb = malloc(CH * NT * sizeof(uint32_t));
    uint64_t off = 0;
    double t0 = now();
    for (uint64_t base = 0; base < N; base += CH) {
        uint64_t cnt = (N - base < CH) ? N - base : CH;
        if (fread(buf, sizeof(uint32_t), cnt * NT, in) != cnt * NT) { fprintf(stderr, "short read\n"); exit(2); }
        uint64_t nt = 0;
        for (uint64_t k = 0; k < cnt; k++) {
            uint64_t s = base + k;                       /* 0-based state index */
            if (s % BLK == 0) { uint64_t o = off + nt; fwrite(&o, sizeof o, 1, fi); }
            const uint32_t *r = buf + k * NT;
            uint32_t hd = 0;
            for (int c = 0; c < 2; c++) {
                int lo = 0, hi = 0;
                while (lo - 1 >= -(H - 1) && r[c * LL + (lo - 1) + H - 1]) lo--;
                while (hi + 1 <= H - 1 && r[c * LL + (hi + 1) + H - 1]) hi++;
                if (!r[c * LL + H - 1]) { fprintf(stderr, "l=0 rejected at state %llu\n", (unsigned long long)s + 1); exit(3); }
                for (int l = -(H - 1); l <= H - 1; l++) {
                    int inside = (l >= lo && l <= hi);
                    if (inside != (r[c * LL + l + H - 1] != 0)) { fprintf(stderr, "non-contiguous range at state %llu\n", (unsigned long long)s + 1); exit(3); }
                }
                hd |= (uint32_t)(lo + H - 1) << (10 * c);
                hd |= (uint32_t)(hi + H - 1) << (10 * c + 5);
                for (int l = lo; l <= hi; l++) tb[nt++] = r[c * LL + l + H - 1];
            }
            hb[k] = hd;
        }
        fwrite(hb, sizeof(uint32_t), cnt, fh);
        fwrite(tb, sizeof(uint32_t), nt, ft);
        off += nt;
    }
    fwrite(&off, sizeof off, 1, fi);   /* sentinel: total number of targets */
    fclose(in); fclose(fh); fclose(ft); fclose(fi);
    printf("h=%d N=%llu targets=%llu (%.2f per state) in %.0fs\n", H, (unsigned long long)N,
           (unsigned long long)off, (double)off / N, now() - t0);
}

/* ---------------- streaming eigen ------------------------------------- */
static uint32_t *HD;          /* headers, 0-based */
static uint64_t *IDX, NB;     /* NB = number of index blocks */
static float *V, *VN;   /* float storage halves memory; sums are in double */
static int tfd;

typedef struct {
    uint32_t *buf[2];
    uint64_t cap;
    uint64_t blk_lo[2], blk_hi[2];   /* block range held in each buffer */
    int ready[2];
    int stop;
    uint64_t next_blk;
    pthread_mutex_t mu;
    pthread_cond_t cv;
} Reader;
static Reader R;
static uint64_t BLOCKS_PER_BUF;

static void *reader_main(void *arg) {
    (void)arg;
    int which = 0;
    for (;;) {
        pthread_mutex_lock(&R.mu);
        while (R.ready[which] && !R.stop) pthread_cond_wait(&R.cv, &R.mu);
        if (R.stop) { pthread_mutex_unlock(&R.mu); return NULL; }
        uint64_t b0 = R.next_blk;
        uint64_t b1 = b0 + BLOCKS_PER_BUF;
        if (b1 > NB) b1 = NB;
        R.next_blk = (b1 >= NB) ? 0 : b1;          /* wrap around for the next iteration */
        pthread_mutex_unlock(&R.mu);
        uint64_t o0 = IDX[b0], o1 = IDX[b1];
        uint64_t bytes = (o1 - o0) * sizeof(uint32_t);
        uint64_t done = 0;
        while (done < bytes) {
            uint64_t want = bytes - done; if (want > (1ULL << 30)) want = 1ULL << 30;   /* macOS caps one read at INT_MAX */
            ssize_t got = pread(tfd, (char *)R.buf[which] + done, want, (off_t)(o0 * sizeof(uint32_t) + done));
            if (got <= 0) { perror("pread"); exit(2); }
            done += (uint64_t)got;
        }
        pthread_mutex_lock(&R.mu);
        R.blk_lo[which] = b0; R.blk_hi[which] = b1; R.ready[which] = 1;
        pthread_cond_broadcast(&R.cv);
        pthread_mutex_unlock(&R.mu);
        which ^= 1;
    }
}

static void iterate(const double *w, double *rmin, double *rmax) {
    static int which = 0;
    double mn = INFINITY, mx = 0;
    uint64_t blocks_done = 0;
    while (blocks_done < NB) {
        pthread_mutex_lock(&R.mu);
        while (!R.ready[which]) pthread_cond_wait(&R.cv, &R.mu);
        pthread_mutex_unlock(&R.mu);
        uint64_t b0 = R.blk_lo[which], b1 = R.blk_hi[which];
        const uint32_t *buf = R.buf[which];
        uint64_t nblk = b1 - b0;
        double *mins = malloc(nblk * sizeof(double)), *maxs = malloc(nblk * sizeof(double));
        dispatch_apply(nblk, dispatch_get_global_queue(QOS_CLASS_USER_INITIATED, 0), ^(size_t bi) {
            uint64_t b = b0 + bi;
            const uint32_t *t = buf + (IDX[b] - IDX[b0]);
            uint64_t s0 = b * BLK, s1 = s0 + BLK;
            if (s1 > N) s1 = N;
            double a = INFINITY, c2 = 0;
            for (uint64_t s = s0; s < s1; s++) {
                uint32_t hd = HD[s];
                double sum = 0;
                for (int c = 0; c < 2; c++) {
                    int lo = (int)((hd >> (10 * c)) & 31) - (H - 1);
                    int hi = (int)((hd >> (10 * c + 5)) & 31) - (H - 1);
                    for (int l = lo; l <= hi; l++) sum += w[l + H - 1] * V[*t++];
                }
                VN[s + 1] = (float)sum;
                double q = sum / (double)V[s + 1];
                if (q < a) a = q;
                if (q > c2) c2 = q;
            }
            mins[bi] = a; maxs[bi] = c2;
        });
        for (uint64_t i = 0; i < nblk; i++) { if (mins[i] < mn) mn = mins[i]; if (maxs[i] > mx) mx = maxs[i]; }
        free(mins); free(maxs);
        blocks_done += nblk;
        pthread_mutex_lock(&R.mu);
        R.ready[which] = 0;
        pthread_cond_broadcast(&R.cv);
        pthread_mutex_unlock(&R.mu);
        which ^= 1;
    }
    *rmin = mn; *rmax = mx;
}

int main(int argc, char **argv) {
    if (argc < 4) return 1;
    H = atoi(argv[2]); LL = 2 * H - 1; NT = 2 * LL;
    const char *dir = argv[3];
    if (!strcmp(argv[1], "convert")) { convert(dir); return 0; }
    char p[1024];
    snprintf(p, sizeof p, "%s/ctrans_h%d.hdr", dir, H);
    FILE *fh = fopen(p, "rb");
    if (!fh || fread(&N, sizeof N, 1, fh) != 1) { fprintf(stderr, "no header file\n"); return 2; }
    HD = malloc(N * sizeof(uint32_t));
    if (fread(HD, sizeof(uint32_t), N, fh) != N) return 2;
    fclose(fh);
    NB = (N + BLK - 1) / BLK;
    IDX = malloc((NB + 1) * sizeof(uint64_t));
    snprintf(p, sizeof p, "%s/ctrans_h%d.idx", dir, H);
    FILE *fi = fopen(p, "rb");
    if (!fi || fread(IDX, sizeof(uint64_t), NB + 1, fi) != NB + 1) { fprintf(stderr, "bad index\n"); return 2; }
    fclose(fi);
    snprintf(p, sizeof p, "%s/ctrans_h%d.tgt", dir, H);
    tfd = open(p, O_RDONLY);
    fcntl(tfd, F_NOCACHE, 1);
    BLOCKS_PER_BUF = 16;                               /* about 64M states per buffer */
    uint64_t maxb = 0;
    for (uint64_t b = 0; b < NB; b += BLOCKS_PER_BUF) {
        uint64_t e = b + BLOCKS_PER_BUF; if (e > NB) e = NB;
        if (IDX[e] - IDX[b] > maxb) maxb = IDX[e] - IDX[b];
    }
    R.cap = maxb;
    R.buf[0] = malloc(maxb * sizeof(uint32_t));
    R.buf[1] = malloc(maxb * sizeof(uint32_t));
    pthread_mutex_init(&R.mu, NULL); pthread_cond_init(&R.cv, NULL);
    V = malloc((N + 1) * sizeof(float)); VN = malloc((N + 1) * sizeof(float));
    V[0] = VN[0] = 0;
    double *tmpd = NULL;
    if (strcmp(argv[4], "-")) {
        FILE *g = fopen(argv[4], "rb");
        int32_t eh; uint64_t en; double ez;
        /* accept either a raw vector or an eigen_h<h>.bin file */
        tmpd = malloc(N * sizeof(double));
        if (fread(&eh, 4, 1, g) == 1 && eh == H && fread(&en, 8, 1, g) == 1 && en == N && fread(&ez, 8, 1, g) == 1) {
            if (fread(tmpd, sizeof(double), N, g) != N) return 2;
        } else {
            fseek(g, 0, SEEK_SET);
            if (fread(tmpd, sizeof(double), N, g) != N) return 2;
        }
        for (uint64_t s = 0; s < N; s++) V[s + 1] = (float)tmpd[s];
        free(tmpd);
        fclose(g);
    } else {
        for (uint64_t s = 1; s <= N; s++) V[s] = 1.0;
    }
    pthread_t th;
    pthread_create(&th, NULL, reader_main, NULL);
    printf("h=%d states=%llu blocks=%llu buffer=%.1f GB\n", H, (unsigned long long)N, (unsigned long long)NB, maxb * 4.0 / 1e9);
    fflush(stdout);
    for (int ai = 5; ai < argc; ai++) {
        double z = atof(argv[ai]);
        double w[64];
        for (int l = -(H - 1); l <= H - 1; l++) {
            int al = l < 0 ? -l : l;
            w[l + H - 1] = (al == H - 1) ? pow(z, H + 1) / (1.0 - z) : pow(z, al + 2);
        }
        double rmin, rmax, t0 = now();
        double best = INFINITY;
        int it = 0, since_best = 0;
        for (;;) {
            iterate(w, &rmin, &rmax);
            it++;
            /* rmax certifies the vector V that was just multiplied; keep the best one */
            if (rmax < best) {
                best = rmax; since_best = 0;
                if (it >= 60) {
                    snprintf(p, sizeof p, "%s/eigen_h%d.bin", dir, H);
                    FILE *fb = fopen(p, "wb");
                    int32_t hh0 = H; uint64_t n0 = N;
                    fwrite(&hh0, 4, 1, fb); fwrite(&n0, 8, 1, fb); fwrite(&z, 8, 1, fb);
                    double *row = malloc((1 << 20) * sizeof(double));
                    for (uint64_t s = 1; s <= N; s += (1 << 20)) {
                        uint64_t c = (N - s + 1 < (1 << 20)) ? N - s + 1 : (1 << 20);
                        for (uint64_t q = 0; q < c; q++) row[q] = V[s + q];
                        fwrite(row, sizeof(double), c, fb);
                    }
                    free(row);
                    fclose(fb);
                    fprintf(stderr, "    saved vector with max ratio %.15f\n", rmax);
                }
            } else since_best++;
            float mx = 0;
            for (uint64_t s = 1; s <= N; s++) if (VN[s] > mx) mx = VN[s];
            for (uint64_t s = 1; s <= N; s++) V[s] = VN[s] / mx;
            fprintf(stderr, "    z=%.6f it=%d ratios [%.12f, %.12f] bound<=%.7f (%.0fs)\n", z, it, rmin, rmax,
                    2.0 * log(2.0 * z) / log(rmax), now() - t0);
            if (rmax / rmin <= 1 + 2e-7 || since_best >= 30 || it >= 3000) break;
            if (0) {   /* checkpointing is done by the best-vector save above */
                snprintf(p, sizeof p, "%s/eigen_h%d.bin", dir, H);
                FILE *f = fopen(p, "wb");
                int32_t hh = H; uint64_t n = N;
                fwrite(&hh, 4, 1, f); fwrite(&n, 8, 1, f); fwrite(&z, 8, 1, f);
                for (uint64_t s = 1; s <= N; s++) { double dv = V[s]; fwrite(&dv, 8, 1, f); }
                fclose(f);
            }
        }
        printf("z=%.6f best max ratio=%.15f bound=%.9f iterations=%d (vector saved at the best)\n", z, best, 2.0 * log(2.0 * z) / log(best), it);
        fflush(stdout);
    }
    pthread_mutex_lock(&R.mu); R.stop = 1; pthread_cond_broadcast(&R.cv); pthread_mutex_unlock(&R.mu);
    return 0;
}
