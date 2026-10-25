/*
 * kbverify.c: kverify.c for kbig's table files, which do not fit in memory.
 *
 * Reads kbig_k<k>_h<h>.bin (records streamed from disk in chunks),
 * kvec_k<k>_h<h>.bin (k, h, N, z, v[0..N-1]), weights_k<k>_h<h>.txt (per weight
 * type an upper and a lower bound on the exact weight at that z, hex doubles)
 * and wtypes_k<k>_h<h>.txt (the weight types the weights were computed for,
 * which must equal the ones stored in the table header).  The stored doubles v
 * are taken as the exact vector.  All terms are nonnegative, so:
 *
 *   with FE_UPWARD in every worker block and the upper weights, the computed
 *   max_s (A_z v)_s / v_s is an upper bound lambda_cert with A_z v <= lambda_cert v;
 *
 *   with FE_DOWNWARD and the lower weights, the computed min_s (A_z v)_s / v_s
 *   is a lower bound mu with A_z v >= mu v, hence mu <= spectral radius of A_z.
 *
 * Each state's sum is formed in the same order and with the same operations as
 * in kverify.c, so on any table both programs print the same numbers.
 *
 * Usage: kbverify <k> <h> <dir>
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <fenv.h>
#include <sys/stat.h>
#include <dispatch/dispatch.h>

#pragma STDC FENV_ACCESS ON

typedef struct { int32_t size, ta, tb, na, nb; } WType;
typedef struct { uint32_t tgt; uint16_t wt; uint16_t cnt; } Rec;

int main(int argc, char **argv) {
    if (argc != 4) { fprintf(stderr, "usage: kbverify k h dir\n"); return 1; }
    int K = atoi(argv[1]), H = atoi(argv[2]);
    const char *dir = argv[3];
    char p1[1024], p2[1024], p3[1024], p4[1024];
    snprintf(p1, sizeof p1, "%s/kbig_k%d_h%d.bin", dir, K, H);
    snprintf(p2, sizeof p2, "%s/kvec_k%d_h%d.bin", dir, K, H);
    snprintf(p3, sizeof p3, "%s/weights_k%d_h%d.txt", dir, K, H);
    snprintf(p4, sizeof p4, "%s/wtypes_k%d_h%d.txt", dir, K, H);
    FILE *ft = fopen(p1, "rb"), *fe = fopen(p2, "rb"), *fw = fopen(p3, "r"), *fy = fopen(p4, "r");
    if (!ft || !fe || !fw || !fy) { perror("open"); return 1; }
    char magic[4];
    int32_t hk[2], hi2[2];
    uint64_t hn[2];
    if (fread(magic, 4, 1, ft) != 1 || memcmp(magic, "KBIG", 4) || fread(hk, sizeof hk, 1, ft) != 1 ||
        fread(hn, sizeof hn, 1, ft) != 1 || fread(hi2, sizeof hi2, 1, ft) != 1 || hk[0] != K || hk[1] != H) {
        fprintf(stderr, "bad table header\n"); return 1;
    }
    uint64_t N = hn[0], R = hn[1];
    int NW = hi2[1];
    WType wt[256];
    if (NW < 1 || NW > 256 || fread(wt, sizeof(WType), NW, ft) != (size_t)NW) { fprintf(stderr, "bad weight types\n"); return 1; }
    uint64_t recbase = 4 + 8 + 16 + 8 + (uint64_t)NW * sizeof(WType);
    struct stat sb;
    if (fstat(fileno(ft), &sb) || (uint64_t)sb.st_size != recbase + 8 * R + 8 * (N + 1) + 8 * N) {
        fprintf(stderr, "table file size does not match its header\n"); return 2;
    }
    /* the weight types in the table header must be the ones the weights file was made for */
    for (int i = 0; i < NW; i++) {
        int idx, s, ta, tb, na, nb;
        if (fscanf(fy, "%d %d %d %d %d %d", &idx, &s, &ta, &tb, &na, &nb) != 6 || idx != i || s != wt[i].size ||
            ta != wt[i].ta || tb != wt[i].tb || na != wt[i].na || nb != wt[i].nb) {
            fprintf(stderr, "weight type %d in the table differs from %s\n", i, p4); return 2;
        }
    }
    { int extra; if (fscanf(fy, "%d", &extra) == 1) { fprintf(stderr, "%s has more weight types than the table\n", p4); return 2; } }
    fclose(fy);
    uint64_t *off = malloc((N + 1) * sizeof(uint64_t));
    if (!off) { fprintf(stderr, "out of memory\n"); return 1; }
    if (fseeko(ft, (off_t)(recbase + 8 * R), SEEK_SET) || fread(off, sizeof(uint64_t), N + 1, ft) != N + 1) {
        fprintf(stderr, "short offsets\n"); return 1;
    }
    int32_t ek[2]; uint64_t en; double z;
    if (fread(ek, sizeof ek, 1, fe) != 1 || fread(&en, sizeof en, 1, fe) != 1 || fread(&z, sizeof z, 1, fe) != 1) return 1;
    if (ek[0] != K || ek[1] != H || en != N) { fprintf(stderr, "vector does not match table\n"); return 1; }
    double *v = malloc(N * sizeof(double));
    if (!v || fread(v, sizeof(double), N, fe) != N) { fprintf(stderr, "short vector\n"); return 1; }
    { char extra; if (fread(&extra, 1, 1, fe) == 1) { fprintf(stderr, "vector file too long\n"); return 1; } }
    fclose(fe);
    double wup[256], wdn[256];
    for (int i = 0; i < NW; i++) {
        char buf[256];
        if (!fgets(buf, sizeof buf, fw)) { fprintf(stderr, "short weights\n"); return 1; }
        char *e;
        wup[i] = strtod(buf, &e);
        wdn[i] = strtod(e, NULL);
        if (!(wup[i] > 0) || !(wdn[i] > 0) || wdn[i] > wup[i]) { fprintf(stderr, "bad weight line %d\n", i); return 1; }
    }
    fclose(fw);
    if (off[0] != 0 || off[N] != R) { fprintf(stderr, "bad offsets\n"); return 2; }
    for (uint64_t s = 0; s < N; s++) {
        if (!(v[s] > 0) || !isfinite(v[s])) { fprintf(stderr, "v[%llu] not positive\n", (unsigned long long)s); return 2; }
        if (off[s + 1] < off[s]) { fprintf(stderr, "offsets not monotone\n"); return 2; }
    }
    /* stream the records: chunks of whole states, at most CAP records each */
    const uint64_t CAP = KBV_CAP;
    Rec *rec = malloc(CAP * sizeof(Rec));
    if (!rec) { fprintf(stderr, "out of memory\n"); return 1; }
    if (fseeko(ft, (off_t)recbase, SEEK_SET)) { perror("fseeko"); return 1; }
    const double *pu = wup, *pd = wdn;
    const int chunks = 512;
    double *maxs = calloc(chunks, sizeof(double)), *mins = calloc(chunks, sizeof(double));
    double lam = 0, mu = INFINITY;
    uint64_t nchunk = 0;
    for (uint64_t s0 = 0; s0 < N;) {
        uint64_t s1 = s0 + 1;
        while (s1 < N && off[s1 + 1] - off[s0] <= CAP) s1++;
        uint64_t r0 = off[s0], nr = off[s1] - r0;
        if (nr > CAP) { fprintf(stderr, "state %llu has too many records\n", (unsigned long long)s0); return 2; }
        if (fread(rec, sizeof(Rec), nr, ft) != nr) { fprintf(stderr, "short records\n"); return 1; }
        for (uint64_t r = 0; r < nr; r++)
            if (rec[r].tgt >= N || rec[r].wt >= NW || rec[r].cnt == 0) { fprintf(stderr, "bad record %llu\n", (unsigned long long)(r0 + r)); return 2; }
        const Rec *rc = rec;
        uint64_t per = (s1 - s0 + chunks - 1) / chunks;
        dispatch_apply(chunks, dispatch_get_global_queue(QOS_CLASS_UTILITY, 0), ^(size_t ci) {
            uint64_t lo = s0 + ci * per, hi = lo + per;
            if (hi > s1) hi = s1;
            fesetround(FE_UPWARD);
            volatile double one = 1.0, tiny = 1e-30;
            if (!(one + tiny > one)) { fprintf(stderr, "upward rounding not in effect\n"); exit(4); }
            double mx = 0;
            for (uint64_t s = lo; s < hi; s++) {
                double sum = 0;
                for (uint64_t r = off[s]; r < off[s + 1]; r++) sum = sum + ((double)rc[r - r0].cnt * pu[rc[r - r0].wt]) * v[rc[r - r0].tgt];
                double q = sum / v[s];
                if (q > mx) mx = q;
            }
            fesetround(FE_DOWNWARD);
            if (!(one - tiny < one) || (one + tiny > one)) { fprintf(stderr, "downward rounding not in effect\n"); exit(4); }
            double mn = INFINITY;
            for (uint64_t s = lo; s < hi; s++) {
                double sum = 0;
                for (uint64_t r = off[s]; r < off[s + 1]; r++) sum = sum + ((double)rc[r - r0].cnt * pd[rc[r - r0].wt]) * v[rc[r - r0].tgt];
                double q = sum / v[s];
                if (q < mn) mn = q;
            }
            fesetround(FE_TONEAREST);
            maxs[ci] = mx; mins[ci] = mn;
        });
        for (int i = 0; i < chunks; i++) { if (maxs[i] > lam) lam = maxs[i]; if (mins[i] < mu) mu = mins[i]; }
        nchunk++;
        s0 = s1;
    }
    fclose(ft);
    printf("k=%d h=%d N=%llu records=%llu z=%a (%.17g)\n", K, H, (unsigned long long)N, (unsigned long long)R, z, z);
    printf("lambda_cert=%a (%.17g)\n", lam, lam);
    printf("rho_lower=%a (%.17g)\n", mu, mu);
    fprintf(stderr, "streamed %llu chunks\n", (unsigned long long)nchunk);
    return 0;
}
