/*
 * kverify.c: rigorous certificate check for the k-letter machine.
 *
 * Reads ktab_k<k>_h<h>.bin, kvec_k<k>_h<h>.bin (k, h, N, z, v[0..N-1]) and
 * weights_k<k>_h<h>.txt (per weight type: an upper and a lower bound on the
 * exact weight at that z, as hex doubles).  The stored doubles v are taken as
 * the exact vector.  All terms are nonnegative, so:
 *
 *   with FE_UPWARD in every worker thread and the upper weights, the computed
 *   max_s (A_z v)_s / v_s is an upper bound lambda_cert with A_z v <= lambda_cert v;
 *
 *   with FE_DOWNWARD and the lower weights, the computed min_s (A_z v)_s / v_s
 *   is a lower bound mu with A_z v >= mu v, hence mu <= spectral radius of A_z.
 *
 * Usage: kverify <k> <h> <dir>
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <fenv.h>
#include <dispatch/dispatch.h>

#pragma STDC FENV_ACCESS ON

typedef struct { int32_t size, ta, tb, na, nb; } WType;
typedef struct { uint32_t tgt; uint16_t wt; uint16_t cnt; } Rec;

int main(int argc, char **argv) {
    if (argc != 4) { fprintf(stderr, "usage: kverify k h dir\n"); return 1; }
    int K = atoi(argv[1]), H = atoi(argv[2]);
    const char *dir = argv[3];
    char p1[1024], p2[1024], p3[1024];
    snprintf(p1, sizeof p1, "%s/ktab_k%d_h%d.bin", dir, K, H);
    snprintf(p2, sizeof p2, "%s/kvec_k%d_h%d.bin", dir, K, H);
    snprintf(p3, sizeof p3, "%s/weights_k%d_h%d.txt", dir, K, H);
    FILE *ft = fopen(p1, "rb"), *fe = fopen(p2, "rb"), *fw = fopen(p3, "r");
    if (!ft || !fe || !fw) { perror("open"); return 1; }
    char magic[4];
    int32_t hk[2], hi2[2];
    uint64_t hn[2];
    if (fread(magic, 4, 1, ft) != 1 || memcmp(magic, "KLUE", 4) || fread(hk, sizeof hk, 1, ft) != 1 ||
        fread(hn, sizeof hn, 1, ft) != 1 || fread(hi2, sizeof hi2, 1, ft) != 1 || hk[0] != K || hk[1] != H) {
        fprintf(stderr, "bad table header\n"); return 1;
    }
    uint64_t N = hn[0], R = hn[1];
    int NW = hi2[1];
    WType wt[256];
    if (NW > 256 || fread(wt, sizeof(WType), NW, ft) != (size_t)NW) { fprintf(stderr, "bad weight types\n"); return 1; }
    uint64_t *off = malloc((N + 1) * sizeof(uint64_t));
    Rec *rec = malloc(R * sizeof(Rec));
    if (!off || !rec) { fprintf(stderr, "out of memory\n"); return 1; }
    if (fread(off, sizeof(uint64_t), N + 1, ft) != N + 1 || fread(rec, sizeof(Rec), R, ft) != R) {
        fprintf(stderr, "short table\n"); return 1;
    }
    fclose(ft);
    int32_t ek[2]; uint64_t en; double z;
    if (fread(ek, sizeof ek, 1, fe) != 1 || fread(&en, sizeof en, 1, fe) != 1 || fread(&z, sizeof z, 1, fe) != 1) return 1;
    if (ek[0] != K || ek[1] != H || en != N) { fprintf(stderr, "vector does not match table\n"); return 1; }
    double *v = malloc(N * sizeof(double));
    if (fread(v, sizeof(double), N, fe) != N) { fprintf(stderr, "short vector\n"); return 1; }
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
    for (uint64_t r = 0; r < R; r++)
        if (rec[r].tgt >= N || rec[r].wt >= NW || rec[r].cnt == 0) { fprintf(stderr, "bad record %llu\n", (unsigned long long)r); return 2; }
    const int chunks = 512;
    double *maxs = calloc(chunks, sizeof(double)), *mins = calloc(chunks, sizeof(double));
    uint64_t per = (N + chunks - 1) / chunks;
    const double *pu = wup, *pd = wdn;
    dispatch_apply(chunks, dispatch_get_global_queue(QOS_CLASS_UTILITY, 0), ^(size_t ci) {
        uint64_t lo = ci * per, hi = lo + per;
        if (hi > N) hi = N;
        fesetround(FE_UPWARD);
        volatile double one = 1.0, tiny = 1e-30;
        if (!(one + tiny > one)) { fprintf(stderr, "upward rounding not in effect\n"); exit(4); }
        double mx = 0;
        for (uint64_t s = lo; s < hi; s++) {
            double sum = 0;
            for (uint64_t r = off[s]; r < off[s + 1]; r++) sum = sum + ((double)rec[r].cnt * pu[rec[r].wt]) * v[rec[r].tgt];
            double q = sum / v[s];
            if (q > mx) mx = q;
        }
        fesetround(FE_DOWNWARD);
        if (!(one - tiny < one) || (one + tiny > one)) { fprintf(stderr, "downward rounding not in effect\n"); exit(4); }
        double mn = INFINITY;
        for (uint64_t s = lo; s < hi; s++) {
            double sum = 0;
            for (uint64_t r = off[s]; r < off[s + 1]; r++) sum = sum + ((double)rec[r].cnt * pd[rec[r].wt]) * v[rec[r].tgt];
            double q = sum / v[s];
            if (q < mn) mn = q;
        }
        fesetround(FE_TONEAREST);
        maxs[ci] = mx; mins[ci] = mn;
    });
    double lam = 0, mu = INFINITY;
    for (int i = 0; i < chunks; i++) { if (maxs[i] > lam) lam = maxs[i]; if (mins[i] < mu) mu = mins[i]; }
    printf("k=%d h=%d N=%llu records=%llu z=%a (%.17g)\n", K, H, (unsigned long long)N, (unsigned long long)R, z, z);
    printf("lambda_cert=%a (%.17g)\n", lam, lam);
    printf("rho_lower=%a (%.17g)\n", mu, mu);
    return 0;
}
