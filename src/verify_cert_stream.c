/*
 * Rigorous check of a certificate for gamma_2 <= 2 log(2z) / log(lambda).
 *
 * Input: the transition table trans_h<h>.bin, a vector file eigen_h<h>.bin
 * (h, N, z, v[1..N]) and a file of weights w_up[k] (k = l + h - 1), each an
 * upper bound on w_{mu,z}(T_{c,l}) for that exact double z.
 *
 * With the rounding mode set upward in every worker thread, the program
 * computes, for every accepting state s, an upper bound on (A_z v)_s / v_s,
 * where v is taken as the exact vector of stored doubles.  All terms are
 * nonnegative, so upward rounding of every product, sum and quotient gives an
 * upper bound.  It prints lambda_cert = max over s, as a decimal and as hex,
 * and checks v_s > 0 for every s.  Then A_z v <= lambda_cert v componentwise,
 * which is the witness Lueker's Theorem 3.3 needs (see (41) in the paper).
 *
 * Usage: verify_cert <trans file> <eigen file> <weights file>
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <fenv.h>
#include <dispatch/dispatch.h>

#pragma STDC FENV_ACCESS ON

int main(int argc, char **argv) {
    if (argc != 4) { fprintf(stderr, "usage: verify_cert trans eigen weights\n"); return 1; }
    FILE *ft = fopen(argv[1], "rb"), *fe = fopen(argv[2], "rb"), *fw = fopen(argv[3], "r");
    if (!ft || !fe || !fw) { perror("open"); return 1; }
    int32_t hdr[2];
    if (fread(hdr, sizeof hdr, 1, ft) != 1) return 1;
    int H = hdr[0];
    uint64_t N = (uint64_t)hdr[1];
    int NT = 2 * (2 * H - 1), L = 2 * H - 1;
    int32_t eh; uint64_t en; double z;
    if (fread(&eh, sizeof eh, 1, fe) != 1 || fread(&en, sizeof en, 1, fe) != 1 ||
        fread(&z, sizeof z, 1, fe) != 1) return 1;
    if (eh != H || en != N) { fprintf(stderr, "eigen file does not match table\n"); return 1; }
    double *v = malloc((N + 1) * sizeof(double));
    v[0] = 0.0;
    if (fread(v + 1, sizeof(double), N, fe) != N) { fprintf(stderr, "short eigen\n"); return 1; }
    double w[64];
    for (int k = 0; k < L; k++) {
        char buf[128];
        if (!fgets(buf, sizeof buf, fw)) { fprintf(stderr, "short weights\n"); return 1; }
        w[k] = strtod(buf, NULL);   /* hex float, exact */
        if (!(w[k] > 0)) { fprintf(stderr, "bad weight\n"); return 1; }
    }
    for (uint64_t s = 1; s <= N; s++) {
        if (!(v[s] > 0) || !isfinite(v[s])) { fprintf(stderr, "v[%llu] not positive\n", (unsigned long long)s); return 2; }
    }
    const uint64_t CH = 1 << 24;
    uint32_t *trbuf = malloc(CH * NT * sizeof(uint32_t));
    const double *wp = w;
    double lam = 0;
    for (uint64_t base = 1; base <= N; base += CH) {
    uint64_t cnt = (N - base + 1 < CH) ? N - base + 1 : CH;
    if (fread(trbuf, sizeof(uint32_t), cnt * NT, ft) != cnt * NT) { fprintf(stderr, "short table\n"); return 1; }
    for (uint64_t i = 0; i < cnt * NT; i++) if (trbuf[i] > N) { fprintf(stderr, "bad state index\n"); return 2; }
    const uint32_t *tr = trbuf;                       /* row of state s is tr + (s - base) * NT */
    const uint64_t N0 = base - 1 + cnt;
    const int chunks = 512;
    double *maxs = calloc(chunks, sizeof(double));
    uint64_t per = (cnt + chunks - 1) / chunks;
    dispatch_apply(chunks, dispatch_get_global_queue(QOS_CLASS_USER_INITIATED, 0), ^(size_t ci) {
        fesetround(FE_UPWARD);
        volatile double one = 1.0, tiny = 1e-30;
        if (!(one + tiny > one)) { fprintf(stderr, "upward rounding not in effect\n"); exit(4); }
        double mx = 0;
        uint64_t lo = base + ci * per, hi = lo + per;
        if (hi > N0 + 1) hi = N0 + 1;
        for (uint64_t s = lo; s < hi; s++) {
            const uint32_t *r = tr + (s - base) * NT;
            double sum = 0;
            for (int k = 0; k < L; k++) sum = sum + wp[k] * v[r[k]];
            for (int k = 0; k < L; k++) sum = sum + wp[k] * v[r[L + k]];
            double q = sum / v[s];
            if (q > mx) mx = q;
        }
        maxs[ci] = mx;
        fesetround(FE_TONEAREST);
    });
    for (int i = 0; i < chunks; i++) if (maxs[i] > lam) lam = maxs[i];
    free(maxs);
    }
    printf("h=%d N=%llu z=%a (%.17g)\n", H, (unsigned long long)N, z, z);
    printf("lambda_cert=%a (%.17g)\n", lam, lam);
    return 0;
}
