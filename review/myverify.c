/* Reviewer tool: independent recomputation of lambda_cert = max_s (A_z v)_s / v_s.
 * Single-threaded (no GCD), streaming the table with F_NOCACHE.
 *  (1) upward rounding set in THIS thread, same summation order as verify_cert -> lam_up
 *  (2) round-to-nearest, then a-priori inflation by (1-u)^-36, u = 2^-53 -> lam_apriori
 *      (each term passes through 1 product and <= 34 sums, plus 1 division)
 * Also reports: min v, max v, argmax state, count of v <= 0 / non-finite, any target > N,
 * and whether any row has a zero l=0 entry.  Weights read as hex floats.
 * usage: myverify <trans> <eigen> <weights> */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <fenv.h>
#include <fcntl.h>
#include <unistd.h>
#pragma STDC FENV_ACCESS ON

int main(int argc, char **argv) {
    int fd = open(argv[1], O_RDONLY); fcntl(fd, F_NOCACHE, 1);
    FILE *fe = fopen(argv[2], "rb"), *fw = fopen(argv[3], "r");
    int32_t hdr[2]; if (pread(fd, hdr, 8, 0) != 8) return 1;
    int H = hdr[0]; uint64_t N = (uint64_t)hdr[1]; int L = 2 * H - 1, NT = 2 * L;
    int32_t eh; uint64_t en; double z;
    if (fread(&eh, 4, 1, fe) != 1 || fread(&en, 8, 1, fe) != 1 || fread(&z, 8, 1, fe) != 1) return 1;
    if (eh != H || en != N) { printf("eigen/table mismatch\n"); return 1; }
    double *v = malloc((N + 1) * 8); v[0] = 0;
    if (fread(v + 1, 8, N, fe) != N) { printf("short eigen\n"); return 1; }
    { int c = fgetc(fe); printf("eigen trailing bytes: %s\n", c == EOF ? "none" : "PRESENT"); }
    double w[64]; char buf[256];
    for (int k = 0; k < L; k++) { if (!fgets(buf, sizeof buf, fw)) return 1; w[k] = strtod(buf, NULL); }
    { int extra = 0; while (fgets(buf, sizeof buf, fw)) extra++; printf("extra weight lines: %d\n", extra); }
    double vmin = INFINITY, vmax = 0; long bad = 0;
    for (uint64_t s = 1; s <= N; s++) { if (!(v[s] > 0) || !isfinite(v[s])) bad++; if (v[s] < vmin) vmin = v[s]; if (v[s] > vmax) vmax = v[s]; }
    printf("h=%d N=%llu z=%a (%.17g) v in [%.6e, %.6e] nonpositive/nonfinite=%ld\n", H, (unsigned long long)N, z, z, vmin, vmax, bad);
    const uint64_t CH = 1 << 20;
    uint32_t *tb = malloc(CH * NT * 4);
    double lam_up = 0, lam_rn = 0; uint64_t arg_up = 0, arg_rn = 0; long badidx = 0, zero_l0 = 0;
    for (uint64_t base = 1; base <= N; base += CH) {
        uint64_t cnt = (N - base + 1 < CH) ? N - base + 1 : CH, bytes = cnt * NT * 4, off = 8 + (base - 1) * NT * 4, d = 0;
        while (d < bytes) { ssize_t r = pread(fd, (char *)tb + d, bytes - d, off + d); if (r <= 0) return 1; d += r; }
        for (uint64_t k = 0; k < cnt * NT; k++) if (tb[k] > N) badidx++;
        for (uint64_t k = 0; k < cnt; k++) { const uint32_t *r = tb + k * NT; if (!r[H - 1] || !r[L + H - 1]) zero_l0++; }
        /* (1) upward, this thread */
        fesetround(FE_UPWARD);
        volatile double one = 1.0, tiny = 1e-300; if (!(one + tiny > one)) { printf("upward not in effect\n"); return 4; }
        for (uint64_t k = 0; k < cnt; k++) {
            const uint32_t *r = tb + k * NT; uint64_t s = base + k; double sum = 0;
            for (int q = 0; q < L; q++) sum = sum + w[q] * v[r[q]];
            for (int q = 0; q < L; q++) sum = sum + w[q] * v[r[L + q]];
            double qq = sum / v[s];
            if (qq > lam_up) { lam_up = qq; arg_up = s; }
        }
        fesetround(FE_TONEAREST);
        /* (2) nearest */
        for (uint64_t k = 0; k < cnt; k++) {
            const uint32_t *r = tb + k * NT; uint64_t s = base + k; double sum = 0;
            for (int q = 0; q < L; q++) sum = sum + w[q] * v[r[q]];
            for (int q = 0; q < L; q++) sum = sum + w[q] * v[r[L + q]];
            double qq = sum / v[s];
            if (qq > lam_rn) { lam_rn = qq; arg_rn = s; }
        }
    }
    printf("targets > N: %ld   rows with l=0 rejected: %ld\n", badidx, zero_l0);
    printf("lam_up      = %a (%.17g) at state %llu\n", lam_up, lam_up, (unsigned long long)arg_up);
    printf("lam_nearest = %a (%.17g) at state %llu\n", lam_rn, lam_rn, (unsigned long long)arg_rn);
    /* a-priori: lam_rn / (1-u)^36, computed upward */
    fesetround(FE_UPWARD);
    double f = 1.0, onemu = 1.0 - 0x1p-53;   /* exact */
    for (int i = 0; i < 36; i++) f = f / onemu;
    double lam_ap = lam_rn * f;
    fesetround(FE_TONEAREST);
    printf("lam_apriori = %a (%.17g)  [factor %a]\n", lam_ap, lam_ap, f);
    return 0;
}
