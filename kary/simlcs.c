/* Monte Carlo sanity check: LCS of two random strings of length n over k letters.
 * Since E[L_n]/n increases to gamma_k (superadditivity), these are estimates from below. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
static uint64_t s = 88172645463325252ULL;
static inline uint64_t rnd(void) { s ^= s << 13; s ^= s >> 7; s ^= s << 17; return s; }
int main(int argc, char **argv) {
    int k = atoi(argv[1]), n = atoi(argv[2]), reps = atoi(argv[3]);
    s ^= (uint64_t)k * 0x9E3779B97F4A7C15ULL;
    unsigned char *a = malloc(n), *b = malloc(n);
    int *prev = calloc(n + 1, sizeof(int)), *cur = calloc(n + 1, sizeof(int));
    double sum = 0;
    for (int r = 0; r < reps; r++) {
        for (int i = 0; i < n; i++) { a[i] = rnd() % k; b[i] = rnd() % k; }
        for (int j = 0; j <= n; j++) prev[j] = 0;
        for (int i = 1; i <= n; i++) {
            cur[0] = 0;
            for (int j = 1; j <= n; j++) {
                int v = a[i - 1] == b[j - 1] ? prev[j - 1] + 1 : (prev[j] > cur[j - 1] ? prev[j] : cur[j - 1]);
                cur[j] = v;
            }
            int *t = prev; prev = cur; cur = t;
        }
        sum += (double)prev[n] / n;
        printf("k=%d n=%d rep %d: %.5f\n", k, n, r, (double)prev[n] / n);
    }
    printf("k=%d n=%d mean of %d: %.5f\n", k, n, reps, sum / reps);
    return 0;
}
