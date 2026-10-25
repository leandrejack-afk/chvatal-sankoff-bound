/* rvb_emu_test.c: the emulated directed operations of rvb_emu.h against the hardware's
 * FE_UPWARD / FE_DOWNWARD on random positive operands (reviewer, 2026-10-25).
 * clang -O2 -ffp-model=strict -march=native -o rvb_emu_test rvb_emu_test.c -lm
 */
#include <fenv.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rvb_emu.h"

#pragma STDC FENV_ACCESS ON

static uint64_t st = 0x9e3779b97f4a7c15ULL;
static uint64_t rnd(void) { st ^= st << 13; st ^= st >> 7; st ^= st << 17; return st; }
static double rdouble(void) {
    /* positive doubles with exponents in [-60, 4], random mantissas; every 4th one short */
    uint64_t m = rnd() & ((1ULL << 52) - 1);
    if ((rnd() & 3) == 0) m &= ~((1ULL << 40) - 1);
    int e = (int)(rnd() % 65) - 60;
    uint64_t bits = ((uint64_t)(e + 1023) << 52) | m;
    double d;
    memcpy(&d, &bits, 8);
    return d;
}

int main(void) {
    long bad = 0, n = 0, moved = 0;
    for (long i = 0; i < 20000000; i++) {
        double a = rdouble(), b = rdouble();
        if (i % 5 == 0) a = (double)(1 + rnd() % 9748);       /* a record count */
        double hu[3], hd[3], eu[3], ed[3];
        fesetround(FE_UPWARD);
        volatile double va = a, vb = b;
        hu[0] = va * vb; hu[1] = va + vb; hu[2] = va / vb;
        fesetround(FE_DOWNWARD);
        hd[0] = va * vb; hd[1] = va + vb; hd[2] = va / vb;
        fesetround(FE_TONEAREST);
        eu[0] = up_mul(a, b); eu[1] = up_add(a, b); eu[2] = up_div(a, b);
        ed[0] = dn_mul(a, b); ed[1] = dn_add(a, b); ed[2] = dn_div(a, b);
        for (int j = 0; j < 3; j++) {
            n += 2;
            if (memcmp(&hu[j], &eu[j], 8) || memcmp(&hd[j], &ed[j], 8)) {
                if (++bad <= 5) printf("MISMATCH op %d a=%a b=%a hw up %a emu %a hw dn %a emu %a\n", j, a, b, hu[j], eu[j], hd[j], ed[j]);
            }
            if (hu[j] != hd[j]) moved++;
        }
    }
    printf("emulated directed rounding vs hardware: %ld operations, %ld mismatches (%ld inexact)\n", n, bad, moved);
    return bad != 0;
}
