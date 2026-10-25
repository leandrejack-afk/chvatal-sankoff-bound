/* rvb_emu.h: directed rounding emulated in round-to-nearest (reviewer, 2026-10-25).
 *
 * The FPU stays in round-to-nearest. Each operation is computed in round-to-nearest,
 * its exact error is recovered with an error-free transformation (fma for products and
 * for the division remainder, Knuth's TwoSum for sums), and the result is moved one ulp
 * up (or down) when the exact value lies above (or below) the rounded one. When the
 * round-to-nearest result r of an exact value x satisfies x <= r, r is already the
 * smallest double >= x (otherwise the double below r would be nearer to x), so this gives
 * exactly IEEE round-toward-+inf (or -inf).
 *
 * Exactness conditions: fma(a, b, -a*b) and fma(-q, b, a) are exact unless the product
 * underflows (|a*b| below about 2^-970); TwoSum is exact barring overflow. Callers count
 * operands below RVB_TINY and refuse to certify if any occur.
 *
 * Needs -ffp-contract=off (no fusing of a*b - p or s - a etc.).
 */
#include <math.h>

#define RVB_TINY 0x1p-900

static inline double up_mul(double a, double b) {
    double p = a * b;
    double e = fma(a, b, -p);
    return e > 0 ? nextafter(p, INFINITY) : p;
}
static inline double dn_mul(double a, double b) {
    double p = a * b;
    double e = fma(a, b, -p);
    return e < 0 ? nextafter(p, -INFINITY) : p;
}
static inline double two_sum_err(double a, double b, double s) {
    double bb = s - a;
    return (a - (s - bb)) + (b - bb);
}
static inline double up_add(double a, double b) {
    double s = a + b;
    return two_sum_err(a, b, s) > 0 ? nextafter(s, INFINITY) : s;
}
static inline double dn_add(double a, double b) {
    double s = a + b;
    return two_sum_err(a, b, s) < 0 ? nextafter(s, -INFINITY) : s;
}
/* b > 0 */
static inline double up_div(double a, double b) {
    double q = a / b;
    double r = fma(-q, b, a);
    return r > 0 ? nextafter(q, INFINITY) : q;
}
static inline double dn_div(double a, double b) {
    double q = a / b;
    double r = fma(-q, b, a);
    return r < 0 ? nextafter(q, -INFINITY) : q;
}
