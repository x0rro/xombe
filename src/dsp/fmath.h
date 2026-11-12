#ifndef XOMBE_FMATH_H
#define XOMBE_FMATH_H

/* exp2, log2, and sincos for the inner loops: libm in float builds, which keep it as the reference, polynomials
   within a few ulp of libm in q15 builds (XOMBE_KERNELS=q15) */

#include <math.h>
#include <stdint.h>
#include <string.h>
#include "xombe_config.h"

#if XOMBE_KERNELS_Q15

/* round to nearest for |x| < 2^22: x + 1.5*2^23 holds the integer in its fraction field, offset by 2^22 */
#define XOMBE_FMATH_RND 12582912.0f

static inline float xombe__fmath_from_bits(uint32_t u)
{
    float f;
    memcpy(&f, &u, sizeof(f));
    return f;
}

static inline int32_t xombe__fmath_rnd_int(float t)
{
    uint32_t u;
    memcpy(&u, &t, sizeof(u));
    return (int32_t)(u & 0x7fffffu) - 0x400000;
}

/* x = k + f, |f| <= 1/2, 2^f by the Cephes exp2f polynomial, 2^k in the exponent field. libm outside
   (-126, 128), where the result may be subnormal, infinite, or NaN */
static inline float xombe__exp2f(float x)
{
    float t, f, p;
    int32_t k;
    if (!(x > -126.0f && x < 128.0f))
        return exp2f(x);
    t = x + XOMBE_FMATH_RND;
    k = xombe__fmath_rnd_int(t);
    f = x - (t - XOMBE_FMATH_RND);
    p = 1.535336188319500e-4f;
    p = p * f + 1.339887440266574e-3f;
    p = p * f + 9.618437357674640e-3f;
    p = p * f + 5.550332471162809e-2f;
    p = p * f + 2.402264791363012e-1f;
    p = p * f + 6.931472028550421e-1f;
    p = p * f + 1.0f;
    if (k > 127) /* x >= 127.5, 2^128 is not a float */
        return p * xombe__fmath_from_bits((uint32_t)(127 + 127) << 23) * 2.0f;
    return p * xombe__fmath_from_bits((uint32_t)(k + 127) << 23);
}

/* x = 2^e * (1 + m), 1 + m in [sqrt(1/2), sqrt 2), ln(1 + m) by the Cephes logf polynomial, times log2(e) in two
   parts (1 and log2(e) - 1) so m and e enter unrounded. libm for zero, subnormals, negatives, infinity, and NaN */
static inline float xombe__log2f(float x)
{
    uint32_t u;
    int32_t e;
    float m, z, y, p;
    if (!(x >= 1.17549435e-38f && x <= 3.40282347e38f))
        return log2f(x);
    memcpy(&u, &x, sizeof(u));
    e = (int32_t)(u >> 23) - 127;
    u = (u & 0x7fffffu) | 0x3f800000u;
    if (u > 0x3fb504f3u) { /* above sqrt 2 */
        u -= 0x00800000u;
        e++;
    }
    m = xombe__fmath_from_bits(u) - 1.0f;
    z = m * m;
    p = 7.0376836292e-2f;
    p = p * m - 1.1514610310e-1f;
    p = p * m + 1.1676998740e-1f;
    p = p * m - 1.2420140846e-1f;
    p = p * m + 1.4249322787e-1f;
    p = p * m - 1.6668057665e-1f;
    p = p * m + 2.0000714765e-1f;
    p = p * m - 2.4999993993e-1f;
    p = p * m + 3.3333331174e-1f;
    y = m * z * p - 0.5f * z;
    return y * 0.44269504088896341f + m * 0.44269504088896341f + y + m + (float)e;
}

/* x = q*pi/2 + r, |r| <= pi/4, with pi/2 in three parts (Cody-Waite), the first two short enough that q*part is
   exact for |q| < 2^11, then the Cephes sinf and cosf polynomials. libm for |x| >= 2^11 */
static inline void xombe__sincosf(float x, float *s, float *c)
{
    float t, r, z, sp, cp;
    int32_t q;
    if (!(fabsf(x) < 2048.0f)) {
        *s = sinf(x);
        *c = cosf(x);
        return;
    }
    t = x * 0.636619772f + XOMBE_FMATH_RND;
    q = xombe__fmath_rnd_int(t);
    t -= XOMBE_FMATH_RND;
    r = x - t * 1.5703125f;
    r -= t * 4.837512969970703125e-4f;
    r -= t * 7.54978995489188216e-8f;
    z = r * r;
    sp = -1.9515295891e-4f;
    sp = sp * z + 8.3321608736e-3f;
    sp = sp * z - 1.6666654611e-1f;
    sp = r + r * z * sp;
    cp = 2.443315711809948e-5f;
    cp = cp * z - 1.388731625493765e-3f;
    cp = cp * z + 4.166664568298827e-2f;
    cp = 1.0f - 0.5f * z + z * z * cp;
    switch (q & 3) {
    case 0: *s = sp; *c = cp; break;
    case 1: *s = cp; *c = -sp; break;
    case 2: *s = -sp; *c = -cp; break;
    default: *s = -cp; *c = sp; break;
    }
}

#else

static inline float xombe__exp2f(float x)
{
    return exp2f(x);
}

static inline float xombe__log2f(float x)
{
    return log2f(x);
}

static inline void xombe__sincosf(float x, float *s, float *c)
{
    *s = sinf(x);
    *c = cosf(x);
}

#endif

#endif
