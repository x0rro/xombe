#ifndef XOMBE_FFT_Q15_H
#define XOMBE_FFT_Q15_H

/* block-floating-point complex FFT on packed Q15 pairs (q15ops.h: re low, im high), decimation in time in
   radix-2^2 passes: stages 2 and 4 with twiddles 1 and -i, then stage pairs (len, 2*len), then a radix-2 stage
   when log2 m is odd
   a bound b on the magnitudes decides halving: a stage halves (SHADD16/SHSUB16 and the exchanged forms) only
   when unhalved outputs (QADD16/QSUB16) could pass 32767. b doubles through an unhalved stage, holds through a
   halved one, and is tightened after each pass from the largest output component, read off the OR of
   v ^ (v << 1). twiddle products are rounded in Q30 by SMLSD/SMLADX */

#include "q15ops.h"

/* largest bound a stage passes unhalved, its outputs reach 2b + 2 with the twiddle's rounding */
#define XOMBE_CFFT_Q15_FREE 16382
/* largest input peak, sqrt 2 * peak <= XOMBE_CFFT_Q15_FREE so the first stage runs unhalved */
#define XOMBE_CFFT_Q15_PEAK 11583

/* ceil(sqrt 2 * v) for 0 <= v < 2^16 */
static inline int32_t xombe__q15_sqrt2_up(int32_t v)
{
    return (int32_t)(((uint32_t)v * 46341u + 32767u) >> 15);
}

#if defined(__GNUC__)
#define XOMBE_Q15_INLINE static inline __attribute__((always_inline))
#else
#define XOMBE_Q15_INLINE static inline
#endif

/* b*w, rounded */
XOMBE_Q15_INLINE uint32_t xombe__q15_cmul(uint32_t b, uint32_t w)
{
    return q15_pack_lsl1((uint32_t)(q15_smlsd(b, w, 0x4000) >> 15), (uint32_t)q15_smladx(b, w, 0x4000));
}

/* a + b, a - b, a + (-i)*b, and a - (-i)*b, halved when h */
XOMBE_Q15_INLINE uint32_t xombe__q15_add(uint32_t a, uint32_t b, int h)
{
    return h ? q15_shadd16(a, b) : q15_qadd16(a, b);
}

XOMBE_Q15_INLINE uint32_t xombe__q15_sub(uint32_t a, uint32_t b, int h)
{
    return h ? q15_shsub16(a, b) : q15_qsub16(a, b);
}

XOMBE_Q15_INLINE uint32_t xombe__q15_addni(uint32_t a, uint32_t b, int h)
{
    return h ? q15_shsax(a, b) : q15_qsax(a, b);
}

XOMBE_Q15_INLINE uint32_t xombe__q15_subni(uint32_t a, uint32_t b, int h)
{
    return h ? q15_shasx(a, b) : q15_qasx(a, b);
}

/* detector bits of v: bit q of each half set where bits q and q - 1 differ */
XOMBE_Q15_INLINE uint32_t xombe__q15_det(uint32_t v)
{
    return v ^ (v << 1);
}

/* stages 2 and 4 (twiddles 1 and -i) on each group of four points, halved per h1 and h2 */
XOMBE_Q15_INLINE uint32_t xombe__cfft_q15_first(uint32_t *x, int m, int h1, int h2)
{
    uint32_t det = 0;
    int i;
    for (i = 0; i < m; i += 4) {
        uint32_t x0 = x[i], x1 = x[i + 1], x2 = x[i + 2], x3 = x[i + 3], y0, y1, y2, y3, z0, z1, z2, z3;
        y0 = xombe__q15_add(x0, x1, h1);
        y1 = xombe__q15_sub(x0, x1, h1);
        y2 = xombe__q15_add(x2, x3, h1);
        y3 = xombe__q15_sub(x2, x3, h1);
        x[i] = z0 = xombe__q15_add(y0, y2, h2);
        x[i + 1] = z1 = xombe__q15_addni(y1, y3, h2);
        x[i + 2] = z2 = xombe__q15_sub(y0, y2, h2);
        x[i + 3] = z3 = xombe__q15_subni(y1, y3, h2);
        det |= xombe__q15_det(z0) | xombe__q15_det(z1) | xombe__q15_det(z2) | xombe__q15_det(z3);
    }
    return det;
}

/* stages len and 2*len as one radix-2^2 pass, halved per h1 and h2: in each block of 2*len, points j, j + q,
   j + 2q, and j + 3q for j < q = len/2, the second stage's odd twiddle the first's times -i */
XOMBE_Q15_INLINE uint32_t xombe__cfft_q15_r22(uint32_t *x, int m, int len, const uint32_t *tw, int h1, int h2)
{
    uint32_t det = 0;
    int q = len >> 1, sa = m / len, sb = sa >> 1, i, j;
    for (i = 0; i < m; i += 2 * len) {
        uint32_t *p = x + i;
        for (j = 0; j < q; j++) {
            uint32_t wa = tw[j * sa], wb = tw[j * sb], y0, y1, y2, y3, z0, z1, z2, z3, t;
            t = xombe__q15_cmul(p[j + q], wa);
            y0 = xombe__q15_add(p[j], t, h1);
            y1 = xombe__q15_sub(p[j], t, h1);
            t = xombe__q15_cmul(p[j + 3 * q], wa);
            y2 = xombe__q15_add(p[j + 2 * q], t, h1);
            y3 = xombe__q15_sub(p[j + 2 * q], t, h1);
            t = xombe__q15_cmul(y2, wb);
            p[j] = z0 = xombe__q15_add(y0, t, h2);
            p[j + 2 * q] = z2 = xombe__q15_sub(y0, t, h2);
            t = xombe__q15_cmul(y3, wb);
            p[j + q] = z1 = xombe__q15_addni(y1, t, h2);
            p[j + 3 * q] = z3 = xombe__q15_subni(y1, t, h2);
            det |= xombe__q15_det(z0) | xombe__q15_det(z1) | xombe__q15_det(z2) | xombe__q15_det(z3);
        }
    }
    return det;
}

/* last stage (len = m) alone */
XOMBE_Q15_INLINE uint32_t xombe__cfft_q15_r2(uint32_t *x, int m, const uint32_t *tw, int h)
{
    uint32_t det = 0;
    int q = m >> 1, j;
    for (j = 0; j < q; j++) {
        uint32_t t = xombe__q15_cmul(x[j + q], tw[j]), a = x[j], u, v;
        x[j] = u = xombe__q15_add(a, t, h);
        x[j + q] = v = xombe__q15_sub(a, t, h);
        det |= xombe__q15_det(u) | xombe__q15_det(v);
    }
    return det;
}

/* whether the next stage halves, updating bound b and halving count s */
static inline int xombe__cfft_q15_halve(int32_t *b, int *s)
{
    if (*b > XOMBE_CFFT_Q15_FREE) {
        *b += 1;
        ++*s;
        return 1;
    }
    *b = 2 * *b + 2;
    return 0;
}

/* tightens b from a pass's detector bits: every component has |v| <= 2^q for the top bit q (bit 0 of each half
   may hold the other half's sign) */
static inline void xombe__cfft_q15_tighten(int32_t *b, uint32_t det)
{
    uint32_t d = (det | det >> 16) & 0xfffeu;
    int32_t bq = xombe__q15_sqrt2_up((int32_t)1 << (d ? 31 - q15_clz32(d) : 0));
    if (*b > bq)
        *b = bq;
}

/* transforms m points x (m a power of two >= 4, bit-reversed order) in place into natural order, twiddles
   tw[j] = exp(-2 pi i j/m) for j < m/2 in Q15. no input component may exceed peak, sqrt 2 * peak <=
   XOMBE_CFFT_Q15_FREE. returns the halving count s, x holds the transform times 2^-s
   each pass decides its stages' halvings from the bound alone and is checked once at its end */
static int xombe__cfft_q15(uint32_t *x, int m, const uint32_t *tw, int32_t peak)
{
    int32_t b = xombe__q15_sqrt2_up(peak);
    int len, s = 0, h1, h2;
    uint32_t det;
    h1 = xombe__cfft_q15_halve(&b, &s);
    h2 = xombe__cfft_q15_halve(&b, &s);
    switch (h1 << 1 | h2) {
    case 0: det = xombe__cfft_q15_first(x, m, 0, 0); break;
    case 1: det = xombe__cfft_q15_first(x, m, 0, 1); break;
    case 2: det = xombe__cfft_q15_first(x, m, 1, 0); break;
    default: det = xombe__cfft_q15_first(x, m, 1, 1); break;
    }
    xombe__cfft_q15_tighten(&b, det);
    for (len = 8; 2 * len <= m; len <<= 2) {
        h1 = xombe__cfft_q15_halve(&b, &s);
        h2 = xombe__cfft_q15_halve(&b, &s);
        switch (h1 << 1 | h2) {
        case 0: det = xombe__cfft_q15_r22(x, m, len, tw, 0, 0); break;
        case 1: det = xombe__cfft_q15_r22(x, m, len, tw, 0, 1); break;
        case 2: det = xombe__cfft_q15_r22(x, m, len, tw, 1, 0); break;
        default: det = xombe__cfft_q15_r22(x, m, len, tw, 1, 1); break;
        }
        xombe__cfft_q15_tighten(&b, det);
    }
    if (len <= m) {
        if (xombe__cfft_q15_halve(&b, &s))
            xombe__cfft_q15_r2(x, m, tw, 1);
        else
            xombe__cfft_q15_r2(x, m, tw, 0);
    }
    return s;
}

#endif
