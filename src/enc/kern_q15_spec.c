/* q15 analysis spectrum and its line fit preparation (kernels.h, XOMBE_KERNELS=q15) */
#include <string.h>
#include "kern_q15.h"
#include "../dsp/fft_q15.h"
#include "../tables/dsp_tables.h"

/* the XOMBE_FFT-point real transform runs as an SP_M-point complex one on sample pairs, then a real-input split */
#define SP_M (XOMBE_FFT / 2)
#define SP_LOG2 9
XOMBE_STATIC_ASSERT(1 << SP_LOG2 == SP_M, sp_log2);
XOMBE_STATIC_ASSERT(XOMBE_MAX_WINDOW < XOMBE_FFT, sp_window_fits);
XOMBE_STATIC_ASSERT(XOMBE_LF_BINS == XOMBE_SPEC_BINS, sp_bins);

/* keeps one copy of the bulky stages when the compiler clones kern_spectrum for constant arguments */
#if defined(__GNUC__)
#define SP_NOINLINE __attribute__((noinline))
#else
#define SP_NOINLINE
#endif
#if defined(__GNUC__) && !defined(__clang__)
#define SP_NOCLONE __attribute__((noclone))
#else
#define SP_NOCLONE
#endif

/* 2^e, -126 <= e <= 127 */
static float sp_pow2(int e)
{
    uint32_t u = (uint32_t)(e + 127) << 23;
    float f;
    memcpy(&f, &u, sizeof f);
    return f;
}

static inline uint32_t sp_abs(int32_t v)
{
    return (uint32_t)(v < 0 ? -v : v);
}

/* bins k and SP_M - k from transform pairs a = x[k], b = x[SP_M - k] and split twiddle w, as r[0..3] = re, im of
   bin k then of bin SP_M - k, in Q14 of the bins. no overflow: a dual product is at most |w| * |a| <= 2^15.5 * 2^15,
   a bin component at most (1 + sqrt 2) * 2^15 * 2^14 */
static inline void sp_bins(uint32_t a, uint32_t b, uint32_t w, int32_t r[4])
{
    uint32_t bx = b >> 16 | b << 16;
    int32_t er = (q15_lo(a) + q15_lo(b)) * 8192, ei = (q15_hi(a) - q15_hi(b)) * 8192;
    /* w * (a - conj b) / i: re = wr*ai + wi*ar + wr*bi - wi*br, im = wr*br + wi*bi - (wr*ar - wi*ai) */
    int32_t pr = (q15_smladx(w, a, 0) >> 2) + (q15_smlsd(w, bx, 0) >> 2);
    int32_t pi = (q15_smuad(w, b) >> 2) - (q15_smlsd(w, a, 0) >> 2);
    r[0] = er + pr;
    r[1] = ei + pi;
    r[2] = er - pr;
    r[3] = pi - ei;
}

/* top halves of |v0| and |v1| as a pair, each |v| taken as v ^ (v >> 31), which is |v| or |v| - 1 */
static inline uint32_t sp_top(int32_t v0, int32_t v1)
{
    return ((uint32_t)(v1 ^ (v1 >> 31)) & 0xffff0000u) | (uint32_t)(v0 ^ (v0 >> 31)) >> 16;
}

/* splits transform x of the pairs in place, returns h with every |component| <= (h + 1) * 2^16. layout after:
   bin 0 re at x[0] (im 0), bins 1..SP_M/2 - 1 re at x[k] and im at x[SP_M - k], x[SP_M/2] a copy of bin
   SP_M/2's im, bins SP_M/2..SP_M re then im at x[SP_M + 2*(k - SP_M/2)]. every write lands on a word already read.
   split twiddle exp(-2 pi i k/XOMBE_FFT) is xombe__fft_tw_q15[k/2] for even k, times SP_W1 for odd k */
#define SP_W1 0xff377fffu /* exp(-2 pi i/XOMBE_FFT) in Q15 */
static SP_NOINLINE uint32_t sp_split(uint32_t *x)
{
    const uint32_t *tw = xombe__fft_tw_q15;
    int32_t r[4];
    uint32_t m = 0;
    int k;
    /* k = 0 pairs x[0] with itself, its im (0) lands on x[SP_M] before k = SP_M/2 writes there */
    for (k = 0; k <= SP_M / 2; k++) {
        uint32_t *u = x + SP_M + 2 * (SP_M / 2 - k), w = tw[k >> 1];
        if (k & 1)
            w = xombe__q15_cmul(w, SP_W1);
        sp_bins(x[k], x[(SP_M - k) & (SP_M - 1)], w, r);
        x[k] = (uint32_t)r[0];
        x[SP_M - k] = (uint32_t)r[1];
        u[0] = (uint32_t)r[2];
        u[1] = (uint32_t)r[3];
        m = q15_umax16(m, q15_umax16(sp_top(r[0], r[1]), sp_top(r[2], r[3])));
    }
    return (m >> 16) > (m & 0xffffu) ? m >> 16 : m & 0xffffu;
}

/* (v0, v1) * 2^p * c / 2^32 rounded, as an int16 pair */
static inline uint32_t sp_q2(int32_t v0, int32_t v1, int p, int32_t c)
{
    return (uint16_t)q15_smmulr((int32_t)((uint32_t)v0 << p), c) |
           (uint32_t)q15_smmulr((int32_t)((uint32_t)v1 << p), c) << 16;
}

static inline void sp_st2(int16_t *d, uint32_t v)
{
    memcpy(d, &v, sizeof v);
}

/* packs centre[j] * w[half + j] * g rounded, |j| <= half, as sample j mod XOMBE_FFT of pairs in bit-reversed order */
static SP_NOINLINE void sp_input(uint32_t *x, const float *centre, const float *w, int half, float g)
{
    int j;
    memset(x, 0, SP_M * sizeof(x[0]));
    for (j = -half; j <= half; j++) {
        int i = (j + XOMBE_FFT) & (XOMBE_FFT - 1);
        int16_t q = (int16_t)round_q(centre[j] * w[half + j] * g);
        memcpy((char *)x + (q15_rbit((uint32_t)(i >> 1)) >> (32 - SP_LOG2)) * 4 + (i & 1) * 2, &q, sizeof q);
    }
}

/* S->v from split output v, |v| <= mx, mx * 2^p in [2^30, 2^31), c = 32767 * 2^32 / (mx * 2^p) */
static SP_NOINLINE void sp_store(xombe_spec *S, const int32_t *v, int p, int32_t c)
{
    int16_t *re = S->v, *im = re + XOMBE_LF_IM;
    int k;
    for (k = 0; k < SP_M / 2; k += 2) {
        sp_st2(re + k, sp_q2(v[k], v[k + 1], p, c));
        sp_st2(im + k, sp_q2(v[SP_M - k], v[SP_M - k - 1], p, c));
    }
    for (k = SP_M / 2; k < SP_M; k += 2) {
        const int32_t *u = v + SP_M + 2 * (k - SP_M / 2);
        sp_st2(re + k, sp_q2(u[0], u[2], p, c));
        sp_st2(im + k, sp_q2(u[1], u[3], p, c));
    }
    sp_st2(re + SP_M, sp_q2(v[2 * SP_M], 0, p, c));
    sp_st2(im + SP_M, 0);
    im[0] = 0; /* the loop took it from x[SP_M] */
}

/* samples block-scaled to XOMBE_CFFT_Q15_PEAK, packed in pairs in bit-reversed order for the block-floating-point
   FFT, split in integer, bins scaled to int16 */
SP_NOINLINE SP_NOCLONE void xombe__kern_spectrum(const float *centre, const float *w, int n, xombe_spec *S)
{
    uint32_t x[2 * SP_M + 2];
    const int32_t *v = (const int32_t *)x;
    float m = 0.0f, g;
    uint32_t mx;
    int32_t c;
    int half = (n - 1) / 2, j, k, s, p;
    for (j = -half; j <= half; j++) {
        float a = fabsf(centre[j] * w[half + j]);
        if (a > m)
            m = a;
    }
    g = q15_gain(m, (float)XOMBE_CFFT_Q15_PEAK);
    S->sc = 0.0f;
    if (g == 0.0f)
        goto zero;
    sp_input(x, centre, w, half, g);
    s = xombe__cfft_q15(x, SP_M, xombe__fft_tw_q15, XOMBE_CFFT_Q15_PEAK);
    /* the bound from the top halves is within 2^-10 of the largest component, else an exact scan */
    mx = sp_split(x);
    if (mx >= 1024) {
        mx = (mx + 1) << 16;
    } else {
        for (mx = 0, k = 0; k < 2 * SP_M + 2; k++)
            if (sp_abs(v[k]) > mx)
                mx = sp_abs(v[k]);
        if (mx == 0)
            goto zero;
    }
    /* |v| * 2^p <= mx * 2^p < 2^31, so |v| * 2^p * c / 2^32 <= 32767 after rounding for either sign. the largest
       component rounds to 32767 for an exact mx */
    p = q15_clz32(mx) - 1;
    c = (int32_t)(((uint64_t)32767 << 32) / (mx << p));
    sp_store(S, v, p, c);
    /* x held the transform of the scaled input times 2^-s, v is Q14 of it */
    S->sc = sp_pow2(18 + s - p) / ((float)c * g);
    return;
zero:
    memset(S->v, 0, sizeof(S->v));
}

/* exact |S[k]|^2 < 2^31, one rounding to float */
void xombe__kern_bin_power(const xombe_spec *S, float sum_w2, float *P)
{
    const int16_t *re = S->v, *im = re + XOMBE_LF_IM;
    float scale = S->sc * S->sc / (XOMBE_FFT * sum_w2);
    int k;
    for (k = 0; k < XOMBE_SPEC_BINS; k++)
        P[k] = (float)((int32_t)re[k] * re[k] + (int32_t)im[k] * im[k]) * scale;
}

/* prefix sums of |S|^2 over kmin..kmax. the line fit's residual cancels in 31-bit sums, so a range whose peak is
   below 16384 is shifted up into e->kq.spec, otherwise the bins are read in place */
void xombe__kern_lf_prepare(xombe_enc *e, const xombe_spec *S, int kmin, int kmax, xombe_lf_spec *sp)
{
    const int16_t *re = S->v, *im = re + XOMBE_LF_IM;
    uint64_t *ep = e->kq.ep, acc = 0;
    uint32_t det = 0, d;
    int32_t f;
    int k, l;
    /* every component fits q + 1 signed bits for the top detector bit q (bit 0 of each half may hold the other
       half's sign), shifted up by 15 - q it stays within int16 */
    for (k = kmin; k <= kmax; k += 2) /* may take in bin kmax + 1 <= XOMBE_SPEC_BINS - 1 */
        det |= xombe__q15_det(q15_ld2(re + k)) | xombe__q15_det(q15_ld2(im + k));
    d = (det | det >> 16) & 0xfffeu;
    l = d ? q15_clz32(d) - 16 : 14;
    f = (int32_t)1 << l;
    ep[kmin] = 0;
    for (k = kmin; k <= kmax; k++) {
        int32_t r = re[k] * f, i = im[k] * f;
        if (l) {
            e->kq.spec.v[k] = (int16_t)r;
            e->kq.spec.v[XOMBE_LF_IM + k] = (int16_t)i;
        }
        acc += (uint32_t)(r * r) + (uint32_t)(i * i);
        ep[k + 1] = acc;
    }
    sp->S = S;
    if (l) {
        e->kq.spec.sc = S->sc * sp_pow2(-l);
        sp->S = &e->kq.spec;
    }
    sp->g = sp->S->sc > 0.0f ? 1.0f / sp->S->sc : 0.0f;
}
