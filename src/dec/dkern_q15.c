/* q15 decoder kernels (dkernels.h, XOMBE_KERNELS=q15): unvoiced transforms through the block-floating-point Q15
   FFT. the wavetable keeps the float one, Q15's ~65 dB on its sparse spectra moves PCM by up to 10 LSB */
#include <float.h>
#include <math.h>
#include <string.h>
#include "dkernels.h"
#include "../tables/dsp_tables.h"
#include "../dsp/fft_q15.h"
#include "../dsp/fmath.h"

#define GRAIN_RAMP XOMBE_GRAIN_RAMP
#define WT_SIZE XOMBE_WT_SIZE
#define UV_FFT XOMBE_UV_FFT
#define UV_BINS XOMBE_UV_BINS
#define QM (UV_FFT / 2) /* complex points of the unvoiced transforms */
#define QM_LOG2 7
#define TW_STEP (XOMBE_TW_N / UV_FFT) /* xombe__fft_tw[k*TW_STEP] = exp(-2 pi i k/UV_FFT) */

XOMBE_STATIC_ASSERT(WT_SIZE <= XOMBE_IRFFT_MAX, wavetable_fits_irfft);
XOMBE_STATIC_ASSERT(XOMBE_TW_N % WT_SIZE == 0, twiddles_cover_wavetable);
XOMBE_STATIC_ASSERT(1 << QM_LOG2 == QM, qm_log2);
XOMBE_STATIC_ASSERT(XOMBE_TW_N % UV_FFT == 0, twiddles_cover_uv_fft);
XOMBE_STATIC_ASSERT(2 * XOMBE_HALF_MAX <= UV_FFT, segment_fits_uv_fft);
XOMBE_STATIC_ASSERT(XOMBE_LMAX <= XOMBE_HALF_MAX, sqrt_table_covers_L);

/* exp(-2 pi i j/QM) for j < QM/2, Q15 rounded and held below 1, as int16 pairs (re, im) */
static const uint32_t tw_q15[QM / 2] = {
    0x00007fffu, 0xf9b87fd9u, 0xf3747f62u, 0xed387e9du, 0xe7077d8au, 0xe0e67c2au,
    0xdad87a7du, 0xd4e17885u, 0xcf047642u, 0xc94673b6u, 0xc3a970e3u, 0xbe326dcau,
    0xb8e36a6eu, 0xb3c066d0u, 0xaecc62f2u, 0xaa0a5ed7u, 0xa57e5a82u, 0xa12955f6u,
    0x9d0e5134u, 0x99304c40u, 0x9592471du, 0x923641ceu, 0x8f1d3c57u, 0x8c4a36bau,
    0x89be30fcu, 0x877b2b1fu, 0x85832528u, 0x83d61f1au, 0x827618f9u, 0x816312c8u,
    0x809e0c8cu, 0x80270648u, 0x80000000u, 0x8027f9b8u, 0x809ef374u, 0x8163ed38u,
    0x8276e707u, 0x83d6e0e6u, 0x8583dad8u, 0x877bd4e1u, 0x89becf04u, 0x8c4ac946u,
    0x8f1dc3a9u, 0x9236be32u, 0x9592b8e3u, 0x9930b3c0u, 0x9d0eaeccu, 0xa129aa0au,
    0xa57ea57eu, 0xaa0aa129u, 0xaecc9d0eu, 0xb3c09930u, 0xb8e39592u, 0xbe329236u,
    0xc3a98f1du, 0xc9468c4au, 0xcf0489beu, 0xd4e1877bu, 0xdad88583u, 0xe0e683d6u,
    0xe7078276u, 0xed388163u, 0xf374809eu, 0xf9b88027u
};

/* nearest integer (ties to even), |v| < 2^22 */
static int32_t round_q(float v)
{
    uint32_t u;
    v += 12582912.0f;
    memcpy(&u, &v, sizeof(u));
    return (int32_t)(u & 0x7fffffu) - 0x400000;
}

static uint32_t pack(int32_t lo, int32_t hi)
{
    return ((uint32_t)lo & 0xffffu) | (uint32_t)hi << 16;
}

static float block_gain(float m)
{
    float g = m > 0.0f ? (float)XOMBE_CFFT_Q15_PEAK / m : 0.0f;
    return g <= FLT_MAX ? g : 0.0f;
}

static uint32_t rev(int k)
{
    return q15_rbit((uint32_t)k) >> (32 - QM_LOG2);
}

/* xombe__irfft for UV_FFT points: float pre-twiddle, then the inverse complex transform as the forward one on
   swapped halves (swap(FFT(swap(z))) is the unscaled inverse) */
static void irfft_q15(const xombe_cpx *in, float *out)
{
    xombe_cpx z[QM];
    uint32_t x[QM];
    float m = 0.0f, g, sc;
    int k, s;
    /* points k and QM - k from the same two bins, w[QM - k] = -conj(w[k]) */
    for (k = 0; k <= QM / 2; k++) {
        xombe_cpx a = in[k], b = in[QM - k], w = xombe__fft_tw[k * TW_STEP];
        float er = 0.5f * (a.re + b.re), ei = 0.5f * (a.im - b.im), dr = 0.5f * (a.re - b.re), di = 0.5f * (a.im + b.im);
        float pr = dr * w.re + di * w.im, pi = di * w.re - dr * w.im;
        z[k].re = er - pi;
        z[k].im = ei + pr;
        m = fmaxf(m, fmaxf(fabsf(z[k].re), fabsf(z[k].im)));
        if (k == 0 || k == QM / 2)
            continue;
        /* bins swapped: er, -ei, -dr, and di against -conj(w), the same pr */
        pi = dr * w.im - di * w.re;
        z[QM - k].re = er - pi;
        z[QM - k].im = pr - ei;
        m = fmaxf(m, fmaxf(fabsf(z[QM - k].re), fabsf(z[QM - k].im)));
    }
    g = block_gain(m);
    if (g == 0.0f) {
        memset(out, 0, sizeof(float) * UV_FFT);
        return;
    }
    for (k = 0; k < QM; k++)
        x[rev(k)] = pack(round_q(z[k].im * g), round_q(z[k].re * g));
    s = xombe__cfft_q15(x, QM, tw_q15, XOMBE_CFFT_Q15_PEAK);
    sc = (float)(1 << s) / (g * (float)QM);
    for (k = 0; k < QM; k++) {
        out[2 * k] = (float)q15_hi(x[k]) * sc;
        out[2 * k + 1] = (float)q15_lo(x[k]) * sc;
    }
}

/* xombe__rfft of UV_FFT samples q[] times 1/g (q block-scaled to XOMBE_CFFT_Q15_PEAK), bins 0..QM */
static void rfft_q15(const int16_t *q, float g, xombe_cpx *X)
{
    uint32_t x[QM];
    float h;
    int k, s;
    for (k = 0; k < QM; k++)
        x[rev(k)] = q15_ld2(q + 2 * k);
    s = xombe__cfft_q15(x, QM, tw_q15, XOMBE_CFFT_Q15_PEAK);
    h = 0.5f * (float)(1 << s) / g;
    for (k = 0; k <= QM / 2; k++) {
        uint32_t a = x[k], b = x[(QM - k) & (QM - 1)];
        float er = (float)(q15_lo(a) + q15_lo(b)), ei = (float)(q15_hi(a) - q15_hi(b));
        float orr = (float)(q15_hi(a) + q15_hi(b)), oi = (float)(q15_lo(b) - q15_lo(a)), pr, pi;
        xombe_cpx w;
        if (k < QM / 2) {
            w = xombe__fft_tw[k * TW_STEP];
        } else {
            w.re = 0.0f;
            w.im = -1.0f;
        }
        pr = w.re * orr - w.im * oi;
        pi = w.re * oi + w.im * orr;
        X[k].re = h * (er + pr);
        X[k].im = h * (ei + pi);
        X[QM - k].re = h * (er - pr);
        X[QM - k].im = h * (pi - ei);
    }
}

/* float wavetable transform: z of QM points (bit-reversed order) transformed in place, forward, radix-2^2 passes
   as in fft_q15.h, twiddles exp(-2 pi i j/QM) = xombe__fft_tw[j*2*TW_STEP] */
static void cfft_wt(xombe_cpx *x)
{
    const int ts = 2 * TW_STEP;
    int len, i, j;
    for (i = 0; i < QM; i += 4) {
        xombe_cpx x0 = x[i], x1 = x[i + 1], x2 = x[i + 2], x3 = x[i + 3];
        float ar = x0.re + x1.re, ai = x0.im + x1.im, br = x0.re - x1.re, bi = x0.im - x1.im;
        float cr = x2.re + x3.re, ci = x2.im + x3.im, dr = x2.re - x3.re, di = x2.im - x3.im;
        x[i].re = ar + cr;
        x[i].im = ai + ci;
        x[i + 2].re = ar - cr;
        x[i + 2].im = ai - ci;
        x[i + 1].re = br + di; /* b + (-i)*d */
        x[i + 1].im = bi - dr;
        x[i + 3].re = br - di;
        x[i + 3].im = bi + dr;
    }
    for (len = 8; 2 * len <= QM; len <<= 2) {
        int q = len >> 1, sa = ts * (QM / len), sb = sa >> 1;
        for (i = 0; i < QM; i += 2 * len) {
            xombe_cpx *p = x + i;
            for (j = 0; j < q; j++) {
                xombe_cpx wa = xombe__fft_tw[j * sa], wb = xombe__fft_tw[j * sb];
                xombe_cpx u0 = p[j], u1 = p[j + q], u2 = p[j + 2 * q], u3 = p[j + 3 * q];
                float t1r = u1.re * wa.re - u1.im * wa.im, t1i = u1.re * wa.im + u1.im * wa.re;
                float t3r = u3.re * wa.re - u3.im * wa.im, t3i = u3.re * wa.im + u3.im * wa.re;
                float y0r = u0.re + t1r, y0i = u0.im + t1i, y1r = u0.re - t1r, y1i = u0.im - t1i;
                float y2r = u2.re + t3r, y2i = u2.im + t3i, y3r = u2.re - t3r, y3i = u2.im - t3i;
                float t2r = y2r * wb.re - y2i * wb.im, t2i = y2r * wb.im + y2i * wb.re;
                float t4r = y3r * wb.re - y3i * wb.im, t4i = y3r * wb.im + y3i * wb.re;
                p[j].re = y0r + t2r;
                p[j].im = y0i + t2i;
                p[j + 2 * q].re = y0r - t2r;
                p[j + 2 * q].im = y0i - t2i;
                p[j + q].re = y1r + t4i; /* y1 + (-i)*t4 */
                p[j + q].im = y1i - t4r;
                p[j + 3 * q].re = y1r - t4i;
                p[j + 3 * q].im = y1i + t4r;
            }
        }
    }
    if (len <= QM) {
        int q = QM >> 1;
        for (j = 0; j < q; j++) {
            xombe_cpx w = xombe__fft_tw[j * ts], a = x[j], b = x[j + q];
            float tr = b.re * w.re - b.im * w.im, ti = b.re * w.im + b.im * w.re;
            x[j].re = a.re + tr;
            x[j].im = a.im + ti;
            x[j + q].re = a.re - tr;
            x[j + q].im = a.im - ti;
        }
    }
}

XOMBE_STATIC_ASSERT(WT_SIZE == UV_FFT && XOMBE_LMAX < QM / 2, wavetable_bins_below_mirror);

/* xombe__irfft of a wavetable spectrum, bins above XOMBE_LMAX zero: pre-twiddle into swapped halves, the forward
   transform above, scaled by 1/QM */
static void irfft_wt(const xombe_cpx *in, float *out)
{
    xombe_cpx z[QM];
    int k;
    memset(z, 0, sizeof(z));
    for (k = 1; k <= XOMBE_LMAX; k++) {
        xombe_cpx a = in[k], w = xombe__fft_tw[k * TW_STEP];
        /* in[QM - k] = 0, so z[k] and z[QM - k] both come from in[k] alone, w[QM - k] = -conj(w) */
        float er = 0.5f * a.re, ei = 0.5f * a.im;
        float pr = er * w.re + ei * w.im, pi = ei * w.re - er * w.im;
        z[rev(k)].re = ei + pr; /* swapped (im, re) */
        z[rev(k)].im = er - pi;
        pi = er * w.im - ei * w.re;
        z[rev(QM - k)].re = pr - ei;
        z[rev(QM - k)].im = er - pi;
    }
    cfft_wt(z);
    for (k = 0; k < QM; k++) {
        out[2 * k] = z[k].im * (1.0f / QM);
        out[2 * k + 1] = z[k].re * (1.0f / QM);
    }
}

void xombe__dk_wavetable(const xombe_set *s, const float *M, const uint8_t *vh, float tab[XOMBE_WT_SIZE])
{
    xombe_cpx spec[WT_SIZE / 2 + 1];
    int l;
    memset(spec, 0, sizeof(spec));
    for (l = 1; l <= s->L && l <= XOMBE_LMAX; l++) {
        float psi, cp, sp;
        if (!vh[l] || M[l] == 0.0f)
            continue;
        psi = s->phi[l] - 2.0f * XOMBE_PI_F * l * XOMBE_EPOCH_REF * s->f0;
        xombe__sincosf(psi, &sp, &cp);
        spec[l].re = WT_SIZE * M[l] * cp;
        spec[l].im = WT_SIZE * M[l] * sp;
    }
    irfft_wt(spec, tab);
}

static float ramp(int j, float cf, float sf)
{
    return 0.5f * (1.0f - (xombe__ramp_cos[j] * cf - xombe__ramp_sin[j] * sf));
}

/* tab at phase ph (cycles in 2^-32 units), linear between entries */
static float tab_at(const float *tab, uint32_t ph)
{
    int i = (int)(ph >> 24);
    float t0 = tab[i], t1 = tab[(i + 1) & (WT_SIZE - 1)];
    return t0 + (float)(ph & 0xffffffu) * (1.0f / 16777216.0f) * (t1 - t0);
}

XOMBE_STATIC_ASSERT(WT_SIZE == 256, wavetable_index_is_top_byte);

/* phase (n - P)*fhat as a 32-bit fraction accumulated from n0, fade-in, flat, and fade-out runs split ahead */
void xombe__dk_add_grain(float *buf, int N, float P, float fhat, float period, const float *tab)
{
    float E = P + period, th;
    int nb = (int)ceilf(P), ne = (int)floorf(E), n0 = nb, n1 = ne + GRAIN_RAMP, nin = (int)ceilf(P + GRAIN_RAMP), n;
    float ai = XOMBE_PI_F * ((float)nb - P) / GRAIN_RAMP, ao = XOMBE_PI_F * (E - (float)ne) / GRAIN_RAMP;
    float cs, ss, ce, se;
    uint32_t ph, inc = (uint32_t)(fhat * 4294967296.0f); /* fhat < 1/16 */
    xombe__sincosf(ai, &ss, &cs);
    xombe__sincosf(ao, &se, &ce);
    if (n0 < 0)
        n0 = 0;
    if (n1 > N + XOMBE_TAIL - 1)
        n1 = N + XOMBE_TAIL - 1;
    th = ((float)n0 - P) * fhat;
    th -= floorf(th);
    ph = ((uint32_t)(th * 16777216.0f) & 0xffffffu) << 8;
    for (n = n0; n <= n1 && n < nin; n++, ph += inc)
        buf[n] += ramp(n - nb, cs, ss) * tab_at(tab, ph);
    for (; n <= n1 && n <= ne; n++, ph += inc)
        buf[n] += tab_at(tab, ph);
    for (; n <= n1; n++, ph += inc)
        buf[n] += ramp(GRAIN_RAMP + ne - n, ce, se) * tab_at(tab, ph);
}

/* noise samples to the FFT's input peak, the fades are at most 1 */
#define NOISE_Q ((float)XOMBE_CFFT_Q15_PEAK / 32768.0f)

static float fade(int n, float rs)
{
    return xombe__sqrt_n[n] * rs;
}

static int round_half_up(float x)
{
    return (int)floorf(x + 0.5f);
}

void xombe__dk_unvoiced(const int16_t *seg, int N, const xombe_set *cur, const float *M, const uint8_t *vh,
                        float avg_f0, float y[XOMBE_UV_FFT])
{
    int16_t q[UV_FFT];
    xombe_cpx X[UV_BINS], Y[UV_BINS];
    float beta = UV_FFT * cur->f0, rs = xombe__rsqrt_half[N - XOMBE_HALF_MIN], sqrt_L = xombe__sqrt_n[cur->L];
    int n, k, l, upper = 0, low;
    int silence = cur->cls == XOMBE_CLASS_SILENCE;
    memset(q, 0, sizeof(q));
    for (n = 0; n < N; n++)
        q[n] = (int16_t)round_q((float)seg[n] * (fade(n, rs) * NOISE_Q));
    for (; n < 2 * N; n++)
        q[n] = (int16_t)round_q((float)seg[n] * (fade(2 * N - n, rs) * NOISE_Q));
    rfft_q15(q, NOISE_Q, X);
    memset(Y, 0, sizeof(Y));
    for (l = 1; l <= cur->L; l++) {
        int lo = round_half_up((l - 0.5f) * beta), hi = round_half_up((l + 0.5f) * beta);
        if (lo > UV_BINS - 1)
            lo = UV_BINS - 1;
        if (hi > UV_BINS - 1)
            hi = UV_BINS - 1;
        if (l == cur->L)
            upper = hi;
        if (hi <= lo)
            continue;
        if (silence) {
            float g = XOMBE_CN_GAIN * sqrt_L * M[l];
            for (k = lo; k < hi; k++) {
                Y[k].re = X[k].re * g;
                Y[k].im = X[k].im * g;
            }
        } else if (!vh[l]) {
            float e = 0.0f;
            float g;
            for (k = lo; k < hi; k++)
                e += X[k].re * X[k].re + X[k].im * X[k].im;
            e /= (hi - lo);
            if (e <= 0.0f)
                continue;
            g = XOMBE_UV_GAIN * (sqrt_L * xombe__sqrt_n[N]) * M[l] / sqrtf(e);
            for (k = lo; k < hi; k++) {
                Y[k].re = X[k].re * g;
                Y[k].im = X[k].im * g;
            }
        }
    }
    if (upper > XOMBE_UV_TOP_BIN)
        upper = XOMBE_UV_TOP_BIN;
    for (k = upper; k < UV_BINS; k++)
        Y[k].re = Y[k].im = 0.0f;
    low = (int)(UV_FFT * avg_f0);
    for (k = 0; k < low && k < UV_BINS; k++)
        Y[k].re = Y[k].im = 0.0f;
    irfft_q15(Y, y);
}
