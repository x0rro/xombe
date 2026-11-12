#include <math.h>
#include <string.h>
#include "../internal.h"
#include "../tables/dsp_tables.h"
#include "../dsp/biquad.h"
#if XOMBE_KERNELS_Q15
#include "kern_q15.h"
#include "../dsp/fft_q15.h"
#endif

#define EXPLAINED_MIN 0.99f
#define DUAL_GATE_SHARE 0.7f
#define DUAL_GATE_SPAN 0.10f
#define DUAL_TWIST_DB 10.0f
#define DUAL_PREFERENCE 0.01f
#define SINGLE_ID_HZ 31.25f /* single-tone IDs are multiples of this frequency */
#define GATE_LO_HZ 100.0f   /* dual-tone gate reference band GATE_LO_HZ..GATE_HI_HZ */
#define GATE_HI_HZ 3900.0f
#define BIN_HZ (XOMBE_FS / XOMBE_TONE_FFT)

XOMBE_STATIC_ASSERT(XOMBE_TONE_FFT >= 2 * XOMBE_TONE_SEG, tone_fft_pads_the_segment);

#if XOMBE_KERNELS_Q15
#define DUAL_SLACK 1e-3f        /* rounding allowance on an explained fraction of 1 */
#define TWIST_RATIO 3.16227766f /* 10^(DUAL_TWIST_DB/20) */
#define GATE_NUM 7u             /* DUAL_GATE_SHARE = GATE_NUM/GATE_DEN */
#define GATE_DEN 10u
typedef uint64_t pow_sum;

/* sin and cos of 2*pi*t for |t| < 2^21, reduced exactly to a quarter turn q/4 plus f, |f| <= 1/8, Taylor
   polynomials to a^9 and a^8 in a = 2*pi*f, truncation under 1e-9 */
static void sincos_turn(float t, float *s, float *c)
{
    int q = round_q(4.0f * t);
    float a = 2.0f * XOMBE_PI_F * (t - 0.25f * (float)q), a2 = a * a, sn, cs;
    sn = a * (1.0f + a2 * (-1.0f / 6 + a2 * (1.0f / 120 + a2 * (-1.0f / 5040 + a2 * (1.0f / 362880)))));
    cs = 1.0f + a2 * (-0.5f + a2 * (1.0f / 24 + a2 * (-1.0f / 720 + a2 * (1.0f / 40320))));
    switch (q & 3) {
    case 0: *s = sn; *c = cs; break;
    case 1: *s = cs; *c = -sn; break;
    case 2: *s = -sn; *c = -cs; break;
    default: *s = -cs; *c = sn; break;
    }
}
#else
typedef float pow_sum;
#endif

/* one sinusoid of the fit */
typedef struct {
    float hz;     /* Hz, w below is radians/sample */
    float rc, rs; /* sum of x[n]*cos(w*n) and of x[n]*sin(w*n) */
    float h[2];   /* cos and sin of w/2 */
    float m[2];   /* cos and sin of (2*len - 1)*w/2 */
    float g[4];   /* own Gram block gram(c, c), set by fit_component */
} component;

#if XOMBE_KERNELS_Q15
/* hz, h, and m of c */
static void angles(int len, float hz, component *c)
{
    float u = hz * (0.5f / XOMBE_FS); /* w/2 in turns */
    c->hz = hz;
    sincos_turn(u, &c->h[1], &c->h[0]);
    sincos_turn((float)(2 * len - 1) * u, &c->m[1], &c->m[0]);
}
#endif

/* correlates x with cos and sin at hz, Goertzel run from x[len - 1] down so its final states need no rotation,
   in Reinsch's form (d = s[n] -/+ s[n - 1]), as plain Goertzel loses about 5 bits near 0 Hz and fs/2 */
static void correlate(const float *x, int len, float hz, component *c)
{
#if XOMBE_KERNELS_Q15
    float hc, hs;
    angles(len, hz, c);
    hc = c->h[0];
    hs = c->h[1];
#else
    float w = 2.0f * XOMBE_PI_F * hz / XOMBE_FS, hc = cosf(0.5f * w), hs = sinf(0.5f * w);
#endif
    float cw = hc * hc - hs * hs, s = 0.0f, d = 0.0f, prev;
    int n = len - 1;
    if (cw >= 0.0f) {
        float l = -4.0f * hs * hs;
        for (; n > 0; n -= 2) {
            d += x[n] + l * s;
            s += d;
            d += x[n - 1] + l * s;
            s += d;
        }
        if (n == 0) {
            d += x[0] + l * s;
            s += d;
        }
        prev = s - d;
    } else {
        float l = 4.0f * hc * hc;
        for (; n > 0; n -= 2) {
            d = x[n] + l * s - d;
            s = d - s;
            d = x[n - 1] + l * s - d;
            s = d - s;
        }
        if (n == 0) {
            d = x[0] + l * s - d;
            s = d - s;
        }
        prev = d - s;
    }
    c->rc = s - cw * prev;
    c->rs = 2.0f * hs * hc * prev;
#if !XOMBE_KERNELS_Q15
    c->hz = hz;
    c->h[0] = hc;
    c->h[1] = hs;
    c->m[0] = cosf((float)(2 * len - 1) * 0.5f * w);
    c->m[1] = sinf((float)(2 * len - 1) * 0.5f * w);
#endif
}

/* sums of cos(th*n) and sin(th*n) over the segment, a Dirichlet kernel, from cos and sin of th/2 (h) and of
   (2*len - 1)*th/2 (m) */
static void dsum(const float h[2], const float m[2], int len, float *c, float *s)
{
    if (fabsf(h[1]) < 1e-6f) {
        *c = (float)len;
        *s = 0.0f;
        return;
    }
    *c = 0.5f + 0.5f * m[1] / h[1];
    *s = 0.5f * (h[0] - m[0]) / h[1];
}

/* Gram block of a and b, rows cos_a, sin_a, columns cos_b, sin_b, row-major, dsum angles by angle addition */
static void gram(const component *a, const component *b, int len, float g[4])
{
    float cd, sd, cs, ss, hd[2], md[2], hs[2], ms[2];
    hd[0] = a->h[0] * b->h[0] + a->h[1] * b->h[1];
    hd[1] = a->h[1] * b->h[0] - a->h[0] * b->h[1];
    md[0] = a->m[0] * b->m[0] + a->m[1] * b->m[1];
    md[1] = a->m[1] * b->m[0] - a->m[0] * b->m[1];
    hs[0] = a->h[0] * b->h[0] - a->h[1] * b->h[1];
    hs[1] = a->h[1] * b->h[0] + a->h[0] * b->h[1];
    ms[0] = a->m[0] * b->m[0] - a->m[1] * b->m[1];
    ms[1] = a->m[1] * b->m[0] + a->m[0] * b->m[1];
    dsum(hd, md, len, &cd, &sd);
    dsum(hs, ms, len, &cs, &ss);
    g[0] = 0.5f * (cd + cs);
    g[1] = 0.5f * (ss - sd);
    g[2] = 0.5f * (ss + sd);
    g[3] = 0.5f * (cd - cs);
}

/* correlate plus own Gram block, so a trial moving one component recomputes only its own */
static void fit_component(const float *x, int len, float hz, component *c)
{
    correlate(x, len, hz, c);
    gram(c, c, len, c->g);
}

#if XOMBE_KERNELS_Q15
/* correlations of x already taken in one search, keyed by frequency (hz 0: empty slot). pairs share
   frequencies (DTMF rows and columns) and revisit trial steps, so on the worst frames over half the
   correlations repeat */
#define MEMO_LOG2 7
#define MEMO_SLOTS (1 << MEMO_LOG2)
#define MEMO_PROBES 8
typedef struct {
    float hz[MEMO_SLOTS], rc[MEMO_SLOTS], rs[MEMO_SLOTS];
} corr_memo;

/* fit_component through the memo, a full probe run computes without storing */
static void memo_component(corr_memo *mo, const float *x, int len, float hz, component *c)
{
    uint32_t bits;
    int i, probe;
    memcpy(&bits, &hz, sizeof(bits));
    i = (int)((bits * 2654435761u) >> (32 - MEMO_LOG2));
    for (probe = 0; probe < MEMO_PROBES; probe++, i = (i + 1) & (MEMO_SLOTS - 1)) {
        if (mo->hz[i] == hz) {
            angles(len, hz, c);
            c->rc = mo->rc[i];
            c->rs = mo->rs[i];
            gram(c, c, len, c->g);
            return;
        }
        if (mo->hz[i] == 0.0f) {
            fit_component(x, len, hz, c);
            mo->hz[i] = hz;
            mo->rc[i] = c->rc;
            mo->rs[i] = c->rs;
            return;
        }
    }
    fit_component(x, len, hz, c);
}
/* pair_fit's memo parameter and argument, absent in float builds */
#define PAIR_COMPONENT(x, len, hz, c) memo_component(mo, x, len, hz, c)
#define PAIR_MEMO corr_memo *mo,
#define PAIR_MEMO_ARG &memo,
#else
#define PAIR_COMPONENT(x, len, hz, c) fit_component(x, len, hz, c)
#define PAIR_MEMO
#define PAIR_MEMO_ARG
#endif

/* solves 2x2 m*p = r, returns 0 if m is singular */
static int solve2(const float m[4], const float r[2], float p[2])
{
    float det = m[0] * m[3] - m[1] * m[2];
    if (fabsf(det) <= 1e-6f * (m[0] * m[0] + m[3] * m[3]))
        return 0;
    p[0] = (m[3] * r[0] - m[1] * r[1]) / det;
    p[1] = (m[0] * r[1] - m[2] * r[0]) / det;
    return 1;
}

static float magnitude(const float p[2])
{
#if XOMBE_KERNELS_Q15
    return sqrtf(p[0] * p[0] + p[1] * p[1]);
#else
    return hypotf(p[0], p[1]);
#endif
}

/* least-squares fit of nf (1 or 2) sinusoids, returns the fraction of energy explained and each amplitude,
   0 and amplitudes 0 for a singular system. two components solved by blocks (Schur complement of the first block) */
static float fit(float energy, const component *c, int nf, int len, float amp[2])
{
    const float *A = c[0].g;
    float ra[2], pa[2], explained;
    amp[0] = amp[1] = 0.0f;
    ra[0] = c[0].rc;
    ra[1] = c[0].rs;
    if (nf == 1) {
        if (!solve2(A, ra, pa))
            return 0.0f;
        explained = pa[0] * ra[0] + pa[1] * ra[1];
    } else {
        const float *B = c[1].g;
        float C[4], Y[4], y[2], S[4], rhs[2], pb[2], col[2], cj[2], rb[2];
        int j;
        rb[0] = c[1].rc;
        rb[1] = c[1].rs;
        gram(&c[0], &c[1], len, C);
        if (!solve2(A, ra, y))
            return 0.0f;
        for (j = 0; j < 2; j++) { /* Y = A^-1*C, column by column */
            cj[0] = C[j];
            cj[1] = C[2 + j];
            if (!solve2(A, cj, col))
                return 0.0f;
            Y[j] = col[0];
            Y[2 + j] = col[1];
        }
        /* S = B - C^T*Y, rhs = rb - C^T*y */
        S[0] = B[0] - (C[0] * Y[0] + C[2] * Y[2]);
        S[1] = B[1] - (C[0] * Y[1] + C[2] * Y[3]);
        S[2] = B[2] - (C[1] * Y[0] + C[3] * Y[2]);
        S[3] = B[3] - (C[1] * Y[1] + C[3] * Y[3]);
        rhs[0] = rb[0] - (C[0] * y[0] + C[2] * y[1]);
        rhs[1] = rb[1] - (C[1] * y[0] + C[3] * y[1]);
        if (!solve2(S, rhs, pb))
            return 0.0f;
        pa[0] = y[0] - (Y[0] * pb[0] + Y[1] * pb[1]);
        pa[1] = y[1] - (Y[2] * pb[0] + Y[3] * pb[1]);
        amp[1] = magnitude(pb);
        explained = pa[0] * ra[0] + pa[1] * ra[1] + pb[0] * rb[0] + pb[1] * rb[1];
    }
    amp[0] = magnitude(pa);
    return explained / energy;
}

static int in_gate(int k, float nominal)
{
    return fabsf((float)k * BIN_HZ - nominal) <= DUAL_GATE_SPAN * nominal;
}

/* bins *lo..*hi of nominal's dual-tone gate, empty when *lo > *hi. the gate is at least a bin wide and k*BIN_HZ
   is exact, so in_gate holds on one run of k containing the nearest bin, each end a step or two from its estimate */
static void gate_run(float nominal, int *lo, int *hi)
{
    int k0 = (int)(nominal / BIN_HZ + 0.5f), k;
    if (k0 > XOMBE_TONE_FFT / 2)
        k0 = XOMBE_TONE_FFT / 2;
    if (!in_gate(k0, nominal)) {
        *lo = 1;
        *hi = 0;
        return;
    }
    k = (int)ceilf((1.0f - DUAL_GATE_SPAN) * nominal / BIN_HZ);
    k = k < 0 ? 0 : k > k0 ? k0 : k;
    while (k > 0 && in_gate(k - 1, nominal))
        k--;
    while (k < k0 && !in_gate(k, nominal))
        k++;
    *lo = k;
    k = (int)floorf((1.0f + DUAL_GATE_SPAN) * nominal / BIN_HZ);
    k = k < k0 ? k0 : k > XOMBE_TONE_FFT / 2 ? XOMBE_TONE_FFT / 2 : k;
    while (k < XOMBE_TONE_FFT / 2 && in_gate(k + 1, nominal))
        k++;
    while (k > k0 && !in_gate(k, nominal))
        k--;
    *hi = k;
}

static int ad_from_level(float level_db_minus_filter)
{
    int ad = (int)floorf(127.0f + level_db_minus_filter / XOMBE_TONE_DB_PER_AD + 0.5f);
    return ad < 0 ? 0 : ad > 127 ? 127 : ad;
}

static float peak_db(float amp)
{
    return 20.0f * log10f(amp / 32768.0f + 1e-30f);
}

#if !XOMBE_KERNELS_Q15
static void spectrum(xombe_enc *e, const float *x, int len, const float *win, float pw[XOMBE_TONE_FFT / 2 + 1])
{
    xombe_cpx *X = e->tone_spec;
    int n, k;
    /* packed two samples per bin for xombe__rfft_packed */
    memset(e->tone_spec, 0, sizeof(e->tone_spec));
    for (n = 0; n < len; n++) {
        float v = x[n] * win[n];
        if (n & 1)
            X[n >> 1].im = v;
        else
            X[n >> 1].re = v;
    }
    xombe__rfft_packed(X, XOMBE_TONE_FFT, xombe__fft_tw, XOMBE_TW_N);
    for (k = 0; k <= XOMBE_TONE_FFT / 2; k++)
        pw[k] = X[k].re * X[k].re + X[k].im * X[k].im;
}
#else
#define SPEC_M (XOMBE_TONE_FFT / 2)
#define SPEC_LOG2 8
XOMBE_STATIC_ASSERT(1 << SPEC_LOG2 == SPEC_M, spec_log2);

/* bin powers of the windowed segment through the SPEC_M point block-floating-point FFT of sample pairs, in
   place over pw: samples block-scaled to XOMBE_CFFT_Q15_PEAK, real-input split in integer, one power unit for
   every bin of the segment */
static void spectrum(xombe_enc *e, const float *x, int len, const float *win, uint32_t pw[XOMBE_TONE_FFT / 2 + 1])
{
    uint32_t *z = pw;
    float m = 0.0f, g;
    int n, k;
    (void)e;
    for (n = 0; n < len; n++) {
        float a = fabsf(x[n] * win[n]);
        m = a > m ? a : m;
    }
    g = q15_gain(m, (float)XOMBE_CFFT_Q15_PEAK);
    memset(z, 0, SPEC_M * sizeof(*z));
    if (g == 0.0f) {
        pw[SPEC_M] = 0;
        return;
    }
    for (n = 0; n < len; n += 2) {
        int32_t lo = round_q(x[n] * win[n] * g), hi = n + 1 < len ? round_q(x[n + 1] * win[n + 1] * g) : 0;
        z[q15_rbit((uint32_t)(n >> 1)) >> (32 - SPEC_LOG2)] = ((uint32_t)lo & 0xffffu) | (uint32_t)hi << 16;
    }
    xombe__cfft_q15(z, SPEC_M, xombe__tone_tw_q15, XOMBE_CFFT_Q15_PEAK);
    /* bins k and SPEC_M - k from Z[k] and Z[SPEC_M - k], 2*X[k] = E + W^k*O with W^k = xombe__fft_tw_q15[k]:
       components of E and O within +/-2^16, so |2*X| <= 2^17.5 and |2*X|^2 >> 4 <= 2^31. the Q15 product
       sums reach 2^32 and take 64 bits */
    for (k = 0; k <= SPEC_M / 2; k++) {
        uint32_t a = z[k], b = z[(SPEC_M - k) & (SPEC_M - 1)];
        int32_t er = q15_lo(a) + q15_lo(b), ei = q15_hi(a) - q15_hi(b);
        int32_t orr = q15_hi(a) + q15_hi(b), oi = q15_lo(b) - q15_lo(a), wr = 0, wi = -32768, pr, pi, r0, i0, r1, i1;
        if (k < SPEC_M / 2) {
            wr = q15_lo(xombe__fft_tw_q15[k]);
            wi = q15_hi(xombe__fft_tw_q15[k]);
        }
        pr = (int32_t)(((int64_t)wr * orr - (int64_t)wi * oi + 16384) >> 15);
        pi = (int32_t)(((int64_t)wr * oi + (int64_t)wi * orr + 16384) >> 15);
        r0 = er + pr;
        i0 = ei + pi;
        r1 = er - pr;
        i1 = pi - ei;
        pw[k] = (uint32_t)(((uint64_t)((int64_t)r0 * r0) + (uint64_t)((int64_t)i0 * i0)) >> 4);
        pw[SPEC_M - k] = (uint32_t)(((uint64_t)((int64_t)r1 * r1) + (uint64_t)((int64_t)i1 * i1)) >> 4);
    }
}
#endif

/* mean-removed segment into x, returns its energy or 0 below the -50 dBFS sinusoid floor. above it, fills pw
   with the zero-padded bin powers of the segment under win */
static float prepare(xombe_enc *e, const int16_t *seg, int len, const float *win, float *x,
                     xombe_tone_pow pw[XOMBE_TONE_FFT / 2 + 1])
{
    static const float min_amp = 32768.0f * 0.0031622777f; /* -50 dBFS peak */
    float mean, energy = 0.0f;
    long sum = 0;
    int n;
    for (n = 0; n < len; n++)
        sum += seg[n];
    mean = (float)sum / len;
    for (n = 0; n < len; n++) {
        x[n] = (float)seg[n] - mean;
        energy += x[n] * x[n];
    }
    if (energy / len < 0.5f * min_amp * min_amp)
        return 0.0f;
    spectrum(e, x, len, win, pw);
    return energy;
}

static pow_sum band_total(const xombe_tone_pow *pw)
{
    pow_sum total = 0;
    int k;
    for (k = (int)ceilf(GATE_LO_HZ / BIN_HZ); k <= (int)floorf(GATE_HI_HZ / BIN_HZ); k++)
        total += pw[k];
    return total;
}

/* at least DUAL_GATE_SHARE of total within DUAL_GATE_SPAN of either nominal frequency */
static int gate_pass(const xombe_tone_pow *pw, const float nominal[2], pow_sum total)
{
    pow_sum near = 0;
    int k, r0lo, r0hi, r1lo, r1hi, klo, khi;
    /* summed in ascending bins over the union of the two runs */
    gate_run(nominal[0], &r0lo, &r0hi);
    gate_run(nominal[1], &r1lo, &r1hi);
    klo = r0lo <= r0hi ? r0lo : r1lo;
    khi = r0lo <= r0hi ? r0hi : r1hi;
    if (r1lo <= r1hi) {
        klo = r1lo < klo ? r1lo : klo;
        khi = r1hi > khi ? r1hi : khi;
    }
    for (k = klo; k <= khi; k++)
        if ((k >= r0lo && k <= r0hi) || (k >= r1lo && k <= r1hi))
            near += pw[k];
#if XOMBE_KERNELS_Q15
    return near * GATE_DEN >= total * GATE_NUM;
#else
    return !(near < DUAL_GATE_SHARE * total);
#endif
}

#define FINE_FFT 8192 /* single-tone peak grid */
#define FINE_HZ (XOMBE_FS / FINE_FFT)
#define FINE_LO 103  /* FINE_FFT bins in GATE_LO_HZ..GATE_HI_HZ */
#define FINE_HI 3993
#define REFINE_SPAN 0.4f /* single-tone refinement, +/- fine bins around the peak */
#define REFINE_ITERS 12
/* explained fraction moves at most ~0.15 per Hz for a 160-sample segment (the fit subspace turns at most
   2*pi*rms(n)/fs radians per Hz), so refining within REFINE_SPAN (0.39 Hz) gains under 0.06 and a fit this far
   below the threshold cannot reach it */
#define REFINE_REACH 0.1f
#define DUAL_STEP 0.0025f /* relative dual search step, steps are 4, 2, and 1 of these */
#define DUAL_STEPS_MAX 10 /* +/-2.5% of nominal in steps */
#define PAIR_TRIALS 24    /* pair_fit trials, 2 passes x 3 steps x 2 components x 2 signs */

float xombe__tone_hpf_gain_db(float hz)
{
    const float p1 = (float)HPF_POLE1, a1 = (float)HPF_A1, a2 = (float)HPF_A2;
    float w = 2.0f * XOMBE_PI_F * hz / XOMBE_FS, cr = cosf(w), ci = -sinf(w), c2r = cosf(2.0f * w),
          c2i = -sinf(2.0f * w);
    float zr = 1.0f - cr, zi = -ci, nr = zr * zr - zi * zi, ni = 2.0f * zr * zi;
    float n3r = nr * zr - ni * zi, n3i = nr * zi + ni * zr;
    float d1r = 1.0f - p1 * cr, d1i = -p1 * ci, d2r = 1.0f - a1 * cr + a2 * c2r, d2i = -a1 * ci + a2 * c2i;
    float num = n3r * n3r + n3i * n3i, den = (d1r * d1r + d1i * d1i) * (d2r * d2r + d2i * d2i);
    return 10.0f * log10f(num / den);
}

/* windowed segment power at FINE_FFT bin k */
static float fine_power(const float *xw, int len, int k)
{
    component c;
    correlate(xw, len, (float)k * FINE_HZ, &c);
    return c.rc * c.rc + c.rs * c.rs;
}

/* largest pw bin in the gate band, in bins, placed by a parabola through log magnitudes: the Hann peak is
   near-Gaussian, so within 0.1 Hz */
static float coarse_peak(const xombe_tone_pow *pw)
{
    int k, kc = -1, klo = (int)ceilf(GATE_LO_HZ / BIN_HZ), khi = (int)floorf(GATE_HI_HZ / BIN_HZ);
    float la, lb, lc, d;
    for (k = klo; k <= khi; k++)
        if (kc < 0 || pw[k] > pw[kc])
            kc = k;
    la = 0.5f * logf((float)pw[kc - 1] + 1e-30f);
    lb = 0.5f * logf((float)pw[kc] + 1e-30f);
    lc = 0.5f * logf((float)pw[kc + 1] + 1e-30f);
    d = la - 2.0f * lb + lc;
    return (float)kc + (d != 0.0f ? 0.5f * (la - lc) / d : 0.0f);
}

/* largest FINE_FFT bin in the gate band plus a parabola through log magnitudes, found by climbing the fine grid
   from the coarse peak */
static float peak_hz(const float *x, const float *win, int len, const xombe_tone_pow *pw)
{
    int kf = (int)floorf(coarse_peak(pw) * (float)(FINE_FFT / XOMBE_TONE_FFT) + 0.5f), n;
    float xw[XOMBE_TONE_SEG], pm, pl, pr, a, b, c, den;
    for (n = 0; n < len; n++)
        xw[n] = x[n] * win[n];
    kf = kf < FINE_LO ? FINE_LO : kf > FINE_HI ? FINE_HI : kf;
    pm = fine_power(xw, len, kf);
    pl = fine_power(xw, len, kf - 1);
    pr = fine_power(xw, len, kf + 1);
    while (pr > pm && kf < FINE_HI) {
        kf++;
        pl = pm;
        pm = pr;
        pr = fine_power(xw, len, kf + 1);
    }
    while (pl > pm && kf > FINE_LO) {
        kf--;
        pr = pm;
        pm = pl;
        pl = fine_power(xw, len, kf - 1);
    }
    a = 0.5f * logf(pl + 1e-30f);
    b = 0.5f * logf(pm + 1e-30f);
    c = 0.5f * logf(pr + 1e-30f);
    den = a - 2.0f * b + c;
    return ((float)kf + (den != 0.0f ? 0.5f * (a - c) / den : 0.0f)) * FINE_HZ;
}

/* ternary search for the best single-tone fit within REFINE_SPAN fine bins of f1 */
static float fit_refined(const float *x, int len, float energy, float f1, float amp[2])
{
    float lo = f1 - REFINE_SPAN * FINE_HZ, hi = f1 + REFINE_SPAN * FINE_HZ, unused[2];
    component c;
    int it;
    for (it = 0; it < REFINE_ITERS; it++) {
        float fa = lo + (hi - lo) / 3.0f, fb = hi - (hi - lo) / 3.0f, ea, eb;
        fit_component(x, len, fa, &c);
        ea = fit(energy, &c, 1, len, unused);
        fit_component(x, len, fb, &c);
        eb = fit(energy, &c, 1, len, unused);
        if (ea > eb)
            hi = fb;
        else
            lo = fa;
    }
    fit_component(x, len, 0.5f * (lo + hi), &c);
    return fit(energy, &c, 1, len, amp);
}

/* coordinate search in relative steps of 4, 2, and 1 DUAL_STEP, each pass trying each step on the first
   component then the second, up then down */
static float pair_fit(PAIR_MEMO const float *x, int len, float energy, const float nominal[2], float hz[2],
                      float amp[2])
{
    static const int steps[3] = {4, 2, 1};
    component comp[2];
    int d[2] = {0, 0}, seen[PAIR_TRIALS + 1], pass, s, c, sign, nseen = 1;
    float best, trial_amp[2] = {0.0f, 0.0f};
    PAIR_COMPONENT(x, len, nominal[0], &comp[0]);
    PAIR_COMPONENT(x, len, nominal[1], &comp[1]);
    best = fit(energy, comp, 2, len, amp);
    seen[0] = 0;
    for (pass = 0; pass < 2; pass++)
        for (s = 0; s < 3; s++)
            for (c = 0; c < 2; c++)
                for (sign = 1; sign >= -1; sign -= 2) {
                    int nd = d[c] + sign * steps[s], key, i;
                    component keep = comp[c];
                    float te;
                    if (nd > DUAL_STEPS_MAX || nd < -DUAL_STEPS_MAX)
                        continue;
                    key = c ? d[0] * 64 + nd : nd * 64 + d[1];
                    for (i = 0; i < nseen && seen[i] != key; i++)
                        ;
                    if (i < nseen) /* a pair already evaluated explains no more than best */
                        continue;
                    seen[nseen++] = key;
                    PAIR_COMPONENT(x, len, nominal[c] * (1.0f + DUAL_STEP * (float)nd), &comp[c]);
                    te = fit(energy, comp, 2, len, trial_amp);
                    if (te > best) {
                        best = te;
                        d[c] = nd;
                        amp[0] = trial_amp[0];
                        amp[1] = trial_amp[1];
                    } else {
                        comp[c] = keep;
                    }
                }
    hz[0] = comp[0].hz;
    hz[1] = comp[1].hz;
    return best;
}

/* components more than DUAL_TWIST_DB apart, amp > 0 */
static int twisted(const float amp[2])
{
#if XOMBE_KERNELS_Q15
    return amp[0] > TWIST_RATIO * amp[1] || amp[1] > TWIST_RATIO * amp[0];
#else
    return fabsf(20.0f * log10f(amp[0] / amp[1])) > DUAL_TWIST_DB;
#endif
}

void xombe__tone_prepare(xombe_enc *e, const int16_t *seg, const xombe_tone_search *ts, xombe_tone_seg *sg)
{
    sg->energy = prepare(e, seg, ts->len, ts->win, sg->x, sg->pw);
}

int xombe__tone_search(const xombe_tone_seg *sg, const xombe_tone_search *ts, xombe_tone_hit *hit)
{
    const float *x = sg->x;
    const xombe_tone_pow *pw = sg->pw;
    pow_sum total;
    float energy = sg->energy, e1 = 0.0f, best = -1.0f, best_level = 0.0f, best_res = 0.0f;
    int i, single = 0, dual = -1, len = ts->len, npairs = ts->npairs;
#if XOMBE_KERNELS_Q15
    corr_memo memo;
    int memo_ready = 0;
#endif
    hit->id = -1;
    if (energy == 0.0f)
        return 0;
    if (ts->singles) {
        float f1 = ts->fine_peak ? peak_hz(x, ts->win, len, pw) : coarse_peak(pw) * BIN_HZ, amp[2];
        int cand = (int)floorf(f1 / SINGLE_ID_HZ + 0.5f);
        component c;
        fit_component(x, len, f1, &c);
        e1 = fit(energy, &c, 1, len, amp);
        if (ts->fine_peak && e1 >= ts->emin - REFINE_REACH)
            e1 = fit_refined(x, len, energy, f1, amp);
        if (e1 >= ts->emin && cand >= 5 && cand <= 122) {
            single = 1;
            hit->id = cand;
            hit->e = e1;
            hit->level = peak_db(amp[0] > 1e-9f ? amp[0] : 1e-9f) - xombe__tone_hpf_gain_db(f1);
            hit->res = energy * (1.0f - e1);
        }
    }
#if XOMBE_KERNELS_Q15
    /* a least-squares fit explains at most the whole energy, so a dual cannot beat e1 + DUAL_PREFERENCE past 1
       plus rounding slack */
    if (single && e1 + DUAL_PREFERENCE > 1.0f + DUAL_SLACK)
        npairs = 0;
#endif
    total = band_total(pw);
    for (i = 0; i < npairs; i++) {
        const xombe_tone_pair *p = &ts->pairs[i];
        float hz[2], amp[2], ex;
        if (!gate_pass(pw, p->hz, total))
            continue;
#if XOMBE_KERNELS_Q15
        if (!memo_ready) {
            memset(memo.hz, 0, sizeof(memo.hz));
            memo_ready = 1;
        }
#endif
        ex = pair_fit(PAIR_MEMO_ARG x, len, energy, p->hz, hz, amp);
        if (ex < ts->emin || amp[0] <= 0.0f || amp[1] <= 0.0f || twisted(amp))
            continue;
        if (ex > best) {
            best = ex;
            dual = p->id;
            best_level = 0.5f * ((peak_db(amp[0]) - xombe__tone_hpf_gain_db(hz[0])) +
                                 (peak_db(amp[1]) - xombe__tone_hpf_gain_db(hz[1])));
            best_res = energy * (1.0f - ex);
        }
    }
    if (dual >= 0 && (!single || best > e1 + DUAL_PREFERENCE)) {
        hit->id = dual;
        hit->e = best;
        hit->level = best_level;
        hit->res = best_res;
    }
    return hit->id >= 0;
}

int xombe__tone_detect(xombe_enc *e, const int16_t *seg, int *id, int *ad)
{
    xombe_tone_pair pairs[XOMBE_TONE_DUALS];
    xombe_tone_search ts;
    xombe_tone_seg sg;
    xombe_tone_hit h;
    int i;
    ts.len = XOMBE_TONE_SEG;
    ts.win = xombe__tone_win;
    ts.emin = EXPLAINED_MIN;
    ts.singles = 1;
    ts.fine_peak = 0;
    ts.pairs = pairs;
    ts.npairs = XOMBE_TONE_DUALS;
    xombe__tone_prepare(e, seg, &ts, &sg);
    if (sg.energy == 0.0f)
        return 0;
    for (i = 0; i < ts.npairs; i++) {
        pairs[i].id = (short)(XOMBE_TONE_DUAL_FIRST + i);
        xombe__tone_freqs(pairs[i].id, &pairs[i].hz[0], &pairs[i].hz[1]);
    }
    if (!xombe__tone_search(&sg, &ts, &h))
        return 0;
    *id = h.id;
    *ad = ad_from_level(h.level);
    return 1;
}

int xombe__tone_decide(xombe_enc *e, int detected, int det_id, int det_ad, int *id, int *ad)
{
    int tone = 0;
    if (detected && det_id == e->tone_prev_id) {
        e->tone_confirmed = 1;
        e->tone_hold_used = 0;
        e->tone_last_id = det_id;
        e->tone_last_ad = det_ad;
        tone = 1;
    } else if (!detected && e->tone_confirmed && !e->tone_hold_used) {
        e->tone_hold_used = 1;
        tone = 1;
    } else {
        e->tone_confirmed = 0;
    }
    e->tone_prev_id = detected ? det_id : -1;
    *id = e->tone_last_id;
    *ad = e->tone_last_ad;
    return tone;
}
