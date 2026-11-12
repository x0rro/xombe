/* q15 pitch kernels (kernels.h, XOMBE_KERNELS=q15): int16 block-scaled input, SMLAD low-pass two outputs per
   pass, low-pass output and its windowed product block-shifted to int16, energy and autocorrelation exact in 64
   bits, lags normalized to r[0] = 32767 by one reciprocal (float scale in pacf_scale). quarter-sample periods make
   the interpolation exact in Q2. also the analysis window and the line fit's window table */
#include "kern_q15.h"
#include "../tables/dsp_tables.h"

#define PE_N (XOMBE_PITCH_WIN + 1)                            /* window padded to even */
#define PE_X (XOMBE_PITCH_WIN + 2 * XOMBE_PITCH_LPF_HALF + 1) /* input samples with one pad, even */
#define PE_P (PE_N / 2 + 1)                                   /* pairs of the window, zero-padded by one */

XOMBE_STATIC_ASSERT(PE_N % 2 == 0 && PE_X % 2 == 0, pe_even);
XOMBE_STATIC_ASSERT(PE_N <= XOMBE_PITCH_WIN + 1, pe_rq_fits_pacf);
XOMBE_STATIC_ASSERT((long)XOMBE_WTAB_OVERSAMPLE * XOMBE_FFT * (XOMBE_MAX_WINDOW - 1) < (1L << 23), window_phase_fits);

/* 2^k, -126 <= k <= 127 */
static float pow2f(int k)
{
    uint32_t b = (uint32_t)(k + 127) << 23;
    float f;
    memcpy(&f, &b, sizeof(f));
    return f;
}

/* v * 2^-sh rounded, held below 2^15, |v| < 2^(15 + sh) */
static int16_t shr_q15(int32_t v, int sh)
{
    if (sh > 0)
        v = ((v >> (sh - 1)) + 1) >> 1;
    else
        v = (int32_t)((uint32_t)v << -sh);
    return (int16_t)(v > 32767 ? 32767 : v);
}

/* Q15 pitch window at i, mirrored */
static int32_t w2q(int i)
{
    return xombe__pitch_w2_q15[i <= XOMBE_PITCH_HALF ? i : XOMBE_PITCH_WIN - 1 - i];
}

static int mag_bits(uint32_t m)
{
    return m ? 32 - q15_clz32(m) : 0;
}

void xombe__kern_pitch_acf(xombe_enc *e, const float *centre, int slot)
{
    int16_t x[PE_X], *rq = e->pacf[slot].q;
    int32_t sv[PE_N];
    uint32_t up[2 * PE_P]; /* up[j] = (u[2j], u[2j + 1]), up[PE_P + j] = (u[2j + 1], u[2j + 2]) */
    float m = 0.0f, gx, cs, cu, b = 0.0f;
    const float *c = centre - XOMBE_PITCH_HALF - XOMBE_PITCH_LPF_HALF;
    int i, t, sh1, sh2, sh = 0;
    uint32_t smax = 0, umax = 0;
    int64_t e0 = 0;
    for (i = 0; i < PE_X - 1; i++)
        if (fabsf(c[i]) > m)
            m = fabsf(c[i]);
    gx = q15_gain(m, 32767.0f);
    for (i = 0; i < PE_X - 1; i++)
        x[i] = (int16_t)round_q(c[i] * gx);
    x[PE_X - 1] = 0;
    /* |sv| < 2^31: sum |taps| <= 65536 */
    for (i = 0; i < XOMBE_PITCH_WIN; i += 2) {
        int32_t sa = 0, sb = 0;
        const int16_t *p = x + i;
        int q;
        for (q = 0; q < XOMBE_PITCH_LPF_PAIRS; q++) {
            uint32_t d = q15_ld2(p + 2 * q);
            sa = q15_smlad(d, xombe__pitch_lpf_q15[0][q], sa);
            sb = q15_smlad(d, xombe__pitch_lpf_q15[1][q], sb);
        }
        if (i + 1 == XOMBE_PITCH_WIN)
            sb = 0;
        sv[i] = sa;
        sv[i + 1] = sb;
        smax |= (uint32_t)(sa < 0 ? -sa : sa) | (uint32_t)(sb < 0 ? -sb : sb);
    }
    /* the OR of the magnitudes has the max's top bit, so s = sv * 2^-sh1 and u = s*w2 * 2^-sh2 peak in
       [2^14, 2^15). s*(s*w2) < 2^45, e0 < 2^54. x holds s, then u */
    sh1 = mag_bits(smax) - 15;
    for (i = 0; i < XOMBE_PITCH_WIN; i++) {
        int32_t s = shr_q15(sv[i], sh1), p = s * w2q(i);
        x[i] = (int16_t)s;
        e0 += (int64_t)s * p;
        umax |= (uint32_t)(p < 0 ? -p : p);
    }
    sh2 = mag_bits(umax) - 15;
    for (i = 0; i < XOMBE_PITCH_WIN; i++)
        x[i] = shr_q15(x[i] * w2q(i), sh2);
    for (; i < 2 * PE_P + 1; i++)
        x[i] = 0;
    for (i = 0; i < PE_P; i++) {
        up[i] = (uint16_t)x[2 * i] | (uint32_t)(uint16_t)x[2 * i + 1] << 16;
        up[PE_P + i] = (uint16_t)x[2 * i + 1] | (uint32_t)(uint16_t)x[2 * i + 2] << 16;
    }
    /* float value of one unit of s and of u */
    cs = gx > 0.0f ? pow2f(sh1) / (gx * xombe__pitch_lpf_q15_gain) : 0.0f;
    cu = cs * pow2f(sh2) / xombe__pitch_w2_q15_gain;
    e->pe0[slot] = i64_to_float(e0) * cs * cs / xombe__pitch_w2_q15_gain;
    e->pacf_scale[slot] = 0.0f;
    /* lags t..t + 3 over the pairs of lag t, zero padding makes the extra products vanish. all four read at one
       pointer, so the loop's accumulators stay in registers */
    for (t = 0; t < XOMBE_PITCH_WIN; t += 4) {
        int64_t a[4] = {0, 0, 0, 0};
        int k;
        const uint32_t *pv = up, *pe = up + (XOMBE_PITCH_WIN - t + 1) / 2, *pb = up + t / 2;
        for (; pv != pe; pv++, pb++) {
            uint32_t v = *pv;
            a[0] = q15_smlald(v, pb[0], a[0]);
            a[2] = q15_smlald(v, pb[1], a[2]);
            a[1] = q15_smlald(v, pb[PE_P], a[1]);
            a[3] = q15_smlald(v, pb[PE_P + 1], a[3]);
        }
        if (t == 0) {
            if (a[0] <= 0) { /* all-zero input */
                memset(rq, 0, PE_N * sizeof(*rq));
                return;
            }
            /* r[0] >= 2^28 (u peaks >= 2^14) shifted below 2^31, |r[t]| <= r[0] (Cauchy-Schwarz),
               b = 32767*2^44/r[0] < 2^31 */
            sh = 33 - q15_clz64((uint64_t)a[0]);
            sh = sh > 0 ? sh : 0;
            e->pacf_scale[slot] = i64_to_float(a[0]) * cu * cu * (1.0f / (32767.0f * 4.0f));
            b = 32767.0f * 17592186044416.0f / (float)(int32_t)(a[0] >> sh);
        }
        for (k = 0; k < 4 && t + k < XOMBE_PITCH_WIN; k++)
            rq[t + k] = (int16_t)(((int64_t)(int32_t)(a[k] >> sh) * (int32_t)b + ((int64_t)1 << 43)) >> 44);
    }
    rq[0] = 32767;
    rq[PE_N - 1] = 0; /* lag XOMBE_PITCH_WIN, past the window */
}

/* reads lags in pairs (r[i], r[i + 1]) at any i */
float xombe__kern_error_at(const xombe_enc *e, int slot, float P)
{
    const int16_t *rq = e->pacf[slot].q;
    float e0 = e->pe0[slot], val;
    /* K = floor((XOMBE_PITCH_WIN - 1)/P) for P = p4/4 */
    int k, p4 = (int)(4.0f * P + 0.5f), K = 4 * (XOMBE_PITCH_WIN - 1) / p4, t4 = p4;
    int32_t sum = 0;
    for (k = 1; k <= K; k++, t4 += p4) {
        uint32_t f = (uint32_t)t4 & 3u;
        sum = q15_smlad(q15_ld2(rq + (t4 >> 2)), (f << 16) | (4u - f), sum);
    }
    /* sum is in quarters, r[0] = 32767 */
    val = (e0 - P * e->pacf_scale[slot] * (float)(2 * sum + 4 * 32767)) / (e0 * (1.0f - P * xombe__pitch_sum_w4));
    return val > 0.0f ? val : 0.0f;
}

/* sin(pi*y), |y| <= 1/2, minimax fit of sin(pi*y)/y, relative error < 2e-7 in float */
static float sinpi(float y)
{
    float y2 = y * y;
    return y * (3.141592741f + y2 * (-5.167709827f + y2 * (2.550069809f + y2 * (-0.598242104f + y2 * 0.0775603801f))));
}

/* sin(pi*y), |y| < 0.032, Taylor to y^5, relative error < 3e-10 */
static float sinpi_small(float y)
{
    float y2 = y * y;
    return y * (3.14159274f + y2 * (-5.16771278f + y2 * 2.55016404f));
}

/* xombe__sin_ratio for 0 <= r < 2m, m < 2^23, im = 1/m */
static float sin_frac(int32_t r, int32_t m, float im)
{
    if (r > m)
        r -= 2 * m;
    if (2 * r > m)
        r = m - r;
    else if (2 * r < -m)
        r = -m - r;
    return sinpi((float)r * im);
}

/* e->win = xombe__hamming(n), n = e->win_len, and the line fit's table of its transform, xombe__hamming_dtft at
   q/g cycles/sample (g = XOMBE_WTAB_OVERSAMPLE*XOMBE_FFT), 0.54*D(x) + 0.23*(D(x - b) + D(x + b)) with
   D(x) = sin(n*x)/sin(x), x = pi*q/g, b = pi/(n - 1). every phase is an exact integer ratio stepped in q:
   x -/+ b = pi*(q*(n - 1) -/+ g)/(g*(n - 1)).
   table in Q15 scaled to peak, mirrored into kq.wt indexed by the integer part s of the signed Q16 position, so it
   steps XOMBE_WTAB_OVERSAMPLE per bin across the centre. s >= 0 holds (wq[s], wq[s + 1]), s < 0 holds
   (wq[-s], wq[-s - 1]), past XOMBE_WTAB_LEN entries it holds the last value as the float fit's table does */
void xombe__kern_window_init(xombe_enc *e)
{
    enum { G = XOMBE_WTAB_OVERSAMPLE * XOMBE_FFT };
    uint32_t *wt = e->kq.wt + XOMBE_LF_WT_OFF;
    int16_t wq[XOMBE_WTAB_LEN];
    int n = e->win_len, i, q, cur, next;
    /* m < 2^23 for n <= XOMBE_MAX_WINDOW. |rm|, |rp| = |q*(n - 1) -/+ g| < 0.032 m (q/g < 0.027, 1/(n - 1) < 0.0051),
       tc, tm, tp = n*q, n*rm, n*rp mod 2g or 2m */
    int32_t m = G * (n - 1), rm = -G, rp = G, tc = 0, tm = 2 * m - n * G, tp = n * G;
    float im = 1.0f / (float)m, ih = 1.0f / (float)(2 * (n - 1)), gain = 0.0f;
    /* cos(2*pi*i/(n - 1)) = sin(pi*(n - 1 - 4i)/(2(n - 1))), symmetric */
    for (i = 0; i <= (n - 1) / 2; i++)
        e->win[i] = e->win[n - 1 - i] = 0.54f - 0.46f * sinpi((float)(n - 1 - 4 * i) * ih);
    for (q = 0; q < XOMBE_WTAB_LEN; q++) {
        float a = sinpi_small((float)q * (1.0f / G)), bm = sinpi_small((float)rm * im);
        float bp = sinpi_small((float)rp * im);
        float A = sin_frac(tc, G, 1.0f / G), Bm = sin_frac(tm, m, im), Bp = sin_frac(tp, m, im), w;
        if (q > 0 && rm != 0) /* one division, a*bm*bp >= 1e-12 */
            w = (0.54f * A * bm * bp + 0.23f * a * (Bm * bp + Bp * bm)) / (a * bm * bp);
        else
            w = 0.54f * (q > 0 ? A / a : (float)n) + 0.23f * ((rm != 0 ? Bm / bm : (float)n) + Bp / bp);
        if (q == 0)
            gain = q15_gain(w, 32767.0f); /* a positive window's transform peaks at 0 */
        wq[q] = (int16_t)round_q(w * gain);
        rm += n - 1;
        rp += n - 1;
        tc += n;
        if (tc >= 2 * G)
            tc -= 2 * G;
        tm += n * (n - 1);
        if (tm >= 2 * m)
            tm -= 2 * m;
        tp += n * (n - 1);
        if (tp >= 2 * m)
            tp -= 2 * m;
    }
    /* wq[q] = entry q, last entry past the end, entry -q is entry q - 1 swapped */
    cur = wq[0];
    for (q = 0; q <= XOMBE_LF_WT_OFF; q++, cur = next) {
        next = wq[q + 1 < XOMBE_WTAB_LEN ? q + 1 : XOMBE_WTAB_LEN - 1];
        wt[q] = (uint16_t)cur | (uint32_t)(uint16_t)next << 16;
        if (q > 0)
            wt[-q] = wt[q - 1] >> 16 | wt[q - 1] << 16;
    }
}
