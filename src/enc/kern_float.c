#include <string.h>
#include "kernels.h"
#include "../tables/dsp_tables.h"

void xombe__kern_spectrum(const float *centre, const float *w, int n, xombe_spec *S)
{
    float buf[XOMBE_FFT];
    int half = (n - 1) / 2, j;
    memset(buf, 0, sizeof(buf));
    for (j = -half; j <= half; j++)
        buf[(j + XOMBE_FFT) % XOMBE_FFT] = centre[j] * w[half + j];
    xombe__rfft(buf, S->b, XOMBE_FFT, xombe__fft_tw, XOMBE_TW_N);
}

void xombe__kern_bin_power(const xombe_spec *S, float sum_w2, float *P)
{
    int k;
    float scale = 1.0f / (XOMBE_FFT * sum_w2);
    for (k = 0; k < XOMBE_SPEC_BINS; k++)
        P[k] = (S->b[k].re * S->b[k].re + S->b[k].im * S->b[k].im) * scale;
}

void xombe__kern_window_init(xombe_enc *e)
{
    (void)e;
}

void xombe__kern_pitch_acf(xombe_enc *e, const float *centre, int slot)
{
    float u[XOMBE_PITCH_WIN], e0 = 0.0f;
    int i, j, t;
    for (i = 0; i < XOMBE_PITCH_WIN; i++) {
        float s = 0.0f;
        for (j = -XOMBE_PITCH_LPF_HALF; j <= XOMBE_PITCH_LPF_HALF; j++)
            s += xombe__pitch_lpf[j + XOMBE_PITCH_LPF_HALF] * centre[i - XOMBE_PITCH_HALF - j];
        e0 += s * s * xombe__pitch_w2[i];
        u[i] = s * xombe__pitch_w2[i];
    }
    for (t = 0; t < XOMBE_PITCH_WIN; t++) {
        float r = 0.0f;
        for (i = 0; i + t < XOMBE_PITCH_WIN; i++)
            r += u[i] * u[i + t];
        e->pacf[slot].f[t] = r;
    }
    e->pe0[slot] = e0;
}

/* P on a quarter-sample grid, so K, k*P, and its integer part are exact */
float xombe__kern_error_at(const xombe_enc *e, int slot, float P)
{
    const float *r = e->pacf[slot].f;
    float e0 = e->pe0[slot], sum = r[0], val;
    int k, K = (int)((float)(XOMBE_PITCH_WIN - 1) / P);
    for (k = 1; k <= K; k++) {
        float t = (float)k * P, f;
        int i = (int)t;
        f = t - (float)i;
        sum += 2.0f * (f > 0.0f ? (1.0f - f) * r[i] + f * r[i + 1] : r[i]);
    }
    val = (e0 - P * sum) / (e0 * (1.0f - P * xombe__pitch_sum_w4));
    return val > 0.0f ? val : 0.0f;
}

void xombe__kern_lf_prepare(xombe_enc *e, const xombe_spec *S, int kmin, int kmax, xombe_lf_spec *sp)
{
    (void)e;
    (void)kmin;
    (void)kmax;
    sp->S = S;
    sp->g = 1.0f;
}

/* window values kept per harmonic for the residual pass, a wider harmonic recomputes the rest */
#define LINE_FIT_KEEP 64

/* window transform nu_bins from its centre, interpolated in e->wtab */
static float window_value(const xombe_enc *e, float nu_bins)
{
    float x = fabsf(nu_bins) * XOMBE_WTAB_OVERSAMPLE;
    int i = (int)x;
    float f = x - (float)i;
    if (i >= XOMBE_WTAB_LEN - 1)
        return e->wtab[XOMBE_WTAB_LEN - 1];
    return (1.0f - f) * e->wtab[i] + f * e->wtab[i + 1];
}

float xombe__kern_line_fit(const xombe_enc *e, const xombe_lf_spec *sp, float period, int kmin, int kmax,
                           float band_resid[8], float band_energy[8])
{
    const xombe_cpx *S = sp->S->b;
    float per_harm = XOMBE_FFT / period, total = 0.0f, wk[LINE_FIT_KEEP];
    int k = kmin < 1 ? 1 : kmin, l = xombe__harmonic_of(k, period), next = l;
    memset(band_resid, 0, 8 * sizeof(float));
    memset(band_energy, 0, 8 * sizeof(float));
    while (k <= kmax) {
        int k0 = k, n = 0, band;
        float sw_re = 0.0f, sw_im = 0.0f, ww = 0.0f, ee = 0.0f, resid = 0.0f, A_re = 0.0f, A_im = 0.0f;
        do {
            float w = window_value(e, (float)k - (float)l * per_harm);
            sw_re += S[k].re * w;
            sw_im += S[k].im * w;
            ww += w * w;
            ee += S[k].re * S[k].re + S[k].im * S[k].im;
            if (n < LINE_FIT_KEEP)
                wk[n] = w;
            n++;
            k++;
        } while (k <= kmax && (next = xombe__harmonic_of(k, period)) == l);
        if (l == 0) {
            l = next;
            continue;
        }
        /* residual summed term by term, ee - |sw|^2/ww cancels on a strongly harmonic band (relative error 1e-3
           against 1e-4 in float) */
        if (ww > 0.0f) {
            A_re = sw_re / ww;
            A_im = sw_im / ww;
        }
        for (n = 0; k0 < k; k0++, n++) {
            float w = n < LINE_FIT_KEEP ? wk[n] : window_value(e, (float)k0 - (float)l * per_harm),
                  dr = S[k0].re - A_re * w, di = S[k0].im - A_im * w;
            resid += dr * dr + di * di;
        }
        band = xombe__band_of(l, period);
        band_resid[band] += resid;
        band_energy[band] += ee;
        total += resid;
        l = next;
    }
    return total;
}
