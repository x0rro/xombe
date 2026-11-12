#include <math.h>
#include <string.h>
#include "dkernels.h"
#include "../tables/dsp_tables.h"

#define GRAIN_RAMP XOMBE_GRAIN_RAMP
#define WT_SIZE XOMBE_WT_SIZE
#define UV_FFT XOMBE_UV_FFT
#define UV_BINS XOMBE_UV_BINS

XOMBE_STATIC_ASSERT(WT_SIZE <= XOMBE_IRFFT_MAX, wavetable_fits_irfft);
XOMBE_STATIC_ASSERT(XOMBE_TW_N % WT_SIZE == 0, twiddles_cover_wavetable);
XOMBE_STATIC_ASSERT(UV_FFT <= XOMBE_IRFFT_MAX, uv_fft_fits_irfft);
XOMBE_STATIC_ASSERT(XOMBE_TW_N % UV_FFT == 0, twiddles_cover_uv_fft);
XOMBE_STATIC_ASSERT(XOMBE_LMAX <= XOMBE_HALF_MAX, sqrt_table_covers_L);

void xombe__dk_wavetable(const xombe_set *s, const float *M, const uint8_t *vh, float tab[XOMBE_WT_SIZE])
{
    xombe_cpx spec[WT_SIZE / 2 + 1];
    int l;
    memset(spec, 0, sizeof(spec));
    for (l = 1; l <= s->L && l <= XOMBE_LMAX; l++) {
        float psi;
        if (!vh[l] || M[l] == 0.0f)
            continue;
        psi = s->phi[l] - 2.0f * XOMBE_PI_F * l * XOMBE_EPOCH_REF * s->f0;
        spec[l].re = WT_SIZE * M[l] * cosf(psi);
        spec[l].im = WT_SIZE * M[l] * sinf(psi);
    }
    xombe__irfft(spec, tab, WT_SIZE, xombe__fft_tw, XOMBE_TW_N);
}

/* 0.5*(1 - cos(pi*(j + frac)/GRAIN_RAMP)) for integer j, cf and sf the cos and sin of pi*frac/GRAIN_RAMP */
static float ramp(int j, float cf, float sf)
{
    return 0.5f * (1.0f - (xombe__ramp_cos[j] * cf - xombe__ramp_sin[j] * sf));
}

void xombe__dk_add_grain(float *buf, int N, float P, float fhat, float period, const float *tab)
{
    float E = P + period;
    int nb = (int)ceilf(P), ne = (int)floorf(E), n0 = nb, n1 = ne + GRAIN_RAMP, n;
    /* fade-in over n - P = (n - nb) + (nb - P), fade-out over E + GRAIN_RAMP - n = (GRAIN_RAMP + ne - n) + (E - ne):
       integer table steps rotated by a fraction fixed for the grain */
    float ai = XOMBE_PI_F * ((float)nb - P) / GRAIN_RAMP, ao = XOMBE_PI_F * (E - (float)ne) / GRAIN_RAMP;
    float cs = cosf(ai), ss = sinf(ai), ce = cosf(ao), se = sinf(ao);
    if (n0 < 0)
        n0 = 0;
    if (n1 > N + XOMBE_TAIL - 1)
        n1 = N + XOMBE_TAIL - 1;
    for (n = n0; n <= n1; n++) {
        float th = (n - P) * fhat, x, fr, w;
        int i;
        th -= floorf(th);
        x = th * WT_SIZE;
        i = (int)x;
        fr = x - i;
        if (n < P + GRAIN_RAMP)
            w = ramp(n - nb, cs, ss);
        else if (n <= E)
            w = 1.0f;
        else
            w = ramp(GRAIN_RAMP + ne - n, ce, se);
        buf[n] += w * (tab[i & (WT_SIZE - 1)] * (1.0f - fr) + tab[(i + 1) & (WT_SIZE - 1)] * fr);
    }
}

/* sqrt(n/N), crossfade gain n samples into a half-frame of N, rs = 1/sqrt(N) */
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
    float buf[UV_FFT];
    xombe_cpx X[UV_BINS], Y[UV_BINS];
    float beta = UV_FFT * cur->f0, rs = xombe__rsqrt_half[N - XOMBE_HALF_MIN], sqrt_L = xombe__sqrt_n[cur->L];
    int n, k, l, upper = 0, low;
    int silence = cur->cls == XOMBE_CLASS_SILENCE;
    memset(buf, 0, sizeof(buf));
    for (n = 0; n < 2 * N; n++) {
        float w = n < N ? fade(n, rs) : fade(2 * N - n, rs);
        buf[n] = w * (float)seg[n];
    }
    xombe__rfft(buf, X, UV_FFT, xombe__fft_tw, XOMBE_TW_N);
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
    xombe__irfft(Y, y, UV_FFT, xombe__fft_tw, XOMBE_TW_N);
}
