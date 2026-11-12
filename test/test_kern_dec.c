/* decoder kernels (src/dec/dkernels.h) and src/dsp/fmath.h against float, argv[1] the vectors directory */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "../src/internal.h"
#include "../src/dec/dkernels.h"
#include "../src/dsp/fmath.h"

#if XOMBE_KERNELS_Q15
/* the float kernels as the reference, renamed */
#define xombe__dk_wavetable ref_wavetable
#define xombe__dk_add_grain ref_add_grain
#define xombe__dk_unvoiced ref_unvoiced
void ref_wavetable(const xombe_set *s, const float *M, const uint8_t *vh, float tab[XOMBE_WT_SIZE]);
void ref_add_grain(float *buf, int N, float P, float fhat, float period, const float *tab);
void ref_unvoiced(const int16_t *seg, int N, const xombe_set *cur, const float *M, const uint8_t *vh,
                  float avg_f0, float y[XOMBE_UV_FFT]);
#include "../src/dec/dkern_float.c"
#undef xombe__dk_wavetable
#undef xombe__dk_add_grain
#undef xombe__dk_unvoiced
#endif

/* fmath.h against libm: exact in float builds, in q15 builds exp2 within EXP2_TOL relative, log2 and sincos
   within LOG2_TOL (absolute, relative above 1) and SINCOS_TOL absolute, sincos over |x| <= 520 (wavetable
   phases reach about 360) */
#define EXP2_TOL 2e-7
#define LOG2_TOL 1.5e-7
#define SINCOS_TOL 2e-7

static double rel_err(double v, double ref)
{
    return fabs(v - ref) / fabs(ref);
}

static double log_err(double v, double ref)
{
    return fabs(v - ref) / fmax(1.0, fabs(ref));
}

static int test_fmath(void)
{
    double e_exp = 0.0, e_log = 0.0, e_sc = 0.0, d;
    int i, fails = 0;
    for (i = -1300000; i <= 1300000; i++) {
        float x = (float)i * 9e-5f, y = exp2f((float)i * 9.6e-5f), z = (float)(i + 1300001) * 1e-6f, s, c;
        xombe__sincosf(4.4f * x, &s, &c);
        d = fmax(fabs(s - sin(4.4f * x)), fabs(c - cos(4.4f * x)));
        e_sc = fmax(e_sc, d);
        e_exp = fmax(e_exp, rel_err(xombe__exp2f(x), exp2(x)));
        e_log = fmax(e_log, fmax(log_err(xombe__log2f(y), log2(y)), log_err(xombe__log2f(z), log2(z))));
        if (!XOMBE_KERNELS_Q15) {
            xombe__sincosf(x, &s, &c);
            if (xombe__exp2f(x) != exp2f(x) || s != sinf(x) || c != cosf(x) || xombe__log2f(y) != log2f(y))
                fails++;
        }
    }
    /* libm's results at the edges */
    if (xombe__exp2f(-200.0f) != exp2f(-200.0f) || !isinf(xombe__exp2f(200.0f)) || isinf(xombe__exp2f(127.9f)) ||
        !isinf(xombe__log2f(0.0f)) || xombe__log2f(1e-40f) != log2f(1e-40f))
        fails++;
    if (XOMBE_KERNELS_Q15 && (e_exp > EXP2_TOL || e_log > LOG2_TOL || e_sc > SINCOS_TOL))
        fails++;
    printf("kern_dec: fmath max error exp2 %.2e (relative), log2 %.2e, sincos %.2e: %s\n", e_exp, e_log, e_sc,
           fails ? "FAIL" : "ok");
    return fails;
}

#if XOMBE_KERNELS_Q15
/* q15 kernels against float on random sets, SNR limits in dB over each output: the wavetable differs only by
   fmath.h, the unvoiced transforms carry the Q15 FFT's 60-70 dB, down to ~52 dB when a few strong unvoiced
   bands hold nearly all the power, the block exponent then set by bins the error spreads over all */
#define WT_SNR 100.0
#define GRAIN_SNR 95.0
#define UV_SNR 50.0
#define CN_SNR 55.0
#define SETS 400

static uint32_t rng = 1;

static float urand(void)
{
    rng = rng * 1664525u + 1013904223u;
    return (float)(rng >> 8) * (1.0f / 16777216.0f);
}

static void snr_add(const float *a, const float *b, int n, double *worst)
{
    double e = 0.0, p = 0.0, r;
    int i;
    for (i = 0; i < n; i++) {
        e += (double)(a[i] - b[i]) * (a[i] - b[i]);
        p += (double)b[i] * b[i];
    }
    r = e > 0.0 ? 10.0 * log10(p / e) : 200.0;
    if (p > 0.0 && r < *worst)
        *worst = r;
}

/* f0 80..400 Hz, log2 amplitudes falling about 0.1 per harmonic plus up to 6 of noise, random phases and
   voicing, every class but tone */
static void random_set(xombe_set *s, float *M, uint8_t *vh)
{
    int l;
    memset(s, 0, sizeof(*s));
    s->f0 = XOMBE_HZ(80.0f + 320.0f * urand());
    s->L = xombe__L_from_f0(s->f0);
    s->cls = urand() < 0.2f ? XOMBE_CLASS_SILENCE : XOMBE_CLASS_VOICE;
    for (l = 0; l <= XOMBE_LMAX; l++) {
        M[l] = l ? exp2f(12.0f - 0.1f * (float)l - 6.0f * urand()) : 0.0f;
        vh[l] = (uint8_t)(urand() < 0.5f);
        s->phi[l] = 400.0f * (urand() - 0.5f);
    }
}

static int test_kernels(void)
{
    double wt = 200.0, uv = 200.0, cn = 200.0, gr = 200.0;
    int i, n, fails;
    for (i = 0; i < SETS; i++) {
        xombe_set s;
        float M[XOMBE_LMAX + 1], a[XOMBE_UV_FFT], b[XOMBE_UV_FFT];
        int16_t seg[2 * XOMBE_HALF_MAX];
        uint8_t vh[XOMBE_LMAX + 1];
        int N = XOMBE_HALF_MIN + (int)(urand() * (XOMBE_HALF_MAX - XOMBE_HALF_MIN + 1));
        random_set(&s, M, vh);
        xombe__dk_wavetable(&s, M, vh, a);
        ref_wavetable(&s, M, vh, b);
        snr_add(a, b, XOMBE_WT_SIZE, &wt);
        {
            float ga[XOMBE_HALF_MAX + XOMBE_TAIL], gb[XOMBE_HALF_MAX + XOMBE_TAIL], fh = s.f0 * (0.9f + 0.2f * urand());
            float P = 1.5f * (float)N * urand() - 0.5f * (float)N;
            memset(ga, 0, sizeof(ga));
            memset(gb, 0, sizeof(gb));
            xombe__dk_add_grain(ga, N, P, fh, 1.0f / s.f0, b);
            ref_add_grain(gb, N, P, fh, 1.0f / s.f0, b);
            snr_add(ga, gb, N + XOMBE_TAIL, &gr);
        }
        for (n = 0; n < 2 * N; n++)
            seg[n] = (int16_t)(urand() * 65535.0f - 32768.0f);
        xombe__dk_unvoiced(seg, N, &s, M, vh, XOMBE_HZ(100.0f), a);
        ref_unvoiced(seg, N, &s, M, vh, XOMBE_HZ(100.0f), b);
        snr_add(a, b, 2 * N, s.cls == XOMBE_CLASS_SILENCE ? &cn : &uv);
    }
    fails = (wt < WT_SNR) + (gr < GRAIN_SNR) + (uv < UV_SNR) + (cn < CN_SNR);
    printf("kern_dec: worst SNR vs float over %d sets, wavetable %.1f dB, grain %.1f dB, unvoiced %.1f dB, comfort "
           "noise %.1f dB: %s\n", SETS, wt, gr, uv, cn, fails ? "FAIL" : "ok");
    return fails;
}
#endif

int main(int argc, char **argv)
{
    int fails = 0;
    (void)argc;
    (void)argv;
    fails += test_fmath();
#if XOMBE_KERNELS_Q15
    fails += test_kernels();
#endif
    return fails != 0;
}
