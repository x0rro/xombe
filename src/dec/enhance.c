#include <math.h>
#include "../internal.h"
#include "../dsp/fmath.h"

#define PHASE_TAIL_SLOPE 0.7212f
#define PHASE_MIRROR 64

static float phase_seq(const float *a, int L, int n)
{
    if (n < 0)
        n = -n;
    if (n > PHASE_MIRROR)
        n = 2 * PHASE_MIRROR - n;
    if (n == 0)
        return 0.0f;
    if (n <= L)
        return a[n];
    return a[L] - PHASE_TAIL_SLOPE * (float)(n - L);
}

#define PHASE_TAPS 19 /* odd m = 1..PHASE_TAPS */

/* d/m, a reciprocal multiply in q15 builds */
static float tap(float d, int m)
{
#if XOMBE_KERNELS_Q15
    static const float inv_m[PHASE_TAPS + 1] = {0.0f, 1.0f, 0.0f, 1.0f / 3, 0.0f, 1.0f / 5, 0.0f, 1.0f / 7, 0.0f,
                                                1.0f / 9, 0.0f, 1.0f / 11, 0.0f, 1.0f / 13, 0.0f, 1.0f / 15, 0.0f,
                                                1.0f / 17, 0.0f, 1.0f / 19};
    return d * inv_m[m];
#else
    return d / (float)m;
#endif
}

void xombe__phases(const float *a, int L, float *phi)
{
    const float k = (float)(2.0 * 0.69314718055994531 / XOMBE_PI); /* 2 ln 2 / pi */
    float t[XOMBE_LMAX + PHASE_TAPS + 1]; /* phase_seq at 0..L + PHASE_TAPS */
    int l, m;
    for (l = 0; l <= L + PHASE_TAPS; l++)
        t[l] = phase_seq(a, L, l);
    phi[0] = 0.0f;
    for (l = 1; l <= XOMBE_LMAX; l++) {
        float s = 0.0f;
        if (l <= L)
            for (m = 1; m <= PHASE_TAPS; m += 2)
                s += tap(t[l + m] - t[l >= m ? l - m : m - l], m);
        phi[l] = k * s;
    }
}

#define ENH_GAIN (0.96f * XOMBE_PI_F)
#define ENH_LOG2_W_MIN (-1.0f)        /* log2 0.5 */
#define ENH_LOG2_W_MAX 0.263034406f   /* log2 1.2 */

/* c[l], s[l] = cos, sin of w0*l for l = 1..L, by rotation in q15 builds (error about l ulp) */
static void harmonic_cs(float w0, int L, float *c, float *s)
{
    int l;
#if XOMBE_KERNELS_Q15
    float c1, s1;
    xombe__sincosf(w0, &s1, &c1);
    c[1] = c1;
    s[1] = s1;
    for (l = 2; l <= L; l++) {
        c[l] = c[l - 1] * c1 - s[l - 1] * s1;
        s[l] = s[l - 1] * c1 + c[l - 1] * s1;
    }
#else
    for (l = 1; l <= L; l++)
        xombe__sincosf(w0 * (float)l, &s[l], &c[l]);
#endif
}

/* spectral enhancement in the log domain: e = a + log2 W + log2 gamma
   W = sqrt(M)*(num/den)^(1/4) above harmonic L/8 (else 1), clamped to 0.5..1.2,
   num/den = 0.96 pi (r0^2 + r1^2 - 2 r0 r1 cos(w0 l))/(w0 r0 (r0^2 - r1^2)), r0, r1 the power and first
   autocorrelation of M = 2^a, gamma restores the power r0
   powers are relative to the strongest harmonic so every sum is scale-free: with p = 2^(2(a - amax)) and
   rho = r1/r0, num/den = 0.96 pi ((1 - rho cos)^2 + (rho sin)^2)/(w0 r0 (1 - rho)(1 + rho)) * 2^(-2 amax) */
void xombe__enhance(const float *a, int L, float f0, float *e)
{
    float w0 = 2.0f * XOMBE_PI_F * f0, amax = a[1], r0 = 0.0f, r1 = 0.0f, sum = 0.0f, rho, lg;
    float lw[XOMBE_LMAX + 1], c[XOMBE_LMAX + 1], s[XOMBE_LMAX + 1];
    int l;
    for (l = 2; l <= L; l++)
        if (a[l] > amax)
            amax = a[l];
    harmonic_cs(w0, L, c, s);
    for (l = 1; l <= L; l++) {
        float p = xombe__exp2f(2.0f * (a[l] - amax));
        r0 += p;
        r1 += p * c[l];
    }
    rho = r1 / r0; /* r0 >= 1, the strongest harmonic contributes 1 */
    for (l = 1; l <= L; l++) {
        float w = 0.0f;
        if (8 * l > L) {
            float u = 1.0f - rho * c[l], v = rho * s[l];
            float den = w0 * r0 * ((1.0f - rho) * (1.0f + rho));
            if (den > 0.0f) {
                w = 0.5f * (a[l] - amax) + 0.25f * xombe__log2f(ENH_GAIN * (u * u + v * v) / den);
                if (w < ENH_LOG2_W_MIN)
                    w = ENH_LOG2_W_MIN;
                if (w > ENH_LOG2_W_MAX)
                    w = ENH_LOG2_W_MAX;
            }
        }
        lw[l] = w;
        sum += xombe__exp2f(2.0f * (w + a[l] - amax));
    }
    lg = 0.5f * (xombe__log2f(r0) - xombe__log2f(sum)); /* sum >= 1/4 from the strongest harmonic, W >= 1/2 */
    for (l = 1; l <= L; l++)
        e[l] = a[l] + lw[l] + lg;
}

#define SHAPE_LOW_MAX XOMBE_HZ(200.0f) /* low roll-off corner: the average pitch, at most this */

/* xombe__hf_tilt through fmath.h, must match it */
static float hf_tilt(float f)
{
    return f > XOMBE_HF_CORNER ? xombe__log2f(f / XOMBE_HF_CORNER) : 0.0f;
}

void xombe__shape(const float *e, float *M, int L, float f0, float avg_f0)
{
    float fc = avg_f0 < SHAPE_LOW_MAX ? avg_f0 : SHAPE_LOW_MAX;
    int l;
    for (l = 1; l <= L; l++) {
        float f = (float)l * f0, g = e[l];
        if (f < fc)
            g += 2.0f * xombe__log2f(f / fc);
        g += hf_tilt(f);
        M[l] = xombe__exp2f(g);
    }
}
