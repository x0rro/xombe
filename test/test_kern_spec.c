/* analysis spectrum kernels (src/enc/kernels.h) against a double-precision DFT, argv[1] the vectors directory */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../src/internal.h"
#include "../src/enc/kernels.h"
#include "../src/tables/dsp_tables.h"

#define PI 3.14159265358979323846
#define HMAX (XOMBE_MAX_WINDOW / 2)

/* worst bin amplitude error relative to the peak bin, and least SNR over all bins. q15: the block-floating-point
   FFT sets both (61 dB, 0.0032 on these frames with a float split), prepared bins against their range's peak */
#if XOMBE_KERNELS_Q15
#define AMP_TOL 5e-3
#define SNR_MIN 58.0
#define LF_TOL 0.05
#else
#define AMP_TOL 1e-5
#define SNR_MIN 100.0
#endif

static int fails;
static double worst_snr = 1e9, worst_amp, sum_snr;
#if XOMBE_KERNELS_Q15
static double worst_lf;
#endif
static int nsnr;

static double urand(uint32_t *s)
{
    *s = *s * 1664525u + 1013904223u;
    return (double)(*s >> 8) / 8388608.0 - 1.0;
}

static void fail(const char *what, const char *msg, double v)
{
    if (fails++ < 10)
        printf("%s: %s %g\n", what, msg, v);
}

#if XOMBE_KERNELS_Q15
/* block layout, then lf_prepare on the whole band and on a high band: prepared peak at least 16384, exact prefix
   sums, and bins / g against the reference relative to the range's peak */
static void check_q15(const char *what, const xombe_spec *S, const double *rr, const double *ri)
{
    static xombe_enc e;
    const int16_t *re = S->v, *im = re + XOMBE_LF_IM;
    int k, j, big = 0, kmax = XOMBE_SPEC_BINS - 2;
    if (im[0] != 0 || im[XOMBE_SPEC_BINS - 1] != 0)
        fail(what, "imaginary DC or Nyquist", 0);
    for (k = 0; k < XOMBE_SPEC_BINS; k++) {
        if (abs(re[k]) > big)
            big = abs(re[k]);
        if (abs(im[k]) > big)
            big = abs(im[k]);
    }
    if (big < 32735 || big > 32767)
        fail(what, "largest component", big);
    for (j = 0; j < 2; j++) {
        xombe_lf_spec sp;
        uint64_t acc = 0;
        double lpeak = 0.0, lf = 0.0;
        const int16_t *pr, *pi;
        int kmin = j ? 300 : 1, top = 0;
        xombe__kern_lf_prepare(&e, S, kmin, kmax, &sp);
        if (sp.g == 0.0f || (sp.S != S && sp.S != &e.kq.spec)) {
            fail(what, "prepared g 0 or a foreign spectrum, from kmin", kmin);
            continue;
        }
        pr = sp.S->v;
        pi = pr + XOMBE_LF_IM;
        if (e.kq.ep[kmin] != 0)
            fail(what, "prefix sum start", (double)e.kq.ep[kmin]);
        for (k = kmin; k <= kmax; k++) {
            acc += (uint64_t)((int32_t)pr[k] * pr[k]) + (uint64_t)((int32_t)pi[k] * pi[k]);
            if (e.kq.ep[k + 1] != acc)
                fail(what, "prefix sum at", k);
            if (abs(pr[k]) > top)
                top = abs(pr[k]);
            if (abs(pi[k]) > top)
                top = abs(pi[k]);
            if (hypot(rr[k], ri[k]) > lpeak)
                lpeak = hypot(rr[k], ri[k]);
        }
        if (top < 16384)
            fail(what, "prepared peak", top);
        for (k = kmin; k <= kmax; k++) {
            double d = hypot(pr[k] / sp.g - rr[k], pi[k] / sp.g - ri[k]) / lpeak;
            if (d > lf)
                lf = d;
        }
        if (lf > worst_lf)
            worst_lf = lf;
        if (lf > LF_TOL)
            fail(what, "prepared bins off, of the range peak, by", lf);
    }
}
#endif

static void check_frame(const char *what, const float *centre, const float *w, int n, double sum_w2)
{
    static double rr[XOMBE_SPEC_BINS], ri[XOMBE_SPEC_BINS];
    static float P[XOMBE_SPEC_BINS];
    static xombe_spec S;
    double peak = 0.0, sig = 0.0, err = 0.0, amp = 0.0, scale = 1.0 / (XOMBE_FFT * sum_w2), snr;
    int half = (n - 1) / 2, j, k;
    for (k = 0; k < XOMBE_SPEC_BINS; k++) {
        double xr = 0.0, xi = 0.0;
        for (j = -half; j <= half; j++) {
            double a = -2.0 * PI * (double)(((j + XOMBE_FFT) * k) % XOMBE_FFT) / XOMBE_FFT;
            double v = (double)centre[j] * w[half + j];
            xr += v * cos(a);
            xi += v * sin(a);
        }
        rr[k] = xr;
        ri[k] = xi;
        if (hypot(xr, xi) > peak)
            peak = hypot(xr, xi);
    }
    xombe__kern_spectrum(centre, w, n, &S);
    xombe__kern_bin_power(&S, (float)sum_w2, P);
    if (peak == 0.0) {
        for (k = 0; k < XOMBE_SPEC_BINS; k++)
            if (P[k] != 0.0f)
                fail(what, "power of a zero frame at bin", k);
        return;
    }
    for (k = 0; k < XOMBE_SPEC_BINS; k++) {
        double d = sqrt(P[k] / scale) - hypot(rr[k], ri[k]);
        sig += rr[k] * rr[k] + ri[k] * ri[k];
        err += d * d;
        if (fabs(d) / peak > amp)
            amp = fabs(d) / peak;
    }
    snr = err > 0.0 ? 10.0 * log10(sig / err) : 999.0;
    if (amp > worst_amp)
        worst_amp = amp;
    sum_snr += snr;
    nsnr++;
    if (snr < worst_snr)
        worst_snr = snr;
    if (amp > AMP_TOL)
        fail(what, "bin amplitude error, of the peak", amp);
    if (snr < SNR_MIN)
        fail(what, "SNR dB", snr);
#if XOMBE_KERNELS_Q15
    check_q15(what, &S, rr, ri);
#endif
}

int main(int argc, char **argv)
{
    static float buf[2 * HMAX + 1], wl[XOMBE_MAX_WINDOW];
    static int16_t pcm[100000];
    float *centre = buf + HMAX;
    const float *ws[2];
    int ns[2] = {XOMBE_BASE_WINDOW, XOMBE_MAX_WINDOW}, i, j, k, frames = 0, len = 0;
    double s2[2] = {0.0, 0.0};
    uint32_t seed = 1;
    char path[1024], what[64];
    FILE *f;
    for (j = 0; j < XOMBE_MAX_WINDOW; j++)
        wl[j] = (float)(0.54 - 0.46 * cos(2.0 * PI * j / (XOMBE_MAX_WINDOW - 1)));
    ws[0] = xombe__base_win;
    ws[1] = wl;
    for (i = 0; i < 2; i++)
        for (j = 0; j < ns[i]; j++)
            s2[i] += (double)ws[i][j] * ws[i][j];
    snprintf(path, sizeof path, "%s/e3/es1_orig.raw", argc > 1 ? argv[1] : "vectors");
    f = fopen(path, "rb");
    if (f) {
        len = (int)fread(pcm, sizeof(pcm[0]), sizeof(pcm) / sizeof(pcm[0]), f);
        fclose(f);
    } else {
        fail(path, "can't read", 0);
    }
    /* speech every 37 ms through es1 (16-bit little-endian), both windows */
    for (k = HMAX; k + HMAX < len; k += 296) {
        for (i = 0; i < 2; i++) {
            int h = (ns[i] - 1) / 2;
            for (j = -h; j <= h; j++) {
                const uint8_t *b = (const uint8_t *)&pcm[k + j];
                centre[j] = (float)(int16_t)(b[0] | b[1] << 8);
            }
            snprintf(what, sizeof what, "es1 at %d, window %d", k, ns[i]);
            check_frame(what, centre, ws[i], ns[i], s2[i]);
            frames++;
        }
    }
    /* harmonics with a falling envelope and noise 70 dB down, loud and quiet, an impulse, and silence */
    for (i = 0; i < 3; i++) {
        double f0 = 0.0125 + 0.02 * i;
        int l;
        for (j = -HMAX; j <= HMAX; j++) {
            double v = 3e-4 * urand(&seed);
            for (l = 1; l * f0 < 0.5; l++)
                v += cos(2.0 * PI * l * f0 * j + l) / l;
            centre[j] = (float)(v * (i == 2 ? 1e-6 : 3000.0));
        }
        snprintf(what, sizeof what, "harmonics f0 %.4f", f0);
        check_frame(what, centre, wl, XOMBE_MAX_WINDOW, s2[1]);
        frames++;
    }
    memset(buf, 0, sizeof buf);
    centre[0] = 1000.0f;
    check_frame("impulse", centre, wl, XOMBE_MAX_WINDOW, s2[1]);
    centre[0] = 0.0f;
    check_frame("silence", centre, wl, XOMBE_MAX_WINDOW, s2[1]);
    frames += 2;
    printf("kern_spec: %d frames, SNR mean %.2f dB least %.1f dB, worst bin amplitude error %.2g of peak", frames,
           sum_snr / nsnr, worst_snr, worst_amp);
#if XOMBE_KERNELS_Q15
    printf(", prepared bins %.2g of the range peak", worst_lf);
#endif
    printf(", %d failures\n", fails);
    return fails != 0;
}
