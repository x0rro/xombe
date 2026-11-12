/* tone detector kernels against float, argv[1] the vectors directory: q15 bin powers against a double DFT, gate
   shares, sincos_turn, and the tone twiddles */
#include "../src/enc/tone_detect.c" /* its spectrum and trig are static */
#include <stdio.h>
#include <stdlib.h>

#if XOMBE_KERNELS_Q15
static int failures, checks;

static void check(int ok, const char *what)
{
    checks++;
    if (!ok) {
        failures++;
        printf("FAIL: %s\n", what);
    }
}

#define PEAK_TOL 2.5e-3 /* largest bin power error, relative to the peak bin (amplitude error near -64 dB) */
#define SHARE_TOL 3e-4  /* gate share error, against DUAL_GATE_SHARE 0.7 */

static struct xombe_enc enc;
static double worst_bin, worst_share;

/* bin powers of x under win, double DFT */
static void ref_power(const float *x, int len, const float *win, double *pw)
{
    int k, n;
    for (k = 0; k <= XOMBE_TONE_FFT / 2; k++) {
        double re = 0.0, im = 0.0;
        for (n = 0; n < len; n++) {
            double v = (double)x[n] * win[n], a = 2.0 * XOMBE_PI * k * n / XOMBE_TONE_FFT;
            re += v * cos(a);
            im -= v * sin(a);
        }
        pw[k] = re * re + im * im;
    }
}

/* q15 powers of one segment against ref_power, scaled by their ratio of sums, and the gate share of each
   AMBE+2 dual against the reference share */
static void one_segment(const int16_t *seg, int len, const float *win, const char *what)
{
    xombe_tone_seg sg;
    double ref[XOMBE_TONE_FFT / 2 + 1], peak = 0.0, num = 0.0, den = 0.0, sc, err = 0.0, rtot = 0.0, serr = 0.0;
    double qtot;
    int k, id;
    char msg[160];
    sg.energy = prepare(&enc, seg, len, win, sg.x, sg.pw);
    if (sg.energy == 0.0f)
        return;
    ref_power(sg.x, len, win, ref);
    for (k = 0; k <= XOMBE_TONE_FFT / 2; k++) {
        peak = ref[k] > peak ? ref[k] : peak;
        num += (double)sg.pw[k];
        den += ref[k];
    }
    sc = num / den;
    for (k = 0; k <= XOMBE_TONE_FFT / 2; k++) {
        double d = fabs((double)sg.pw[k] - sc * ref[k]) / (sc * peak);
        err = d > err ? d : err;
    }
    worst_bin = err > worst_bin ? err : worst_bin;
    snprintf(msg, sizeof(msg), "%s, %d samples: bin power error %.2e of peak", what, len, err);
    check(err <= PEAK_TOL, msg);
    for (k = (int)ceilf(GATE_LO_HZ / BIN_HZ); k <= (int)floorf(GATE_HI_HZ / BIN_HZ); k++)
        rtot += ref[k];
    qtot = (double)band_total(sg.pw);
    for (id = XOMBE_TONE_DUAL_FIRST; id < XOMBE_TONE_DUAL_FIRST + XOMBE_TONE_DUALS; id++) {
        float nominal[2];
        int r0lo, r0hi, r1lo, r1hi;
        double near = 0.0, qnear = 0.0;
        xombe__tone_freqs(id, &nominal[0], &nominal[1]);
        gate_run(nominal[0], &r0lo, &r0hi);
        gate_run(nominal[1], &r1lo, &r1hi);
        for (k = 0; k <= XOMBE_TONE_FFT / 2; k++)
            if ((k >= r0lo && k <= r0hi) || (k >= r1lo && k <= r1hi)) {
                near += ref[k];
                qnear += (double)sg.pw[k];
            }
        if (rtot > 0.0 && qtot > 0.0) {
            double d = fabs(qnear / qtot - near / rtot);
            serr = d > serr ? d : serr;
        }
    }
    worst_share = serr > worst_share ? serr : worst_share;
    snprintf(msg, sizeof(msg), "%s, %d samples: gate share error %.2e", what, len, serr);
    check(serr <= SHARE_TOL, msg);
}

static void spectra(const char *vec)
{
    static int16_t seg[XOMBE_TONE_SEG];
    static const double hz[] = {300.0, 697.0, 1209.0, 1633.0, 2000.0, 3000.0, 3800.0};
    static const double dbfs[] = {-3.0, -20.0, -45.0};
    char path[4096], what[64];
    unsigned lcg = 1u;
    int i, j, n, segs = 0;
    FILE *f;
    for (i = 0; i < (int)(sizeof(hz) / sizeof(hz[0])); i++)
        for (j = 0; j < 3; j++) {
            double a = 32768.0 * pow(10.0, dbfs[j] / 20.0) / 1.5;
            for (n = 0; n < XOMBE_TONE_SEG; n++) {
                lcg = lcg * 1664525u + 1013904223u;
                seg[n] = (int16_t)floor(a * cos(2.0 * XOMBE_PI * hz[i] * n / XOMBE_FS + 0.3) +
                                        0.5 * a * cos(2.0 * XOMBE_PI * 1.37 * hz[i] * n / XOMBE_FS) +
                                        (double)(lcg >> 28) - 7.5 + 0.5);
            }
            snprintf(what, sizeof(what), "%.0f Hz %.0f dBFS", hz[i], dbfs[j]);
            one_segment(seg, XOMBE_TONE_SEG, xombe__tone_win, what);
            one_segment(seg, 160, xombe__tone_win + 40, what); /* shorter segment, any window */
        }
    snprintf(path, sizeof(path), "%s/e3/es1_orig.raw", vec);
    f = fopen(path, "rb");
    if (!f) {
        perror(path);
        exit(2);
    }
    while (segs < 120 && fread(seg, sizeof(int16_t), XOMBE_TONE_SEG, f) == XOMBE_TONE_SEG) {
        snprintf(what, sizeof(what), "es1 segment %d", segs++);
        one_segment(seg, XOMBE_TONE_SEG, xombe__tone_win, what);
    }
    fclose(f);
    printf("kern_tone: worst bin power error %.2e of peak, gate share error %.2e\n", worst_bin, worst_share);
}

static void trig(void)
{
    double worst = 0.0;
    int i;
    for (i = -400000; i <= 400000; i++) {
        float t = (float)i * 3.0e-4f, s, c; /* to 120 turns, as far as (2*len - 1)*w/2 reaches */
        double ds = sin(2.0 * XOMBE_PI * (double)t), dc = cos(2.0 * XOMBE_PI * (double)t), e;
        sincos_turn(t, &s, &c);
        e = fmax(fabs(s - ds), fabs(c - dc));
        worst = e > worst ? e : worst;
    }
    printf("kern_tone: sincos_turn worst error %.2e\n", worst);
    check(worst < 3e-7, "sincos_turn within 3e-7");
}

static void twiddles(void)
{
    int j, ok = 1;
    for (j = 0; j < XOMBE_TONE_FFT / 4; j++)
        ok &= xombe__tone_tw_q15[j] == xombe__fft_tw_q15[2 * j];
    check(ok, "xombe__tone_tw_q15 is every second xombe__fft_tw_q15");
}
#endif

int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr, "usage: test_kern_tone VEC\n");
        return 2;
    }
#if XOMBE_KERNELS_Q15
    twiddles();
    trig();
    spectra(argv[1]);
    printf("kern_tone: %d/%d checks passed\n", checks - failures, checks);
    return failures != 0;
#else
    (void)argv;
    printf("kern_tone: float build, no checks\n");
    return 0;
#endif
}
