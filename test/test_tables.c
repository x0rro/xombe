#include <math.h>
#include <stdio.h>
#include "../src/internal.h"
#include "../src/tables/tables.h"
#include "../src/tables/dsp_tables.h"
#include "ref.h"

/* generated tables (scripts/gen_dsp_tables.py) against their closed forms, spectral blocks against the kernel
   table's range, and the harmonic interpolation's end rules */

static int failures, checks;

static void check(int ok, const char *what, int i)
{
    checks++;
    if (!ok) {
        failures++;
        fprintf(stderr, "FAIL %s (item %d)\n", what, i);
    }
}

static int rounds_to(float t, double v)
{
    return fabs((double)t - v) <= 6e-8 * (fabs(v) > 1.0 ? fabs(v) : 1.0);
}

int main(void)
{
    int J, m, i, L, k, N, n;
    float a[XOMBE_LMAX + 1];

    for (m = 0; m < 8; m++)
        for (i = 0; i < 8; i++)
            check(rounds_to(xombe__dct8[m][i], cos(XOMBE_PI * m * (i + 0.5) / 8.0)), "dct8", 8 * m + i);
    for (J = 1; J <= XOMBE_DCT_BLK_J_MAX; J++) {
        const float *kern = xombe__dct_blk + xombe__dct_blk_off[J];
        check(xombe__dct_blk_off[J] + XOMBE_DCT_BLK_ROWS * J <= XOMBE_DCT_BLK_SIZE, "dct_blk extent", J);
        for (m = 0; m < XOMBE_DCT_BLK_ROWS; m++)
            for (i = 0; i < J; i++)
                check(rounds_to(kern[m * J + i], cos(XOMBE_PI * m * (i + 0.5) / J)), "dct_blk", 100 * J + 10 * m + i);
    }
    for (L = 9; L <= XOMBE_LMAX; L++) {
        int sum = 0;
        for (k = 0; k < 4; k++) {
            check(xombe__blk[L][k] >= 1 && xombe__blk[L][k] <= XOMBE_DCT_BLK_J_MAX, "block length in kernel range",
                  4 * L + k);
            sum += xombe__blk[L][k];
        }
        check(sum == L, "blocks cover L", L);
        check(rounds_to(xombe__half_log2[L], 0.5 * log2((double)L)), "half_log2", L);
    }
    for (i = 0; i < XOMBE_GRAIN_RAMP; i++) {
        check(rounds_to(xombe__ramp_cos[i], cos(XOMBE_PI * i / XOMBE_GRAIN_RAMP)), "ramp_cos", i);
        check(rounds_to(xombe__ramp_sin[i], sin(XOMBE_PI * i / XOMBE_GRAIN_RAMP)), "ramp_sin", i);
    }
    for (n = 0; n <= XOMBE_HALF_MAX; n++)
        check(rounds_to(xombe__sqrt_n[n], sqrt((double)n)), "sqrt_n", n);
    for (N = XOMBE_HALF_MIN; N <= XOMBE_HALF_MAX; N++)
        check(rounds_to(xombe__rsqrt_half[N - XOMBE_HALF_MIN], 1.0 / sqrt((double)N)), "rsqrt_half", N);
    for (i = 0; i < XOMBE_B0_COUNT - 1; i++) {
        check(rounds_to(xombe__w0_mid[i], sqrt((double)xombe__w0[i] * xombe__w0[i + 1])), "w0_mid", i);
        check(xombe__w0_mid[i] < xombe__w0[i] && xombe__w0_mid[i] > xombe__w0[i + 1], "w0_mid between", i);
    }

    for (k = 0; k < XOMBE_TW_N / 2; k++) {
        double ang = -2.0 * XOMBE_PI * k / XOMBE_TW_N;
        check(rounds_to(xombe__fft_tw[k].re, cos(ang)) && rounds_to(xombe__fft_tw[k].im, sin(ang)), "fft_tw", k);
    }

#if XOMBE_HAVE_ENCODER
    {
        /* init-time tables against their closed forms */
        double s = 0.0, s2 = 0.0, s4 = 0.0, h[2 * XOMBE_PITCH_LPF_HALF + 1];
        float lo, hi;
        int j;
        for (i = 0; i < XOMBE_BASE_WINDOW; i++)
            check(rounds_to(xombe__base_win[i], ref_hamming(i, XOMBE_BASE_WINDOW)), "base_win", i);
        for (j = -XOMBE_PITCH_LPF_HALF; j <= XOMBE_PITCH_LPF_HALF; j++) {
            double x = 0.25 * j, sinc = j == 0 ? 1.0 : sin(XOMBE_PI * x) / (XOMBE_PI * x);
            h[j + XOMBE_PITCH_LPF_HALF] = 0.25 * sinc * (0.54 + 0.46 * cos(XOMBE_PI * j / XOMBE_PITCH_LPF_HALF));
            s += h[j + XOMBE_PITCH_LPF_HALF];
        }
        for (j = 0; j < 2 * XOMBE_PITCH_LPF_HALF + 1; j++)
            check(rounds_to(xombe__pitch_lpf[j], h[j] / s), "pitch_lpf", j);
        for (j = 0; j < XOMBE_PITCH_WIN; j++)
            s2 += ref_hamming(j, XOMBE_PITCH_WIN) * ref_hamming(j, XOMBE_PITCH_WIN);
        for (j = 0; j < XOMBE_PITCH_WIN; j++) {
            double w2 = ref_hamming(j, XOMBE_PITCH_WIN) * ref_hamming(j, XOMBE_PITCH_WIN) / s2;
            check(rounds_to(xombe__pitch_w2[j], w2), "pitch_w2", j);
            s4 += w2 * w2;
        }
        check(fabs(xombe__pitch_sum_w4 - s4) <= 6e-8 * s4, "pitch_sum_w4", 0);
        for (i = 0; i < XOMBE_TONE_SEG; i++)
            check(rounds_to(xombe__tone_win[i], 0.5 - 0.5 * cos(2.0 * XOMBE_PI * i / (XOMBE_TONE_SEG - 1))),
                  "tone_win", i);
        for (i = 0; i < XOMBE_TONE_IDS; i++) { /* filter gain at each tone component, float within 2.2e-5 dB */
            if (!xombe__tone_freqs(i, &lo, &hi))
                continue;
            check(fabs(xombe__tone_hpf_gain_db(lo) - ref_hpf_gain_db(lo)) <= 1e-4, "tone_hpf_gain_db", i);
            check(fabs(xombe__tone_hpf_gain_db(hi) - ref_hpf_gain_db(hi)) <= 1e-4, "tone_hpf_gain_db", i);
        }
    }
    /* window transform table reaches the widest harmonic fit: half the spacing at the highest pitch, widened
       by the voicing metric's 2% jitter, plus the interpolation's next point */
    check(0.5 * XOMBE_FFT * xombe__w0[0] * 1.02 * XOMBE_WTAB_OVERSAMPLE + 1.0 < XOMBE_WTAB_LEN - 1,
          "window transform covers the harmonic fits", XOMBE_WTAB_LEN);
    /* every analysis window, its energy, and its tabulated transform against double sums of the window */
    for (n = XOMBE_BASE_WINDOW; n <= XOMBE_MAX_WINDOW; n += 2) {
        static float w[XOMBE_MAX_WINDOW];
        double s2 = 0.0, w0 = 0.0, worst = 0.0;
        int h = (n - 1) / 2, q, j;
        xombe__hamming(w, n);
        for (i = 0; i < n; i++) {
            check(fabs(w[i] - ref_hamming(i, n)) <= 1.2e-7, "hamming", 1000 * n + i);
            s2 += ref_hamming(i, n) * ref_hamming(i, n);
            w0 += ref_hamming(i, n);
        }
        check(fabs(xombe__hamming_sum2(n) - s2) <= 2e-7 * s2, "hamming_sum2", n); /* three roundings */
        for (q = 0; q < XOMBE_WTAB_LEN; q++) {
            double nu = 2.0 * XOMBE_PI * q / ((double)XOMBE_WTAB_OVERSAMPLE * XOMBE_FFT), v = 0.0;
            for (j = -h; j <= h; j++)
                v += ref_hamming(j + h, n) * cos(nu * j);
            v = fabs(xombe__hamming_dtft(n, q, (long)XOMBE_WTAB_OVERSAMPLE * XOMBE_FFT) - v);
            if (v > worst)
                worst = v;
        }
        check(worst <= 1e-6 * w0, "hamming_dtft", n);
    }
#endif

    for (i = 0; i <= XOMBE_LMAX; i++)
        a[i] = (float)i;
    check(xombe__harmonic_lerp(a, 0, 0.7f) == 1.0f, "lerp holds a[1] below harmonic 1", 0);
    check(xombe__harmonic_lerp(a, -3, 0.0f) == 1.0f, "lerp holds a[1] for negative positions", -3);
    check(xombe__harmonic_lerp(a, XOMBE_LMAX, 0.5f) == (float)XOMBE_LMAX, "lerp holds a[LMAX] above", XOMBE_LMAX);
    check(fabsf(xombe__harmonic_lerp(a, 10, 0.25f) - 10.25f) < 1e-5f, "lerp interior", 10);

    printf("tables: %d/%d checks pass\n", checks - failures, checks);
    return failures ? 1 : 0;
}
