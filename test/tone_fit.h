#ifndef XOMBE_TEST_TONE_FIT_H
#define XOMBE_TEST_TONE_FIT_H

/* least-squares sinusoid fits of decoded tones, FIT_LEN samples */

#include <math.h>
#include <string.h>
#include "ref.h"

#define FIT_LEN 1600

static inline double sinusoid_fit(const double *x, int n, const double *freq, int nf, double *amp)
{
    double basis[4][FIT_LEN], g[16], r[4], energy = 0.0, explained = 0.0;
    int c, i, j, m = 2 * nf;
    for (c = 0; c < nf; c++) {
        double w = 2.0 * XOMBE_PI * freq[c] / 8000.0, cr = cos(w), ci = sin(w), pr = 1.0, pi = 0.0;
        for (i = 0; i < n; i++) {
            double t;
            basis[2 * c][i] = pr;
            basis[2 * c + 1][i] = pi;
            t = pr * cr - pi * ci;
            pi = pr * ci + pi * cr;
            pr = t;
        }
    }
    for (i = 0; i < m; i++) {
        r[i] = 0.0;
        for (j = 0; j < m; j++) {
            double s = 0.0;
            int k;
            for (k = 0; k < n; k++)
                s += basis[i][k] * basis[j][k];
            g[i * m + j] = s;
        }
        for (j = 0; j < n; j++)
            r[i] += basis[i][j] * x[j];
    }
    {
        double rr[4];
        memcpy(rr, r, sizeof(rr));
        if (!ref_solve(g, r, m))
            return 1e300;
        for (i = 0; i < m; i++)
            explained += r[i] * rr[i];
    }
    for (i = 0; i < n; i++)
        energy += x[i] * x[i];
    for (c = 0; c < nf; c++)
        amp[c] = hypot(r[2 * c], r[2 * c + 1]);
    return energy - explained;
}

static inline void search_freqs(const double *x, double *freq, int nf, double span)
{
    int c, round;
    double amp[2];
    for (round = 0; round < 2; round++) {
        for (c = 0; c < nf; c++) {
            double centre = freq[c], step = centre * span / 100.0, best = freq[c], bestr = 1e300;
            int pass, k;
            for (pass = 0; pass < 4; pass++) {
                for (k = -100; k <= 100; k++) {
                    double r;
                    freq[c] = centre + k * step;
                    r = sinusoid_fit(x, FIT_LEN, freq, nf, amp);
                    if (r < bestr) {
                        bestr = r;
                        best = freq[c];
                    }
                }
                centre = best;
                step /= 50.0;
            }
            freq[c] = best;
        }
    }
}

static inline double expected_db(int ad, double f)
{
    return 20.0 * log10(2.0 * exp2(xombe__tone_log2_amp(ad)) / 32768.0) + ref_hpf_gain_db(f);
}

#endif
