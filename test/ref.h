#ifndef XOMBE_TEST_REF_H
#define XOMBE_TEST_REF_H

/* double-precision references for the generated tables (scripts/gen_dsp_tables.py) and the tests */

#include <math.h>
#include "../src/internal.h"
#include "../src/dsp/biquad.h"

static inline double ref_hamming(int i, int n)
{
    return 0.54 - 0.46 * cos(2.0 * XOMBE_PI * i / (n - 1));
}

/* solves a*x = b (n x n, row-major) in place by Gaussian elimination with partial pivoting, x replaces b
   returns 0 if a is singular */
static inline int ref_solve(double *a, double *b, int n)
{
    int i, j, k;
    for (i = 0; i < n; i++) {
        int piv = i;
        for (j = i + 1; j < n; j++)
            if (fabs(a[j * n + i]) > fabs(a[piv * n + i]))
                piv = j;
        if (fabs(a[piv * n + i]) < 1e-12)
            return 0;
        for (k = 0; k < n; k++) {
            double t = a[i * n + k];
            a[i * n + k] = a[piv * n + k];
            a[piv * n + k] = t;
        }
        {
            double t = b[i];
            b[i] = b[piv];
            b[piv] = t;
        }
        for (j = i + 1; j < n; j++) {
            double m = a[j * n + i] / a[i * n + i];
            for (k = i; k < n; k++)
                a[j * n + k] -= m * a[i * n + k];
            b[j] -= m * b[i];
        }
    }
    for (i = n - 1; i >= 0; i--) {
        for (k = i + 1; k < n; k++)
            b[i] -= a[i * n + k] * b[k];
        b[i] /= a[i * n + i];
    }
    return 1;
}

/* codec high-pass power gain at f_hz, dB */
static inline double ref_hpf_gain_db(double f_hz)
{
    double w = 2.0 * XOMBE_PI * f_hz / XOMBE_FS;
    double cr = cos(w), ci = -sin(w), c2r = cos(2 * w), c2i = -sin(2 * w);
    double zr = 1.0 - cr, zi = -ci;
    double nr = zr * zr - zi * zi, ni = 2.0 * zr * zi, n3r = nr * zr - ni * zi, n3i = nr * zi + ni * zr;
    double d1r = 1.0 - HPF_POLE1 * cr, d1i = -HPF_POLE1 * ci;
    double d2r = 1.0 - HPF_A1 * cr + HPF_A2 * c2r, d2i = -HPF_A1 * ci + HPF_A2 * c2i;
    double num = n3r * n3r + n3i * n3i;
    double den = (d1r * d1r + d1i * d1i) * (d2r * d2r + d2i * d2i);
    return 10.0 * log10(num / den);
}

#endif
