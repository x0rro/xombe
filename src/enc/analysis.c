#include <math.h>
#include <string.h>
#include "../internal.h"
#include "../tables/dsp_tables.h"
#include "kernels.h"

void xombe__enc_set_window(xombe_enc *e, int n)
{
    if (e->win_len == n)
        return;
    e->win_len = n;
    e->win_sum2 = xombe__hamming_sum2(n);
#if !XOMBE_KERNELS_Q15 /* q15: the kernel builds e->win and its own table, e->wtab stays unused */
    {
        int q;
        xombe__hamming(e->win, n);
        for (q = 0; q < XOMBE_WTAB_LEN; q++)
            e->wtab[q] = xombe__hamming_dtft(n, q, (long)XOMBE_WTAB_OVERSAMPLE * XOMBE_FFT);
    }
#endif
    xombe__kern_window_init(e);
}

void xombe__enc_integrate(const float *P, int L, float f0, float *a)
{
    float step = f0 * XOMBE_FFT;
    int l, k;
    a[0] = 0.0f;
    for (l = 1; l <= L; l++) {
        float lo = ((float)l - 0.5f) * step, hi = ((float)l + 0.5f) * step, m2 = 0.0f;
        int k0 = (int)floorf(lo + 0.5f), k1 = (int)ceilf(hi + 0.5f);
        if (k0 < 0)
            k0 = 0;
        if (k1 > XOMBE_SPEC_BINS - 1)
            k1 = XOMBE_SPEC_BINS - 1;
        /* bin overlap with [lo, hi], bounds may include bins that do not overlap */
        for (k = k0; k <= k1; k++) {
            float c = fminf((float)k + 0.5f, hi) - fmaxf((float)k - 0.5f, lo);
            if (c > 0.0f)
                m2 += (c > 1.0f ? 1.0f : c) * P[k];
        }
        a[l] = m2 > 0.0f ? 0.5f * log2f(m2) : XOMBE_AMP_MIN;
        if (a[l] < XOMBE_AMP_MIN)
            a[l] = XOMBE_AMP_MIN;
    }
}

void xombe__enc_hf_precomp(float *a, int L, float f0)
{
    int l;
    for (l = 1; l <= L; l++)
        a[l] -= 0.5f * xombe__hf_tilt((float)l * f0);
}
