#include <math.h>
#include "../internal.h"

#define LAMBDA_FAST 0.98f
#define SMOOTH_GAP 2.0f
#define SMOOTH_MARGIN 5.5f

void xombe__smooth(xombe_dec_core *d, float *a, int L, float f0)
{
    float sum = 0.0f, r, tau;
    int l;
    for (l = 1; l <= L; l++)
        sum += a[l] * a[l];
    r = sqrtf(sum / (float)L);
    if (r > d->lambda - 1.0f)
        d->lambda = LAMBDA_FAST * d->lambda + (1.0f - LAMBDA_FAST) * r;
    else
        d->lambda -= d->rules->lambda_decay;
    if (d->eps_t > 3 || d->eps_1 > 1) {
        float ratio = f0 / d->pi.f0;
        tau = d->lambda + SMOOTH_MARGIN - 0.5f * (float)d->eps_1;
        for (l = 1; l <= L; l++) {
            float x = (float)l * ratio;
            int i = (int)x;
            float ps = xombe__harmonic_lerp(d->pi.a, i, x - (float)i);
            float g = 6.0f * (float)l * f0;
            float e;
            if (g > 2.0f)
                g = 2.0f;
            e = a[l] - (tau - g);
            if (a[l] - ps > SMOOTH_GAP && e > 0.0f)
                a[l] += (ps - a[l]) * (e < 0.5f ? e : 0.5f);
            a[l] = xombe__clamp_amp(a[l]);
        }
    }
    for (l = L + 1; l <= XOMBE_LMAX; l++)
        a[l] = a[L];
}
