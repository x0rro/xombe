#include <string.h>
#include "kernels.h"

#define METRIC_TOP XOMBE_HZ(3950.0f) /* metric fits harmonics below this */
#define JITTER_STEPS 9
#define JITTER_STEP 0.005f
#define SOFT_OFFSET (-0.661f)
#define WEIGHT_FLOOR 1e-3f
#define FIRST_HALF_WEIGHT 0.5f
#define B1_SEARCH_MAX 17

#define METRIC_KMAX ((int)(METRIC_TOP * XOMBE_FFT) - 1)

void xombe__voicing_prepare(xombe_enc *e, const xombe_spec *S, xombe_lf_spec *sp)
{
    xombe__kern_lf_prepare(e, S, 1, METRIC_KMAX, sp);
}

void xombe__voicing_metric_sp(xombe_enc *e, const xombe_lf_spec *sp, float fhat, float m[8], float w[8])
{
    float br[8], be[8], be0[8];
    int kmax = METRIC_KMAX, j, s;
    for (j = 0; j < 8; j++)
        m[j] = 1e30f;
    for (s = 0; s < JITTER_STEPS; s++) {
        float eps = (float)(s - JITTER_STEPS / 2) * JITTER_STEP;
        xombe__kern_line_fit(e, sp, 1.0f / (fhat * (1.0f + eps)), 1, kmax, br, be);
        if (s == JITTER_STEPS / 2)
            memcpy(be0, be, sizeof(be0));
        for (j = 0; j < 8; j++) {
            float mj = be[j] > 0.0f ? br[j] / be[j] : 1.0f;
            if (mj < m[j])
                m[j] = mj;
        }
    }
    for (j = 0; j < 8; j++) {
        float wj = be0[j] / (XOMBE_FFT * e->win_sum2);
        w[j] = wj > WEIGHT_FLOOR ? wj : WEIGHT_FLOOR;
    }
}

void xombe__voicing_metric(xombe_enc *e, const xombe_spec *S, float fhat, float m[8], float w[8])
{
    xombe_lf_spec sp;
    xombe__voicing_prepare(e, S, &sp);
    xombe__voicing_metric_sp(e, &sp, fhat, m, w);
}

float xombe__voicing_map(float sm)
{
    float soft = sm > 0.0f ? SOFT_OFFSET - 0.5f * log2f(sm) : 1.0f;
    return soft < 0.0f ? 0.0f : soft > 1.0f ? 1.0f : soft;
}

void xombe__voicing_soft_sp(xombe_enc *e, const xombe_lf_spec *sp, float fhat, float scale, float a[8], float w[8])
{
    float m[8];
    int j;
    xombe__voicing_metric_sp(e, sp, fhat, m, w);
    for (j = 0; j < 8; j++)
        a[j] = xombe__voicing_map(scale * m[j]);
}

void xombe__voicing_soft(xombe_enc *e, const xombe_spec *S, float fhat, float scale, float a[8], float w[8])
{
    xombe_lf_spec sp;
    xombe__voicing_prepare(e, S, &sp);
    xombe__voicing_soft_sp(e, &sp, fhat, scale, a, w);
}

int xombe__choose_b1(const float a2[8], const float w2[8], const float a1[8], const float w1[8],
                     const uint8_t v2_prev[8])
{
    float best = 0.0f;
    int b1, j, choice = XOMBE_UNVOICED_B1;
    for (b1 = 0; b1 <= B1_SEARCH_MAX; b1++) {
        uint8_t v2[8], v1[8];
        float c = 0.0f;
        xombe__v2_from_b1(b1, v2);
        xombe__f3_first_half(v2_prev, v2, v1);
        for (j = 0; j < 8; j++) {
            float d2 = v2[j] ? 1.0f - a2[j] : a2[j];
            float d1 = v1[j] ? 1.0f - a1[j] : a1[j];
            c += w2[j] * d2 * d2 + FIRST_HALF_WEIGHT * w1[j] * d1 * d1;
        }
        if (b1 == 0 || c < best) {
            best = c;
            choice = b1;
        }
    }
    return choice;
}
