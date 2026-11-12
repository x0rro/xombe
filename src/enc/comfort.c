#include <string.h>
#include "../internal.h"

#define CN_FAST_UPDATES 16
#define CN_CLIP_RATIO 4.0f
#define CN_CLIP_MIN_POWER 1.0f
#define CN_LEVEL_OFFSET 0.21f

void xombe__cn_half(xombe_enc *e, const float *P, int noise, int confident_run)
{
    int k;
    if (!e->cn_ready) {
        if (noise) {
            memcpy(e->cn_power, P, sizeof(e->cn_power));
            e->cn_ready = 1;
            e->cn_updates = 0;
        }
        return;
    }
    if (confident_run >= XOMBE_CONFIDENT_RUN) {
        float alpha = e->cn_updates < CN_FAST_UPDATES ? 0.5f : 0.1f;
        for (k = 0; k < XOMBE_SPEC_BINS; k++) {
            float ceiling = CN_CLIP_RATIO * (e->cn_power[k] > CN_CLIP_MIN_POWER ? e->cn_power[k] : CN_CLIP_MIN_POWER);
            float p = P[k] < ceiling ? P[k] : ceiling;
            e->cn_power[k] = (1.0f - alpha) * e->cn_power[k] + alpha * p;
        }
        e->cn_updates++;
    }
}

void xombe__cn_targets(const xombe_enc *e, float *a)
{
    int l;
    xombe__enc_integrate(e->cn_power, XOMBE_SILENCE_L, XOMBE_SILENCE_F0, a);
    for (l = 1; l <= XOMBE_SILENCE_L; l++)
        a[l] = xombe__clamp_amp(a[l] + CN_LEVEL_OFFSET);
}
