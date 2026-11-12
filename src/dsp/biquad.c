#include "../internal.h"
#include "biquad.h"

void xombe__hpf_init(xombe_hpf *f, float gain)
{
    f->x1 = f->s1 = f->s2 = f->y1 = f->y2 = 0.0f;
    f->gain = gain;
}

float xombe__hpf_step(xombe_hpf *f, float x)
{
    float s = x - f->x1 + (float)HPF_POLE1 * f->s1;
    float y = s - 2.0f * f->s1 + f->s2 + (float)HPF_A1 * f->y1 - (float)HPF_A2 * f->y2;
    f->x1 = x;
    f->s2 = f->s1;
    f->s1 = s;
    f->y2 = f->y1;
    f->y1 = y;
    return y * f->gain;
}
