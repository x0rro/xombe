#ifndef XOMBE_KERN_Q15_H
#define XOMBE_KERN_Q15_H

/* helpers shared by the q15 kernels (kern_q15_*.c) */

#include <float.h>
#include "kernels.h"
#include "../dsp/q15ops.h"

static inline int round_q(float v)
{
    return (int)(v + (v < 0.0f ? -0.5f : 0.5f));
}

/* block scale taking peak magnitude m to peak, 0 when m is 0 or subnormal enough to overflow the scale (block
   counts as zero) */
static inline float q15_gain(float m, float peak)
{
    float g = m > 0.0f ? peak / m : 0.0f;
    return g <= FLT_MAX ? g : 0.0f;
}

/* |v| < 2^63 to float without the library call, top 31 bits through a 32-bit conversion */
static inline float i64_to_float(int64_t v)
{
    int64_t a = v < 0 ? -v : v;
    int sh = (a >> 31) ? 33 - q15_clz64((uint64_t)a) : 0;
    float f = (float)(int32_t)(a >> sh) * (float)(1u << sh);
    return v < 0 ? -f : f;
}

#endif
