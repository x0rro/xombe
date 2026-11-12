#ifndef XOMBE_DKERNELS_H
#define XOMBE_DKERNELS_H

/* decoder inner loops, dkern_float.c or dkern_q15.c (XOMBE_KERNELS=q15) */

#include "../internal.h"

#define XOMBE_EPOCH_REF 16 /* sample of a grain the wavetable's phases refer to */

#define XOMBE_UV_FFT 256
#define XOMBE_UV_BINS (XOMBE_UV_FFT / 2 + 1)
/* sqrt(12)/2^16 takes uniform 16-bit noise to unit variance. unvoiced bands are normalized by their energy, so only
   the comfort-noise gain carries this scale */
#define XOMBE_UV_NOISE_SCALE (3.4641016151377544f / 65536.0f)
#define XOMBE_UV_GAIN 1.449f
#define XOMBE_CN_GAIN (1.147f * XOMBE_UV_NOISE_SCALE)
#define XOMBE_UV_TOP_BIN 121

/* one period of set s's voiced harmonics (vh[l], amplitudes M[l]) in XOMBE_WT_SIZE samples */
void xombe__dk_wavetable(const xombe_set *s, const float *M, const uint8_t *vh, float tab[XOMBE_WT_SIZE]);
/* adds one pitch period (period = 1/fhat samples) of tab starting at P to buf[0..N + XOMBE_TAIL), faded in over
   the first XOMBE_GRAIN_RAMP samples and out over XOMBE_GRAIN_RAMP samples after its end */
void xombe__dk_add_grain(float *buf, int N, float P, float fhat, float period, const float *tab);
/* y[0..2N): noise seg[0..2N) under sqrt(n/N) fades in and out, its spectrum shaped to set cur's unvoiced bands
   (every band for silence) with gains M[l], cleared below avg_f0 and above band L and XOMBE_UV_TOP_BIN */
void xombe__dk_unvoiced(const int16_t *seg, int N, const xombe_set *cur, const float *M, const uint8_t *vh,
                        float avg_f0, float y[XOMBE_UV_FFT]);

#endif
