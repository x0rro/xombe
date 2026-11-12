#ifndef XOMBE_KERNELS_H
#define XOMBE_KERNELS_H

/* encoder inner loops, kern_float.c or kern_q15_*.c (XOMBE_KERNELS=q15) */

#include <math.h>
#include "../internal.h"

/* analysis spectrum S (xombe_spec, layout per backend in kern_state.h): XOMBE_SPEC_BINS bins of the
   XOMBE_FFT-point transform of centre[-h..h] under window w of n = 2h + 1 points, centred on sample 0 */
void xombe__kern_spectrum(const float *centre, const float *w, int n, xombe_spec *S);
/* P[k] = |S[k]|^2 / (XOMBE_FFT * sum_w2), XOMBE_SPEC_BINS float bins */
void xombe__kern_bin_power(const xombe_spec *S, float sum_w2, float *P);

/* called by xombe__enc_set_window for a new e->win_len, after it rebuilds e->win and e->wtab (float), or to build
   e->win itself (q15, which leaves e->wtab unused) */
void xombe__kern_window_init(xombe_enc *e);

/* autocorrelation of the low-passed, windowed signal around centre into slot's pacf, energy into pe0[slot] */
void xombe__kern_pitch_acf(xombe_enc *e, const float *centre, int slot);
/* P on the quarter-sample grid, slot not silent */
float xombe__kern_error_at(const xombe_enc *e, int slot, float P);

/* spectrum prepared for line fits, may keep data in e so valid only until the next prepare */
typedef struct {
    const xombe_spec *S;
    float g;
} xombe_lf_spec;
/* kmin >= 1, later fits must stay within kmin..kmax */
void xombe__kern_lf_prepare(xombe_enc *e, const xombe_spec *S, int kmin, int kmax, xombe_lf_spec *sp);
/* returns the residual of the least-squares fit of the window transform, one complex amplitude per harmonic,
   to bins kmin..kmax, with per-band residual and bin energy.
   q15: XOMBE_FFT / (2*XOMBE_LF_PAD_BINS - 4) <= period <= XOMBE_FFT, kmax <= XOMBE_FFT / 2 */
float xombe__kern_line_fit(const xombe_enc *e, const xombe_lf_spec *sp, float period, int kmin, int kmax,
                           float band_resid[8], float band_energy[8]);

/* voicing metric on a spectrum prepared once, to try several pitches, valid until the next prepare */
void xombe__voicing_prepare(xombe_enc *e, const xombe_spec *S, xombe_lf_spec *sp);
void xombe__voicing_metric_sp(xombe_enc *e, const xombe_lf_spec *sp, float fhat, float m[8], float w[8]);
void xombe__voicing_soft_sp(xombe_enc *e, const xombe_lf_spec *sp, float fhat, float scale, float a[8], float w[8]);

/* nearest harmonic of bin k, midway goes up. computed from period, not 1/period, so 1/8-sample grid periods
   are exact and midway bins assign the same in any precision */
static inline int xombe__harmonic_of(int k, float period)
{
    return (int)floorf((float)k * period / XOMBE_FFT + 0.5f);
}

static inline int xombe__band_of(int l, float period)
{
    int band = (int)floorf((float)l * XOMBE_BANDS_PER_UNIT / period);
    return band > 7 ? 7 : band;
}

#endif
