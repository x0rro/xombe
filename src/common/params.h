#ifndef XOMBE_PARAMS_H
#define XOMBE_PARAMS_H

#include "../state.h"

/* spectral-amplitude tables and constants of one codec, AMBE+2 (xombe__amp_ambe2) or D-STAR (xombe__amp_dstar)
     t = block-DCT spectral shape of b3..b8 (xombe__spectral_shape_ts), harmonics 1..L
     G = min(dg[b2] + grho*G_prev, gmax)
     a = clamp((t - mean t) + rho*(pt - mean pt) + G - 1/2 log2 L, amin, amax), pt = prediction from the previous frame
   ok_* mark entries the closed-loop quantizer may pick (1) and unfitted ones it skips (0), NULL means all */
typedef struct {
    const float *dg;
    int n_dg;
    const float (*prba24)[3];
    int n_prba24;
    const float (*prba58)[4];
    int n_prba58;
    const float *hoc[4]; /* b5..b8: rows of 4 floats, the first min(J, 6) - 2 used */
    int n_hoc[4];
    const unsigned char (*blk)[4]; /* block lengths J per L, [57][4] */
    const unsigned char *ok_dg, *ok_prba24, *ok_prba58, *ok_hoc[4];
    float rho, grho, gmax, amin, amax;
} xombe_amp_tables;

extern const xombe_amp_tables xombe__amp_ambe2;

typedef struct {
    int L;
    int pulse;
    float f0;
    float gamma;
    float a[XOMBE_LMAX + 1];
} xombe_params;

/* log amplitude of harmonics with no component (tone sets' other harmonics), synthesizes as 0 */
#define XOMBE_AMP_OFF (-1000.0f)

float xombe__clamp_amp(float a);
/* decoder's high-frequency boost: log2(f/XOMBE_HF_CORNER) above the corner, else 0 */
float xombe__hf_tilt(float f);
int xombe__L_from_f0(float f0);
/* row of tab (rows x stride floats) nearest t[0..used) in squared error, first on ties, rows with mask[j] 0
   skipped (mask NULL: none), 0 if every row is */
int xombe__nearest_row(const float *tab, int rows, int stride, int used, const float *t, const unsigned char *mask);
void xombe__pred_init(xombe_pred *pi);
/* (1 - f)*a[i] + f*a[i + 1] over log amplitudes a[1..XOMBE_LMAX], holding a[1] below harmonic 1 and
   a[XOMBE_LMAX] from harmonic XOMBE_LMAX up */
float xombe__harmonic_lerp(const float a[XOMBE_LMAX + 1], int i, float f);
void xombe__predict(const xombe_pred *pi, int L, float pt[XOMBE_LMAX + 1]);
/* prediction pt[1..L] from a previous frame of Lp harmonics, ap[1..XOMBE_LMAX] extended past Lp */
void xombe__predict_from(int Lp, const float ap[XOMBE_LMAX + 1], int L, float pt[XOMBE_LMAX + 1]);
/* HOC codebook of block k (b5..b8) and its row count, rows of 4 floats, bound to the AMBE+2 tables */
const float *xombe__hoc_table(int k, int *rows);
/* xombe__spectral_shape_ts bound to the AMBE+2 tables */
void xombe__spectral_shape(const int b[9], int L, float t[XOMBE_LMAX + 1]);
void xombe__spectral_shape_ts(const xombe_amp_tables *ts, const int b[9], int L, float t[XOMBE_LMAX + 1]);
/* xombe__assemble_ts bound to the AMBE+2 tables */
void xombe__assemble(const float t[XOMBE_LMAX + 1], const float pt[XOMBE_LMAX + 1], int L, float gamma,
                     float a[XOMBE_LMAX + 1]);
/* a[1..L] by the xombe_amp_tables model on set ts, a[L + 1..XOMBE_LMAX] = a[L] */
void xombe__assemble_ts(const xombe_amp_tables *ts, const float t[XOMBE_LMAX + 1], const float pt[XOMBE_LMAX + 1],
                        int L, float gamma, float a[XOMBE_LMAX + 1]);
/* G = min(dg[b2] + grho*Gprev, gmax) on set ts */
float xombe__gain_ts(const xombe_amp_tables *ts, int b2, float Gprev);
/* AMBE+2 gain, as xombe__gain_ts but a silence frame takes off XOMBE_SILENCE_GAIN_OFFSET before the cap */
float xombe__gain(int b2, float prev_gamma, int silence);
void xombe__dequant(const int b[9], const xombe_pred *pi, int silence, xombe_params *out);
void xombe__pred_store_voice(xombe_pred *pi, const xombe_params *p, const uint8_t v2[8]);
void xombe__pred_store_other(xombe_pred *pi, int cls);

#endif
