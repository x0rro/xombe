#include <math.h>
#include <string.h>
#include "params.h"
#include "consts.h"
#include "../tables/tables.h"
#include "../tables/dsp_tables.h"

const xombe_amp_tables xombe__amp_ambe2 = {
    xombe__dg, 32, xombe__prba24, 512, xombe__prba58, 128,
    {&xombe__hoc5[0][0], &xombe__hoc6[0][0], &xombe__hoc7[0][0], &xombe__hoc8[0][0]}, {32, 16, 16, 8},
    xombe__blk, 0, 0, 0, {0, 0, 0, 0},
    XOMBE_RHO, 0.5f, XOMBE_GAIN_MAX, XOMBE_AMP_MIN, XOMBE_AMP_MAX};

float xombe__clamp_amp(float a)
{
    if (a < XOMBE_AMP_MIN)
        return XOMBE_AMP_MIN;
    if (a > XOMBE_AMP_MAX)
        return XOMBE_AMP_MAX;
    return a;
}

float xombe__hf_tilt(float f)
{
    return f > XOMBE_HF_CORNER ? log2f(f / XOMBE_HF_CORNER) : 0.0f;
}

/* L = floor(0.4627/f0): harmonics below about 3.7 kHz, clamped to 9..XOMBE_LMAX */
#define L_SCALE 0.4627f
#define L_MIN 9

int xombe__L_from_f0(float f0)
{
    int L = (int)floorf(L_SCALE / f0);
    if (L < L_MIN)
        L = L_MIN;
    if (L > XOMBE_LMAX)
        L = XOMBE_LMAX;
    return L;
}

int xombe__nearest_row(const float *tab, int rows, int stride, int used, const float *t, const unsigned char *mask)
{
    float best = 0.0f;
    int j, m, bi = 0, found = 0;
    for (j = 0; j < rows; j++) {
        float d = 0.0f;
        if (mask && !mask[j])
            continue;
        for (m = 0; m < used; m++) {
            float e = tab[j * stride + m] - t[m];
            d += e * e;
        }
        if (!found || d < best) {
            found = 1;
            best = d;
            bi = j;
        }
    }
    return bi;
}

void xombe__pred_init(xombe_pred *pi)
{
    memset(pi, 0, sizeof(*pi));
    pi->cls = XOMBE_CLASS_VOICE;
    pi->L = XOMBE_SILENCE_L;
    pi->f0 = XOMBE_SILENCE_F0;
}

float xombe__harmonic_lerp(const float a[XOMBE_LMAX + 1], int i, float f)
{
    if (i < 1)
        return a[1];
    if (i >= XOMBE_LMAX)
        return a[XOMBE_LMAX];
    return (1.0f - f) * a[i] + f * a[i + 1];
}

void xombe__predict(const xombe_pred *pi, int L, float pt[XOMBE_LMAX + 1])
{
    xombe__predict_from(pi->L, pi->a, L, pt);
}

void xombe__predict_from(int Lp, const float ap[XOMBE_LMAX + 1], int L, float pt[XOMBE_LMAX + 1])
{
    int l;
    pt[0] = 0.0f;
    for (l = 1; l <= L; l++) {
        /* harmonic l maps to l*Lp/L of the previous frame, split exactly into i + num % L / L */
        int num = l * Lp;
        int i = num / L;
        float f = (float)(num % L) / (float)L;
        /* i reaches XOMBE_LMAX only for l = L with Lp = XOMBE_LMAX, so f = 0 and a[XOMBE_LMAX] is in bounds */
        pt[l] = xombe__harmonic_lerp(ap, i, f);
    }
}

const float *xombe__hoc_table(int k, int *rows)
{
    *rows = xombe__amp_ambe2.n_hoc[k];
    return xombe__amp_ambe2.hoc[k];
}

void xombe__spectral_shape(const int b[9], int L, float t[XOMBE_LMAX + 1])
{
    xombe__spectral_shape_ts(&xombe__amp_ambe2, b, L, t);
}

void xombe__spectral_shape_ts(const xombe_amp_tables *ts, const int b[9], int L, float t[XOMBE_LMAX + 1])
{
    float g[9], r[8];
    int i, m, k, pos = 1;
    g[1] = 0.0f;
    for (m = 0; m < 3; m++)
        g[2 + m] = ts->prba24[b[3]][m];
    for (m = 0; m < 4; m++)
        g[5 + m] = ts->prba58[b[4]][m];
    for (i = 0; i < 8; i++) {
        float s = 0.0f;
        for (m = 1; m <= 8; m++)
            s += (m == 1 ? 1.0f : 2.0f) * g[m] * xombe__dct8[m - 1][i];
        r[i] = s;
    }
    for (k = 0; k < 4; k++) {
        int J = ts->blk[L][k];
        int top = J < 6 ? J : 6;
        float c[7];
        const float *hoc = ts->hoc[k] + 4 * b[5 + k];
        const float *kern = xombe__dct_blk + xombe__dct_blk_off[J];
        c[1] = (r[2 * k] + r[2 * k + 1]) * 0.5f;
        c[2] = (r[2 * k] - r[2 * k + 1]) * 0.35355339f; /* 1/(2*sqrt 2) */

        for (m = 3; m <= 6; m++)
            c[m] = m <= top ? hoc[m - 3] : 0.0f;
        for (i = 0; i < J; i++) {
            float s = 0.0f;
            for (m = 2; m <= top; m++)
                s += c[m] * kern[(m - 1) * J + i];
            t[pos++] = c[1] + 2.0f * s;
        }
    }
    t[0] = 0.0f;
}

void xombe__assemble(const float t[XOMBE_LMAX + 1], const float pt[XOMBE_LMAX + 1], int L, float gamma,
                     float a[XOMBE_LMAX + 1])
{
    xombe__assemble_ts(&xombe__amp_ambe2, t, pt, L, gamma, a);
}

void xombe__assemble_ts(const xombe_amp_tables *ts, const float t[XOMBE_LMAX + 1], const float pt[XOMBE_LMAX + 1],
                        int L, float gamma, float a[XOMBE_LMAX + 1])
{
    float mt = 0.0f, mp = 0.0f, off;
    int l;
    for (l = 1; l <= L; l++) {
        mt += t[l];
        mp += pt[l];
    }
    mt /= (float)L;
    mp /= (float)L;
    off = gamma - xombe__half_log2[L];
    a[0] = 0.0f;
    for (l = 1; l <= L; l++) {
        float v = (t[l] - mt) + ts->rho * (pt[l] - mp) + off;
        a[l] = v < ts->amin ? ts->amin : v > ts->amax ? ts->amax : v;
    }
    for (l = L + 1; l <= XOMBE_LMAX; l++)
        a[l] = a[L];
}

float xombe__gain_ts(const xombe_amp_tables *ts, int b2, float Gprev)
{
    float G = ts->dg[b2] + ts->grho * Gprev;
    return G < ts->gmax ? G : ts->gmax;
}

float xombe__gain(int b2, float prev_gamma, int silence)
{
    const xombe_amp_tables *ts = &xombe__amp_ambe2;
    float g = ts->dg[b2] + ts->grho * prev_gamma - (silence ? XOMBE_SILENCE_GAIN_OFFSET : 0.0f);
    return g < ts->gmax ? g : ts->gmax;
}

void xombe__dequant(const int b[9], const xombe_pred *pi, int silence, xombe_params *out)
{
    float t[XOMBE_LMAX + 1], pt[XOMBE_LMAX + 1];
    out->pulse = 0;
    if (silence) {
        out->f0 = XOMBE_SILENCE_F0;
        out->L = XOMBE_SILENCE_L;
    } else if (b[1] >= XOMBE_PULSE_B1_MIN) {
        out->pulse = 1;
        out->f0 = XOMBE_PULSE_F0;
        out->L = XOMBE_PULSE_L;
    } else {
        out->f0 = xombe__w0[b[0]];
        out->L = xombe__lt[b[0]];
    }
    out->gamma = xombe__gain(b[2], pi->gamma, silence);
    xombe__spectral_shape(b, out->L, t);
    xombe__predict(pi, out->L, pt);
    xombe__assemble(t, pt, out->L, out->gamma, out->a);
}

void xombe__pred_store_voice(xombe_pred *pi, const xombe_params *p, const uint8_t v2[8])
{
    int l;
    pi->cls = XOMBE_CLASS_VOICE;
    pi->L = p->L;
    pi->f0 = p->f0;
    pi->gamma = p->gamma;
    for (l = 0; l <= XOMBE_LMAX; l++)
        pi->a[l] = l <= p->L ? p->a[l] : p->a[p->L];
    memcpy(pi->v2, v2, 8);
}

void xombe__pred_store_other(xombe_pred *pi, int cls)
{
    pi->cls = cls;
    memset(pi->v2, 0, 8);
}
