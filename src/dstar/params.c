#include <math.h>
#include <string.h>
#include "dstar_internal.h"

const xombe_amp_tables xombe__amp_dstar = {
    xombe__ds_dg, 64, xombe__ds_prba24, 512, xombe__ds_prba58, 128,
    {&xombe__ds_hoc5[0][0], &xombe__ds_hoc6[0][0], &xombe__ds_hoc7[0][0], &xombe__ds_hoc8[0][0]}, {16, 16, 16, 8},
    xombe__blk, xombe__ds_ok_dg, xombe__ds_ok_prba24, xombe__ds_ok_prba58,
    {xombe__ds_ok_hoc5, xombe__ds_ok_hoc6, xombe__ds_ok_hoc7, xombe__ds_ok_hoc8},
    XOMBE_DS_RHO, XOMBE_DS_GRHO, XOMBE_DS_GAIN_MAX, XOMBE_AMP_MIN, XOMBE_DS_AMP_MAX};

#define DS_RESET_L 15

float xombe__ds_f0(int b0)
{
    return xombe__ds_w0[b0 & (XOMBE_DS_B0_COUNT - 1)];
}

int xombe__ds_b0_from_f0(float f0)
{
    int b = (int)floorf((log2f(f0) - XOMBE_DS_P0) / XOMBE_DS_P1);
    return b < 0 ? 0 : b > XOMBE_DS_B0_VOICE_MAX ? XOMBE_DS_B0_VOICE_MAX : b;
}

int xombe__ds_L(int b0)
{
    return xombe__ds_lt[b0 & (XOMBE_DS_B0_COUNT - 1)];
}

void xombe__ds_harm_voicing(const uint8_t V[8], float f0, int L, uint8_t vl[XOMBE_LMAX + 1])
{
    int l;
    vl[0] = 0;
    for (l = 1; l <= L; l++) {
        /* scripts/gen_dstar_tables.py checks this float index matches the double one for every voice b0 and
           harmonic */
        int j = (int)((float)l * f0 * XOMBE_BANDS_PER_UNIT);
        vl[l] = V[j < 7 ? j : 7];
    }
}

void xombe__ds_amps(const int b[9], int L, int Lp, const float ap[XOMBE_LMAX + 1], float Gp,
                    float Tv[XOMBE_LMAX + 1], float *G, float a[XOMBE_LMAX + 1])
{
    const xombe_amp_tables *ts = &xombe__amp_dstar;
    float t[XOMBE_LMAX + 1], pt[XOMBE_LMAX + 1], g;
    xombe__spectral_shape_ts(ts, b, L, t);
    xombe__predict_from(Lp, ap, L, pt);
    g = xombe__gain_ts(ts, b[2], Gp);
    xombe__assemble_ts(ts, t, pt, L, g, a);
    *G = g;
    if (Tv)
        memcpy(Tv, t, sizeof(t));
}

void xombe__ds_mirror_init(xombe_ds_mirror *m)
{
    memset(m, 0, sizeof(*m));
    m->L = DS_RESET_L;
}

void xombe__ds_mirror_predict(const xombe_ds_mirror *m, int L, float pt[XOMBE_LMAX + 1])
{
    xombe__predict_from(m->L, m->a, L, pt);
}

void xombe__ds_dequant(const int b[9], const xombe_ds_mirror *m, xombe_ds_params *p)
{
    memcpy(p->b, b, sizeof(p->b));
    p->cls = xombe__ds_classify(b[0]);
    p->T = p->tid = p->AD = -1;
    if (p->cls != XOMBE_DS_VOICE) {
        p->L = 0;
        return;
    }
    p->f0 = xombe__ds_f0(b[0]);
    p->L = xombe__ds_L(b[0]);
    memcpy(p->V, xombe__ds_vuv[b[1] & 15], 8);
    xombe__ds_harm_voicing(p->V, p->f0, p->L, p->vl);
    xombe__ds_amps(b, p->L, m->L, m->a, m->G, p->Tv, &p->G, p->a);
}

int xombe__ds_mirror_commit(xombe_ds_mirror *m, const int b[9], xombe_ds_params *p)
{
    xombe_ds_params tmp;
    if (!p)
        p = &tmp;
    xombe__ds_dequant(b, m, p);
    if (p->cls != XOMBE_DS_VOICE)
        return 0;
    m->L = p->L;
    m->G = p->G;
    memcpy(m->a, p->a, sizeof(m->a));
    return 1;
}

int xombe__ds_params_from_frame(xombe_ds_mirror *m, const xombe_ds_frame *f, xombe_ds_params *p)
{
    int voice = xombe__ds_mirror_commit(m, f->b, p);
    p->T = f->T;
    p->tid = f->tid;
    p->AD = f->AD;
    return voice;
}
