#include <math.h>
#include <string.h>
#include "dstar_enc.h"
#include "../tables/dstar_enc_tables.h"
#include "../common/prof.h"

#define DS_GAIN_DB 0.2f
#define DS_GAIN_LOG2 (DS_GAIN_DB / 6.02059991f) /* dB -> log2 amplitude */

XOMBE_STATIC_ASSERT(sizeof(xombe_ds_enc) <= XOMBE_DS_ENC_STATE, ds_enc_state_fits_XOMBE_DS_ENC_STATE);

void xombe__ds_enc_reset(xombe_ds_enc *s)
{
    s->s_vn = XOMBE_DS_S_VN;
    xombe__ds_mirror_init(&s->mirror);
}

void xombe__ds_enc_init(xombe_enc *e)
{
    xombe__ds_enc_reset(XOMBE_DS_ENC(e));
    xombe__ds_tone_init(&XOMBE_DS_ENC(e)->tones);
    XOMBE_DS_ENC(e)->ws.cap = XOMBE_CQ_CAP;
}

void xombe__ds_m_noise(int N, int b0, float mn[8])
{
    int i = (N - XOMBE_DS_MN_N_MIN) / 2, lo, hi, j;
    const unsigned short *r;
    if (i < 0)
        i = 0;
    if (i > XOMBE_DS_MN_WINDOWS - 1)
        i = XOMBE_DS_MN_WINDOWS - 1;
    lo = xombe__ds_mn_lo[i];
    hi = lo + xombe__ds_mn_cnt[i] - 1;
    if (b0 < lo)
        b0 = lo;
    if (b0 > hi)
        b0 = hi;
    r = xombe__ds_mn[xombe__ds_mn_off[i] + b0 - lo];
    for (j = 0; j < 8; j++)
        mn[j] = (float)r[j] * (1.0f / XOMBE_DS_MN_SCALE);
}

int xombe__ds_b1_search(const float a[8], const float w[8])
{
    float best = 0.0f;
    int c, j, bi = 0;
    for (c = 0; c < 16; c++) {
        const unsigned char *V = xombe__ds_vuv[c];
        float cost = 0.0f;
        if (c == XOMBE_DS_B1_SKIP)
            continue;
        for (j = 0; j < 8; j++) {
            float d = V[j] ? 1.0f - a[j] : a[j];
            cost += w[j] * d * d;
        }
        if (c == 0 || cost < best) {
            best = cost;
            bi = c;
        }
    }
    return bi;
}

int xombe__ds_soft(const float m[8], const float mn[8], float s_vn, int speech2, float a[8])
{
    int j, any = 0;
    for (j = 0; j < 8; j++) {
        a[j] = speech2 ? xombe__voicing_map(s_vn * (m[j] / mn[j])) : 0.0f;
        any |= a[j] > 0.5f;
    }
    return !any;
}

void xombe__ds_enc_voicing(xombe_enc *e, xombe_fe_frame *a, xombe_be_ns *ns)
{
    xombe_ds_enc *s = XOMBE_DS_ENC(e);
    xombe_ds_enc_frame *fr = &s->fr;
    memset(fr, 0, sizeof(*fr));
    fr->no_band_voiced = 1;
    fr->cls = XOMBE_DS_CLS_NULL;
    if (e->cfg.tones) {
        fr->det = s->tones.ndet;
        fr->trk = s->tones.ntrk;
        fr->tone = s->tones.tone;
    }
    if (fr->tone.on) {
        fr->cls = fr->tone.T < 0 ? XOMBE_DS_CLS_KNOX : XOMBE_DS_CLS_TONE;
    } else if (a->active || !e->cfg.dstar_dtx) {
        fr->cls = XOMBE_DS_CLS_VOICE;
        fr->voice = 1;
        XOMBE_PF(8, xombe__fe_voice(e, a));
        fr->b0 = xombe__ds_b0_from_f0(a->fstar);
        fr->f0 = xombe__ds_f0(fr->b0);
        fr->L = xombe__ds_L(fr->b0);
        XOMBE_PF(9, xombe__voicing_metric_sp(e, xombe__fe_s2_prepared(e, a), fr->f0, fr->m, fr->w));
        xombe__ds_m_noise(a->nw, fr->b0, fr->mn);
        fr->no_band_voiced = xombe__ds_soft(fr->m, fr->mn, s->s_vn, a->speech2, fr->a);
        fr->b1 = xombe__ds_b1_search(fr->a, fr->w);
        memcpy(fr->V, xombe__ds_vuv[fr->b1], 8);
    }
    ns->no_band_voiced[0] = ns->no_band_voiced[1] = fr->no_band_voiced;
    ns->mem = &s->ns_mem;
    if (e->cfg.ns)
        xombe__ns_bin_voicing(fr->V, 0, ns->vbin);
}

void xombe__ds_enc_targets(const xombe_spec *S, float win_sum2, float f0, int L, const float *G, const uint8_t V[8],
                           float at[XOMBE_LMAX + 1])
{
    float P[XOMBE_SPEC_BINS];
    int l;
    xombe__kern_bin_power(S, win_sum2, P);
    xombe__enc_integrate(P, L, f0, at);
    xombe__enc_hf_precomp(at, L, f0);
    if (G) {
        uint8_t vh[XOMBE_LMAX + 1];
        float da[XOMBE_LMAX + 1];
        xombe__ds_harm_voicing(V, f0, L, vh);
        xombe__ns_harmonic_offsets(G, L, f0, vh, da);
        for (l = 1; l <= L; l++)
            at[l] += da[l];
    }
    for (l = 1; l <= L; l++) {
        at[l] += DS_GAIN_LOG2;
        if (at[l] < XOMBE_AMP_MIN)
            at[l] = XOMBE_AMP_MIN;
    }
}

void xombe__ds_enc_amps(xombe_ds_enc *s, int b0, int b1, const float at[XOMBE_LMAX + 1], int b[9])
{
    float w[XOMBE_LMAX + 1], pt[XOMBE_LMAX + 1];
    xombe_cq_result res;
    int L = xombe__ds_L(b0), k;
    xombe__cq_amp_weights(at, L, XOMBE_CQ_QWEIGHT, XOMBE_CQ_WEIGHT_FLOOR, w);
    xombe__ds_mirror_predict(&s->mirror, L, pt);
    XOMBE_PF(13, xombe__cq_quantize(&s->ws, &xombe__amp_dstar, at, L, pt, s->mirror.G, w, 1, XOMBE_CQ_N_PAIRS,
                                    XOMBE_CQ_N_GAIN, XOMBE_CQ_ITERS, XOMBE_CQ_CLIP_PASSES, &res));
    XOMBE_PF_ADD(15, L);
    b[0] = b0;
    b[1] = b1;
    for (k = 2; k < 9; k++)
        b[k] = res.b[k];
    xombe__ds_mirror_commit(&s->mirror, b, NULL);
}

void xombe__ds_enc_frame(xombe_enc *e, const xombe_fe_frame *a, const xombe_be_ns *ns, uint8_t out[7])
{
    static const uint8_t null_bits[7] = XOMBE_DS_NULL_BITS;
    xombe_ds_enc *s = XOMBE_DS_ENC(e);
    const xombe_ds_enc_frame *fr = &s->fr;
    float at[XOMBE_LMAX + 1];
    int b[9];
    if (fr->cls == XOMBE_DS_CLS_KNOX) { /* voice frame from the tone table, no VAD gate or suppression */
        int k;
        for (k = 0; s->tones.knox[k].id != -fr->tone.T; k++)
            ;
        xombe__ds_knox_amps(s, k, fr->tone.level, b);
        s->fr.b0 = b[0];
        s->fr.b1 = b[1];
        xombe__ds_fields_bits(b, out);
        return;
    }
    if (!fr->voice) { /* tone and NULL frames leave the mirror alone */
        if (fr->cls == XOMBE_DS_CLS_TONE) {
            uint8_t d[XOMBE_DS_DATA_BITS];
            xombe__ds_tone_bits(fr->tone.T, xombe__ds_tone_q(fr->tone.level), d);
            xombe__ds_pack(d, out);
        } else {
            memcpy(out, null_bits, 7);
        }
        return;
    }
    XOMBE_PF(14, xombe__ds_enc_targets(&a->S2, a->win_sum2, fr->f0, fr->L, e->cfg.ns ? ns->G2 : NULL, fr->V,
                                       at));
    xombe__ds_enc_amps(s, fr->b0, fr->b1, at, b);
    xombe__ds_fields_bits(b, out);
}
