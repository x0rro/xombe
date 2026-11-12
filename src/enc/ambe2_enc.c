#include <math.h>
#include <string.h>
#include "frontend.h"
#include "../tables/dsp_tables.h"
#include "../common/cquant.h"
#if XOMBE_HAVE_DSTAR
#include "../dstar/dstar_enc.h"

/* D-STAR back end's workspace: the two emits never overlap, and xombe__cq_prepare rebuilds its maps when the
   table set changes, so D-STAR output is unaffected */
#define A2_CQ_WS(e) (&XOMBE_DS_ENC(e)->ws)
#else
#define A2_CQ_WS(e) ((xombe_cq_ws *)(void *)(e)->cq.bytes)
XOMBE_STATIC_ASSERT(sizeof(xombe_cq_ws) <= XOMBE_ENC_CQ_SIZE, cq_ws_fits_XOMBE_ENC_CQ_SIZE);
#endif

/* AMBE+2 back end, phases B and D of frontend.h */

/* pitch index nearest fhat on a log scale, ties to the higher pitch (xombe__w0 descends) */
static int nearest_pitch_index(float fhat)
{
    int b = 0;
    while (b < XOMBE_B0_COUNT - 1 && fhat < xombe__w0_mid[b])
        b++;
    return b;
}

static int any_above_half(const float a[8])
{
    int j;
    for (j = 0; j < 8; j++)
        if (a[j] > 0.5f)
            return 1;
    return 0;
}

static void choose_pitch_and_voicing(xombe_enc *e, xombe_fe_frame *a, xombe_a2_frame *f)
{
    float w1[8], w2[8], f0, fh1;
    int *b = f->b;
    xombe__fe_voice(e, a);
    b[0] = nearest_pitch_index(a->fstar);
    f0 = xombe__w0[b[0]];
    fh1 = e->pi.cls == XOMBE_CLASS_VOICE ? sqrtf(e->pi.f0 * f0) : f0;
    xombe__voicing_soft(e, &a->S1, fh1, e->voicing_scale, f->soft1, w1);
    xombe__voicing_soft_sp(e, xombe__fe_s2_prepared(e, a), f0, e->voicing_scale, f->soft2, w2);
    if (!a->speech1)
        memset(f->soft1, 0, sizeof(f->soft1));
    if (!a->speech2)
        memset(f->soft2, 0, sizeof(f->soft2));
    b[1] = xombe__choose_b1(f->soft2, w2, f->soft1, w1, e->pi.v2);
    xombe__v2_from_b1(b[1], f->v2);
}

/* voice frames only, silence's gain offset and fixed L stay open loop. the workspace starts zeroed
   (xombe_enc_init), cap is set on every call */
static void quantize_closed(xombe_enc *e, int L, const float *target, int b[9])
{
    xombe_cq_ws *ws = A2_CQ_WS(e);
    float w[XOMBE_LMAX + 1];
    xombe_cq_result res;
    int k;
    ws->cap = XOMBE_CQ_CAP;
    xombe__cq_amp_weights(target, L, XOMBE_CQ_QWEIGHT, XOMBE_CQ_WEIGHT_FLOOR, w);
    xombe__cq_search_ambe2(ws, &e->pi, target, L, w, XOMBE_CQ_N_PAIRS, XOMBE_CQ_N_GAIN, &res, NULL, NULL);
    for (k = 2; k < 9; k++)
        b[k] = res.b[k];
}

static void emit_quantized(xombe_enc *e, int silence, int L, float target[XOMBE_LMAX + 1], int b[9],
                           const uint8_t v2[8], uint8_t d[49])
{
    int l;
    for (l = 1; l <= L; l++)
        target[l] = xombe__clamp_amp(target[l]);
    XOMBE_TAP(e, XOMBE_STAGE_TARGET, silence ? XOMBE_SILENCE_F0 : xombe__w0[b[0]], target + 1, L);
    if (!silence && e->cfg.dmr_quant == XOMBE_ENC_QUANT_CLOSED)
        quantize_closed(e, L, target, b);
    else
        xombe__quantize(&e->pi, silence, L, target, b);
    memset(d, 0, 49);
    xombe__put_fields(b, d);
    if (silence) {
        xombe__pred_store_other(&e->pi, XOMBE_CLASS_SILENCE);
    } else {
        xombe_params rec;
        xombe__dequant(b, &e->pi, 0, &rec);
        xombe__pred_store_voice(&e->pi, &rec, v2);
    }
}

static void emit_voice(xombe_enc *e, const xombe_fe_frame *a, const xombe_a2_frame *f, const xombe_be_ns *ns,
                       int b[9], uint8_t d[49])
{
    float f0 = xombe__w0[b[0]], Pw[XOMBE_SPEC_BINS], target[XOMBE_LMAX + 1];
    int L = xombe__lt[b[0]];
    xombe__kern_bin_power(&a->S2, a->win_sum2, Pw);
    xombe__enc_integrate(Pw, L, f0, target);
    xombe__enc_hf_precomp(target, L, f0);
    if (e->cfg.ns)
        xombe__ns_apply_voice(ns->G2, target, L, f0, f->v2);
    emit_quantized(e, 0, L, target, b, f->v2, d);
}

static void emit_silence(xombe_enc *e, uint8_t d[49])
{
    static const uint8_t no_voicing[8];
    float target[XOMBE_LMAX + 1];
    int b[9];
    xombe__cn_targets(e, target);
    if (e->cfg.ns)
        xombe__ns_apply_silence(e, target, XOMBE_SILENCE_L, XOMBE_SILENCE_F0);
    b[0] = XOMBE_SILENCE_B0;
    b[1] = XOMBE_SILENCE_B1;
    emit_quantized(e, 1, XOMBE_SILENCE_L, target, b, no_voicing, d);
}

void xombe__a2_voicing(xombe_enc *e, xombe_fe_frame *a, xombe_a2_frame *f, xombe_be_ns *ns)
{
    f->cls = a->tone ? XOMBE_CLASS_TONE : (a->active || !e->cfg.dmr_dtx) ? XOMBE_CLASS_VOICE : XOMBE_CLASS_SILENCE;
    memset(f->v2, 0, sizeof(f->v2));
    memset(f->soft1, 0, sizeof(f->soft1));
    memset(f->soft2, 0, sizeof(f->soft2));
    if (f->cls == XOMBE_CLASS_VOICE)
        choose_pitch_and_voicing(e, a, f);
    ns->no_band_voiced[0] = !any_above_half(f->soft1);
    ns->no_band_voiced[1] = !any_above_half(f->soft2);
    ns->mem = &e->ns_mem;
    if (e->cfg.ns)
        xombe__ns_bin_voicing(f->v2, 1, ns->vbin);
}

void xombe__a2_emit(xombe_enc *e, const xombe_fe_frame *a, xombe_a2_frame *f, const xombe_be_ns *ns,
                    uint8_t bits[7])
{
    uint8_t d[49];
    if (f->cls == XOMBE_CLASS_TONE) {
        xombe__tone_frame_bits(a->tone_id, a->tone_ad, d);
        xombe__pred_store_other(&e->pi, XOMBE_CLASS_TONE);
    } else if (f->cls == XOMBE_CLASS_VOICE) {
        emit_voice(e, a, f, ns, f->b, d);
    } else {
        emit_silence(e, d);
    }
    xombe__pack49(d, bits);
}
