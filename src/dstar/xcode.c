#include <math.h>
#include <string.h>
#include "../internal.h"
#include "../common/tone_table.h"
#include "xc_state.h"
#include "trace.h"

/* AMBE+2 rate-33 -> D-STAR parameter transcoder */

#define XC_VTHR 0.5f      /* band voiced if its energy-weighted voiced fraction is above this */
#define XC_KNOX_REST 4.0f /* transcoder's other Knox harmonics this far below the tones, log2 */
#define XC_PULSE_FIRST 18 /* AMBE+2 pulse frames: b1 18..31 */
#define XC_PULSE_LAST 31

/* AMBE+2 pulse frame b1 18..31 -> D-STAR b0 as DVSI converts it, 30 -> 80 inferred as the one unused value of
   the 64..116 step-4 grid */
static const uint8_t pulse_b0[XC_PULSE_LAST - XC_PULSE_FIRST + 1] = {88, 92, 96, 100, 64, 68, 112, 116,
                                                                     72, 76, 104, 108, 80, 84};

void xombe__xc_fwd_init(xombe_xc_fwd *x, xombe_cq_ws *ws)
{
    static const uint8_t null_frame[XOMBE_DS_FRAME_BYTES] = XOMBE_DS_NULL_FRAME;
    memset(x, 0, sizeof(*x));
    x->ws = ws;
    memset(x->ws, 0, sizeof(*x->ws));
    x->ws->cap = XOMBE_CQ_CAP;
    xombe__dec_core_reset(&x->d);
    xombe__fec_init_dmr(&x->fec);
    xombe__ds_mirror_init(&x->mirror);
    memcpy(x->last, null_frame, XOMBE_DS_FRAME_BYTES);
}

/* energy-weighted voiced fraction per 1 kHz band of harmonics 1..L, harmonic voicing with AMBE+2's +250 Hz OR */
static void band_soft(const float *a, int L, float f0, const uint8_t V[8], float s[4])
{
    float num[4] = {0.0f, 0.0f, 0.0f, 0.0f}, den[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    float f = f0 * XOMBE_FS;
    int l, j;
    for (l = 1; l <= L; l++) {
        float fr = (float)l * f, w;
        int j1, j2;
        j = (int)(fr / 1000.0f);
        if (j > 3)
            continue;
        j1 = (int)(fr / 500.0f);
        j2 = (int)((fr + 250.0f) / 500.0f);
        w = exp2f(2.0f * a[l]);
        num[j] += (V[j1 > 7 ? 7 : j1] || V[j2 > 7 ? 7 : j2]) ? w : 0.0f;
        den[j] += w;
    }
    for (j = 0; j < 4; j++)
        s[j] = den[j] > 0.0f ? num[j] / den[j] : 0.0f;
}

/* the code nearest the per-band decisions s > XC_VTHR (bit 3 = 0..1 kHz), never XOMBE_DS_B1_SKIP: least sum of
   |s - XC_VTHR| over the bands it flips, first on ties */
static int b1_from_soft(const float s[4])
{
    float best = 0.0f;
    int c, j, want = 0, bc = -1;
    for (j = 0; j < 4; j++)
        want |= (s[j] > XC_VTHR) << (3 - j);
    for (c = 0; c < 16; c++) {
        float cost = 0.0f;
        if (c == XOMBE_DS_B1_SKIP)
            continue;
        for (j = 0; j < 4; j++)
            if (((c ^ want) >> (3 - j)) & 1)
                cost += fabsf(s[j] - XC_VTHR);
        if (bc < 0 || cost < best) {
            best = cost;
            bc = c;
        }
    }
    return bc;
}

/* a[1..L] at harmonics of f0 -> D-STAR voice frame: resampled linearly in harmonic number to the harmonics of
   b0 (nearest pitch of f0 if b0 < 0), quantized against the mirror and committed */
static void requant(xombe_xc_fwd *x, float f0, const float *a, int L, int b1, int b0, int b[9])
{
    float tgt[XOMBE_LMAX + 1], w[XOMBE_LMAX + 1], pt[XOMBE_LMAX + 1], f0d;
    xombe_cq_result res;
    int Ld, l, k;
    if (b0 < 0)
        b0 = xombe__ds_b0_from_f0(f0);
    f0d = xombe__ds_f0(b0);
    Ld = xombe__ds_L(b0);
    tgt[0] = 0.0f;
    for (l = 1; l <= Ld; l++) {
        float h = (float)l * f0d / f0, fr;
        int i;
        h = h < 1.0f ? 1.0f : h > (float)L ? (float)L : h;
        i = (int)h;
        fr = h - (float)i;
        tgt[l] = (1.0f - fr) * a[i] + fr * a[i < L ? i + 1 : L];
    }
    xombe__cq_amp_weights(tgt, Ld, XOMBE_CQ_QWEIGHT, XOMBE_CQ_WEIGHT_FLOOR, w);
    xombe__ds_mirror_predict(&x->mirror, Ld, pt);
    xombe__cq_quantize(x->ws, &xombe__amp_dstar, tgt, Ld, pt, x->mirror.G, w, 1, XOMBE_CQ_N_PAIRS, XOMBE_CQ_N_GAIN,
                       XOMBE_CQ_ITERS, XOMBE_CQ_CLIP_PASSES, &res);
    b[0] = b0;
    b[1] = b1;
    for (k = 2; k < 9; k++)
        b[k] = res.b[k];
    xombe__ds_mirror_commit(&x->mirror, b, NULL);
    if (x->trace) {
        x->trace->tL = Ld;
        x->trace->tb0 = b0;
        memcpy(x->trace->tgt, tgt, sizeof(tgt));
    }
}

/* Knox tone as a voice frame at the table f0: the two tone harmonics at the tone level
   log2 min(16384 * 10^(0.03555 (AD - 127)), 2^14), the others XC_KNOX_REST below, all bands voiced.
   returns 0 if the tone table has no pair for id */
static int knox(xombe_xc_fwd *x, int id, int AD, int b[9])
{
    float f0hz, f0, lm, a[XOMBE_LMAX + 1];
    int l1, l2, L, l;
    xombe__tone_harmonics(id, &f0hz, &l1, &l2);
    if (l1 <= 0 || l2 <= 0)
        return 0;
    f0 = f0hz / XOMBE_FS;
    lm = XOMBE_TONE_AD_LOG10 * (float)(AD - 127) * XOMBE_LOG2_10;
    lm = (lm > 0.0f ? 0.0f : lm) + 14.0f;
    L = xombe__L_from_f0(f0);
    for (l = 1; l <= L; l++)
        a[l] = lm - XC_KNOX_REST;
    if (l1 <= L)
        a[l1] = lm;
    if (l2 <= L)
        a[l2] = lm;
    requant(x, f0, a, L, XOMBE_DS_B1_VOICED, -1, b);
    return 1;
}

/* silence and pulse pitches 250.60 Hz and 64.346 Hz as the 2^-19 fractions consts.h rounds, the resampling
   positions l * f0d / f0 otherwise move by up to 1.2e-5 relative, enough to change quantizer picks */
static float exact_f0(float f0)
{
    if (f0 == XOMBE_SILENCE_F0)
        return 16423.0f / 524288.0f;
    if (f0 == XOMBE_PULSE_F0)
        return 4217.0f / 524288.0f;
    return f0;
}

int xombe__xc_fwd_frame(xombe_xc_fwd *x, const uint8_t in[9], unsigned flags, uint8_t out[9])
{
    xombe_dec_core *d = &x->d;
    const xombe_decision *dc = &d->decision;
    xombe_xc_trace *tr = x->trace;
    const xombe_set *s = &d->p2;
    uint8_t bits[7];
    int b[9], id = -1, ad = -1, voice_b = 0, cls;
    memset(bits, 0, sizeof(bits));
    xombe__dec_start_air(d, &x->fec, in, flags, bits);
    xombe__dec_skip_synth(d);
    cls = s->cls;
    if (cls == XOMBE_CLASS_TONE && dc->conceal == XOMBE_CONCEAL_NONE) {
        uint8_t d49[49];
        xombe__unpack49(bits, d49);
        xombe__tone_frame_parse(d49, &id, &ad);
    }
    if (tr) {
        uint8_t d49[49];
        memset(tr, 0, sizeof(*tr));
        tr->cls = dc->cls;
        tr->conceal = dc->conceal;
        tr->cls_out = cls;
        tr->tid = tr->AD = tr->tb0 = -1;
        if (dc->conceal == XOMBE_CONCEAL_NONE) {
            xombe__unpack49(bits, d49);
            if (cls != XOMBE_CLASS_TONE) {
                tr->has_rb = tr->has_v1 = 1;
                xombe__get_fields(d49, tr->rb);
                memcpy(tr->V1, d->p1.vb, 8);
                tr->G = d->pi.gamma;
            }
        }
    }
    if (cls == XOMBE_CLASS_TONE && dc->conceal != XOMBE_CONCEAL_NONE) {
        /* the AMBE+2 decoder replays the tone set: resend the last frame, which the D-STAR decoder decodes anew */
        if (x->have_last_b)
            xombe__ds_mirror_commit(&x->mirror, x->last_b, NULL);
        memcpy(out, x->last, XOMBE_DS_FRAME_BYTES);
        x->out = XOMBE_XC_OUT_RESEND;
        if (tr && x->have_last_b) {
            tr->has_b = 1;
            memcpy(tr->b, x->last_b, sizeof(tr->b));
        }
        return 1;
    }
    if (cls == XOMBE_CLASS_TONE) {
        int T = xombe__ds_tone_code(id);
        if (tr) {
            tr->tid = id;
            tr->AD = ad;
        }
        if (T >= 0) {
            xombe__ds_tone_frame(T, xombe__ds_tone_q(XOMBE_DS_AD_LEVEL(ad)), out);
            x->out = XOMBE_XC_OUT_TONE;
        } else if (id >= XOMBE_TONE_KNOX_FIRST && id < XOMBE_TONE_KNOX_FIRST + XOMBE_TONE_KNOX_COUNT &&
                   knox(x, id, ad, b)) {
            voice_b = 1;
            x->out = XOMBE_XC_OUT_KNOX;
        } else {
            memcpy(out, x->last, XOMBE_DS_FRAME_BYTES);
            x->out = XOMBE_XC_OUT_RESEND;
            x->have_last_b = 0;
            return 1;
        }
    } else {
        /* voice, or silence as an all-unvoiced voice frame from its comfort-noise amplitudes */
        float soft[4] = {0.0f, 0.0f, 0.0f, 0.0f};
        int pulse = 0, b0 = -1;
        if (cls == XOMBE_CLASS_VOICE && dc->conceal == XOMBE_CONCEAL_NONE) {
            uint8_t d49[49];
            int rb[9];
            xombe__unpack49(bits, d49);
            xombe__get_fields(d49, rb);
            pulse = rb[1] >= XC_PULSE_FIRST && rb[1] <= XC_PULSE_LAST;
            if (pulse)
                b0 = pulse_b0[rb[1] - XC_PULSE_FIRST];
        }
        if (cls == XOMBE_CLASS_VOICE && !pulse)
            band_soft(s->a, s->L, s->f0, s->vb, soft);
        if (tr) {
            tr->has_par = 1;
            tr->f0 = s->f0;
            tr->L = s->L;
            memcpy(tr->a, s->a, sizeof(tr->a));
            memcpy(tr->V2, s->vb, 8);
            tr->pulse = pulse;
            tr->has_soft = cls == XOMBE_CLASS_VOICE && !pulse;
            memcpy(tr->soft, soft, sizeof(soft));
        }
        requant(x, exact_f0(s->f0), s->a, s->L, b1_from_soft(soft), b0, b);
        voice_b = 1;
        x->out = XOMBE_XC_OUT_VOICE;
    }
    if (voice_b)
        xombe__ds_encode_fields(b, out);
    memcpy(x->last, out, XOMBE_DS_FRAME_BYTES);
    x->have_last_b = voice_b;
    if (voice_b)
        memcpy(x->last_b, b, sizeof(b));
    if (tr && voice_b) {
        tr->has_b = 1;
        memcpy(tr->b, b, sizeof(b));
    }
    return 1;
}
