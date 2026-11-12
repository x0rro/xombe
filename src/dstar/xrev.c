#include <string.h>
#include "../internal.h"
#include "../common/bits.h"
#include "../common/tone_table.h"
#include "../tables/dsp_tables.h"
#include "xc_state.h"
#include "trace.h"

/* D-STAR -> AMBE+2 rate-33 parameter transcoder */

#define XR_OCTAVE_F0 (400.0f / 8000.0f) /* pitches above take the AMBE+2 index of half their f0 */
#define XR_PULSE_F0 (4217.0f / 524288.0f) /* the pulse pitch as the 2^-19 fraction XOMBE_PULSE_F0 rounds */

/* all-unvoiced D-STAR b0 -> AMBE+2 pulse code b1 18..31, DVSI's conversion where seen, else the pulse whose
   AMBE+2 -> D-STAR b0 is nearest, first on ties. 0: a plain all-unvoiced frame */
static const uint8_t pulse_code[XOMBE_DS_B0_VOICE_MAX + 1] = {
    22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22,
    22, 22, 22, 22, 22, 22, 22, 22, 22, 0,  22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22,
    22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 0,  22, 0,  23, 23, 0,  26, 26, 0,  0,
    27, 0,  28, 27, 30, 30, 30, 0,  0,  31, 31, 0,  18, 18, 20, 19, 19, 19, 20, 23, 21, 20, 0,  21, 21,
    0,  20, 20, 20, 28, 28, 27, 28, 29, 31, 27, 28, 24, 24, 24, 25, 25, 25, 19, 0,  25, 25, 25, 25};

/* D-STAR b1 -> AMBE+2 b1 by the previous D-STAR frame: DVSI's conversion where seen, else the code nearest in
   Hamming distance to the D-STAR bands each split in two. rows: previous b1 <= 0 or not a voice frame, 15,
   anything else */
static const uint8_t vuv_map[3][16] = {{16, 16, 14, 5, 16, 0, 4, 0, 12, 8, 6, 6, 10, 7, 9, 0},
                                       {16, 16, 14, 5, 16, 0, 4, 0, 12, 8, 6, 6, 9, 7, 4, 0},
                                       {16, 16, 14, 5, 16, 0, 4, 0, 12, 8, 6, 6, 8, 7, 0, 0}};

static const uint8_t std_silence[7] = XOMBE_SILENCE_FRAME;

void xombe__xc_rev_init(xombe_xc_rev *x, xombe_cq_ws *ws)
{
    memset(x, 0, sizeof(*x));
    xombe__dec_core_reset(&x->d);
    xombe__dsd_logic_init(&x->s, &x->d);
    xombe__fec_init_dstar(&x->fec);
    xombe__pred_init(&x->pi);
    x->ws = ws;
    memset(x->ws, 0, sizeof(*x->ws));
    x->ws->cap = XOMBE_CQ_CAP;
    x->prev_b1 = -1;
}

/* nearest AMBE+2 pitch index in log f0, the lower index on ties */
static int a2_b0(float f0)
{
    int b = 0;
    if (f0 > XR_OCTAVE_F0)
        f0 *= 0.5f;
    while (b < XOMBE_B0_COUNT - 1 && f0 < xombe__w0_mid[b])
        b++;
    return b;
}

/* a[1..L] at harmonics of f0 -> AMBE+2 frame b0, b1 at L2 harmonics of f2: resampled linearly in harmonic
   number, quantized against the prediction memory and stored */
static void requant(xombe_xc_rev *x, float f0, const float *a, int L, int b0, int b1, int L2, float f2,
                    const uint8_t v2[8], uint8_t d49[49])
{
    float tgt[XOMBE_LMAX + 1], w[XOMBE_LMAX + 1];
    xombe_cq_result res;
    xombe_params rec;
    int b[9], l, k;
    tgt[0] = 0.0f;
    for (l = 1; l <= L2; l++) {
        float h = (float)l * f2 / f0, fr;
        int i;
        h = h < 1.0f ? 1.0f : h > (float)L ? (float)L : h;
        i = (int)h;
        fr = h - (float)i;
        tgt[l] = (1.0f - fr) * a[i] + fr * a[i < L ? i + 1 : L];
    }
    xombe__cq_amp_weights(tgt, L2, XOMBE_CQ_QWEIGHT, XOMBE_CQ_WEIGHT_FLOOR, w);
    xombe__cq_search_ambe2(x->ws, &x->pi, tgt, L2, w, XOMBE_CQ_N_PAIRS, XOMBE_CQ_N_GAIN, &res, NULL, NULL);
    b[0] = b0;
    b[1] = b1;
    for (k = 2; k < 9; k++)
        b[k] = res.b[k];
    xombe__dequant(b, &x->pi, 0, &rec);
    xombe__pred_store_voice(&x->pi, &rec, v2);
    memset(d49, 0, 49);
    xombe__put_fields(b, d49);
    if (x->trace) {
        x->trace->tL = L2;
        memcpy(x->trace->tgt, tgt, sizeof(tgt));
        x->trace->has_b = 1;
        memcpy(x->trace->b, b, sizeof(b));
    }
}

/* the voice set (f0, a[1..L]) with fields db of its D-STAR frame -> AMBE+2 voice or pulse frame */
static int voice(xombe_xc_rev *x, const int db[9], float f0, const float *a, int L, uint8_t d49[49])
{
    static const uint8_t unvoiced[8];
    int code = db[1] == 0 ? pulse_code[db[0]] : 0, b0, b1;
    if (code) {
        /* pulse positions by code, DVSI's choice for 26..31 varies (24 or 103) */
        b0 = code <= 23 ? 103 : code <= 25 ? 30 : 24;
        requant(x, f0, a, L, b0, code, XOMBE_PULSE_L, XR_PULSE_F0, unvoiced, d49);
        return XOMBE_XC_OUT_PULSE;
    }
    b1 = vuv_map[x->prev_b1 <= 0 ? 0 : x->prev_b1 == 15 ? 1 : 2][db[1]];
    /* DVSI's first-half codebook choice by the previous output's bands */
    if (b1 == 0 && memchr(x->pi.v2, 0, 8))
        b1 = 1;
    else if (b1 == XOMBE_UNVOICED_B1 && memchr(x->pi.v2, 1, 8))
        b1 = XOMBE_UNVOICED_B1 + 1;
    b0 = a2_b0(f0);
    requant(x, f0, a, L, b0, b1, xombe__lt[b0], xombe__w0[b0], xombe__vuv[b1], d49);
    return XOMBE_XC_OUT_VOICE;
}

static int tone(xombe_xc_rev *x, int id, int ad, uint8_t d49[49])
{
    xombe__tone_frame_bits(id, ad, d49);
    xombe__pred_store_other(&x->pi, XOMBE_CLASS_TONE);
    return XOMBE_XC_OUT_TONE;
}

static int silence(xombe_xc_rev *x, uint8_t d49[49])
{
    xombe__unpack49(std_silence, d49);
    xombe__pred_store_other(&x->pi, XOMBE_CLASS_SILENCE);
    return XOMBE_XC_OUT_SILENCE;
}

int xombe__xc_rev_frame(xombe_xc_rev *x, const uint8_t in[9], unsigned flags, uint8_t out[9])
{
    xombe_dsd_logic *s = &x->s;
    const xombe_set *p = &x->d.p2;
    xombe_xr_trace *tr = x->trace;
    xombe_ds_frame f;
    xombe_fec_stats st;
    uint8_t d49[49], bits[7];
    int lost = (flags & XOMBE_FLAG_LOST) != 0, ic, o;
    if (lost) {
        xombe__dsd_start_frame(s, &x->d, NULL, 0, 0, 0);
    } else {
        uint8_t data[XOMBE_DS_DATA_BITS];
        unsigned fl = xombe_fec_decode(&x->fec, in, bits, &st);
        xombe__ds_unpack(bits, data);
        xombe__ds_frame_from_data(data, &f);
        xombe__dsd_start_frame(s, &x->d, &f, fl, st.eps_t, st.eps1);
        x->d.decision.e0 = st.eps0;
        x->d.decision.e1 = st.eps1;
    }
    xombe__dec_skip_synth(&x->d);
    if (s->log == XOMBE_DSD_LOG_MUTE)
        ic = XOMBE_XR_IN_MUTE;
    else if (s->log == XOMBE_DSD_LOG_NULL)
        ic = XOMBE_XR_IN_NULL;
    else if (x->d.decision.conceal != XOMBE_CONCEAL_NONE)
        ic = XOMBE_XR_IN_REPEAT;
    else
        ic = s->log == XOMBE_DSD_LOG_TONE ? XOMBE_XR_IN_TONE : XOMBE_XR_IN_VOICE;
    if (tr) {
        memset(tr, 0, sizeof(*tr));
        tr->in_cls = ic;
        tr->scls = p->cls;
        tr->tid = tr->AD = -1;
        if (!lost) {
            tr->has_db = 1;
            memcpy(tr->db, f.b, sizeof(f.b));
        }
        if (p->cls == XOMBE_CLASS_VOICE || p->cls == XOMBE_CLASS_SILENCE) {
            tr->has_par = 1;
            tr->f0 = p->f0;
            tr->L = p->L;
            memcpy(tr->a, p->a, sizeof(tr->a));
        }
    }
    if (ic == XOMBE_XR_IN_VOICE) {
        x->have_last = 1;
        memcpy(x->last_db, f.b, sizeof(f.b));
        o = voice(x, f.b, p->f0, p->a, p->L, d49);
        x->prev_b1 = f.b[1];
    } else if (ic == XOMBE_XR_IN_TONE) {
        x->have_tone = 1;
        x->tone_id = f.tid;
        x->tone_ad = f.AD;
        x->prev_b1 = -1;
        o = tone(x, f.tid, f.AD, d49);
        if (tr) {
            tr->tid = f.tid;
            tr->AD = f.AD;
        }
    } else if (ic == XOMBE_XR_IN_REPEAT && p->cls == XOMBE_CLASS_VOICE && x->have_last) {
        o = voice(x, x->last_db, p->f0, p->a, p->L, d49);
    } else if (ic == XOMBE_XR_IN_REPEAT && p->cls == XOMBE_CLASS_TONE && x->have_tone) {
        o = tone(x, x->tone_id, x->tone_ad, d49);
        if (tr)
            tr->tid = x->tone_id;
    } else {
        if (ic != XOMBE_XR_IN_REPEAT)
            x->prev_b1 = -1;
        o = silence(x, d49);
    }
    xombe__pack49(d49, bits);
    xombe__fec_encode_dmr(bits, out);
    x->out = o;
    if (tr)
        memcpy(tr->d49, d49, 49);
    return 1;
}
