#include <string.h>
#include "../internal.h"
#include "dsd_state.h"

/* D-STAR decoder on the AMBE+2 set machinery and synthesis (dec/decoder.c), differing in frame logic here,
   parameters (params.c), and ds_rules */

/* xombe__ds_harm_voicing leaves vh above L alone */
static void harm_voicing(const uint8_t vb[8], float f0, int L, uint8_t vh[XOMBE_LMAX + 1])
{
    memset(vh, 0, XOMBE_LMAX + 1);
    xombe__ds_harm_voicing(vb, f0, L, vh);
}

static const struct xombe_dec_rules ds_rules = {harm_voicing, XOMBE_DS_LAMBDA_DECAY};

void xombe__dsd_logic_init(xombe_dsd_logic *s, xombe_dec_core *d)
{
    memset(s, 0, sizeof(*s));
    d->rules = &ds_rules;
}

static void voice(xombe_dec_core *d, const int b[9])
{
    xombe_params p;
    xombe_set c;
    uint8_t v1[8];
    const uint8_t *v2 = xombe__ds_vuv[b[1] & 15];
    memset(&p, 0, sizeof(p));
    p.f0 = xombe__ds_f0(b[0]);
    p.L = xombe__ds_L(b[0]);
    xombe__ds_amps(b, p.L, d->pi.L, d->pi.a, d->pi.gamma, NULL, &p.gamma, p.a);
    XOMBE_TAP(d, XOMBE_STAGE_DEQUANT, p.f0, p.a + 1, p.L);
    xombe__f3_first_half(d->pi.v2, v2, v1);
    xombe__smooth(d, p.a, p.L, p.f0);
    XOMBE_TAP(d, XOMBE_STAGE_SMOOTH, p.f0, p.a + 1, p.L);
    xombe__dec_set_from_params(&c, XOMBE_CLASS_VOICE, p.f0, p.L, p.a, v2);
    xombe__dec_play(d, &c, v1, &p, v2);
}

/* concealed and null frames: the second half is a fresh copy without the repeat flag, so a silence set there
   plays as the mute set */
static void play_pair(xombe_dec_core *d, const xombe_set *s)
{
    d->p1 = *s;
    d->p1.repeat = 1;
    d->p2 = *s;
    d->p2.repeat = 0;
}

static void conceal(xombe_dsd_logic *s, xombe_dec_core *d, int c)
{
    play_pair(d, c == XOMBE_CONCEAL_MUTE ? &d->mute : &d->played);
    s->log = c == XOMBE_CONCEAL_MUTE ? XOMBE_DSD_LOG_MUTE : XOMBE_DSD_LOG_REPEAT;
}

/* the rest of the frame's decision once c, what plays, is known */
static void decided(xombe_dec_core *d, int c)
{
    d->decision.conceal = (enum xombe_conceal)c;
    d->decision.run = d->run;
    d->decision.mute_src = c == XOMBE_CONCEAL_MUTE ? d->mute_src : -1;
    d->prev_concealed = c != XOMBE_CONCEAL_NONE;
}

void xombe__dsd_start_frame(xombe_dsd_logic *s, xombe_dec_core *d, const xombe_ds_frame *f, unsigned flags,
                            int eps_t, int eps1)
{
    xombe_decision *dc = &d->decision;
    xombe_ds_flags ef;
    int c = XOMBE_CONCEAL_NONE;
    d->frame++;
    dc->e0 = dc->e1 = -1;
    if (!f) {
        dc->cls = XOMBE_CLASS_ERASURE;
        dc->bf = dc->hb = 0;
        s->eT = s->e1 = 0;
        c = xombe__conceal_decide(d, 1, 0);
        conceal(s, d, c);
        decided(d, c);
        return;
    }
    ef.bf = (flags & XOMBE_FLAG_BF) != 0;
    ef.hb = (flags & XOMBE_FLAG_HB) != 0;
    ef.eT = eps_t < 6 ? 0 : eps_t;
    ef.e1 = eps_t < 6 ? 0 : eps1;
    dc->bf = ef.bf;
    dc->hb = ef.hb;
    s->eT = d->eps_t = ef.eT;
    s->e1 = d->eps_1 = ef.e1;
    if (f->cls == XOMBE_DS_NULL) {
        /* fixed silence set, the memory keeps its amplitudes but its class becomes silence, so the next voice
           frame starts at half amplitude. bf and hb are ignored */
        xombe_set z;
        dc->cls = XOMBE_CLASS_SILENCE;
        xombe__conceal_decide(d, 0, 0);
        xombe__dec_silence_set(&z);
        play_pair(d, &z);
        xombe__pred_store_other(&d->pi, XOMBE_CLASS_SILENCE);
        s->log = XOMBE_DSD_LOG_NULL;
    } else if (f->cls == XOMBE_DS_VOICE && ef.bf) {
        /* bad voice frame updates the memory as if good (unless hb, then nothing changes) and plays a repeat of
           the last played set, never a mute */
        xombe_set last = d->played;
        dc->cls = XOMBE_CLASS_VOICE;
        xombe__conceal_decide(d, 1, 0); /* counts as concealed */
        c = XOMBE_CONCEAL_REPEAT;
        if (ef.hb) {
            s->log = XOMBE_DSD_LOG_MUTE;
        } else {
            voice(d, f->b);
            s->log = XOMBE_DSD_LOG_VOICE;
        }
        play_pair(d, &last);
    } else if (f->cls == XOMBE_DS_VOICE) {
        dc->cls = XOMBE_CLASS_VOICE;
        c = xombe__conceal_decide(d, 0, ef.hb);
        if (c != XOMBE_CONCEAL_NONE) {
            conceal(s, d, c);
        } else {
            voice(d, f->b);
            s->log = XOMBE_DSD_LOG_VOICE;
        }
    } else {
        /* played as the AMBE+2 tone frame of its ID and level, an unknown tone index is a bad frame */
        dc->cls = XOMBE_CLASS_TONE;
        c = xombe__conceal_decide(d, ef.bf || f->tid < 0, ef.hb);
        if (c != XOMBE_CONCEAL_NONE) {
            conceal(s, d, c);
        } else {
            static const uint8_t none[8];
            xombe_set t;
            xombe__tone_set(&t, f->tid, f->AD);
            xombe__dec_play(d, &t, none, NULL, none);
            s->log = XOMBE_DSD_LOG_TONE;
        }
    }
    decided(d, c);
}

static void start_bits(xombe_dec *d, const uint8_t bits[7], unsigned flags)
{
    xombe_ds_frame f;
    uint8_t data[XOMBE_DS_DATA_BITS];
    if (flags & XOMBE_FLAG_LOST) {
        xombe__dsd_start_frame(&d->ds, &d->c, NULL, flags, d->err_t, d->err_1);
        return;
    }
    xombe__ds_unpack(bits, data);
    xombe__ds_frame_from_data(data, &f);
    xombe__dsd_start_frame(&d->ds, &d->c, &f, flags, d->err_t, d->err_1);
}

const struct xombe_dec_ops xombe__dsd_ops = {start_bits};
