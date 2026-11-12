#include <math.h>
#include <string.h>
#include "frontend.h"
#include "../tables/dsp_tables.h"
#include "../common/prof.h"

/* each call appends one frame to the XOMBE_ENC_BUF-sample history. the encoded frame starts at frame_base(e)
   <= FRAME_BASE, pitch analysis runs on the frame before the newest (so its one-frame look-ahead is complete) */
#define FRAME XOMBE_FRAME
#define FRAME_BASE (XOMBE_ENC_BUF - (XOMBE_ENC_DELAY + 1 + XOMBE_ENC_LOOKAHEAD) * FRAME)
#define PITCH_FRAME_BASE (XOMBE_ENC_BUF - 2 * FRAME)
#define HALF1_CENTRE 40
#define HALF2_CENTRE 120
#define BASE_WINDOW XOMBE_BASE_WINDOW

XOMBE_STATIC_ASSERT(FRAME_BASE + HALF1_CENTRE - (XOMBE_MAX_WINDOW - 1) / 2 >= 0, analysis_window_in_buffer);
XOMBE_STATIC_ASSERT(FRAME_BASE + HALF2_CENTRE + (XOMBE_MAX_WINDOW - 1) / 2 < XOMBE_ENC_BUF, analysis_window_end);
XOMBE_STATIC_ASSERT(FRAME_BASE == XOMBE_ENC_BUF - XOMBE_ENC_RAW, raw_history_starts_at_frame);
XOMBE_STATIC_ASSERT(HALF1_CENTRE - (BASE_WINDOW - 1) / 2 >= -XOMBE_ENC_RAW_LEAD, vad_flat_window_in_raw_lead);
XOMBE_STATIC_ASSERT(XOMBE_TONE_SEG <= XOMBE_ENC_RAW, tone_segment_in_buffer);
XOMBE_STATIC_ASSERT(PITCH_FRAME_BASE + HALF2_CENTRE - XOMBE_PITCH_HALF - XOMBE_PITCH_LPF_HALF >= 0,
                    pitch_window_in_buffer);
XOMBE_STATIC_ASSERT(PITCH_FRAME_BASE + HALF2_CENTRE + XOMBE_PITCH_HALF + XOMBE_PITCH_LPF_HALF < XOMBE_ENC_BUF,
                    pitch_window_end);
/* one twiddle table serves both transform sizes */
XOMBE_STATIC_ASSERT(XOMBE_TW_N % XOMBE_FFT == 0, twiddles_cover_analysis_fft);
XOMBE_STATIC_ASSERT(XOMBE_TW_N % XOMBE_TONE_FFT == 0, twiddles_cover_tone_fft);

#define HANGOVER_FRAMES 2

/* returns 1 if active, counting hangover */
static int detect_activity(xombe_enc *e, xombe_fe_frame *f)
{
    float tap[2];
    const int h = (BASE_WINDOW - 1) / 2;
    xombe__kern_spectrum(f->xh + HALF1_CENTRE, xombe__base_win, BASE_WINDOW, &f->S1);
    xombe__kern_spectrum(f->xh + HALF2_CENTRE, xombe__base_win, BASE_WINDOW, &f->S2);
    xombe__kern_bin_power(&f->S1, e->base_win_sum2, f->P1);
    xombe__kern_bin_power(&f->S2, e->base_win_sum2, f->P2);

    xombe__vad_half(e, f->P1, xombe__vad_flat(f->raw + HALF1_CENTRE - h, BASE_WINDOW), &f->speech1, &f->vad1);
    e->confident_run = (!f->speech1 && f->vad1 < XOMBE_CONFIDENT) ? e->confident_run + 1 : 0;
    f->run1 = e->confident_run;
    xombe__vad_half(e, f->P2, xombe__vad_flat(f->raw + HALF2_CENTRE - h, BASE_WINDOW), &f->speech2, &f->vad2);
    e->confident_run = (!f->speech2 && f->vad2 < XOMBE_CONFIDENT) ? e->confident_run + 1 : 0;
    f->run2 = e->confident_run;
    tap[0] = (float)f->speech1;
    tap[1] = f->vad1;
    XOMBE_TAP(e, XOMBE_STAGE_VAD, 0.0f, tap, 2);
    tap[0] = (float)f->speech2;
    tap[1] = f->vad2;
    XOMBE_TAP(e, XOMBE_STAGE_VAD, 0.0f, tap, 2);
    if (f->speech1 || f->speech2) {
        e->hangover = HANGOVER_FRAMES;
        return 1;
    }
    if (e->hangover > 0) {
        e->hangover--;
        return 1;
    }
    return 0;
}

static int frame_base(const xombe_enc *e)
{
    return XOMBE_ENC_BUF - (XOMBE_ENC_DELAY + 1 + e->lookahead) * FRAME;
}

/* pitch analyses sit in a ring of XOMBE_ENC_DELAY slots, the tracked frame is the oldest, whose slot the next
   call overwrites. its period goes to e->period[0], the encoded frame is e->lookahead older */
int xombe__fe_push(xombe_enc *e, const int16_t pcm[XOMBE_FRAME])
{
    int k = e->calls;
    int i;
    memmove(e->raw, e->raw + FRAME, sizeof(e->raw[0]) * (XOMBE_ENC_RAW_LEAD + XOMBE_ENC_RAW - FRAME));
    memmove(e->xh, e->xh + FRAME, sizeof(float) * (XOMBE_ENC_BUF - FRAME));
    for (i = 0; i < FRAME; i++) {
        e->raw[XOMBE_ENC_RAW_LEAD + XOMBE_ENC_RAW - FRAME + i] = pcm[i];
        e->xh[XOMBE_ENC_BUF - FRAME + i] = xombe__hpf_step(&e->hpf, (float)pcm[i]);
    }
    if (k < XOMBE_ENC_DELAY + e->lookahead)
        e->calls++;
    if (k >= 1) {
        xombe__pitch_errors(e, e->xh + PITCH_FRAME_BASE + HALF2_CENTRE, e->pitch_slot);
        e->pitch_slot = (e->pitch_slot + 1) % XOMBE_ENC_DELAY;
    }
    if (k >= XOMBE_ENC_DELAY) {
        int cur = e->pitch_slot;
#if XOMBE_ENC_LOOKAHEAD
        for (i = e->lookahead; i > 0; i--)
            e->period[i] = e->period[i - 1];
#endif
        e->period[0] = xombe__pitch_track(e, cur, (cur + 1) % XOMBE_ENC_DELAY, (cur + 2) % XOMBE_ENC_DELAY);
    }
    return k >= XOMBE_ENC_DELAY + e->lookahead;
}

const int16_t *xombe__fe_raw(const xombe_enc *e, int back)
{
    return e->raw + XOMBE_ENC_RAW_LEAD + XOMBE_ENC_RAW - (back + 1) * FRAME;
}

void xombe__fe_analyze(xombe_enc *e, xombe_fe_frame *a)
{
    a->xh = e->xh + frame_base(e);
    a->raw = xombe__fe_raw(e, XOMBE_ENC_DELAY + e->lookahead);
    a->active = detect_activity(e, a);
    a->tone_det = a->tone_det_id = a->tone_det_ad = 0;
    a->tone = a->tone_id = a->tone_ad = 0;
    if (e->cfg.tones && (e->cfg.outputs & XOMBE_ENC_DMR)) {
        XOMBE_PF(10, a->tone_det = xombe__tone_detect(e, a->raw, &a->tone_det_id, &a->tone_det_ad));
        a->tone = xombe__tone_decide(e, a->tone_det, a->tone_det_id, a->tone_det_ad, &a->tone_id, &a->tone_ad);
    }
    a->period = e->period[e->lookahead];
    a->have_voice = 0;
    a->have_lf2 = 0;
    a->nw = BASE_WINDOW;
    a->fstar = 0.0f;
    a->win_sum2 = e->base_win_sum2;
}

void xombe__fe_voice(xombe_enc *e, xombe_fe_frame *a)
{
    if (a->have_voice)
        return;
    a->nw = xombe__pitch_window_len(a->period);
    xombe__enc_set_window(e, a->nw);
    if (a->nw != BASE_WINDOW) {
        xombe__kern_spectrum(a->xh + HALF1_CENTRE, e->win, a->nw, &a->S1);
        xombe__kern_spectrum(a->xh + HALF2_CENTRE, e->win, a->nw, &a->S2);
    }
    a->fstar = 1.0f / xombe__pitch_refine(e, &a->S2, a->period);
    a->win_sum2 = e->win_sum2;
    a->have_voice = 1;
    a->have_lf2 = 0;
}

const xombe_lf_spec *xombe__fe_s2_prepared(xombe_enc *e, xombe_fe_frame *a)
{
    if (!a->have_lf2) {
        xombe__voicing_prepare(e, &a->S2, &a->lf2);
        a->have_lf2 = 1;
    }
    return &a->lf2;
}

void xombe__fe_background(xombe_enc *e, const xombe_fe_frame *a, xombe_be_ns *be[XOMBE_BE_COUNT])
{
    if (e->cfg.ns) {
        float X1[XOMBE_NS_BINS], X2[XOMBE_NS_BINS];
        const xombe_be_ns *primary = be[XOMBE_BE_AMBE2] ? be[XOMBE_BE_AMBE2] : be[XOMBE_BE_DSTAR];
        int i;
        xombe__ns_bins(a->P1, X1);
        xombe__ns_bins(a->P2, X2);
        xombe__ns_update(e, X1, a->speech1, a->run1, primary->no_band_voiced[0]);
        for (i = 0; i < XOMBE_BE_COUNT; i++)
            if (be[i])
                xombe__ns_gain(e, be[i]->mem, X1, a->speech1, a->vad1, be[i]->vbin, be[i]->G1);
        xombe__ns_update(e, X2, a->speech2, a->run2, primary->no_band_voiced[1]);
        for (i = 0; i < XOMBE_BE_COUNT; i++)
            if (be[i])
                xombe__ns_gain(e, be[i]->mem, X2, a->speech2, a->vad2, be[i]->vbin, be[i]->G2);
    }
    xombe__cn_half(e, a->P1, !a->speech1, a->run1);
    xombe__cn_half(e, a->P2, !a->speech2, a->run2);
}
