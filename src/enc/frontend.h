#ifndef XOMBE_FRONTEND_H
#define XOMBE_FRONTEND_H

/* encoder front end, shared by the back ends (AMBE+2 ambe2_enc.c, D-STAR dstar/dstar_enc.c) enabled in
   e->cfg.outputs. encoder.c runs each frame in four phases:
   A  front end once per frame: xombe__fe_push, xombe__fe_analyze, and on demand xombe__fe_voice (costly,
      later calls free)
   B  voicing per back end, filling its xombe_be_ns
   C  shared background tracking (xombe__fe_background): one noise update per half, no_band_voiced taken from
      the primary back end (AMBE+2 if enabled, else D-STAR), then per back end gains, then comfort noise
   D  emit per back end with its own gains, a forced tone (encoder.c) replaces the class here, B and C saw the
      frame as unforced
   a disabled back end does no B or D work and keeps no gain memory, so AMBE+2 output does not depend on
   whether D-STAR is enabled.
   after xombe__fe_voice a back end may also read nw, fstar, win_sum2, and S1/S2 under nw, and e->win, e->wtab,
   and the kernels' window tables stay set to nw until the next frame */

#include "kernels.h"

/* phase A record of the frame being encoded */
typedef struct {
    const float *xh;    /* high-passed input at frame start */
    const int16_t *raw; /* unfiltered input at frame start, XOMBE_ENC_RAW_LEAD samples of history before it */
    /* half-frame spectra centred on samples 40 and 120, under XOMBE_BASE_WINDOW until xombe__fe_voice, then
       under the nw-point pitch-adaptive window */
    xombe_spec S1, S2;
    float P1[XOMBE_SPEC_BINS], P2[XOMBE_SPEC_BINS]; /* bin powers under the base window */
    int speech1, speech2;  /* per-half VAD decisions */
    float vad1, vad2;      /* per-half VAD metrics, log2 units */
    int run1, run2;        /* confident-background run after each half */
    int active;            /* either half speech, or within hangover */
    /* AMBE+2 tone detector result (all 0 when cfg.tones or AMBE+2 is off), then the decision after confirmation
       and hold: tone 1 with the ID and amplitude to send */
    int tone_det, tone_det_id, tone_det_ad;
    int tone, tone_id, tone_ad;
    float period; /* tracked pitch period, samples */
    /* valid when have_voice */
    int have_voice;
    int nw;         /* pitch-adaptive window length, e->win, e->wtab, and the kernels' window follow it */
    float fstar;    /* refined fundamental, cycles/sample */
    float win_sum2; /* window sum of squares */
    /* S2 under nw prepared for the voicing metric, valid when have_lf2 */
    int have_lf2;
    xombe_lf_spec lf2;
} xombe_fe_frame;

/* back end indexes into the xombe__fe_background array */
#define XOMBE_BE_AMBE2 0
#define XOMBE_BE_DSTAR 1
#define XOMBE_BE_COUNT 2

/* one back end's part in background tracking, phase B -> C -> D */
typedef struct {
    /* in, from phase B */
    int no_band_voiced[2]; /* per half: no band clearly voiced, noise update input when primary */
    uint8_t vbin[XOMBE_NS_BINS]; /* voicing at each suppression bin centre, both halves,
                                    needed only when cfg.ns */
    xombe_ns_mem *mem;     /* back end's own gain smoothing memory */
    /* out, from phase C when cfg.ns: per-half suppression gains */
    float G1[XOMBE_NS_BINS], G2[XOMBE_NS_BINS];
} xombe_be_ns;

/* returns 1 once a frame is ready to encode, after XOMBE_ENC_DELAY + e->lookahead calls */
int xombe__fe_push(xombe_enc *e, const int16_t pcm[XOMBE_FRAME]);
/* unfiltered input from the start of the frame `back` frames before the newest (0 = newest, the encoded frame is
   XOMBE_ENC_DELAY + e->lookahead back), XOMBE_ENC_RAW_LEAD samples of history before the oldest */
const int16_t *xombe__fe_raw(const xombe_enc *e, int back);
/* call once for each push that returned 1 */
void xombe__fe_analyze(xombe_enc *e, xombe_fe_frame *a);
/* sets a->have_voice, no-op when already done */
void xombe__fe_voice(xombe_enc *e, xombe_fe_frame *a);
/* S2 under nw prepared for the voicing metric, shared by the back ends. kernels keep one prepared spectrum, so
   call right before the S2 line fits */
const xombe_lf_spec *xombe__fe_s2_prepared(xombe_enc *e, xombe_fe_frame *a);
/* phase C, be[XOMBE_BE_*] NULL for disabled back ends, at least one enabled */
void xombe__fe_background(xombe_enc *e, const xombe_fe_frame *a, xombe_be_ns *be[XOMBE_BE_COUNT]);

/* AMBE+2 back end (ambe2_enc.c) */
typedef struct {
    int cls;                  /* XOMBE_CLASS_VOICE, SILENCE or TONE */
    int b[9];                 /* voice frames: b0 and b1 after phase B, all after phase D */
    float soft1[8], soft2[8]; /* soft band voicing per half, 0 unless voice */
    uint8_t v2[8];            /* second-half band voicing, 0 unless voice */
} xombe_a2_frame;

void xombe__a2_voicing(xombe_enc *e, xombe_fe_frame *a, xombe_a2_frame *f, xombe_be_ns *ns);
void xombe__a2_emit(xombe_enc *e, const xombe_fe_frame *a, xombe_a2_frame *f, const xombe_be_ns *ns,
                    uint8_t bits[7]);

#endif
