#ifndef XOMBE_DSTAR_ENC_H
#define XOMBE_DSTAR_ENC_H

/* D-STAR back end, phases B and D of enc/frontend.h (voice, NULL, tone, and Knox frames). encoder.c calls:
     xombe__ds_tone_push    every call with cfg.tones, right after the input is taken, even before the delay
                            has passed
     xombe__ds_enc_voicing  phase B, after xombe__fe_analyze and the AMBE+2 phase B: class (tone or Knox from
                            the tone pipeline, else NULL when cfg.dstar_dtx and inactive, else voice), b0, and
                            band voicing. tone and Knox frames count as unvoiced like NULL frames. voice frames
                            share the prepared S2 with AMBE+2. fills ns->no_band_voiced[2] (no a_j above 1/2),
                            ns->mem, and with cfg.ns ns->vbin (band flags without the +250 Hz OR)
     xombe__fe_background   phase C, D-STAR gains in ns->G1/G2 when cfg.ns, the frame uses G2
     xombe__ds_enc_force_tone  phase D for a forced tone, before xombe__ds_enc_frame: tone frame at the AMBE+2
                            level ad, a Knox ID a Knox frame
     xombe__ds_enc_frame    phase D, writes the voice bits
   state: e->ds, zeroed by xombe_enc_init before xombe__ds_enc_init, XOMBE_DS_ENC(e) is the typed view */

#include "../enc/frontend.h"
#include "../common/cquant.h"
#include "dstar_internal.h"

#define XOMBE_DS_S_VN 0.35f /* voicing scale on the noise-normalised metric */

/* tone pipeline (dstar_tone.c). after the call appending input frame t, the detector has run on t - 1 (segment
   ends 20 samples into t), the run rule has decided t - 3 (reads two detections ahead), and t - 4, the encoded
   frame, has its final tone (tone-in-noise reads the next frame's decision) */
#define XOMBE_DS_TONE_RING 4  /* indexed by frame mod 4 */
#define XOMBE_DS_NZ_FLOOR 50  /* noise floor is the quietest of the last N + 1 frames */

typedef struct {
    int on;      /* 0: no tone */
    int T;       /* D-STAR tone index, or -ID for a Knox tone */
    float level; /* real-valued AMBE+2 level, 127 + level term / 0.711 + TONE_ADOFF */
} xombe_ds_tone;

typedef struct {
    int n_frames;
    xombe_ds_tone detected[XOMBE_DS_TONE_RING];
    int run[XOMBE_DS_TONE_RING];                 /* run length ending at each detection, capped at 3 */
    xombe_ds_tone confirmed[XOMBE_DS_TONE_RING]; /* detections the run rule keeps, and held tones */
    xombe_ds_tone held;
    int hold_left;
    float level_db[XOMBE_DS_NZ_FLOOR + 1];       /* ring of frame levels */
    int n_levels, level_pos;
    xombe_ds_tone ndet, ntrk, tone;         /* encoded frame's detection, run-rule tone, and final tone */
    int nknox;                              /* Knox IDs in the tone table, with their pairs */
    xombe_tone_pair knox[XOMBE_TONE_KNOX_COUNT];
    float knox_f0[XOMBE_TONE_KNOX_COUNT];   /* table f0 in cycles/sample, and harmonics */
    uint8_t knox_l[XOMBE_TONE_KNOX_COUNT][2];
} xombe_ds_tones;

#define XOMBE_DS_CLS_VOICE 0
#define XOMBE_DS_CLS_NULL 1
#define XOMBE_DS_CLS_TONE 2
#define XOMBE_DS_CLS_KNOX 3

/* phase B result */
typedef struct {
    int cls;            /* XOMBE_DS_CLS_* */
    xombe_ds_tone det, trk, tone; /* detection, run rule, final */
    int voice;          /* 1 voice, 0 NULL (dtx) */
    int b0, b1, L;
    float f0;
    float m[8], mn[8];  /* band metric and its white-noise reference */
    float a[8], w[8];   /* soft voicing after the VAD gate, band weights */
    uint8_t V[8];       /* xombe__ds_vuv[b1], 0 on NULL frames */
    int no_band_voiced; /* no a_j above 1/2, 1 on NULL frames */
} xombe_ds_enc_frame;

typedef struct {
    float s_vn;
    xombe_ns_mem ns_mem;    /* D-STAR gain smoothing */
    xombe_ds_mirror mirror; /* decoder's prediction memory */
    xombe_ds_tones tones;
    xombe_ds_enc_frame fr;  /* phase B -> D */
    xombe_cq_ws ws;         /* only the last L's maps live between frames */
} xombe_ds_enc;

#define XOMBE_DS_ENC(e) ((xombe_ds_enc *)(void *)(e)->ds.bytes)

void xombe__ds_enc_init(xombe_enc *e);
/* every call with cfg.tones, even when the front end returns no frame */
void xombe__ds_tone_push(xombe_enc *e);
void xombe__ds_tone_init(xombe_ds_tones *s);
/* Knox voice frame fields for entry k at level AD, committed to the mirror */
void xombe__ds_knox_amps(xombe_ds_enc *s, int k, float AD, int b[9]);
void xombe__ds_enc_voicing(xombe_enc *e, xombe_fe_frame *a, xombe_be_ns *ns);
void xombe__ds_enc_force_tone(xombe_enc *e, int id, int ad);
void xombe__ds_enc_frame(xombe_enc *e, const xombe_fe_frame *a, const xombe_be_ns *ns, uint8_t out[7]);

/* pieces run by test/dstar_dump_enc_be.c on recorded inputs */
void xombe__ds_enc_reset(xombe_ds_enc *s);          /* s_vn default, mirror reset */
void xombe__ds_m_noise(int N, int b0, float mn[8]); /* tables/dstar_enc_tables.h */
int xombe__ds_b1_search(const float a[8], const float w[8]);
/* noise normalisation and VAD gate, returns no_band_voiced */
int xombe__ds_soft(const float m[8], const float mn[8], float s_vn, int speech2, float a[8]);
/* targets at[1..L], G NULL for no noise suppression */
void xombe__ds_enc_targets(const xombe_spec *S, float win_sum2, float f0, int L, const float *G, const uint8_t V[8],
                           float at[XOMBE_LMAX + 1]);
/* weights, closed-loop quantization against the mirror, and commit */
void xombe__ds_enc_amps(xombe_ds_enc *s, int b0, int b1, const float at[XOMBE_LMAX + 1], int b[9]);

#endif
