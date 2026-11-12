#ifndef XOMBE_STATE_H
#define XOMBE_STATE_H

/* decoder and encoder state behind the opaque handles of xombe.h */

#include <stdint.h>
#include "xombe.h"

#define XOMBE_LMAX 56
#define XOMBE_TAIL 84

#define XOMBE_WT_SIZE 256 /* voiced wavetable length, a power of 2 */

typedef struct {
    float re, im;
} xombe_cpx;

/* per-stage taps for tests (test/trace_dump.c), each passes n values and their f0 */
enum xombe_stage {
    XOMBE_STAGE_DEQUANT,  /* a[1..L] after dequantization */
    XOMBE_STAGE_SMOOTH,   /* a[1..L] after smoothing */
    XOMBE_STAGE_SET,      /* a[1..L] of the set a half-frame plays (silence and voice sets) */
    XOMBE_STAGE_ENHANCED, /* linear M[1..L] after enhancement (silence sets: 2^a) */
    XOMBE_STAGE_PHASE,    /* phi[1..L] (voice sets) */
    XOMBE_STAGE_SHAPED,   /* M[1..L] as synthesized, after the tilt (tone sets: tone amplitudes) */
    XOMBE_STAGE_SYNTH,    /* a half-frame's N samples before the output high-pass */
    XOMBE_STAGE_TARGET,   /* encoder amplitude targets a[1..L] as quantized */
    XOMBE_STAGE_VAD,      /* encoder VAD per half: decision (1 = speech) and metric v in log2 units, f0 = 0 */
    XOMBE_STAGE_COUNT
};
typedef void (*xombe_tap_fn)(void *ctx, enum xombe_stage stage, float f0, const float *v, int n);
#define XOMBE_TAP(o, stage, f0, v, n)                                                                          \
    do {                                                                                                       \
        if ((o)->tap)                                                                                          \
            (o)->tap((o)->tap_ctx, (stage), (f0), (v), (n));                                                   \
    } while (0)

typedef struct {
    int cls;
    int L;
    float f0;
    float gamma;
    float a[XOMBE_LMAX + 1];
    uint8_t v2[8];
} xombe_pred;

typedef struct {
    int cls;
    int repeat;
    int prepared;
    int L;
    float f0;
    float a[XOMBE_LMAX + 1];
    float e[XOMBE_LMAX + 1]; /* log2 amplitudes before the tilt: enhanced a (voice), a (silence), or tone
                                harmonics (tone, XOMBE_AMP_OFF elsewhere) */
    float phi[XOMBE_LMAX + 1];
    uint8_t vb[8];
} xombe_set;

typedef struct {
    float x1, s1, s2, y1, y2;
    float gain;
} xombe_hpf;

struct xombe_dec_rules;
struct xombe_dec_ops;

/* frame logic up to the sets a frame plays, all the transcoders keep of a decoder */
typedef struct {
    const struct xombe_dec_rules *rules; /* AMBE+2's, or D-STAR's (xombe__dsd_logic_init) */
    xombe_pred pi;
    float lambda;
    float avg_f0; /* synthesis reads and updates it, here it fills the core's padding */
    int eps_t, eps_1;
    int prev_concealed;
    unsigned hist;
    int run;
    int64_t frame;
    int64_t mute_src;
    xombe_decision decision;
    xombe_set p1, p2, played, mute;
    xombe_tap_fn tap;
    void *tap_ctx;
} xombe_dec_core;

#if XOMBE_HAVE_DSTAR
/* D-STAR frame logic state besides the core (dstar/decoder.c) */
typedef struct {
    int log;    /* XOMBE_DSD_LOG_* of the last frame */
    int eT, e1; /* error rule smoothing inputs for the last frame */
} xombe_dsd_logic;
#endif

struct xombe_dec {
    xombe_dec_core c;
    const struct xombe_dec_ops *ops; /* the codec's frame logic */
    int second_half_next; /* 1 between the two half calls of a frame */
    int err_t, err_1;     /* error counts of the next *_bits frame, xombe_dec_set_errors or the FEC */
    float tail[XOMBE_TAIL];
    float pitch_phase, epoch;
    float wprev[XOMBE_WT_SIZE];
    float prev_f0;
    int prev_voiced;
    float uv_overlap[XOMBE_TAIL];
    int16_t noise_hist[XOMBE_TAIL];
    uint16_t lcg;
    uint32_t prng;
    int noise_mode;
    xombe_hpf hpf;
    xombe_fec fec; /* for the frame calls */
#if XOMBE_HAVE_DSTAR
    xombe_dsd_logic ds;
#endif
};

#if XOMBE_HAVE_ENCODER
#define XOMBE_FRAME 160
#define XOMBE_ENC_BUF 1024 /* input history in samples, newest last */
/* frames of look-ahead beyond XOMBE_ENC_DELAY for the D-STAR tone rules (dstar/dstar_enc.h), with the D-STAR
   output on the encoder codes a frame this many frames older */
#if XOMBE_HAVE_DSTAR
#define XOMBE_ENC_LOOKAHEAD 1
#else
#define XOMBE_ENC_LOOKAHEAD 0
#endif
/* unfiltered history from the frame being encoded through the newest */
#define XOMBE_ENC_RAW ((XOMBE_ENC_DELAY + 1 + XOMBE_ENC_LOOKAHEAD) * XOMBE_FRAME)
/* unfiltered samples kept before that frame: the VAD's digital-silence test reads 59, D-STAR's tone-in-noise
   window 80 */
#if XOMBE_HAVE_DSTAR
#define XOMBE_ENC_RAW_LEAD 80
#else
#define XOMBE_ENC_RAW_LEAD 64
#endif
#define XOMBE_FFT 1024
#define XOMBE_SPEC_BINS (XOMBE_FFT / 2 + 1)
#define XOMBE_TONE_FFT 512 /* uses the analysis FFT's twiddles */
#define XOMBE_TONE_IDS 164
#define XOMBE_TONE_SEG 240
/* analysis window widens to XOMBE_MAX_WINDOW for long periods, its transform tabulated to XOMBE_WTAB_BINS bins:
   fits read at most half a harmonic spacing away, the widest about 52.2 bins (xombe__w0[0] with +2% jitter) */
#define XOMBE_BASE_WINDOW 199
#define XOMBE_MAX_WINDOW 321
#define XOMBE_WTAB_OVERSAMPLE 16
#define XOMBE_WTAB_BINS 27
#define XOMBE_WTAB_LEN (XOMBE_WTAB_BINS * XOMBE_WTAB_OVERSAMPLE + 1)
/* pitch grid: 21, 21.5, ..., 123 samples */
#define XOMBE_PITCH_LPF_HALF 10
#define XOMBE_PITCH_HALF 150
#define XOMBE_PITCH_WIN (2 * XOMBE_PITCH_HALF + 1)
#define XOMBE_PITCH_GRID 205
#define XOMBE_VAD_HIST 300 /* half-frames of band levels for the stationarity rule */
#define XOMBE_NS_BINS 40

#include "enc/kern_state.h"

#if XOMBE_HAVE_DSTAR
/* D-STAR back end state bytes (dstar/dstar_enc.h checks its struct fits): XOMBE_ENC_DSTAR_SIZE less what the
   front end's longer input history (XOMBE_ENC_LOOKAHEAD) takes, so growing one grows the other */
#define XOMBE_ENC_LOOKAHEAD_STATE 512
#define XOMBE_DS_ENC_STATE (XOMBE_ENC_DSTAR_SIZE - XOMBE_ENC_LOOKAHEAD_STATE)
#endif

/* a back end's noise-suppression gain smoothing (xombe__ns_gain): last half's gains before the floor */
typedef struct {
    float prev_gain[XOMBE_NS_BINS];
    int have_prev;
} xombe_ns_mem;

struct xombe_enc {
    xombe_enc_cfg cfg;
    float ns_strength, voicing_scale; /* cfg's, as fractions */
    int calls;      /* counted up to the delay */
    int lookahead;  /* XOMBE_ENC_LOOKAHEAD with the D-STAR output, else 0 */
    int pitch_slot; /* perr/pacf slot the next pitch analysis writes */
    int16_t raw[XOMBE_ENC_RAW_LEAD + XOMBE_ENC_RAW];
    float xh[XOMBE_ENC_BUF];
    xombe_hpf hpf;
    float perr[XOMBE_ENC_DELAY][XOMBE_PITCH_GRID];
    /* pitch autocorrelation per slot, float lags or int16 lags (padded by one) scaled by pacf_scale */
    union {
        float f[XOMBE_PITCH_WIN];
        int16_t q[XOMBE_PITCH_WIN + 1];
    } pacf[XOMBE_ENC_DELAY];
    float pacf_scale[XOMBE_ENC_DELAY];
    float pe0[XOMBE_ENC_DELAY];
    int psilent[XOMBE_ENC_DELAY];
    float p_prev, e_prev1, e_prev2;
    float period[XOMBE_ENC_LOOKAHEAD + 1]; /* tracked periods, newest first, the frame encoded takes [lookahead] */
    float base_win_sum2;
    float win[XOMBE_MAX_WINDOW];
    int win_len;
    float win_sum2;
    float wtab[XOMBE_WTAB_LEN];
    float vad_nl[8], vad_smooth[8];
    int16_t vad_hist[XOMBE_VAD_HIST][8];  /* log2 of the smoothed band powers, Q10 (vad.c) */
    int16_t vad_hist_fl[XOMBE_VAD_HIST];  /* log2 of each half's floor, Q10 */
    int64_t vad_hist_sum[8];              /* history's band powers, dequantized to Q20 and summed */
    int16_t vad_above[8];                 /* halves in a row with band SNR > 0 dB, up to XOMBE_VAD_HIST */
    int vad_hist_n, vad_hist_pos, vad_ready, hangover;
    int confident_run;
    float ns_noise[XOMBE_NS_BINS], ns_floor; /* shared background estimate (xombe__ns_update) */
    int ns_ready;
    xombe_ns_mem ns_mem;                     /* AMBE+2 back end's gain smoothing */
    float cn_power[XOMBE_SPEC_BINS];
    int cn_ready, cn_updates;
    xombe_cpx tone_spec[XOMBE_TONE_FFT / 2 + 1]; /* scratch, not carried between frames */
    int tone_prev_id, tone_confirmed, tone_hold_used, tone_last_id, tone_last_ad;
    uint8_t tone_req[XOMBE_ENC_DELAY + XOMBE_ENC_LOOKAHEAD + 1][2]; /* forced tone ID (0: none) and AD of each
                                                                      input frame in the delay, newest first */
    xombe_pred pi; /* AMBE+2 back end's prediction memory */
    xombe_tap_fn tap;
    void *tap_ctx;
    XOMBE_KERN_STATE
#if !XOMBE_HAVE_DSTAR
    /* AMBE+2 closed-loop quantizer workspace, opaque here (enc/ambe2_enc.c: A2_CQ_WS(e)) */
    union {
        int64_t align_i;
        void *align_p;
        float align_f;
        unsigned char bytes[XOMBE_ENC_CQ_SIZE];
    } cq;
#endif
#if XOMBE_HAVE_DSTAR
    /* D-STAR back end state, opaque here (dstar/dstar_enc.h: XOMBE_DS_ENC(e)) */
    union {
        int64_t align_i;
        void *align_p;
        float align_f;
        unsigned char bytes[XOMBE_DS_ENC_STATE];
    } ds;
#endif
};
#endif

#endif
