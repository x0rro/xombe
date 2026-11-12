#ifndef XOMBE_INTERNAL_H
#define XOMBE_INTERNAL_H

#include <stdint.h>
#include "xombe.h"
#include "state.h"
#include "common/consts.h"
#include "common/bits.h"
#include "common/params.h"
#include "common/voicing.h"
#include "common/tone_table.h"

/* C99 has no _Static_assert, a negative array size fails the build */
#define XOMBE_STATIC_ASSERT(cond, name) typedef char xombe_static_assert_##name[(cond) ? 1 : -1]

/* largest xombe__irfft size, it keeps n/2 points of scratch on the stack */
#define XOMBE_IRFFT_MAX 256

/* n real samples in, n/2 + 1 bins out */
void xombe__rfft(const float *in, xombe_cpx *out, int n, const xombe_cpx *tw, int tw_n_real);
/* in place: io[k] holds samples 2k and 2k + 1 for k < n/2 on entry, bin k for k <= n/2 on return */
void xombe__rfft_packed(xombe_cpx *io, int n, const xombe_cpx *tw, int tw_n_real);
/* n <= XOMBE_IRFFT_MAX */
void xombe__irfft(const xombe_cpx *in, float *out, int n, const xombe_cpx *tw, int tw_n_real);

float xombe__sin_ratio(long r, long m);
void xombe__hamming(float *w, int n);
float xombe__hamming_sum2(int n);
float xombe__hamming_dtft(int n, long q, long g);

void xombe__hpf_init(xombe_hpf *f, float gain);
float xombe__hpf_step(xombe_hpf *f, float x);

uint16_t xombe__lcg_next(uint16_t u);
int16_t xombe__lcg_draw(uint16_t *u); /* next value as a signed 16-bit sample */

uint32_t xombe__golay23_encode(uint32_t data12);
/* corrects up to 3 errors in a 23-bit word, returns the count, *codeword gets the word (data in bits 22..11)
   an uncorrectable word also returns 3 and comes back unchanged, the FEC statistics depend on that count */
int xombe__golay23_decode(uint32_t received, uint32_t *codeword);
uint32_t xombe__pn_mask(uint32_t d0_12);
void xombe__interleave(const uint8_t w[72], uint8_t air_bits[72]);
void xombe__deinterleave(const uint8_t air_bits[72], uint8_t w[72]);
/* DMR error rule and running error rate from st->eps0 and eps1 */
void xombe__fec_stats_update(xombe_fec *f, xombe_fec_stats *st);

/* installs a stage tap (state.h), fn NULL removes it */
void xombe__dec_set_tap(xombe_dec *d, xombe_tap_fn fn, void *ctx);

/* where the D-STAR decoder's rules differ from AMBE+2's in the shared decoder, set at reset */
struct xombe_dec_rules {
    /* vh[0..XOMBE_LMAX] from the band voicing, 0 above L */
    void (*harm_voicing)(const uint8_t vb[8], float f0, int L, uint8_t vh[XOMBE_LMAX + 1]);
    float lambda_decay; /* smoothing level decay per frame */
};
/* a codec's frame logic, set by its init: start_bits begins a frame from voice bits (unread if lost) and the
   error counts in err_t and err_1 */
struct xombe_dec_ops {
    void (*start_bits)(xombe_dec *d, const uint8_t bits[7], unsigned flags);
};
void xombe__smooth(xombe_dec_core *d, float *a, int L, float f0);
void xombe__phases(const float *a, int L, float *phi);
/* e[1..L]: log2 amplitudes a[1..L] after spectral enhancement */
void xombe__enhance(const float *a, int L, float f0, float *e);
/* M[1..L] = 2^(e + spectral tilt): low harmonics rolled off below the average pitch, high ones boosted */
void xombe__shape(const float *e, float *M, int L, float f0, float avg_f0);
void xombe__tone_set(xombe_set *s, int id, int ad);
void xombe__synth_voiced(xombe_dec *d, const xombe_set *cur, const float *M, const uint8_t *vh, int cur_voiced,
                         int N, float *out);
void xombe__synth_unvoiced(xombe_dec *d, const xombe_set *cur, const float *M, const uint8_t *vh, int N,
                           float *out);
void xombe__noise_init(xombe_dec *d);
int xombe__conceal_decide(xombe_dec_core *d, int bad, int hb);

/* frame-independent decoder machinery, shared by the AMBE+2 parser (dec/decoder.c) and dstar/decoder.c
   a parser either conceals (xombe__dec_conceal_sets: mute or last played set, first half flagged repeat) or plays
   the frame's set (xombe__dec_play: picks the first-half set from the previous class, makes cur the second half,
   stores voice parameters p and band flags v2 in the prediction memory, v1 is a voice frame's first-half
   voicing), then xombe__dec_synth synthesizes each half. the frame logic takes the core, synthesis the whole
   decoder */
void xombe__dec_core_reset(xombe_dec_core *d);
void xombe__dec_reset(xombe_dec *d);
void xombe__dec_silence_set(xombe_set *s); /* initial silence set: f0 250 Hz, L 14, every amplitude 0.1 */
void xombe__dec_set_from_params(xombe_set *s, int cls, float f0, int L, const float *a, const uint8_t *vb);
void xombe__dec_play(xombe_dec_core *d, const xombe_set *cur, const uint8_t v1[8], const xombe_params *p,
                     const uint8_t v2[8]);
void xombe__dec_conceal_sets(xombe_dec_core *d, int conceal);
/* half-call order check: -1 if second_half is not the half due, else advances it and returns n clamped to
   XOMBE_HALF_MIN..XOMBE_HALF_MAX */
int xombe__dec_half_begin(xombe_dec *d, int second_half, int n);
void xombe__dec_synth(xombe_dec *d, int second_half, int N, int16_t *pcm);
/* AMBE+2 frame logic without synthesis: xombe__dec_start_frame classifies, conceals, or decodes 7 bytes (unread if
   lost) into p1 and p2, then xombe__dec_skip_synth updates the mute and played sets as the two half-frame calls
   would */
void xombe__dec_start_frame(xombe_dec_core *d, const uint8_t bytes[7], unsigned flags);
void xombe__dec_skip_synth(xombe_dec_core *d);
/* xombe__dec_start_frame of a DMR frame through fec, flags 0 or XOMBE_FLAG_LOST, bits gets the decoded voice
   bits (zeros if lost) */
void xombe__dec_start_air(xombe_dec_core *d, xombe_fec *fec, const uint8_t air[9], unsigned flags, uint8_t bits[7]);

#if XOMBE_HAVE_ENCODER
/* VAD threshold, metric in log2 power units
   confident noise half: metric below XOMBE_CONFIDENT, noise suppression and comfort noise
   adapt from the XOMBE_CONFIDENT_RUN-th in a row */
#define XOMBE_VAD_THRESHOLD XOMBE_DB(4.0f)
#define XOMBE_CONFIDENT XOMBE_DB(2.0f)
#define XOMBE_CONFIDENT_RUN 4

void xombe__enc_set_tap(xombe_enc *e, xombe_tap_fn fn, void *ctx);

void xombe__enc_set_window(xombe_enc *e, int n);
void xombe__enc_integrate(const float *P, int L, float f0, float *a);
void xombe__enc_hf_precomp(float *a, int L, float f0);

void xombe__pitch_init(xombe_enc *e);
void xombe__pitch_errors(xombe_enc *e, const float *centre, int slot);
float xombe__pitch_track(xombe_enc *e, int cur, int next1, int next2);
float xombe__pitch_refine(xombe_enc *e, const xombe_spec *S, float period);
int xombe__pitch_window_len(float period);

/* voicing metric per 500 Hz band (residual over energy of the harmonic line fit, min over +/-2% pitch jitter)
   and band weights, soft value of a scaled metric sm = clip(-0.661 - log2(sm)/2, 0, 1), and both for one scale */
void xombe__voicing_metric(xombe_enc *e, const xombe_spec *S, float fhat, float m[8], float w[8]);
float xombe__voicing_map(float sm);
void xombe__voicing_soft(xombe_enc *e, const xombe_spec *S, float fhat, float scale, float a[8], float w[8]);
int xombe__choose_b1(const float a2[8], const float w2[8], const float a1[8], const float w1[8],
                     const uint8_t v2_prev[8]);

int xombe__vad_flat(const int16_t *x, int n);
void xombe__vad_half(xombe_enc *e, const float *P, int flat, int *speech, float *v);

void xombe__ns_bins(const float *P, float Xm[XOMBE_NS_BINS]);
float xombe__ns_shape(float f);
/* noise suppression per half-frame (enc/frontend.h phase C)
   xombe__ns_update: shared noise estimate (from the confident_run-th confident noise half, or slowly in halves
   with no band voiced) and suppression floor (up in speech, down otherwise)
   xombe__ns_gain: one back end's gains G from that estimate, its own smoothing memory m, and vbin, the voicing it
   codes at each bin's centre
   xombe__ns_bin_voicing: vbin from 500 Hz band flags, ambe2_rule voices a bin if its band or the band 250 Hz
   above is (xombe__voiced_at), else the band alone */
void xombe__ns_update(xombe_enc *e, const float Xm[XOMBE_NS_BINS], int speech, int confident_run,
                      int no_band_voiced);
void xombe__ns_bin_voicing(const uint8_t vb[8], int ambe2_rule, uint8_t vbin[XOMBE_NS_BINS]);
void xombe__ns_gain(const xombe_enc *e, xombe_ns_mem *m, const float Xm[XOMBE_NS_BINS], int speech, float v,
                    const uint8_t vbin[XOMBE_NS_BINS], float G[XOMBE_NS_BINS]);
/* log2 gain offsets da[1..L] of the harmonics of f0 from gains G, voicing vh[1..L] (voiced harmonics below
   800 Hz take part of the suppression), xombe__ns_apply_voice adds them to a[1..L] with AMBE+2 band voicing v2,
   floored at XOMBE_AMP_MIN */
void xombe__ns_harmonic_offsets(const float G[XOMBE_NS_BINS], int L, float f0, const uint8_t *vh, float *da);
void xombe__ns_apply_voice(const float G[XOMBE_NS_BINS], float *a, int L, float f0, const uint8_t v2[8]);
void xombe__ns_apply_silence(const xombe_enc *e, float *a, int L, float f0);

void xombe__cn_half(xombe_enc *e, const float *P, int noise, int confident_run);
void xombe__cn_targets(const xombe_enc *e, float *a);

/* AMBE+2 tone detector on XOMBE_TONE_SEG samples from seg, gives tone ID and level code
   (xombe__tone_search with singles, the 36 dual tones, emin 0.99, coarse peak) */
int xombe__tone_detect(xombe_enc *e, const int16_t *seg, int *id, int *ad);
/* least-squares tone detector over a segment length and tone set, the level at each component's own frequency */
typedef struct {
    short id;
    float hz[2]; /* in the fit's order */
} xombe_tone_pair;
typedef struct {
    int len;          /* segment length, at most XOMBE_TONE_SEG */
    const float *win; /* its Hann window */
    float emin;       /* explained-fraction threshold */
    int singles;      /* single tones (IDs 5..122) are candidates */
    int fine_peak;    /* single-tone peak on an 8192-point grid refined in the fit (ternary search, +/-0.4 bin),
                         else the parabola on the XOMBE_TONE_FFT grid */
    const xombe_tone_pair *pairs;
    int npairs;
} xombe_tone_search;
typedef struct {
    int id;      /* -1: none */
    float e;     /* fit's explained fraction */
    float level; /* dB re full scale less the decoder filter's gain, AD = 127 + level/0.711 */
    float res;   /* fit's residual energy */
} xombe_tone_hit;
#if XOMBE_KERNELS_Q15
typedef uint32_t xombe_tone_pow; /* bin power in a per-segment unit, only ratios within a segment mean anything */
#else
typedef float xombe_tone_pow;
#endif
/* segment as the search reads it: x less its mean, pw bin powers of x under the window, zero-padded */
typedef struct {
    float x[XOMBE_TONE_SEG];
    xombe_tone_pow pw[XOMBE_TONE_FFT / 2 + 1];
    float energy; /* of x, 0 below the detector's floor */
} xombe_tone_seg;
/* sg for ts->len samples from seg under ts->win, which any search on sg must share */
void xombe__tone_prepare(xombe_enc *e, const int16_t *seg, const xombe_tone_search *ts, xombe_tone_seg *sg);
int xombe__tone_search(const xombe_tone_seg *sg, const xombe_tone_search *ts, xombe_tone_hit *hit);
/* decoder output filter's gain at hz, dB */
float xombe__tone_hpf_gain_db(float hz);
int xombe__tone_decide(xombe_enc *e, int detected, int det_id, int det_ad, int *id, int *ad);

void xombe__quantize(const xombe_pred *pi, int silence, int L, const float *target, int b[9]);
#endif

#endif
