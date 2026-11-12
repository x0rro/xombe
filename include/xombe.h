#ifndef XOMBE_H
#define XOMBE_H

#include <stddef.h>
#include <stdint.h>
#include "xombe_config.h"

/* DMR (AMBE+2 half rate) and D-STAR (AMBE 2400) decoders, encoder, channel FEC, and transcoders between the two.
   every call that handles frames takes the codec as an enum xombe_codec, at init or as an argument, and treats
   both codecs the same way

   build options, recorded in xombe_config.h as 0/1 and guarding the declarations below:
     XOMBE_HAVE_ENCODER  the encoder (CMake XOMBE_ENCODER, default on)
     XOMBE_HAVE_DSTAR    the D-STAR codec and the transcoders (XOMBE_DSTAR, default on). without it an init for
                         XOMBE_DSTAR returns NULL
     XOMBE_KERNELS_Q15   fixed-point kernels (XOMBE_KERNELS=q15), larger encoder state
   the decoder and the FEC are always built

   frames: a frame is the 72-bit channel frame as sent, XOMBE_FRAME_BYTES bytes MSB first, FEC included: DMR's
   rate-33 frame or D-STAR's voice frame. the voice bits without FEC are XOMBE_BITS_BYTES bytes MSB first, DMR
   d0..d48, D-STAR d0..d47 and a zero bit. the *_bits calls and xombe_fec_* take them, everything else frames

   conventions, for every object here:
   - state is opaque and caller-allocated, the library never allocates and has no mutable globals. pass init
     the *_SIZE bytes, aligned for any scalar type: the matching *_mem union (static, automatic, or a struct
     member) or malloc(*_size()). init returns the handle (the same address), or NULL if mem is NULL or the
     settings are invalid. there is no deinit: drop or reuse the memory, init again to reset
   - *_size() is the exact size, at most *_SIZE. sizes can change between versions, compare xombe_version()
     with XOMBE_VERSION when the library may come from elsewhere
   - no thread safety inside an instance: one thread at a time per instance. distinct instances are independent
   - handles and required buffers must be valid, they are not checked. buffers are caller-owned and only
     accessed during the call
   - PCM frames are 20 ms: 160 samples of 8 kHz signed 16-bit PCM
   - *_frame calls return 1 when the output holds the frame for an input, 0 when it holds something else (the
     encoder's silence frame while its delay fills, a decoder's pending half-frame), -1 for a call the object
     cannot take (nothing written or consumed)
   - decoders and transcoders take the XOMBE_FLAG_* flags and report each frame through *_decision_get
   - unknown flag bits are reserved, pass 0
   - the calls that take a codec at init are inline here and pick a per-codec function (xombe__*, not part of the
     API), so a program that names one codec as a constant links none of the other's code */

#ifdef __cplusplus
extern "C" {
#endif

#define XOMBE_VERSION_MAJOR 0
#define XOMBE_VERSION_MINOR 1
#define XOMBE_VERSION_PATCH 0
#define XOMBE_VERSION (XOMBE_VERSION_MAJOR * 10000 + XOMBE_VERSION_MINOR * 100 + XOMBE_VERSION_PATCH)

/* the XOMBE_VERSION the library was built with */
XOMBE_API int xombe_version(void);

enum xombe_codec {
    XOMBE_DMR = 0,  /* AMBE+2 half rate, 49 voice bits, DMR rate-33 FEC */
    XOMBE_DSTAR = 1 /* AMBE 2400, 48 voice bits, D-STAR FEC */
};

#define XOMBE_FRAME_BYTES 9
#define XOMBE_BITS_BYTES 7

/* flags of every call that decodes or transcodes a frame
   LOST: the frame was not received, its input is not read (may be NULL). concealed as a bad frame (repeat, mute
   after three concealed frames in a row), the error counts and the running error rate are left alone
   BF, HB: the *_bits calls only, the error rule's flags as xombe_fec_decode returns them. calls that take frames
   compute them */
#define XOMBE_FLAG_LOST 0x1u
#define XOMBE_FLAG_BF 0x4u /* bad frame: conceal by repeating */
#define XOMBE_FLAG_HB 0x8u /* high error rate: mute */

#define XOMBE_HALF_MIN 76
#define XOMBE_HALF_MAX 84

enum xombe_class {
    XOMBE_CLASS_VOICE = 0,
    XOMBE_CLASS_ERASURE = 1,
    XOMBE_CLASS_SILENCE = 2,
    XOMBE_CLASS_TONE = 3
};

enum xombe_conceal {
    XOMBE_CONCEAL_NONE = 0,
    XOMBE_CONCEAL_REPEAT = 1,
    XOMBE_CONCEAL_MUTE = 2
};

/* decoder noise source: production (the default), or the conformance generator, which reproduces the test vectors */
enum xombe_noise_mode {
    XOMBE_NOISE_PRODUCTION = 0,
    XOMBE_NOISE_CONFORMANCE = 1
};

/* a decoder's account of the last frame. before the first frame all 0 but mute_src = e0 = e1 = -1 */
typedef struct {
    enum xombe_class cls;       /* ERASURE for a DMR erasure frame or a lost one, SILENCE for a DMR silence frame
                                   or a D-STAR null frame */
    enum xombe_conceal conceal; /* what played in place of the frame: NONE, REPEAT, or MUTE */
    int run;                    /* concealed frames in a row up to this one, 0 if it was not concealed */
    int64_t mute_src;           /* with MUTE: frame (0 = first) the muted set came from, -1 for the initial set */
    int e0, e1;                 /* as xombe_fec_stats eps0 and eps1, -1 if lost or from a *_bits call */
    int bf, hb;                 /* error rule: bad frame, high error rate (*_bits calls: the flags passed), 0 if
                                   lost */
} xombe_decision;

/* channel FEC: voice bits <-> frames, stateless but for the running error rate
   eps0, eps1: Golay corrections in the two protected words, DMR 0..4 and 0..3, D-STAR 0..4 each, 4 =
   uncorrectable. eps_t: their sum. bf, hb: the codec's error rule
   eps_r: running error rate 0.95*eps_r + 0.05*eps_t/n per frame, n 47 (DMR) or 48 (D-STAR), Q24 (1 << 24 = 1.0) */
typedef struct {
    int eps0, eps1, eps_t;
    int bf, hb;
    uint32_t eps_r;
} xombe_fec_stats;

typedef struct {
    int codec;
    uint32_t eps_r;
    float er;
} xombe_fec; /* contents private, set by xombe_fec_init */

XOMBE_API int xombe__fec_init_dmr(xombe_fec *f);
XOMBE_API int xombe__fec_init_dstar(xombe_fec *f);
XOMBE_API void xombe__fec_encode_dmr(const uint8_t bits[XOMBE_BITS_BYTES], uint8_t frame[XOMBE_FRAME_BYTES]);
XOMBE_API int xombe__fec_encode_dstar(const uint8_t bits[XOMBE_BITS_BYTES], uint8_t frame[XOMBE_FRAME_BYTES]);
/* returns 0, -1 for a codec not built (f unchanged) */
static inline int xombe_fec_init(xombe_fec *f, enum xombe_codec codec)
{
    return codec == XOMBE_DMR ? xombe__fec_init_dmr(f) : codec == XOMBE_DSTAR ? xombe__fec_init_dstar(f) : -1;
}
/* returns 0, -1 with nothing written for a codec not built */
static inline int xombe_fec_encode(enum xombe_codec codec, const uint8_t bits[XOMBE_BITS_BYTES],
                                   uint8_t frame[XOMBE_FRAME_BYTES])
{
    if (codec == XOMBE_DMR) {
        xombe__fec_encode_dmr(bits, frame);
        return 0;
    }
    return codec == XOMBE_DSTAR ? xombe__fec_encode_dstar(bits, frame) : -1;
}
/* fills st (may be NULL), returns the frame's XOMBE_FLAG_* word */
XOMBE_API unsigned xombe_fec_decode(xombe_fec *f, const uint8_t frame[XOMBE_FRAME_BYTES],
                                    uint8_t bits[XOMBE_BITS_BYTES], xombe_fec_stats *st);

/* decoder: frames or voice bits in, PCM out, no delay. FEC decoded here, classes, concealment, then synthesis
   DMR: an erasure frame or a bad one repeats the last played set, a high running error rate mutes
   D-STAR: classes voice, tone, or null frame (b0 124 or 125, and the DTX fill b0 127, ff x 9). a bad voice frame
   repeats the last played set, a high running error rate mutes */
typedef struct xombe_dec xombe_dec;
#define XOMBE_DEC_SIZE 5248
typedef union {
    int64_t align_i;
    void *align_p;
    unsigned char bytes[XOMBE_DEC_SIZE];
} xombe_dec_mem;

XOMBE_API size_t xombe_dec_size(void);
XOMBE_API xombe_dec *xombe__dec_init_dmr(void *mem);
XOMBE_API xombe_dec *xombe__dec_init_dstar(void *mem);
static inline xombe_dec *xombe_dec_init(void *mem, enum xombe_codec codec)
{
    return codec == XOMBE_DMR ? xombe__dec_init_dmr(mem) : codec == XOMBE_DSTAR ? xombe__dec_init_dstar(mem) : NULL;
}
/* noise source and its seed (0: default), restarts the source. init sets XOMBE_NOISE_PRODUCTION, seed 0.
   call before the first frame, returns 0, -1 for an unknown mode */
XOMBE_API int xombe_dec_set_noise(xombe_dec *d, enum xombe_noise_mode mode, uint32_t seed);
/* flags: 0 or XOMBE_FLAG_LOST. writes 160 samples, returns 1. after an unpaired first half it only ends that
   frame and returns 0: pcm[0..79] zeros, pcm[80..159] the pending second half, in unread, so pass it again */
XOMBE_API int xombe_dec_frame(xombe_dec *d, const uint8_t in[XOMBE_FRAME_BYTES], unsigned flags, int16_t pcm[160]);
/* one half-frame of n samples, n clamped to 76..84 (pcm holds XOMBE_HALF_MAX whatever n), returns the count
   written. calls alternate first half (second_half 0, reads in and flags) and second half (in unread, may be
   NULL), a call out of that order writes nothing and returns -1. two halves of 80 equal xombe_dec_frame */
XOMBE_API int xombe_dec_half(xombe_dec *d, const uint8_t in[XOMBE_FRAME_BYTES], unsigned flags, int second_half,
                             int n, int16_t *pcm);
/* xombe_dec_frame and xombe_dec_half for voice bits, FEC done by the caller. flags: XOMBE_FLAG_LOST, or
   XOMBE_FLAG_BF | XOMBE_FLAG_HB as xombe_fec_decode returns them */
XOMBE_API int xombe_dec_frame_bits(xombe_dec *d, const uint8_t bits[XOMBE_BITS_BYTES], unsigned flags,
                                   int16_t pcm[160]);
XOMBE_API int xombe_dec_half_bits(xombe_dec *d, const uint8_t bits[XOMBE_BITS_BYTES], unsigned flags,
                                  int second_half, int n, int16_t *pcm);
/* *_bits calls: FEC error counts for adaptive smoothing (xombe_fec_stats eps_t, eps1), kept until set again.
   the frame calls set them from their own FEC */
XOMBE_API void xombe_dec_set_errors(xombe_dec *d, int eps_t, int eps1);
XOMBE_API void xombe_dec_decision_get(const xombe_dec *d, xombe_decision *out);

/* tone frame of AMBE+2 tone ID id (5..122, 128..163) at level ad (0..127), as xombe_enc_tone. D-STAR sends the
   tone at AD max(ad, 8), as the encoder and the transcoder do. returns 0, -1 with nothing written for an id or ad
   out of range, a codec not built, or a D-STAR id without a tone code (Knox IDs 144..159) */
XOMBE_API int xombe__tone_frame_dmr(int id, int ad, uint8_t out[XOMBE_FRAME_BYTES]);
XOMBE_API int xombe__tone_frame_dstar(int id, int ad, uint8_t out[XOMBE_FRAME_BYTES]);
static inline int xombe_tone_frame(enum xombe_codec codec, int id, int ad, uint8_t out[XOMBE_FRAME_BYTES])
{
    return codec == XOMBE_DMR     ? xombe__tone_frame_dmr(id, ad, out)
           : codec == XOMBE_DSTAR ? xombe__tone_frame_dstar(id, ad, out)
                                  : -1;
}

/* encoder outputs, xombe_enc_cfg.outputs. XOMBE_ENC_DSTAR needs XOMBE_HAVE_DSTAR, init refuses it otherwise */
#define XOMBE_ENC_DMR 0x1u
#define XOMBE_ENC_DSTAR 0x2u

/* DMR amplitude quantizers, xombe_enc_cfg.dmr_quant */
#define XOMBE_ENC_QUANT_CLOSED 0 /* weighted closed-loop search, better quality */
#define XOMBE_ENC_QUANT_OPEN 1   /* open loop, faster */

#if XOMBE_HAVE_ENCODER
/* encoder delay, frames: DMR alone, and with the D-STAR output on. D-STAR's tone rules look one frame further
   ahead than the pitch tracker (a tone-in-noise frame needs the next frame's tone decision). both outputs share
   the delay: with both on the DMR frames are those of a DMR-only encoder, one call later */
#define XOMBE_ENC_DELAY 3
#define XOMBE_ENC_DSTAR_DELAY (XOMBE_ENC_DELAY + 1)

typedef struct xombe_enc xombe_enc;
/* state size: larger with q15 kernels, and XOMBE_ENC_DSTAR_SIZE more with D-STAR, whether or not used. without
   D-STAR, XOMBE_ENC_CQ_SIZE more for the closed-loop quantizer's workspace (with D-STAR it shares the D-STAR's) */
#if XOMBE_HAVE_DSTAR
#define XOMBE_ENC_DSTAR_SIZE 21504
#define XOMBE_ENC_CQ_SIZE 0
#else
#define XOMBE_ENC_DSTAR_SIZE 0
#define XOMBE_ENC_CQ_SIZE 19968
#endif
#if XOMBE_KERNELS_Q15
#define XOMBE_ENC_SIZE (39424 + XOMBE_ENC_DSTAR_SIZE + XOMBE_ENC_CQ_SIZE)
#else
#define XOMBE_ENC_SIZE (25088 + XOMBE_ENC_DSTAR_SIZE + XOMBE_ENC_CQ_SIZE)
#endif
typedef union {
    int64_t align_i;
    void *align_p;
    unsigned char bytes[XOMBE_ENC_SIZE];
} xombe_enc_mem;

/* xombe_enc_cfg_default's values in parentheses, flags are on when nonzero */
typedef struct {
    int ns;                     /* noise suppression (on) */
    int ns_strength_permille;   /* noise suppression strength, 0..1000 (500) */
    int tones;                  /* tone detection, both outputs (on) */
    int voicing_scale_permille; /* voicing decision scale, > 0 (350, 600 is less heavily voiced) */
    unsigned outputs;           /* XOMBE_ENC_DMR and/or XOMBE_ENC_DSTAR (DMR, 0 means DMR) */
    int dmr_dtx;                /* DMR silence frames for background (on) */
    int dstar_dtx;              /* D-STAR null frames for background (off) */
    int dmr_quant;              /* DMR amplitude quantizer, XOMBE_ENC_QUANT_* (CLOSED) */
} xombe_enc_cfg;

/* tone sent in place of one input frame's speech, xombe_enc_frame */
typedef struct {
    int id; /* AMBE+2 tone ID: 5..122 single tones at 31.25 Hz * id, 128..163 dual tones (DTMF 128..143, Knox
               144..159, call progress 160..163) */
    int ad; /* level 0..127, each component's amplitude 16384 * 10^(0.03555 * (ad - 127)), 0.711 dB per step */
} xombe_enc_tone;

/* encoder: PCM in, frames of either or both codecs out from one analysis, xombe_enc_delay(e) frames late */
XOMBE_API void xombe_enc_cfg_default(xombe_enc_cfg *cfg);
XOMBE_API size_t xombe_enc_size(void);
/* cfg NULL: the defaults, cfg is copied. returns NULL if mem is NULL or cfg is invalid: ns_strength_permille
   outside 0..1000, voicing_scale_permille <= 0, an unknown outputs bit, XOMBE_ENC_DSTAR without
   XOMBE_HAVE_DSTAR, or an unknown dmr_quant */
XOMBE_API xombe_enc *xombe_enc_init(void *mem, const xombe_enc_cfg *cfg);
/* XOMBE_ENC_DELAY, or XOMBE_ENC_DSTAR_DELAY when the D-STAR output is on */
XOMBE_API int xombe_enc_delay(const xombe_enc *e);
/* returns 1 when the frames hold the input xombe_enc_delay(e) calls earlier, else 0. pcm NULL: 160 zeros, at the
   end of input call it xombe_enc_delay(e) times so every input frame has its output. tone: this frame sent as a
   tone, even with cfg.tones 0, NULL for none. the tone frame is the output xombe_enc_delay(e) calls later, pcm is
   still analysed. D-STAR gets the tone at AD max(ad, 8), or for a Knox ID a voice frame from the tone table
   dmr and dstar may be NULL. an output the encoder does not produce, or any output before the delay has elapsed,
   is written as that codec's silence frame (D-STAR: the null frame). with both on the D-STAR frames differ
   slightly from a D-STAR-only encoder's: the DMR voicing then drives the shared noise tracking
   returns -1 with nothing consumed for a tone id or ad out of range */
XOMBE_API int xombe_enc_frame(xombe_enc *e, const int16_t *pcm, const xombe_enc_tone *tone,
                              uint8_t dmr[XOMBE_FRAME_BYTES], uint8_t dstar[XOMBE_FRAME_BYTES]);
/* xombe_enc_frame writing voice bits */
XOMBE_API int xombe_enc_frame_bits(xombe_enc *e, const int16_t *pcm, const xombe_enc_tone *tone,
                                   uint8_t dmr[XOMBE_BITS_BYTES], uint8_t dstar[XOMBE_BITS_BYTES]);
#endif

#if XOMBE_HAVE_DSTAR
/* transcoder between DMR and D-STAR frames, in the parameter domain: the source decoder's FEC, error rule,
   classes, concealment, and dequantization, then the parameters requantized against the target decoder's
   prediction memory. one frame out per frame in, no delay, nothing to flush. from D-STAR needs
   XOMBE_HAVE_ENCODER, init returns NULL without it
   DMR -> D-STAR
     voice     nearest D-STAR pitch, per-kHz voicing, amplitudes resampled to the D-STAR harmonics
     silence   all-unvoiced voice frame from the comfort-noise amplitudes
     tone      D-STAR tone frame, same tone and level. Knox tones (IDs 144..159, no D-STAR code) become voice
               frames at the tone table's f0
     concealed what the DMR decoder plays (its repeat or mute set), for a concealed tone the last output frame
               again
   D-STAR -> DMR
     voice     nearest DMR pitch, voicing from the previous and current D-STAR b1, amplitudes resampled to the
               DMR harmonics
     pulse     an all-unvoiced D-STAR frame whose b0 maps to a pulse: DMR pulse frame (b1 18..31)
     tone      DMR tone frame, same tone ID and level (Knox tones arrive as voice frames and stay voice)
     silence   null frames and mutes: the standard DMR silence frame
     concealed the D-STAR decoder's repeat set, requantized with the last good frame's voicing and pitch */
#define XOMBE_XC_DELAY 0

typedef struct xombe_xc xombe_xc;
#define XOMBE_XC_SIZE 23168
#define XOMBE_XC_SHARED_SIZE 3712 /* without the workspace, for xombe_xc_init_shared */
typedef union {
    int64_t align_i;
    void *align_p;
    unsigned char bytes[XOMBE_XC_SIZE];
} xombe_xc_mem;
typedef union {
    int64_t align_i;
    void *align_p;
    unsigned char bytes[XOMBE_XC_SHARED_SIZE];
} xombe_xc_shared_mem;

/* a transcoder's quantizer workspace: scratch only, nothing in it lives between calls, so transcoders that never
   run at the same time can share one (xombe_xc_init_shared). XOMBE_XC_SIZE includes its own */
#define XOMBE_QWS_SIZE 19968
typedef union {
    int64_t align_i;
    void *align_p;
    unsigned char bytes[XOMBE_QWS_SIZE];
} xombe_qws_mem;

enum xombe_xc_out {
    XOMBE_XC_OUT_VOICE = 0,   /* voice frame (to D-STAR: also for DMR silence) */
    XOMBE_XC_OUT_PULSE = 1,   /* to DMR: pulse frame, b1 18..31 */
    XOMBE_XC_OUT_TONE = 2,    /* tone frame */
    XOMBE_XC_OUT_KNOX = 3,    /* to D-STAR: voice frame for a Knox tone */
    XOMBE_XC_OUT_SILENCE = 4, /* to DMR: the standard silence frame */
    XOMBE_XC_OUT_RESEND = 5   /* to D-STAR: the previous output again (initially the D-STAR null frame) */
};

typedef struct {
    xombe_decision in;    /* the source decoder's account of the input frame */
    enum xombe_xc_out out;
} xombe_xc_decision;

XOMBE_API size_t xombe_xc_size(void);
XOMBE_API xombe_xc *xombe__xc_init_dmr(void *mem);
XOMBE_API xombe_xc *xombe__xc_init_dstar(void *mem);
XOMBE_API xombe_xc *xombe__xc_init_shared_dmr(void *mem, void *ws);
XOMBE_API xombe_xc *xombe__xc_init_shared_dstar(void *mem, void *ws);
/* from: the input frames' codec, the output is the other one */
static inline xombe_xc *xombe_xc_init(void *mem, enum xombe_codec from)
{
    return from == XOMBE_DMR ? xombe__xc_init_dmr(mem) : from == XOMBE_DSTAR ? xombe__xc_init_dstar(mem) : NULL;
}
/* mem: XOMBE_XC_SHARED_SIZE bytes (an xombe_xc_shared_mem). ws: XOMBE_QWS_SIZE bytes (an xombe_qws_mem),
   zeroed here, shared only between transcoders called from one thread. NULL if either is NULL */
static inline xombe_xc *xombe_xc_init_shared(void *mem, enum xombe_codec from, void *ws)
{
    return from == XOMBE_DMR     ? xombe__xc_init_shared_dmr(mem, ws)
           : from == XOMBE_DSTAR ? xombe__xc_init_shared_dstar(mem, ws)
                                 : NULL;
}
/* flags: 0 or XOMBE_FLAG_LOST. returns 1, out holds this input's frame */
XOMBE_API int xombe_xc_frame(xombe_xc *x, const uint8_t in[XOMBE_FRAME_BYTES], unsigned flags,
                             uint8_t out[XOMBE_FRAME_BYTES]);
XOMBE_API void xombe_xc_decision_get(const xombe_xc *x, xombe_xc_decision *out);
#endif

#ifdef __cplusplus
}
#endif

#endif
