#ifndef XOMBE_DSTAR_INTERNAL_H
#define XOMBE_DSTAR_INTERNAL_H

/* D-STAR (AMBE 2400) core
   frame: 72 bits (9 bytes, MSB first), three interleaved 24-bit words:
   A = Golay(24,12)(da), B = Golay(24,12)(db) ^ pn24(da), C = 24 uncoded bits
   data: d[0..47] = da, db, C (MSB first)
   fields: b0..b8 are 7, 4, 6, 9, 7, 4, 4, 4, 3 bits, scattered over d */

#include <stdint.h>
#include "../internal.h"
#include "../tables/tables.h"
#include "../tables/dstar_tables.h"

#define XOMBE_DS_FRAME_BYTES 9
#define XOMBE_DS_DATA_BITS 48
#define XOMBE_DS_GOLAY_BAD 4 /* error count of an uncorrectable Golay(24,12) word */

/* special b0 codes, b0 0..123 are voice */
#define XOMBE_DS_B0_VOICE_MAX 123
#define XOMBE_DS_B0_SILENCE 124 /* standard silence frame, XOMBE_DS_NULL_FRAME */
#define XOMBE_DS_B0_NULL125 125 /* plays as 124 */
#define XOMBE_DS_B0_TONE 126
#define XOMBE_DS_B0_FILL 127    /* DTX fill (ff x 9), plays as 124 */
/* fields 124 0 4 165 41 10 7 7 5, both Golay words clean */
#define XOMBE_DS_NULL_FRAME {0x9E, 0x8D, 0x32, 0x88, 0x26, 0x1A, 0x3F, 0x61, 0xE8}
#define XOMBE_DS_NULL_BITS {0xF8, 0x54, 0xAA, 0x9D, 0xE8, 0x0A, 0x00} /* its voice bits */
/* null frames play f0 250 Hz, L 14, all unvoiced, M 0.1 */
#define XOMBE_DS_NULL_F0 (250.0f / 8000.0f)
#define XOMBE_DS_NULL_L 14
#define XOMBE_DS_NULL_M 0.1f

#define XOMBE_DS_B1_SKIP 4    /* plays as 12, never sent */
#define XOMBE_DS_B1_VOICED 15 /* all bands voiced */

#define XOMBE_DS_GAIN_MAX 13.0f
/* D-STAR clamps a to [-15, 15 - 1e-9], which rounds to 15.0f (AMBE+2: 15 - 2^-11) */
#define XOMBE_DS_AMP_MAX 15.0f
#define XOMBE_DS_LAMBDA_DECAY (20.0f / 2048.0f) /* smoothing level decay per frame (AMBE+2: 0.0098) */

/* every 7-bit b0 is voice, tone, or null, so xombe__ds_classify never returns XOMBE_DS_ERASURE, it is for
   callers that lose a frame */
enum {
    XOMBE_DS_VOICE = 0,
    XOMBE_DS_TONE,
    XOMBE_DS_NULL,
    XOMBE_DS_ERASURE
};

/* everything the receive side learns from 9 bytes */
typedef struct {
    uint32_t A, B, C;  /* deinterleaved 24-bit words as received */
    uint32_t pn;       /* pn24(da) */
    int da, db;        /* corrected 12-bit data of A and of B ^ pn */
    int ea, eb;        /* corrected errors, XOMBE_DS_GOLAY_BAD = uncorrectable */
    uint8_t d[XOMBE_DS_DATA_BITS]; /* one bit per byte */
    int b[9];
    int cls;           /* XOMBE_DS_VOICE, TONE, or NULL */
    int T, AD;         /* tone frames: 8-bit tone index and AMBE+2 level 8..127, else -1 */
    int tid;           /* AMBE+2 tone ID, -1 for invalid T or not a tone frame */
} xombe_ds_frame;

/* Golay(24,12) = Golay(23,12) (data in bits 23..12) || even parity. decode returns the error count 0..3, or 4
   when uncorrectable (three corrections and a parity mismatch), and the data bits as the 23-bit decoder
   corrected them */
uint32_t xombe__ds_golay24_encode(uint32_t data12);
int xombe__ds_golay24_decode(uint32_t cw24, int *data12);
/* B's 24-bit whitening sequence seeded from A's data, AMBE+2's LCG, 24 steps */
uint32_t xombe__ds_pn24(int da);
/* A = bits 0, 6, ..., 66, 1, 7, ..., 67, B = A + 2, C = A + 4, MSB first */
void xombe__ds_deinterleave(const uint8_t in[9], uint32_t *A, uint32_t *B, uint32_t *C);
void xombe__ds_interleave(uint32_t A, uint32_t B, uint32_t C, uint8_t out[9]);
void xombe__ds_get_fields(const uint8_t d[48], int b[9]);
void xombe__ds_put_fields(const int b[9], uint8_t d[48]); /* each field masked to its width */
int xombe__ds_classify(int b0);
/* full receive path: FEC, data bits, fields, class, and tone fields */
void xombe__ds_fec_decode(const uint8_t in[9], xombe_ds_frame *f);
/* the part after FEC: d, fields, class, and tone fields (A, B, C, pn, da, db, ea, eb untouched) */
void xombe__ds_frame_from_data(const uint8_t d[48], xombe_ds_frame *f);
/* data bits <-> voice bits (7 bytes, MSB first, the last bit zero) */
void xombe__ds_pack(const uint8_t d[48], uint8_t bits[7]);
void xombe__ds_unpack(const uint8_t bits[7], uint8_t d[48]);
/* xombe_fec_encode and xombe_fec_decode for D-STAR, st not NULL */
void xombe__ds_fec_encode_bits(const uint8_t bits[7], uint8_t frame[9]);
unsigned xombe__ds_fec_decode_bits(xombe_fec *fec, const uint8_t frame[9], uint8_t bits[7], xombe_fec_stats *st);
/* xombe_tone_frame for D-STAR */
int xombe__ds_tone_frame_id(int id, int ad, uint8_t out[9]);
void xombe__ds_encode_bits(const uint8_t d[48], uint8_t out[9]);
/* always writes out, a field out of range clamped to its width. returns 0, or -1 if any was clamped */
int xombe__ds_encode_fields(const int b[9], uint8_t out[9]);
/* xombe__ds_encode_fields writing voice bits */
int xombe__ds_fields_bits(const int b[9], uint8_t bits[7]);
/* tone frames carry the level in half steps q = floor(2*level) clamped to 16..254: F = q/2 - 8 in bits
   12..16, 43, 44 and q & 1 in bit 17, so the decoder plays AD = q/2 + (q & 1), rounded half up, capped at 127 */
int xombe__ds_tone_q(float level);
/* AMBE+2 level ad as a tone level: q = 2*ad - 1, which plays as AD ad */
#define XOMBE_DS_AD_LEVEL(ad) ((float)(ad) - 0.5f)
/* D-STAR tone index of an AMBE+2 tone ID, -1 if it has none (Knox) */
int xombe__ds_tone_code(int id);
void xombe__ds_tone_bits(int T, int q, uint8_t d[48]);
void xombe__ds_tone_frame(int T, int q, uint8_t out[9]);
/* level AD = 4*bits(12..16) + bit 17 + bits(43,44) + 8, capped at 127 */
int xombe__ds_tone_index(const uint8_t d[48]);
int xombe__ds_tone_level(const uint8_t d[48]);

/* fitted error rule, stateful over the stream. per frame, with eT = ea + eb:
   er = 0.95*er + 0.05*eT/48, bf (repeat) iff ea >= 4, hb (mute) iff er > 0.096,
   smoothing inputs eT and e1 = eb, both 0 when eT < 6. null frames ignore bf and hb but still advance er */
typedef struct {
    float er;
} xombe_ds_err;

typedef struct {
    int bf, hb, eT, e1;
} xombe_ds_flags;

void xombe__ds_err_init(xombe_ds_err *s);
void xombe__ds_err_update(xombe_ds_err *s, int ea, int eb, xombe_ds_flags *out);

/* D-STAR amplitude tables: rho 0.8, gain recursion 0, gmax 13, clamp [-15, XOMBE_DS_AMP_MAX], fitted codebooks,
   and their ok masks */
extern const xombe_amp_tables xombe__amp_dstar;

/* voice frames fill everything up to a, others only cls, b, and the tone fields */
typedef struct {
    int cls;
    int b[9];
    float f0;                       /* cycles/sample */
    int L;
    uint8_t V[8];                   /* 500 Hz band flags, xombe__ds_vuv[b1] */
    uint8_t vl[XOMBE_LMAX + 1];     /* per-harmonic voicing vl[1..L] */
    float Tv[XOMBE_LMAX + 1];       /* spectral shape before mean removal, Tv[1..L] */
    float G;
    float a[XOMBE_LMAX + 1];        /* log2 amplitudes a[1..L], a[L + 1..XOMBE_LMAX] = a[L] */
    int T, tid, AD;                 /* tone frames, else -1 */
} xombe_ds_params;

/* dequantizer prediction memory, a extended to a[L + 1..XOMBE_LMAX] = a[L], reset L 15, a 0, G 0, copyable by
   assignment */
typedef struct {
    int L;
    float G;
    float a[XOMBE_LMAX + 1];
} xombe_ds_mirror;

/* any b0 0..127 */
float xombe__ds_f0(int b0);
int xombe__ds_L(int b0);
int xombe__ds_b0_from_f0(float f0); /* nearest pitch index in log f0, 0..123 */
/* harmonic l takes band min(floor(l*f0*16), 7) of V, without AMBE+2's +250 Hz OR, vl[0] = 0 */
void xombe__ds_harm_voicing(const uint8_t V[8], float f0, int L, uint8_t vl[XOMBE_LMAX + 1]);

/* Tv, G, and a of fields b at L harmonics against an explicit memory (Lp, ap extended, Gp), Tv may be NULL.
   the decoder's memory is smoothed on error frames, so it calls this, smooths a, then stores its own memory,
   the mirror never smooths */
void xombe__ds_amps(const int b[9], int L, int Lp, const float ap[XOMBE_LMAX + 1], float Gp,
                    float Tv[XOMBE_LMAX + 1], float *G, float a[XOMBE_LMAX + 1]);

void xombe__ds_mirror_init(xombe_ds_mirror *m);
/* prediction pt[1..L] for a frame of L harmonics */
void xombe__ds_mirror_predict(const xombe_ds_mirror *m, int L, float pt[XOMBE_LMAX + 1]);
/* decodes b against the memory without changing it. b0 > 123 gives only the class, tone fields need the data
   bits (xombe__ds_params_from_frame) */
void xombe__ds_dequant(const int b[9], const xombe_ds_mirror *m, xombe_ds_params *p);
/* decodes b and stores L, a, and G for a voice frame. returns 1 for voice, 0 otherwise (memory unchanged),
   p may be NULL */
int xombe__ds_mirror_commit(xombe_ds_mirror *m, const int b[9], xombe_ds_params *p);
/* received frame's parameters against the mirror, committing voice frames */
int xombe__ds_params_from_frame(xombe_ds_mirror *m, const xombe_ds_frame *f, xombe_ds_params *p);

#endif
