#ifndef XOMBE_TONE_TABLE_H
#define XOMBE_TONE_TABLE_H

#include <stdint.h>

/* tone level ad: amplitude 16384 * 10^(XOMBE_TONE_AD_LOG10 * (ad - 127)) */
#define XOMBE_TONE_AD_LOG10 0.03555f
#define XOMBE_TONE_DB_PER_AD 0.7110f
#define XOMBE_TONE_DUAL_FIRST 128
#define XOMBE_TONE_DUALS 36
#define XOMBE_TONE_KNOX_FIRST 144
#define XOMBE_TONE_KNOX_COUNT 16

int xombe__tone_id_accepted(int id);
/* id accepted and level code ad 0..127, what a public tone request may name */
int xombe__tone_valid(int id, int ad);
int xombe__tone_freqs(int id, float *f_low, float *f_high);
void xombe__tone_harmonics(int id, float *f0_hz, int *l1, int *l2);
/* log2 synthesized amplitude of each tone component for level code ad */
float xombe__tone_log2_amp(int ad);
/* tone frame: class 63, AD in d6..d11 + d44, the ID four times from d12
   returns 0 if the copies disagree too much or the ID is not accepted */
int xombe__tone_frame_parse(const uint8_t d[49], int *id, int *ad);
void xombe__tone_frame_bits(int id, int ad, uint8_t d[49]);

#endif
