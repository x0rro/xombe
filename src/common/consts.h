#ifndef XOMBE_CONSTS_H
#define XOMBE_CONSTS_H

#define XOMBE_FS 8000.0f
/* Hz -> cycles/sample, the unit of every frequency (f0 = 1/period) but the tone tables and tone detector */
#define XOMBE_HZ(f) ((f) / XOMBE_FS)
#define XOMBE_DB(db) ((db) / 3.01029996f) /* dB -> log2 power units */
/* 8 voicing bands of 500 Hz, 16 = 8000/500 per cycle/sample, so the band of f is floor(16*f), exact in binary */
#define XOMBE_BANDS_PER_UNIT 16.0f
/* decoder boosts harmonics above this by f/XOMBE_HF_CORNER, encoder targets pre-compensate half of that */
#define XOMBE_HF_CORNER XOMBE_HZ(2400.0f)
#define XOMBE_RHO 0.65f
#define XOMBE_AMP_MIN (-15.0f)
#define XOMBE_AMP_MAX (15.0f - 1.0f / 2048.0f)
#define XOMBE_GAIN_MAX 13.0f
#define XOMBE_SILENCE_F0 0.031324f
#define XOMBE_SILENCE_L 15
#define XOMBE_SILENCE_GAIN_OFFSET 2.0f
#define XOMBE_PULSE_F0 (64.346f / 8000.0f)
#define XOMBE_PULSE_L 56
#define XOMBE_PULSE_B1_MIN 18
#define XOMBE_SILENCE_B0 124
#define XOMBE_SILENCE_B1 16
/* AMBE+2 silence frame, 7 bytes as xombe_dec_frame takes them */
#define XOMBE_SILENCE_FRAME {0xF8, 0x01, 0xA9, 0x9F, 0x8C, 0xE0, 0x80}
#define XOMBE_UNVOICED_B1 16
#define XOMBE_PI 3.14159265358979323846
#define XOMBE_PI_F 3.14159265f
#define XOMBE_LOG2_10 3.32192809f

#endif
