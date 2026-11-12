#include <math.h>
#include "tone_table.h"
#include "bits.h"
#include "../state.h"
#if XOMBE_HAVE_DSTAR
#include "../dstar/dstar_internal.h"
#endif

/* single tones synthesize as harmonic l1 of an f0 at most this */
#define SINGLE_F0_MAX_HZ 250.0f
#define DUAL_FIRST XOMBE_TONE_DUAL_FIRST
#define DUAL_IDS XOMBE_TONE_DUALS

/* dual tones 128..163 by id - 128, low and high frequencies in Hz
   DTMF 128..143 and call progress 160..162: DVSI AMBE-3000R Vocoder Chip Users Manual, Table 103
   Knox 144..159 and call progress 163: US 7,970,606 B2, Table 1 ("MBE tone parameters") */
static const short dual_low[DUAL_IDS] = {941, 697, 697, 697, 770, 770, 770, 852, 852, 852, 697, 770,
                                         852, 941, 941, 941, 820, 606, 606, 606, 672, 672, 672, 743,
                                         743, 743, 606, 672, 743, 820, 820, 820, 350, 440, 480, 350};
static const short dual_high[DUAL_IDS] = {1336, 1209, 1336, 1477, 1209, 1336, 1477, 1209, 1336, 1477, 1633, 1633,
                                          1633, 1633, 1209, 1477, 1162, 1052, 1162, 1297, 1052, 1162, 1297, 1052,
                                          1162, 1297, 1430, 1430, 1430, 1430, 1052, 1297, 440, 480, 620, 490};

/* f0 (Hz) and harmonic numbers l1 < l2 synthesizing each dual tone, US 7,970,606 B2, Table 1 */
static const float dual_f0[DUAL_IDS] = {
    /* 128..143, DTMF */
    78.50f, 173.48f, 70.0f, 87.0f, 109.95f, 191.68f, 70.17f, 71.06f, 121.58f, 212.0f, 116.41f, 96.15f,
    71.0f, 234.26f, 134.38f, 134.35f,
    /* 144..159, Knox */
    68.33f, 150.89f, 67.82f, 86.50f, 95.79f, 166.92f, 67.70f, 74.74f, 105.90f, 92.78f, 101.55f, 84.02f,
    67.83f, 102.30f, 117.0f, 117.49f,
    /* 160..163, call progress */
    87.78f, 70.83f, 122.0f, 70.0f};
static const unsigned char dual_l[DUAL_IDS][2] = {
    {12, 17}, {4, 7}, {10, 19}, {8, 17}, {7, 11}, {4, 7}, {11, 21}, {12, 17},
    {7, 11}, {4, 7}, {6, 14}, {8, 17}, {12, 23}, {4, 7}, {7, 9}, {7, 11},
    {12, 17}, {4, 7}, {9, 17}, {7, 15}, {7, 11}, {4, 7}, {10, 19}, {10, 14},
    {7, 11}, {8, 14}, {6, 14}, {8, 17}, {11, 21}, {8, 14}, {7, 9}, {7, 11},
    {4, 5}, {6, 7}, {4, 5}, {5, 7}};

int xombe__tone_id_accepted(int id)
{
    return (id >= 5 && id <= 122) || (id >= DUAL_FIRST && id < DUAL_FIRST + DUAL_IDS);
}

int xombe__tone_valid(int id, int ad)
{
    return xombe__tone_id_accepted(id) && ad >= 0 && ad <= 127;
}

int xombe__tone_freqs(int id, float *f_low, float *f_high)
{
    if (id >= 5 && id <= 122) {
        *f_low = *f_high = 31.25f * (float)id;
        return 1;
    }
    if (id >= DUAL_FIRST && id < DUAL_FIRST + DUAL_IDS) {
        *f_low = dual_low[id - DUAL_FIRST];
        *f_high = dual_high[id - DUAL_FIRST];
        return 2;
    }
    return 0;
}

void xombe__tone_harmonics(int id, float *f0_hz, int *l1, int *l2)
{
    float lo, hi;
    int kind = xombe__tone_freqs(id, &lo, &hi);
    *l1 = *l2 = 0;
    *f0_hz = SINGLE_F0_MAX_HZ;
    if (kind == 1) {
        *l1 = (int)ceilf(lo / SINGLE_F0_MAX_HZ);
        *f0_hz = lo / (float)*l1;
    } else if (kind == 2) {
        *f0_hz = dual_f0[id - DUAL_FIRST];
        *l1 = dual_l[id - DUAL_FIRST][0];
        *l2 = dual_l[id - DUAL_FIRST][1];
    }
}

/* 16384*10^(0.03555*(ad - 127)), in log2 */
#define TONE_LOG2_FULL 14.0f
#define TONE_LOG2_PER_AD 0.118094543f /* 0.03555*log2 10 */

float xombe__tone_log2_amp(int ad)
{
    return TONE_LOG2_FULL + TONE_LOG2_PER_AD * (float)(ad - 127);
}

int xombe__tone_frame_parse(const uint8_t d[49], int *id, int *ad)
{
    int c[4], i, dist = 0;
    for (i = 0; i < 4; i++)
        c[i] = xombe__get_bits(d, 12 + 8 * i, 8);
    *id = c[0];
    *ad = (xombe__get_bits(d, 6, 6) << 1) | d[44];
    for (i = 1; i < 4; i++)
        dist += xombe__popcount32((uint32_t)(c[0] ^ c[i]));
    if ((c[1] >> 4) != (c[0] >> 4) || dist > 5 || !xombe__tone_id_accepted(c[0]))
        return 0;
    return 1;
}

int xombe__tone_frame_dmr(int id, int ad, uint8_t out[9])
{
    uint8_t d[49], bits[7];
    if (!xombe__tone_valid(id, ad))
        return -1;
    xombe__tone_frame_bits(id, ad, d);
    xombe__pack49(d, bits);
    xombe__fec_encode_dmr(bits, out);
    return 0;
}

int xombe__tone_frame_dstar(int id, int ad, uint8_t out[9])
{
#if XOMBE_HAVE_DSTAR
    return xombe__ds_tone_frame_id(id, ad, out);
#else
    (void)id;
    (void)ad;
    (void)out;
    return -1;
#endif
}

void xombe__tone_frame_bits(int id, int ad, uint8_t d[49])
{
    int i;
    for (i = 0; i < 49; i++)
        d[i] = 0;
    xombe__put_bits(d, 0, 6, 63);
    xombe__put_bits(d, 6, 6, ad >> 1);
    for (i = 0; i < 4; i++)
        xombe__put_bits(d, 12 + 8 * i, 8, id);
    d[44] = (uint8_t)(ad & 1);
}
