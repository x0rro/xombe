#include "bits.h"
#include "xombe.h"

static const uint8_t b0_bits[] = {0, 1, 2, 3, 37, 38, 39};
static const uint8_t b1_bits[] = {4, 5, 6, 7, 35};
static const uint8_t b2_bits[] = {8, 9, 10, 11, 36};
static const uint8_t b3_bits[] = {12, 13, 14, 15, 16, 17, 18, 19, 43};
static const uint8_t b4_bits[] = {20, 21, 22, 23, 40, 41, 42};
static const uint8_t b5_bits[] = {24, 25, 26, 27, 44};
static const uint8_t b6_bits[] = {28, 29, 30, 45};
static const uint8_t b7_bits[] = {31, 32, 33, 46};
static const uint8_t b8_bits[] = {34, 47, 48};

static const uint8_t *const field_bits[9] = {b0_bits, b1_bits, b2_bits, b3_bits, b4_bits,
                                             b5_bits, b6_bits, b7_bits, b8_bits};
static const uint8_t field_width[9] = {7, 5, 5, 9, 7, 5, 4, 4, 3};

void xombe__unpack_bits(const uint8_t *bytes, int count, uint8_t *bits)
{
    int i;
    for (i = 0; i < count; i++)
        bits[i] = (uint8_t)((bytes[i >> 3] >> (7 - (i & 7))) & 1);
}

void xombe__pack_bits(const uint8_t *bits, int count, uint8_t *bytes)
{
    int i;
    for (i = 0; i < (count + 7) / 8; i++)
        bytes[i] = 0;
    for (i = 0; i < count; i++)
        if (bits[i])
            bytes[i >> 3] |= (uint8_t)(0x80 >> (i & 7));
}

void xombe__unpack49(const uint8_t bytes[7], uint8_t d[49])
{
    xombe__unpack_bits(bytes, 49, d);
}

void xombe__pack49(const uint8_t d[49], uint8_t bytes[7])
{
    xombe__pack_bits(d, 49, bytes);
}

int xombe__popcount32(uint32_t v)
{
    int n = 0;
    for (; v; v &= v - 1)
        n++;
    return n;
}

int xombe__get_bits(const uint8_t *d, int first, int count)
{
    int v = 0, i;
    for (i = 0; i < count; i++)
        v = (v << 1) | d[first + i];
    return v;
}

void xombe__put_bits(uint8_t *d, int first, int count, int value)
{
    int i;
    for (i = 0; i < count; i++)
        d[first + i] = (uint8_t)((value >> (count - 1 - i)) & 1);
}

int xombe__frame_class(const uint8_t d[49])
{
    int c = xombe__get_bits(d, 0, 6);
    if (c < 60)
        return XOMBE_CLASS_VOICE;
    if (c < 62)
        return XOMBE_CLASS_ERASURE;
    if (c == 62)
        return XOMBE_CLASS_SILENCE;
    return XOMBE_CLASS_TONE;
}

void xombe__get_fields(const uint8_t d[49], int b[9])
{
    int f, i;
    for (f = 0; f < 9; f++) {
        int v = 0;
        for (i = 0; i < field_width[f]; i++)
            v = (v << 1) | d[field_bits[f][i]];
        b[f] = v;
    }
}

void xombe__put_fields(const int b[9], uint8_t d[49])
{
    int f, i;
    for (f = 0; f < 9; f++) {
        int w = field_width[f];
        for (i = 0; i < w; i++)
            d[field_bits[f][i]] = (uint8_t)((b[f] >> (w - 1 - i)) & 1);
    }
}
