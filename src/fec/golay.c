#include "../internal.h"

#define GOLAY_POLY 0xC75u

static uint32_t golay_remainder(uint32_t v)
{
    int bit;
    for (bit = 22; bit >= 11; bit--)
        if (v & (1u << bit))
            v ^= GOLAY_POLY << (bit - 11);
    return v;
}

uint32_t xombe__golay23_encode(uint32_t data12)
{
    uint32_t shifted = (data12 & 0xFFFu) << 11;
    return shifted | golay_remainder(shifted);
}

int xombe__golay23_decode(uint32_t received, uint32_t *codeword)
{
    uint32_t syn = golay_remainder(received & 0x7FFFFFu), data_syn[12];
    int i, j, k;
    if (xombe__popcount32(syn) <= 3) {
        *codeword = (received ^ syn) & 0x7FFFFFu;
        return xombe__popcount32(syn);
    }
    for (i = 0; i < 12; i++)
        data_syn[i] = golay_remainder(1u << (11 + i));
    for (i = 0; i < 12; i++) {
        uint32_t s = syn ^ data_syn[i];
        if (xombe__popcount32(s) <= 2) {
            *codeword = (received ^ (1u << (11 + i)) ^ s) & 0x7FFFFFu;
            return 1 + xombe__popcount32(s);
        }
    }
    for (i = 0; i < 12; i++)
        for (j = i + 1; j < 12; j++) {
            uint32_t s = syn ^ data_syn[i] ^ data_syn[j];
            if (xombe__popcount32(s) <= 1) {
                *codeword = (received ^ (1u << (11 + i)) ^ (1u << (11 + j)) ^ s) & 0x7FFFFFu;
                return 2 + xombe__popcount32(s);
            }
        }
    for (i = 0; i < 12; i++)
        for (j = i + 1; j < 12; j++)
            for (k = j + 1; k < 12; k++)
                if ((syn ^ data_syn[i] ^ data_syn[j] ^ data_syn[k]) == 0) {
                    *codeword = (received ^ (1u << (11 + i)) ^ (1u << (11 + j)) ^ (1u << (11 + k))) & 0x7FFFFFu;
                    return 3;
                }
    *codeword = received & 0x7FFFFFu;
    return 3;
}
