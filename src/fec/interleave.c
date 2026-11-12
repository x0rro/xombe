#include "../internal.h"

void xombe__interleave(const uint8_t w[72], uint8_t air_bits[72])
{
    int i, k;
    for (i = 0; i < 18; i++)
        for (k = 0; k < 4; k++)
            air_bits[4 * i + k] = w[18 * k + i];
}

void xombe__deinterleave(const uint8_t air_bits[72], uint8_t w[72])
{
    int i, k;
    for (i = 0; i < 18; i++)
        for (k = 0; k < 4; k++)
            w[18 * k + i] = air_bits[4 * i + k];
}
