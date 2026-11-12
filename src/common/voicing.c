#include "voicing.h"
#include "consts.h"
#include "../state.h"
#include "../tables/tables.h"

void xombe__v2_from_b1(int b1, uint8_t v2[8])
{
    int j;
    for (j = 0; j < 8; j++)
        v2[j] = xombe__vuv[b1][j];
}

void xombe__f3_first_half(const uint8_t v2_prev[8], const uint8_t v2[8], uint8_t v1[8])
{
    int j;
    for (j = 0; j < 8; j++)
        v1[j] = (uint8_t)(v2_prev[j] && v2[j]);
}

static int band_of(float f)
{
    int j = (int)(f * XOMBE_BANDS_PER_UNIT);
    return j > 7 ? 7 : j;
}

int xombe__voiced_at(const uint8_t vb[8], float f)
{
    return vb[band_of(f)] || vb[band_of(f + 0.5f / XOMBE_BANDS_PER_UNIT)];
}

void xombe__harmonic_voicing(const uint8_t vb[8], float f0, int L, uint8_t vh[XOMBE_LMAX + 1])
{
    int l;
    vh[0] = 0;
    for (l = 1; l <= XOMBE_LMAX; l++) {
        float f = (float)l * f0;
        vh[l] = (uint8_t)(l <= L && xombe__voiced_at(vb, f));
    }
}

int xombe__count_voiced(const uint8_t vb[8])
{
    int j, n = 0;
    for (j = 0; j < 8; j++)
        n += vb[j] != 0;
    return n;
}
