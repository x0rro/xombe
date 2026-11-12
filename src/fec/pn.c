#include "../internal.h"

uint32_t xombe__pn_mask(uint32_t d0_12)
{
    uint16_t u = (uint16_t)(16u * (d0_12 & 0xFFFu));
    uint32_t mask = 0;
    int n;
    for (n = 0; n < 23; n++) {
        u = xombe__lcg_next(u);
        mask = (mask << 1) | ((u >> 15) & 1u);
    }
    return mask;
}
