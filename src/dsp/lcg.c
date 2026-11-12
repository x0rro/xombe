#include "../internal.h"

uint16_t xombe__lcg_next(uint16_t u)
{
    return (uint16_t)((173u * u + 13849u) & 0xFFFFu);
}

int16_t xombe__lcg_draw(uint16_t *u)
{
    *u = xombe__lcg_next(*u);
    return (int16_t)(*u >= 32768 ? (int)*u - 65536 : (int)*u);
}
