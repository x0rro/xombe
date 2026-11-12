#ifndef XOMBE_VOICING_H
#define XOMBE_VOICING_H

#include <stdint.h>
#include "../state.h"

void xombe__v2_from_b1(int b1, uint8_t v2[8]);
void xombe__f3_first_half(const uint8_t v2_prev[8], const uint8_t v2[8], uint8_t v1[8]);
/* 1 if the band holding f, or the one half a band higher, is voiced */
int xombe__voiced_at(const uint8_t vb[8], float f);
void xombe__harmonic_voicing(const uint8_t vb[8], float f0, int L, uint8_t vh[XOMBE_LMAX + 1]);
int xombe__count_voiced(const uint8_t vb[8]);

#endif
