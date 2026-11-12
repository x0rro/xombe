#ifndef XOMBE_BITS_H
#define XOMBE_BITS_H

#include <stdint.h>

/* one bit per byte <-> packed MSB first */
void xombe__unpack_bits(const uint8_t *bytes, int count, uint8_t *bits);
void xombe__pack_bits(const uint8_t *bits, int count, uint8_t *bytes);
void xombe__unpack49(const uint8_t bytes[7], uint8_t d[49]);
void xombe__pack49(const uint8_t d[49], uint8_t bytes[7]);
int xombe__popcount32(uint32_t v);
int xombe__frame_class(const uint8_t d[49]);
void xombe__get_fields(const uint8_t d[49], int b[9]);
void xombe__put_fields(const int b[9], uint8_t d[49]);
/* count <= 24 bits starting at d[first], MSB first */
int xombe__get_bits(const uint8_t *d, int first, int count);
void xombe__put_bits(uint8_t *d, int first, int count, int value);

#endif
