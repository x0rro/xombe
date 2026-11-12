#ifndef XOMBE_Q15OPS_H
#define XOMBE_Q15OPS_H

/* fixed-point kernel DSP operations on packed int16 pairs (low half first) with ACLE semantics: Arm DSP
   extension instructions where the target has them (__ARM_FEATURE_DSP), the same arithmetic in C elsewhere
   a pair is a uint32_t of the two halves' bits, 32-bit results wrap modulo 2^32 as the instructions do (they
   only set the Q flag) */

#include <stdint.h>
#include <string.h>

static inline int32_t q15_lo(uint32_t a)
{
    return (int16_t)(uint16_t)(a & 0xffffu);
}

static inline int32_t q15_hi(uint32_t a)
{
    return (int16_t)(uint16_t)(a >> 16);
}

#if defined(__ARM_FEATURE_DSP)
#include <arm_acle.h>

/* lo(a)*lo(b) + hi(a)*hi(b) */
static inline int32_t q15_smuad(uint32_t a, uint32_t b)
{
    return __smuad((int16x2_t)a, (int16x2_t)b);
}

/* acc + lo(a)*lo(b) + hi(a)*hi(b) */
static inline int32_t q15_smlad(uint32_t a, uint32_t b, int32_t acc)
{
    return __smlad((int16x2_t)a, (int16x2_t)b, acc);
}

/* acc + lo(a)*lo(b) + hi(a)*hi(b) in 64 bits */
static inline int64_t q15_smlald(uint32_t a, uint32_t b, int64_t acc)
{
    return __smlald((int16x2_t)a, (int16x2_t)b, acc);
}

/* pair (lo(a), lo(b << 2)), PKHBT with a shifted top operand */
static inline uint32_t q15_pack_lsl2(uint32_t a, uint32_t b)
{
    uint32_t r;
    __asm__("pkhbt %0, %1, %2, lsl #2" : "=r"(r) : "r"(a), "r"(b));
    return r;
}

/* x != 0 */
static inline int q15_clz32(uint32_t x)
{
    return (int)__clz(x);
}

/* lo(a)*lo(b) - hi(a)*hi(b) + acc */
static inline int32_t q15_smlsd(uint32_t a, uint32_t b, int32_t acc)
{
    return __smlsd((int16x2_t)a, (int16x2_t)b, acc);
}

/* lo(a)*hi(b) + hi(a)*lo(b) + acc */
static inline int32_t q15_smladx(uint32_t a, uint32_t b, int32_t acc)
{
    return __smladx((int16x2_t)a, (int16x2_t)b, acc);
}

/* halves' sums and differences, saturated (QADD16, QSUB16) or halved and rounded down (SHADD16, SHSUB16) */
static inline uint32_t q15_qadd16(uint32_t a, uint32_t b)
{
    return (uint32_t)__qadd16((int16x2_t)a, (int16x2_t)b);
}

static inline uint32_t q15_qsub16(uint32_t a, uint32_t b)
{
    return (uint32_t)__qsub16((int16x2_t)a, (int16x2_t)b);
}

static inline uint32_t q15_shadd16(uint32_t a, uint32_t b)
{
    return (uint32_t)__shadd16((int16x2_t)a, (int16x2_t)b);
}

static inline uint32_t q15_shsub16(uint32_t a, uint32_t b)
{
    return (uint32_t)__shsub16((int16x2_t)a, (int16x2_t)b);
}

/* exchanged forms a + (-i)*b and a - (-i)*b for complex pairs (re, im): (lo(a) + hi(b), hi(a) - lo(b))
   (QSAX) and (lo(a) - hi(b), hi(a) + lo(b)) (QASX), saturated or halved */
static inline uint32_t q15_qsax(uint32_t a, uint32_t b)
{
    return (uint32_t)__qsax((int16x2_t)a, (int16x2_t)b);
}

static inline uint32_t q15_qasx(uint32_t a, uint32_t b)
{
    return (uint32_t)__qasx((int16x2_t)a, (int16x2_t)b);
}

static inline uint32_t q15_shsax(uint32_t a, uint32_t b)
{
    return (uint32_t)__shsax((int16x2_t)a, (int16x2_t)b);
}

static inline uint32_t q15_shasx(uint32_t a, uint32_t b)
{
    return (uint32_t)__shasx((int16x2_t)a, (int16x2_t)b);
}

/* pair (lo(a), lo(b << 1)) */
static inline uint32_t q15_pack_lsl1(uint32_t a, uint32_t b)
{
    uint32_t r;
    __asm__("pkhbt %0, %1, %2, lsl #1" : "=r"(r) : "r"(a), "r"(b));
    return r;
}

static inline uint32_t q15_rbit(uint32_t x)
{
    return __rbit(x);
}

/* per-half unsigned max (USUB16 then SEL on its GE bits) */
static inline uint32_t q15_umax16(uint32_t a, uint32_t b)
{
    uint32_t r;
    __asm__("usub16 %0, %1, %2\n\tsel %0, %1, %2" : "=&r"(r) : "r"(a), "r"(b) : "cc");
    return r;
}

/* a*b / 2^32 rounded to nearest, halves up */
static inline int32_t q15_smmulr(int32_t a, int32_t b)
{
    int32_t r;
    __asm__("smmulr %0, %1, %2" : "=r"(r) : "r"(a), "r"(b));
    return r;
}

#else

/* each product fits in int32, sums wrap in uint32 as the instructions do */
static inline int32_t q15_smuad(uint32_t a, uint32_t b)
{
    return (int32_t)((uint32_t)(q15_lo(a) * q15_lo(b)) + (uint32_t)(q15_hi(a) * q15_hi(b)));
}

static inline int32_t q15_smlad(uint32_t a, uint32_t b, int32_t acc)
{
    return (int32_t)((uint32_t)acc + (uint32_t)(q15_lo(a) * q15_lo(b)) + (uint32_t)(q15_hi(a) * q15_hi(b)));
}

static inline int64_t q15_smlald(uint32_t a, uint32_t b, int64_t acc)
{
    return (int64_t)((uint64_t)acc + (uint64_t)(int64_t)(q15_lo(a) * q15_lo(b)) +
                     (uint64_t)(int64_t)(q15_hi(a) * q15_hi(b)));
}

static inline uint32_t q15_pack_lsl2(uint32_t a, uint32_t b)
{
    return (a & 0xffffu) | ((b << 2) & 0xffff0000u);
}

static inline int q15_clz32(uint32_t x)
{
    int n = 0;
    if (!(x & 0xffff0000u)) {
        n += 16;
        x <<= 16;
    }
    if (!(x & 0xff000000u)) {
        n += 8;
        x <<= 8;
    }
    if (!(x & 0xf0000000u)) {
        n += 4;
        x <<= 4;
    }
    if (!(x & 0xc0000000u)) {
        n += 2;
        x <<= 2;
    }
    return x & 0x80000000u ? n : n + 1;
}

static inline int32_t q15_smlsd(uint32_t a, uint32_t b, int32_t acc)
{
    return (int32_t)((uint32_t)acc + (uint32_t)(q15_lo(a) * q15_lo(b)) - (uint32_t)(q15_hi(a) * q15_hi(b)));
}

static inline int32_t q15_smladx(uint32_t a, uint32_t b, int32_t acc)
{
    return (int32_t)((uint32_t)acc + (uint32_t)(q15_lo(a) * q15_hi(b)) + (uint32_t)(q15_hi(a) * q15_lo(b)));
}

static inline uint32_t q15_pair(int32_t lo, int32_t hi)
{
    return ((uint32_t)lo & 0xffffu) | (uint32_t)hi << 16;
}

static inline int32_t q15_sat(int32_t v)
{
    return v > 32767 ? 32767 : v < -32768 ? -32768 : v;
}

/* right shift of a negative value is implementation-defined in C, these halve with a floor */
static inline int32_t q15_half(int32_t v)
{
    return (v - (v & 1)) / 2;
}

static inline uint32_t q15_qadd16(uint32_t a, uint32_t b)
{
    return q15_pair(q15_sat(q15_lo(a) + q15_lo(b)), q15_sat(q15_hi(a) + q15_hi(b)));
}

static inline uint32_t q15_qsub16(uint32_t a, uint32_t b)
{
    return q15_pair(q15_sat(q15_lo(a) - q15_lo(b)), q15_sat(q15_hi(a) - q15_hi(b)));
}

static inline uint32_t q15_shadd16(uint32_t a, uint32_t b)
{
    return q15_pair(q15_half(q15_lo(a) + q15_lo(b)), q15_half(q15_hi(a) + q15_hi(b)));
}

static inline uint32_t q15_shsub16(uint32_t a, uint32_t b)
{
    return q15_pair(q15_half(q15_lo(a) - q15_lo(b)), q15_half(q15_hi(a) - q15_hi(b)));
}

static inline uint32_t q15_qsax(uint32_t a, uint32_t b)
{
    return q15_pair(q15_sat(q15_lo(a) + q15_hi(b)), q15_sat(q15_hi(a) - q15_lo(b)));
}

static inline uint32_t q15_qasx(uint32_t a, uint32_t b)
{
    return q15_pair(q15_sat(q15_lo(a) - q15_hi(b)), q15_sat(q15_hi(a) + q15_lo(b)));
}

static inline uint32_t q15_shsax(uint32_t a, uint32_t b)
{
    return q15_pair(q15_half(q15_lo(a) + q15_hi(b)), q15_half(q15_hi(a) - q15_lo(b)));
}

static inline uint32_t q15_shasx(uint32_t a, uint32_t b)
{
    return q15_pair(q15_half(q15_lo(a) - q15_hi(b)), q15_half(q15_hi(a) + q15_lo(b)));
}

static inline uint32_t q15_pack_lsl1(uint32_t a, uint32_t b)
{
    return (a & 0xffffu) | ((b << 1) & 0xffff0000u);
}

static inline uint32_t q15_rbit(uint32_t x)
{
    x = (x >> 1 & 0x55555555u) | (x & 0x55555555u) << 1;
    x = (x >> 2 & 0x33333333u) | (x & 0x33333333u) << 2;
    x = (x >> 4 & 0x0f0f0f0fu) | (x & 0x0f0f0f0fu) << 4;
    x = (x >> 8 & 0x00ff00ffu) | (x & 0x00ff00ffu) << 8;
    return x >> 16 | x << 16;
}

static inline uint32_t q15_umax16(uint32_t a, uint32_t b)
{
    uint32_t l = (a & 0xffffu) > (b & 0xffffu) ? a & 0xffffu : b & 0xffffu;
    return ((a >> 16) > (b >> 16) ? a & 0xffff0000u : b & 0xffff0000u) | l;
}

static inline int32_t q15_smmulr(int32_t a, int32_t b)
{
    return (int32_t)(((int64_t)a * b + 0x80000000) >> 32);
}

#endif

/* 64 for x = 0 */
static inline int q15_clz64(uint64_t x)
{
    uint32_t h = (uint32_t)(x >> 32), l = (uint32_t)x;
    if (h)
        return q15_clz32(h);
    return l ? 32 + q15_clz32(l) : 64;
}

/* p need only be 2-byte aligned */
static inline uint32_t q15_ld2(const int16_t *p)
{
    uint32_t v;
    memcpy(&v, p, 4);
    return v;
}

#endif
