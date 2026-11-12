/* src/dsp/q15ops.h operations against 64-bit reference arithmetic on edge values and a pseudo-random sweep,
   checks the C fallback on a host and the intrinsics on an Arm target with the DSP extension */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../src/dsp/fft_q15.h"

static int fails;

static void check(const char *what, int64_t got, int64_t want, uint32_t a, uint32_t b)
{
    if (got != want && fails++ < 10)
        printf("%s(%08lx, %08lx): %lld, want %lld\n", what, (unsigned long)a, (unsigned long)b, (long long)got,
               (long long)want);
}

static int64_t lo(uint32_t a)
{
    return (int64_t)(a & 0xffffu) - (a & 0x8000u ? 65536 : 0);
}

static int64_t hi(uint32_t a)
{
    return lo(a >> 16);
}

static int64_t wrap32(int64_t v)
{
    v %= 4294967296LL;
    if (v < 0)
        v += 4294967296LL;
    return v >= 2147483648LL ? v - 4294967296LL : v;
}

static int ref_clz64(uint64_t x)
{
    int n = 64;
    while (x) {
        n--;
        x >>= 1;
    }
    return n;
}

static int64_t sat16(int64_t v)
{
    return v > 32767 ? 32767 : v < -32768 ? -32768 : v;
}

static int64_t floor_half(int64_t v)
{
    return (int64_t)floor((double)v / 2.0);
}

static int64_t floor_div32(int64_t v)
{
    return v >= 0 ? v / 4294967296LL : -((-v + 4294967295LL) / 4294967296LL);
}

static int64_t pair(int64_t l, int64_t h)
{
    return (int64_t)(((uint32_t)l & 0xffffu) | ((uint32_t)h & 0xffffu) << 16);
}

static uint32_t ref_rbit(uint32_t x)
{
    uint32_t r = 0;
    int i;
    for (i = 0; i < 32; i++)
        r |= (x >> i & 1u) << (31 - i);
    return r;
}

/* FFT against a double-precision DFT of its quantized input: SNR in dB of the 512-point transform of 1024 real
   samples, packed as xombe__kern_spectrum packs them */
#define FM 512
#define PI 3.14159265358979323846
static double fft_snr(const double *in)
{
    static uint32_t tw[FM / 2], x[FM];
    static double zr[FM], zi[FM];
    double peak = 0.0, g, sig = 0.0, err = 0.0;
    int k, n, s;
    for (k = 0; k < FM / 2; k++) {
        double c = floor(cos(2.0 * PI * k / FM) * 32768.0 + 0.5);
        double si = floor(-sin(2.0 * PI * k / FM) * 32768.0 + 0.5);
        tw[k] = (uint32_t)pair(c > 32767 ? 32767 : (int64_t)c, si > 32767 ? 32767 : (int64_t)si);
    }
    for (n = 0; n < 2 * FM; n++)
        if (fabs(in[n]) > peak)
            peak = fabs(in[n]);
    g = XOMBE_CFFT_Q15_PEAK / peak;
    for (n = 0; n < FM; n++) {
        zr[n] = floor(in[2 * n] * g + 0.5);
        zi[n] = floor(in[2 * n + 1] * g + 0.5);
        x[q15_rbit((uint32_t)n) >> 23] = (uint32_t)pair((int64_t)zr[n], (int64_t)zi[n]);
    }
    s = xombe__cfft_q15(x, FM, tw, XOMBE_CFFT_Q15_PEAK);
    for (k = 0; k < FM; k++) {
        double xr = 0.0, xi = 0.0, er, ei;
        for (n = 0; n < FM; n++) {
            double a = -2.0 * PI * (double)((k * n) % FM) / FM, c = cos(a), si = sin(a);
            xr += zr[n] * c - zi[n] * si;
            xi += zr[n] * si + zi[n] * c;
        }
        er = ldexp((double)lo(x[k]), s) - xr;
        ei = ldexp((double)hi(x[k]), s) - xi;
        sig += xr * xr + xi * xi;
        err += er * er + ei * ei;
    }
    return 10.0 * log10(sig / err);
}

static double urand(uint32_t *s)
{
    *s = *s * 1664525u + 1013904223u;
    return (double)(*s >> 8) / 8388608.0 - 1.0;
}

/* noise, and voiced-like frames: harmonics of f0 with a falling envelope under the 199-point Hamming window,
   centred on sample 0 as xombe__kern_spectrum centres it, with noise 60 dB down */
static void fft_tests(void)
{
    static double in[2 * FM];
    static const double f0s[] = {0.0125, 0.03, 0.06}; /* cycles/sample, 100, 240, 480 Hz */
    uint32_t s = 7;
    double snr, worst = 1e9;
    int n, i, l;
    for (n = 0; n < 2 * FM; n++)
        in[n] = urand(&s);
    snr = fft_snr(in);
    printf("fft_q15 noise: %.1f dB\n", snr);
    if (snr < 62.0)
        fails++;
    for (i = 0; i < 3; i++) {
        memset(in, 0, sizeof in);
        for (n = -99; n <= 99; n++) {
            double v = 0.0, w = 0.54 + 0.46 * cos(PI * n / 99.0);
            for (l = 1; l * f0s[i] < 0.5; l++)
                v += cos(2.0 * PI * l * f0s[i] * n + l) / (double)l;
            in[(n + 2 * FM) % (2 * FM)] = w * (v + 1e-3 * urand(&s));
        }
        snr = fft_snr(in);
        if (snr < worst)
            worst = snr;
        printf("fft_q15 harmonic f0 %.4f: %.1f dB\n", f0s[i], snr);
    }
    if (worst < 57.0)
        fails++;
}

static void one(uint32_t a, uint32_t b, uint32_t c, uint64_t d)
{
    int64_t dot = lo(a) * lo(b) + hi(a) * hi(b);
    int32_t acc = (int32_t)wrap32((int64_t)c);
    int64_t acc64 = (int64_t)(d >> 1) - (int64_t)(d >> 2); /* keeps the sum within int64 */
    uint32_t pk = (a & 0xffffu) | ((b << 2) & 0xffff0000u);
    check("smuad", q15_smuad(a, b), wrap32(dot), a, b);
    check("smlad", q15_smlad(a, b, acc), wrap32(acc + dot), a, b);
    check("smlald", q15_smlald(a, b, acc64), acc64 + dot, a, b);
    check("pack_lsl2", q15_pack_lsl2(a, b), pk, a, b);
    check("smlsd", q15_smlsd(a, b, acc), wrap32(acc + lo(a) * lo(b) - hi(a) * hi(b)), a, b);
    check("smladx", q15_smladx(a, b, acc), wrap32(acc + lo(a) * hi(b) + hi(a) * lo(b)), a, b);
    check("qadd16", q15_qadd16(a, b), pair(sat16(lo(a) + lo(b)), sat16(hi(a) + hi(b))), a, b);
    check("qsub16", q15_qsub16(a, b), pair(sat16(lo(a) - lo(b)), sat16(hi(a) - hi(b))), a, b);
    check("shadd16", q15_shadd16(a, b), pair(floor_half(lo(a) + lo(b)), floor_half(hi(a) + hi(b))), a, b);
    check("shsub16", q15_shsub16(a, b), pair(floor_half(lo(a) - lo(b)), floor_half(hi(a) - hi(b))), a, b);
    check("qsax", q15_qsax(a, b), pair(sat16(lo(a) + hi(b)), sat16(hi(a) - lo(b))), a, b);
    check("qasx", q15_qasx(a, b), pair(sat16(lo(a) - hi(b)), sat16(hi(a) + lo(b))), a, b);
    check("shsax", q15_shsax(a, b), pair(floor_half(lo(a) + hi(b)), floor_half(hi(a) - lo(b))), a, b);
    check("shasx", q15_shasx(a, b), pair(floor_half(lo(a) - hi(b)), floor_half(hi(a) + lo(b))), a, b);
    check("pack_lsl1", q15_pack_lsl1(a, b), (a & 0xffffu) | ((b << 1) & 0xffff0000u), a, b);
    check("rbit", q15_rbit(a), ref_rbit(a), a, 0);
    check("smmulr", q15_smmulr((int32_t)a, (int32_t)b),
          floor_div32((int64_t)(int32_t)a * (int32_t)b + 0x80000000LL), a, b);
    check("umax16", q15_umax16(a, b),
          ((a >> 16) > (b >> 16) ? a >> 16 : b >> 16) << 16 | (((a & 0xffffu) > (b & 0xffffu) ? a : b) & 0xffffu),
          a, b);
    check("clz64", q15_clz64(d), ref_clz64(d), (uint32_t)(d >> 32), (uint32_t)d);
    check("clz64", q15_clz64(d >> (a & 63)), ref_clz64(d >> (a & 63)), a, (uint32_t)d);
}

int main(void)
{
    static const uint32_t edge[] = {0, 1, 0x7fff, 0x8000, 0xffff, 0x7fff7fff, 0x80008000, 0x80007fff,
                                    0x7fff8000, 0xffffffff, 0x00010000, 0x12345678};
    uint32_t s = 1;
    int i, j;
    for (i = 0; i < (int)(sizeof edge / sizeof edge[0]); i++)
        for (j = 0; j < (int)(sizeof edge / sizeof edge[0]); j++) {
            one(edge[i], edge[j], 0x7fffffff, 0);
            one(edge[i], edge[j], 0x80000000u, ~0ull);
        }
    for (i = 0; i < 200000; i++) {
        uint32_t a, b, c;
        uint64_t d;
        s = s * 1664525u + 1013904223u;
        a = s;
        s = s * 1664525u + 1013904223u;
        b = s;
        s = s * 1664525u + 1013904223u;
        c = s;
        s = s * 1664525u + 1013904223u;
        d = (uint64_t)c << 32 | s;
        one(a, b, c, d);
    }
    {
        const int16_t v[3] = {-2, 3, 5};
        uint32_t p = q15_ld2(v + 1);
        check("ld2", lo(p) * 100 + hi(p), 305, p, 0);
    }
    fft_tests();
    printf("q15ops (%s): %s\n",
#if defined(__ARM_FEATURE_DSP)
           "DSP extension",
#else
           "C fallback",
#endif
           fails ? "FAIL" : "ok");
    return fails != 0;
}
