#include "../internal.h"
#if XOMBE_HAVE_DSTAR
#include "../dstar/dstar_internal.h"
#endif

/* running error rate in Q24: eps_r <- 0.95*eps_r + 0.05*eps_t/47 = (893*eps_r + 2^24*eps_t)/940, rounded once
   HB above 0.096, from ceil(0.096*2^24) */
#define EPS_R_ONE (1ull << 24)
#define EPS_R_HB_MIN 1610613u

static uint32_t take(const uint8_t *bits, int count)
{
    return (uint32_t)xombe__get_bits(bits, 0, count);
}

static void put(uint8_t *bits, int count, uint32_t v)
{
    xombe__put_bits(bits, 0, count, (int)v);
}

static int parity(uint32_t v)
{
    return xombe__popcount32(v) & 1;
}

static void dmr_encode(const uint8_t bits[7], uint8_t air[9])
{
    uint8_t d[49], w[72], a[72];
    uint32_t data0, data1, c0, c1;
    int i;
    xombe__unpack49(bits, d);
    data0 = take(d, 12);
    data1 = take(d + 12, 12);
    c0 = xombe__golay23_encode(data0);
    c0 = (c0 << 1) | (uint32_t)parity(c0);
    c1 = xombe__golay23_encode(data1) ^ xombe__pn_mask(data0);
    put(w, 24, c0);
    put(w + 24, 23, c1);
    for (i = 0; i < 25; i++)
        w[47 + i] = d[24 + i];
    xombe__interleave(w, a);
    xombe__pack_bits(a, 72, air);
}

void xombe__fec_encode_dmr(const uint8_t bits[7], uint8_t frame[9])
{
    dmr_encode(bits, frame);
}

int xombe__fec_encode_dstar(const uint8_t bits[7], uint8_t frame[9])
{
#if XOMBE_HAVE_DSTAR
    xombe__ds_fec_encode_bits(bits, frame);
    return 0;
#else
    (void)bits;
    (void)frame;
    return -1;
#endif
}

void xombe__fec_stats_update(xombe_fec *f, xombe_fec_stats *st)
{
    st->eps_t = st->eps0 + st->eps1;
    st->bf = st->eps0 >= 4 || (st->eps0 >= 2 && st->eps_t >= 6);
    f->eps_r = (uint32_t)((893u * (uint64_t)f->eps_r + EPS_R_ONE * (uint64_t)st->eps_t + 470u) / 940u);
    st->eps_r = f->eps_r;
    st->hb = f->eps_r >= EPS_R_HB_MIN;
}

static unsigned dmr_decode(xombe_fec *f, const uint8_t air[9], uint8_t bits[7], xombe_fec_stats *st)
{
    uint8_t a[72], w[72], d[49];
    uint32_t c0, c1, cw0, cw1;
    int i, w0, received_parity;
    xombe__unpack_bits(air, 72, a);
    xombe__deinterleave(a, w);
    c0 = take(w, 24);
    received_parity = (int)(c0 & 1u);
    w0 = xombe__golay23_decode(c0 >> 1, &cw0);
    if ((parity(cw0) ^ received_parity) == 0)
        st->eps0 = w0;
    else
        st->eps0 = w0 <= 2 ? w0 + 1 : 4;
    c1 = take(w + 24, 23) ^ xombe__pn_mask(cw0 >> 11);
    st->eps1 = xombe__golay23_decode(c1, &cw1);
    put(d, 12, cw0 >> 11);
    put(d + 12, 12, cw1 >> 11);
    for (i = 0; i < 25; i++)
        d[24 + i] = w[47 + i];
    xombe__pack49(d, bits);
    xombe__fec_stats_update(f, st);
    return (st->bf ? XOMBE_FLAG_BF : 0u) | (st->hb ? XOMBE_FLAG_HB : 0u);
}

static void reset(xombe_fec *f, int codec)
{
    f->codec = codec;
    f->eps_r = 0;
    f->er = 0.0f;
}

int xombe__fec_init_dmr(xombe_fec *f)
{
    reset(f, XOMBE_DMR);
    return 0;
}

int xombe__fec_init_dstar(xombe_fec *f)
{
#if XOMBE_HAVE_DSTAR
    reset(f, XOMBE_DSTAR);
    return 0;
#else
    (void)f;
    return -1;
#endif
}

/* both codecs' decoders behind one call: the decoder's frame calls reach it, its *_bits calls do not */
unsigned xombe_fec_decode(xombe_fec *f, const uint8_t frame[9], uint8_t bits[7], xombe_fec_stats *st)
{
    xombe_fec_stats local;
    if (!st)
        st = &local;
#if XOMBE_HAVE_DSTAR
    if (f->codec == XOMBE_DSTAR)
        return xombe__ds_fec_decode_bits(f, frame, bits, st);
#endif
    return dmr_decode(f, frame, bits, st);
}
