#include <math.h>
#include <string.h>
#include "dstar_internal.h"

/* data bits of b0..b8, MSB first */
static const uint8_t field_width[9] = {7, 4, 6, 9, 7, 4, 4, 4, 3};
static const uint8_t field_bits[48] = {
    0, 1, 2, 3, 4, 5, 47,              /* b0 */
    37, 38, 39, 40,                    /* b1 */
    6, 7, 8, 9, 41, 42,                /* b2 */
    10, 11, 12, 13, 14, 15, 16, 43, 44, /* b3 */
    17, 18, 19, 20, 21, 45, 46,        /* b4 */
    22, 23, 24, 25,                    /* b5 */
    26, 27, 28, 29,                    /* b6 */
    30, 31, 32, 33,                    /* b7 */
    34, 35, 36,                        /* b8 */
};
/* tone frames: index T at these data bits, 7-bit level F at the next, half step at bit 17 */
static const uint8_t tone_t_bits[8] = {6, 7, 8, 9, 41, 42, 10, 11};
static const uint8_t tone_f_bits[7] = {12, 13, 14, 15, 16, 43, 44};
#define TONE_HALF_BIT 17

uint32_t xombe__ds_golay24_encode(uint32_t data12)
{
    uint32_t c23 = xombe__golay23_encode(data12 & 0xFFFu);
    return (c23 << 1) | (uint32_t)(xombe__popcount32(c23) & 1);
}

int xombe__ds_golay24_decode(uint32_t cw24, int *data12)
{
    uint32_t c23;
    int w = xombe__golay23_decode((cw24 >> 1) & 0x7FFFFFu, &c23);
    int parity_ok = ((xombe__popcount32(c23) + (int)(cw24 & 1u)) & 1) == 0;
    *data12 = (int)(c23 >> 11);
    if (w == 3 && !parity_ok)
        return XOMBE_DS_GOLAY_BAD;
    return w + !parity_ok;
}

uint32_t xombe__ds_pn24(int da)
{
    uint16_t u = (uint16_t)(16u * ((unsigned)da & 0xFFFu));
    uint32_t out = 0;
    int n;
    for (n = 0; n < 24; n++) {
        u = xombe__lcg_next(u);
        out = (out << 1) | ((u >> 15) & 1u);
    }
    return out;
}

/* air bit of A's word bit n (n = 0 the MSB), B and C are 2 and 4 later */
static int pos_a(int n)
{
    return n < 12 ? 6 * n : 6 * (n - 12) + 1;
}

void xombe__ds_deinterleave(const uint8_t in[9], uint32_t *A, uint32_t *B, uint32_t *C)
{
    uint32_t w[3] = {0, 0, 0};
    int n, k;
    for (n = 0; n < 24; n++)
        for (k = 0; k < 3; k++) {
            int p = pos_a(n) + 2 * k;
            w[k] = (w[k] << 1) | (uint32_t)((in[p >> 3] >> (7 - (p & 7))) & 1);
        }
    *A = w[0];
    *B = w[1];
    *C = w[2];
}

void xombe__ds_interleave(uint32_t A, uint32_t B, uint32_t C, uint8_t out[9])
{
    uint32_t w[3];
    int n, k;
    w[0] = A;
    w[1] = B;
    w[2] = C;
    memset(out, 0, 9);
    for (n = 0; n < 24; n++)
        for (k = 0; k < 3; k++) {
            int p = pos_a(n) + 2 * k;
            out[p >> 3] |= (uint8_t)(((w[k] >> (23 - n)) & 1u) << (7 - (p & 7)));
        }
}

void xombe__ds_get_fields(const uint8_t d[48], int b[9])
{
    int k, j, i = 0;
    for (k = 0; k < 9; k++) {
        int v = 0;
        for (j = 0; j < field_width[k]; j++)
            v = (v << 1) | (d[field_bits[i++]] & 1);
        b[k] = v;
    }
}

void xombe__ds_put_fields(const int b[9], uint8_t d[48])
{
    int k, j, i = 0;
    for (k = 0; k < 9; k++)
        for (j = 0; j < field_width[k]; j++)
            d[field_bits[i++]] = (uint8_t)((b[k] >> (field_width[k] - 1 - j)) & 1);
}

int xombe__ds_classify(int b0)
{
    if (b0 <= XOMBE_DS_B0_VOICE_MAX)
        return XOMBE_DS_VOICE;
    if (b0 == XOMBE_DS_B0_TONE)
        return XOMBE_DS_TONE;
    return XOMBE_DS_NULL;
}

static int bits_at(const uint8_t d[48], const uint8_t *idx, int n)
{
    int v = 0, i;
    for (i = 0; i < n; i++)
        v = (v << 1) | (d[idx[i]] & 1);
    return v;
}

int xombe__ds_tone_index(const uint8_t d[48])
{
    return bits_at(d, tone_t_bits, 8);
}

int xombe__ds_tone_level(const uint8_t d[48])
{
    /* 4*bits(12..16) + bits(43, 44) is the 7-bit F */
    int ad = bits_at(d, tone_f_bits, 7) + d[TONE_HALF_BIT] + 8;
    return ad < 127 ? ad : 127;
}

void xombe__ds_fec_decode(const uint8_t in[9], xombe_ds_frame *f)
{
    uint8_t d[XOMBE_DS_DATA_BITS];
    int i;
    xombe__ds_deinterleave(in, &f->A, &f->B, &f->C);
    f->ea = xombe__ds_golay24_decode(f->A, &f->da);
    f->pn = xombe__ds_pn24(f->da);
    f->eb = xombe__ds_golay24_decode(f->B ^ f->pn, &f->db);
    for (i = 0; i < 12; i++) {
        d[i] = (uint8_t)((f->da >> (11 - i)) & 1);
        d[12 + i] = (uint8_t)((f->db >> (11 - i)) & 1);
    }
    for (i = 0; i < 24; i++)
        d[24 + i] = (uint8_t)((f->C >> (23 - i)) & 1u);
    xombe__ds_frame_from_data(d, f);
}

void xombe__ds_frame_from_data(const uint8_t d[48], xombe_ds_frame *f)
{
    memcpy(f->d, d, XOMBE_DS_DATA_BITS);
    xombe__ds_get_fields(f->d, f->b);
    f->cls = xombe__ds_classify(f->b[0]);
    f->T = f->AD = f->tid = -1;
    if (f->cls == XOMBE_DS_TONE) {
        f->T = xombe__ds_tone_index(f->d);
        f->AD = xombe__ds_tone_level(f->d);
        if (xombe__ds_tone_id[f->T] != XOMBE_DS_TONE_NONE)
            f->tid = xombe__ds_tone_id[f->T];
    }
}

void xombe__ds_encode_bits(const uint8_t d[48], uint8_t out[9])
{
    uint32_t da = 0, db = 0, c = 0;
    int i;
    for (i = 0; i < 12; i++) {
        da = (da << 1) | (d[i] & 1u);
        db = (db << 1) | (d[12 + i] & 1u);
    }
    for (i = 24; i < 48; i++)
        c = (c << 1) | (d[i] & 1u);
    xombe__ds_interleave(xombe__ds_golay24_encode(da), xombe__ds_golay24_encode(db) ^ xombe__ds_pn24((int)da), c,
                         out);
}

/* d of fields b, each clamped to its width, returns 0, or -1 if any was clamped */
static int clamped_fields(const int b[9], uint8_t d[48])
{
    int c[9], k, r = 0;
    for (k = 0; k < 9; k++) {
        int top = (1 << field_width[k]) - 1;
        c[k] = b[k] < 0 ? 0 : b[k] > top ? top : b[k];
        r |= c[k] != b[k];
    }
    xombe__ds_put_fields(c, d);
    return r ? -1 : 0;
}

int xombe__ds_encode_fields(const int b[9], uint8_t out[9])
{
    uint8_t d[48];
    int r = clamped_fields(b, d);
    xombe__ds_encode_bits(d, out);
    return r;
}

int xombe__ds_fields_bits(const int b[9], uint8_t bits[7])
{
    uint8_t d[48];
    int r = clamped_fields(b, d);
    xombe__ds_pack(d, bits);
    return r;
}

int xombe__ds_tone_q(float level)
{
    /* 2*level is exact in float, so no epsilon against double rounding */
    float q = floorf(2.0f * level);
    if (!(q >= 16.0f)) /* also NaN */
        return 16;
    return q > 254.0f ? 254 : (int)q;
}

int xombe__ds_tone_code(int id)
{
    int T;
    for (T = 0; T < 256; T++)
        if (xombe__ds_tone_id[T] == id)
            return T;
    return -1;
}

void xombe__ds_pack(const uint8_t d[48], uint8_t bits[7])
{
    memset(bits, 0, 7);
    xombe__pack_bits(d, XOMBE_DS_DATA_BITS, bits);
}

void xombe__ds_unpack(const uint8_t bits[7], uint8_t d[48])
{
    xombe__unpack_bits(bits, XOMBE_DS_DATA_BITS, d);
}

void xombe__ds_fec_encode_bits(const uint8_t bits[7], uint8_t frame[9])
{
    uint8_t d[XOMBE_DS_DATA_BITS];
    xombe__ds_unpack(bits, d);
    xombe__ds_encode_bits(d, frame);
}

unsigned xombe__ds_fec_decode_bits(xombe_fec *fec, const uint8_t frame[9], uint8_t bits[7], xombe_fec_stats *st)
{
    uint32_t A, B, C;
    int da, db, ea, eb;
    xombe_ds_err er;
    xombe_ds_flags fl;
    xombe__ds_deinterleave(frame, &A, &B, &C);
    ea = xombe__ds_golay24_decode(A, &da);
    eb = xombe__ds_golay24_decode(B ^ xombe__ds_pn24(da), &db);
    da &= 0xFFF;
    db &= 0xFFF;
    bits[0] = (uint8_t)(da >> 4);
    bits[1] = (uint8_t)(((da & 15) << 4) | (db >> 8));
    bits[2] = (uint8_t)db;
    bits[3] = (uint8_t)(C >> 16);
    bits[4] = (uint8_t)(C >> 8);
    bits[5] = (uint8_t)C;
    bits[6] = 0;
    er.er = fec->er;
    xombe__ds_err_update(&er, ea, eb, &fl);
    fec->er = er.er;
    st->eps0 = ea;
    st->eps1 = eb;
    st->eps_t = ea + eb;
    st->bf = fl.bf;
    st->hb = fl.hb;
    st->eps_r = (uint32_t)(er.er * 16777216.0f + 0.5f);
    return (fl.bf ? XOMBE_FLAG_BF : 0u) | (fl.hb ? XOMBE_FLAG_HB : 0u);
}

int xombe__ds_tone_frame_id(int id, int ad, uint8_t out[9])
{
    int T = xombe__tone_valid(id, ad) ? xombe__ds_tone_code(id) : -1;
    if (T < 0)
        return -1;
    xombe__ds_tone_frame(T, xombe__ds_tone_q(XOMBE_DS_AD_LEVEL(ad)), out);
    return 0;
}

void xombe__ds_tone_bits(int T, int q, uint8_t d[48])
{
    int j, F;
    if (q < 16)
        q = 16;
    if (q > 254)
        q = 254;
    F = q / 2 - 8;
    memset(d, 0, 48);
    for (j = 0; j < 6; j++)
        d[j] = 1; /* b0 = 126: bits 0..5 set, bit 47 clear */
    for (j = 0; j < 8; j++)
        d[tone_t_bits[j]] = (uint8_t)((T >> (7 - j)) & 1);
    for (j = 0; j < 7; j++)
        d[tone_f_bits[j]] = (uint8_t)((F >> (6 - j)) & 1);
    d[TONE_HALF_BIT] = (uint8_t)(q & 1);
}

void xombe__ds_tone_frame(int T, int q, uint8_t out[9])
{
    uint8_t d[48];
    xombe__ds_tone_bits(T, q, d);
    xombe__ds_encode_bits(d, out);
}

void xombe__ds_err_init(xombe_ds_err *s)
{
    s->er = 0.0f;
}

void xombe__ds_err_update(xombe_ds_err *s, int ea, int eb, xombe_ds_flags *out)
{
    int eT = ea + eb;
    s->er = 0.95f * s->er + 0.05f * (float)eT / 48.0f;
    out->bf = ea >= 4;
    out->hb = s->er > 0.096f;
    if (eT < 6)
        eT = eb = 0;
    out->eT = eT;
    out->e1 = eb;
}
