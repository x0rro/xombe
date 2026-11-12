#include <string.h>
#include "frontend.h"
#include "../common/prof.h"
#include "../dstar/dstar_internal.h"
#if XOMBE_HAVE_DSTAR
#include "../dstar/dstar_enc.h"
#endif

#define INPUT_HPF_GAIN 0.9427f
#if XOMBE_HAVE_DSTAR
#define ALL_OUTPUTS (XOMBE_ENC_DMR | XOMBE_ENC_DSTAR)
#else
#define ALL_OUTPUTS XOMBE_ENC_DMR
#endif

static const uint8_t silence_frame_bytes[7] = XOMBE_SILENCE_FRAME;
static const uint8_t null_bits[7] = XOMBE_DS_NULL_BITS;
#if !XOMBE_HAVE_DSTAR
static const uint8_t null_frame_bytes[XOMBE_DS_FRAME_BYTES] = XOMBE_DS_NULL_FRAME;
#endif

void xombe_enc_cfg_default(xombe_enc_cfg *cfg)
{
    cfg->ns = 1;
    cfg->ns_strength_permille = 500;
    cfg->dmr_dtx = 1;
    cfg->tones = 1;
    cfg->voicing_scale_permille = 350;
    cfg->outputs = XOMBE_ENC_DMR;
    cfg->dstar_dtx = 0;
    cfg->dmr_quant = XOMBE_ENC_QUANT_CLOSED;
}

XOMBE_STATIC_ASSERT(sizeof(struct xombe_enc) <= XOMBE_ENC_SIZE, encoder_state_fits_XOMBE_ENC_SIZE);
#if XOMBE_HAVE_DSTAR
XOMBE_STATIC_ASSERT(XOMBE_ENC_DSTAR_DELAY == XOMBE_ENC_DELAY + XOMBE_ENC_LOOKAHEAD, dstar_delay_is_the_lookahead);
#endif

int xombe_enc_delay(const xombe_enc *e)
{
    return XOMBE_ENC_DELAY + e->lookahead;
}

size_t xombe_enc_size(void)
{
    return sizeof(struct xombe_enc);
}

static int cfg_valid(const xombe_enc_cfg *cfg)
{
    return cfg->ns_strength_permille >= 0 && cfg->ns_strength_permille <= 1000 && cfg->voicing_scale_permille > 0 &&
           (cfg->outputs & ~ALL_OUTPUTS) == 0 &&
           (cfg->dmr_quant == XOMBE_ENC_QUANT_CLOSED || cfg->dmr_quant == XOMBE_ENC_QUANT_OPEN);
}

xombe_enc *xombe_enc_init(void *mem, const xombe_enc_cfg *cfg)
{
    xombe_enc *e = mem;
    xombe_enc_cfg c;
    if (cfg)
        c = *cfg;
    else
        xombe_enc_cfg_default(&c);
    if (!e || !cfg_valid(&c))
        return NULL;
    c.ns = c.ns != 0;
    c.dmr_dtx = c.dmr_dtx != 0;
    c.tones = c.tones != 0;
    c.dstar_dtx = c.dstar_dtx != 0;
    if (!c.outputs)
        c.outputs = XOMBE_ENC_DMR;
    memset(e, 0, sizeof(*e));
    e->cfg = c;
    e->ns_strength = (float)c.ns_strength_permille / 1000.0f;
    e->voicing_scale = (float)c.voicing_scale_permille / 1000.0f;
    xombe__hpf_init(&e->hpf, INPUT_HPF_GAIN);
    xombe__pitch_init(e);
    e->base_win_sum2 = xombe__hamming_sum2(XOMBE_BASE_WINDOW);
    xombe__enc_set_window(e, XOMBE_BASE_WINDOW);
    xombe__pred_init(&e->pi);
    e->tone_prev_id = -1;
#if XOMBE_HAVE_DSTAR
    if (c.outputs & XOMBE_ENC_DSTAR)
        e->lookahead = XOMBE_ENC_LOOKAHEAD;
    xombe__ds_enc_init(e);
#endif
    return e;
}

void xombe__enc_set_tap(xombe_enc *e, xombe_tap_fn fn, void *ctx)
{
    e->tap = fn;
    e->tap_ctx = ctx;
}

#ifdef XOMBE_PROF_DWT
unsigned long xombe__pf[16];
#endif

static int tone_valid(const xombe_enc_tone *t)
{
    return !t || xombe__tone_valid(t->id, t->ad);
}

/* ambe2 and dstar may be NULL, a disabled output gets its silence frame. a forced tone replaces only phase D's
   class, so analysis and background tracking see the frame as unforced */
static int encode(xombe_enc *e, const int16_t pcm[XOMBE_FRAME], const xombe_enc_tone *tone, uint8_t *ambe2,
                  uint8_t *dstar)
{
    xombe_fe_frame a;
    xombe_a2_frame f2;
    xombe_be_ns ns[XOMBE_BE_COUNT], *be[XOMBE_BE_COUNT];
    const uint8_t *req;
    int i, ready;
    memmove(e->tone_req[1], e->tone_req[0], sizeof(e->tone_req) - sizeof(e->tone_req[0]));
    e->tone_req[0][0] = (uint8_t)(tone ? tone->id : 0);
    e->tone_req[0][1] = (uint8_t)(tone ? tone->ad : 0);
    req = e->tone_req[XOMBE_ENC_DELAY + e->lookahead];
    XOMBE_PF(0, ready = xombe__fe_push(e, pcm));
#if XOMBE_HAVE_DSTAR
    if ((e->cfg.outputs & XOMBE_ENC_DSTAR) && e->cfg.tones)
        XOMBE_PF(7, xombe__ds_tone_push(e));
#endif
    if (!ready) {
        if (ambe2)
            memcpy(ambe2, silence_frame_bytes, 7);
        if (dstar)
            memcpy(dstar, null_bits, 7);
        return 0;
    }
    for (i = 0; i < XOMBE_BE_COUNT; i++)
        be[i] = NULL;
    XOMBE_PF(1, xombe__fe_analyze(e, &a));
    if (e->cfg.outputs & XOMBE_ENC_DMR) {
        be[XOMBE_BE_AMBE2] = &ns[XOMBE_BE_AMBE2];
        XOMBE_PF(2, xombe__a2_voicing(e, &a, &f2, be[XOMBE_BE_AMBE2]));
    }
#if XOMBE_HAVE_DSTAR
    if (e->cfg.outputs & XOMBE_ENC_DSTAR) {
        be[XOMBE_BE_DSTAR] = &ns[XOMBE_BE_DSTAR];
        XOMBE_PF(3, xombe__ds_enc_voicing(e, &a, be[XOMBE_BE_DSTAR]));
    }
#endif
    XOMBE_PF(4, xombe__fe_background(e, &a, be));
    if (be[XOMBE_BE_AMBE2]) {
        uint8_t tmp[7];
        if (req[0]) {
            f2.cls = XOMBE_CLASS_TONE;
            a.tone_id = req[0];
            a.tone_ad = req[1];
        }
        XOMBE_PF(5, xombe__a2_emit(e, &a, &f2, be[XOMBE_BE_AMBE2], ambe2 ? ambe2 : tmp));
    } else if (ambe2) {
        memcpy(ambe2, silence_frame_bytes, 7);
    }
#if XOMBE_HAVE_DSTAR
    if (be[XOMBE_BE_DSTAR]) {
        uint8_t tmp[7];
        if (req[0])
            xombe__ds_enc_force_tone(e, req[0], req[1]);
        XOMBE_PF(6, xombe__ds_enc_frame(e, &a, be[XOMBE_BE_DSTAR], dstar ? dstar : tmp));
    } else if (dstar) {
        memcpy(dstar, null_bits, 7);
    }
#else
    if (dstar)
        memcpy(dstar, null_bits, 7);
#endif
    return 1;
}

static const int16_t zeros[XOMBE_FRAME];

int xombe_enc_frame(xombe_enc *e, const int16_t *pcm, const xombe_enc_tone *tone, uint8_t dmr[XOMBE_FRAME_BYTES],
                    uint8_t dstar[XOMBE_FRAME_BYTES])
{
    uint8_t a[7], d[7];
    int r = xombe_enc_frame_bits(e, pcm, tone, dmr ? a : NULL, dstar ? d : NULL);
    if (r < 0)
        return r;
    if (dmr)
        xombe__fec_encode_dmr(a, dmr);
#if XOMBE_HAVE_DSTAR
    if (dstar)
        xombe__ds_fec_encode_bits(d, dstar);
#else
    if (dstar)
        memcpy(dstar, null_frame_bytes, XOMBE_DS_FRAME_BYTES);
#endif
    return r;
}

int xombe_enc_frame_bits(xombe_enc *e, const int16_t *pcm, const xombe_enc_tone *tone,
                         uint8_t dmr[XOMBE_BITS_BYTES], uint8_t dstar[XOMBE_BITS_BYTES])
{
    if (!tone_valid(tone))
        return -1;
    return encode(e, pcm ? pcm : zeros, tone, dmr, dstar);
}
