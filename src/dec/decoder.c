#include <math.h>
#include <string.h>
#include "../internal.h"
#include "../dsp/fmath.h"
#if XOMBE_HAVE_DSTAR
#include "../dstar/dsd_state.h"
#endif

#define AVG_F0_INIT XOMBE_HZ(111.63f)
#define LAMBDA_INIT 7.0f
#define INIT_SET_F0 XOMBE_HZ(250.0f)
#define INIT_SET_L 14
#define INIT_SET_LOG2_AMP (-3.321928095f) /* log2(0.1) */
#define DEFAULT_PRNG_SEED 0x12345678u

static const struct xombe_dec_rules ambe2_rules = {xombe__harmonic_voicing, 0.0098f};

void xombe__dec_silence_set(xombe_set *s)
{
    int l;
    memset(s, 0, sizeof(*s));
    s->cls = XOMBE_CLASS_SILENCE;
    s->prepared = 1;
    s->L = INIT_SET_L;
    s->f0 = INIT_SET_F0;
    s->e[0] = XOMBE_AMP_OFF;
    for (l = 1; l <= XOMBE_LMAX; l++)
        s->a[l] = s->e[l] = INIT_SET_LOG2_AMP;
}

XOMBE_STATIC_ASSERT(sizeof(struct xombe_dec) <= XOMBE_DEC_SIZE, decoder_state_fits_XOMBE_DEC_SIZE);

size_t xombe_dec_size(void)
{
    return sizeof(struct xombe_dec);
}

void xombe__dec_core_reset(xombe_dec_core *d)
{
    memset(d, 0, sizeof(*d));
    d->rules = &ambe2_rules;
    xombe__pred_init(&d->pi);
    d->lambda = LAMBDA_INIT;
    d->avg_f0 = AVG_F0_INIT;
    d->frame = -1;
    d->mute_src = -1;
    d->decision.mute_src = -1;
    d->decision.e0 = d->decision.e1 = -1;
    xombe__dec_silence_set(&d->mute);
    xombe__dec_silence_set(&d->played);
    d->p1 = d->p2 = d->played;
}

void xombe__dec_reset(xombe_dec *d)
{
    memset(d, 0, sizeof(*d));
    xombe__dec_core_reset(&d->c);
    d->prev_f0 = INIT_SET_F0;
    d->noise_mode = XOMBE_NOISE_PRODUCTION;
    d->prng = DEFAULT_PRNG_SEED;
    xombe__noise_init(d);
    xombe__hpf_init(&d->hpf, 1.0f);
    xombe__fec_init_dmr(&d->fec);
}

static void dmr_start_bits(xombe_dec *d, const uint8_t bits[7], unsigned flags)
{
    d->c.eps_t = d->err_t;
    d->c.eps_1 = d->err_1;
    xombe__dec_start_frame(&d->c, bits, flags);
}

static const struct xombe_dec_ops dmr_ops = {dmr_start_bits};

xombe_dec *xombe__dec_init_dmr(void *mem)
{
    xombe_dec *d = mem;
    if (!d)
        return NULL;
    xombe__dec_reset(d);
    d->ops = &dmr_ops;
    return d;
}

xombe_dec *xombe__dec_init_dstar(void *mem)
{
#if XOMBE_HAVE_DSTAR
    xombe_dec *d = mem;
    if (!d)
        return NULL;
    xombe__dec_reset(d);
    xombe__fec_init_dstar(&d->fec);
    xombe__dsd_logic_init(&d->ds, &d->c);
    d->ops = &xombe__dsd_ops;
    return d;
#else
    (void)mem;
    return NULL;
#endif
}

int xombe_dec_set_noise(xombe_dec *d, enum xombe_noise_mode mode, uint32_t seed)
{
    if (mode != XOMBE_NOISE_CONFORMANCE && mode != XOMBE_NOISE_PRODUCTION)
        return -1;
    d->noise_mode = mode;
    d->prng = seed ? seed : DEFAULT_PRNG_SEED;
    xombe__noise_init(d);
    return 0;
}

void xombe__dec_set_tap(xombe_dec *d, xombe_tap_fn fn, void *ctx)
{
    d->c.tap = fn;
    d->c.tap_ctx = ctx;
}

void xombe_dec_set_errors(xombe_dec *d, int eps_t, int eps1)
{
    d->err_t = eps_t;
    d->err_1 = eps1;
}

void xombe_dec_decision_get(const xombe_dec *d, xombe_decision *out)
{
    *out = d->c.decision;
}

/* log amplitude of a at harmonic position x >= 0 */
static float harmonic_interp(const float *a, float x)
{
    int i = (int)x;
    return xombe__harmonic_lerp(a, i, x - (float)i);
}

void xombe__dec_set_from_params(xombe_set *s, int cls, float f0, int L, const float *a, const uint8_t *vb)
{
    memset(s, 0, sizeof(*s));
    s->cls = cls;
    s->f0 = f0;
    s->L = L;
    memcpy(s->a, a, sizeof(s->a));
    if (vb)
        memcpy(s->vb, vb, 8);
}

static void first_half_voice_after_voice(const xombe_dec_core *d, const xombe_set *c, const uint8_t *v1, xombe_set *p1)
{
    float f1, rc, rp, a1[XOMBE_LMAX + 1];
    int l, L1;
    if (!xombe__count_voiced(d->pi.v2))
        f1 = c->f0;
    else if (!xombe__count_voiced(c->vb))
        f1 = d->pi.f0;
    else
        f1 = sqrtf(d->pi.f0 * c->f0);
    L1 = xombe__L_from_f0(f1);
    rc = f1 / c->f0;
    rp = f1 / d->pi.f0;
    a1[0] = 0.0f;
    for (l = 1; l <= L1; l++) {
        float sc = harmonic_interp(c->a, (float)l * rc);
        float sp = harmonic_interp(d->pi.a, (float)l * rp);
        a1[l] = xombe__clamp_amp(0.5f * (sc + sp));
    }
    for (l = L1 + 1; l <= XOMBE_LMAX; l++)
        a1[l] = a1[L1];
    xombe__dec_set_from_params(p1, XOMBE_CLASS_VOICE, f1, L1, a1, v1);
}

static void decode_normal(xombe_dec_core *d, const uint8_t bits[49], int cls)
{
    xombe_set c;
    xombe_params p;
    uint8_t v2[8], v1[8];
    int b[9];
    memset(v2, 0, sizeof(v2));
    memset(v1, 0, sizeof(v1));
    xombe__get_fields(bits, b);
    if (cls == XOMBE_CLASS_VOICE) {
        xombe__dequant(b, &d->pi, 0, &p);
        XOMBE_TAP(d, XOMBE_STAGE_DEQUANT, p.f0, p.a + 1, p.L);
        xombe__v2_from_b1(b[1], v2);
        xombe__f3_first_half(d->pi.v2, v2, v1);
        xombe__smooth(d, p.a, p.L, p.f0);
        XOMBE_TAP(d, XOMBE_STAGE_SMOOTH, p.f0, p.a + 1, p.L);
        xombe__dec_set_from_params(&c, XOMBE_CLASS_VOICE, p.f0, p.L, p.a, v2);
    } else if (cls == XOMBE_CLASS_SILENCE) {
        xombe__dequant(b, &d->pi, 1, &p);
        XOMBE_TAP(d, XOMBE_STAGE_DEQUANT, p.f0, p.a + 1, p.L);
        xombe__smooth(d, p.a, p.L, p.f0);
        XOMBE_TAP(d, XOMBE_STAGE_SMOOTH, p.f0, p.a + 1, p.L);
        xombe__dec_set_from_params(&c, XOMBE_CLASS_SILENCE, p.f0, p.L, p.a, NULL);
    } else {
        int id, ad;
        xombe__tone_frame_parse(bits, &id, &ad);
        xombe__tone_set(&c, id, ad);
    }
    xombe__dec_play(d, &c, v1, &p, v2);
}

void xombe__dec_play(xombe_dec_core *d, const xombe_set *cur, const uint8_t v1[8], const xombe_params *p,
                     const uint8_t v2[8])
{
    const xombe_set c = *cur;
    int l, last = d->pi.cls;
    if (c.cls == XOMBE_CLASS_TONE || last == XOMBE_CLASS_TONE) {
        d->p1 = c;
    } else if (last == XOMBE_CLASS_SILENCE) {
        d->p1 = c;
        memcpy(d->p1.vb, v1, 8);
        if (c.cls == XOMBE_CLASS_VOICE)
            for (l = 0; l <= XOMBE_LMAX; l++)
                d->p1.a[l] = 0.5f * c.a[l];
    } else if (c.cls == XOMBE_CLASS_SILENCE) {
        if (!d->prev_concealed)
            xombe__dec_set_from_params(&d->p1, XOMBE_CLASS_VOICE, d->pi.f0, d->pi.L, d->pi.a, d->pi.v2);
        else
            d->p1 = c;
    } else {
        first_half_voice_after_voice(d, &c, v1, &d->p1);
    }
    d->p2 = c;

    if (c.cls == XOMBE_CLASS_VOICE)
        xombe__pred_store_voice(&d->pi, p, v2);
    else
        xombe__pred_store_other(&d->pi, c.cls);
}

void xombe__dec_conceal_sets(xombe_dec_core *d, int conceal)
{
    d->p1 = conceal == XOMBE_CONCEAL_MUTE ? d->mute : d->played;
    d->p1.repeat = 1;
    d->p2 = d->p1;
}

void xombe__dec_start_frame(xombe_dec_core *d, const uint8_t bytes[7], unsigned flags)
{
    uint8_t bits[49];
    int cls, bad, hb, conceal, tone_ok = 1;
    d->frame++;
    if (flags & XOMBE_FLAG_LOST) {
        memset(bits, 0, sizeof(bits));
        cls = XOMBE_CLASS_ERASURE;
        flags = 0;
    } else {
        xombe__unpack49(bytes, bits);
        cls = xombe__frame_class(bits);
    }
    if (cls == XOMBE_CLASS_TONE) {
        int id, ad;
        tone_ok = xombe__tone_frame_parse(bits, &id, &ad);
    }
    hb = (flags & XOMBE_FLAG_HB) != 0;
    bad = (flags & XOMBE_FLAG_BF) || cls == XOMBE_CLASS_ERASURE || !tone_ok;
    conceal = xombe__conceal_decide(d, bad, hb);
    d->decision.cls = cls;
    d->decision.conceal = conceal;
    d->decision.run = d->run;
    d->decision.mute_src = conceal == XOMBE_CONCEAL_MUTE ? d->mute_src : -1;
    d->decision.e0 = d->decision.e1 = -1;
    d->decision.bf = (flags & XOMBE_FLAG_BF) != 0;
    d->decision.hb = hb;
    if (conceal != XOMBE_CONCEAL_NONE)
        xombe__dec_conceal_sets(d, conceal);
    else
        decode_normal(d, bits, cls);
    d->prev_concealed = conceal != XOMBE_CONCEAL_NONE;
}

static void prepare_set(xombe_set *s)
{
    int l;
    if (s->prepared || s->repeat)
        return;
    for (l = 0; l <= XOMBE_LMAX; l++)
        s->e[l] = l >= 1 && l <= s->L ? s->a[l] : XOMBE_AMP_OFF;
    if (s->cls == XOMBE_CLASS_VOICE) {
        xombe__phases(s->a, s->L, s->phi);
        xombe__enhance(s->a, s->L, s->f0, s->e);
    } else {
        memset(s->phi, 0, sizeof(s->phi));
    }
    s->prepared = 1;
}

static void synth_half(xombe_dec *d, xombe_set *s, int first_half, int N, int16_t *pcm)
{
    float M[XOMBE_LMAX + 1], x[XOMBE_HALF_MAX];
    uint8_t vh[XOMBE_LMAX + 1];
    int n, l, tone = s->cls == XOMBE_CLASS_TONE, voiced;
    prepare_set(s);
    if (!tone && d->c.tap) {
        float lin[XOMBE_LMAX + 1];
        for (l = 1; l <= s->L; l++)
            lin[l] = xombe__exp2f(s->e[l]);
        XOMBE_TAP(&d->c, XOMBE_STAGE_SET, s->f0, s->a + 1, s->L);
        XOMBE_TAP(&d->c, XOMBE_STAGE_ENHANCED, s->f0, lin + 1, s->L);
    }
    if (s->cls == XOMBE_CLASS_VOICE)
        XOMBE_TAP(&d->c, XOMBE_STAGE_PHASE, s->f0, s->phi + 1, s->L);
    if (s->cls == XOMBE_CLASS_SILENCE && !s->repeat) {
        d->c.mute = *s;
        d->c.mute_src = d->c.frame;
    }
    if (tone) {
        for (l = 0; l <= XOMBE_LMAX; l++)
            vh[l] = (uint8_t)(l >= 1 && l <= s->L);
    } else {
        d->c.rules->harm_voicing(s->vb, s->f0, s->L, vh);
    }
    memset(M, 0, sizeof(M));
    if (tone) {
        for (l = 1; l <= s->L; l++)
            if (s->e[l] > XOMBE_AMP_OFF)
                M[l] = xombe__exp2f(s->e[l]);
    } else {
        if (first_half && xombe__count_voiced(s->vb) >= 4)
            d->c.avg_f0 = 0.9f * d->c.avg_f0 + 0.1f * s->f0;
        xombe__shape(s->e, M, s->L, s->f0, d->c.avg_f0);
    }
    XOMBE_TAP(&d->c, XOMBE_STAGE_SHAPED, s->f0, M + 1, s->L);
    memset(x, 0, sizeof(x));
    xombe__synth_unvoiced(d, s, M, vh, N, x);
    voiced = tone || xombe__count_voiced(s->vb) > 0;
    xombe__synth_voiced(d, s, M, vh, voiced, N, x);
    XOMBE_TAP(&d->c, XOMBE_STAGE_SYNTH, s->f0, x, N);
    for (n = 0; n < N; n++) {
        float y = floorf(xombe__hpf_step(&d->hpf, x[n]) + 0.5f);
        if (y > 32767.0f)
            y = 32767.0f;
        if (y < -32768.0f)
            y = -32768.0f;
        pcm[n] = (int16_t)y;
    }
}

void xombe__dec_synth(xombe_dec *d, int second_half, int N, int16_t *pcm)
{
    if (!second_half) {
        synth_half(d, &d->c.p1, 1, N, pcm);
    } else {
        synth_half(d, &d->c.p2, 0, N, pcm);
        d->c.played = d->c.p2;
    }
}

void xombe__dec_skip_synth(xombe_dec_core *d)
{
    const xombe_set *h[2];
    int i;
    h[0] = &d->p1;
    h[1] = &d->p2;
    for (i = 0; i < 2; i++)
        if (h[i]->cls == XOMBE_CLASS_SILENCE && !h[i]->repeat) {
            d->mute = *h[i];
            d->mute_src = d->frame;
        }
    d->played = d->p2;
}

int xombe__dec_half_begin(xombe_dec *d, int second_half, int n)
{
    if ((second_half != 0) != d->second_half_next)
        return -1;
    d->second_half_next = !second_half;
    return n < XOMBE_HALF_MIN ? XOMBE_HALF_MIN : n > XOMBE_HALF_MAX ? XOMBE_HALF_MAX : n;
}

/* the codec's start_bits of a frame through the decoder's FEC */
static void start_frame(xombe_dec *d, const uint8_t in[9], unsigned flags)
{
    uint8_t bits[7];
    xombe_fec_stats st;
    memset(bits, 0, sizeof(bits));
    if (flags & XOMBE_FLAG_LOST) {
        d->ops->start_bits(d, bits, XOMBE_FLAG_LOST);
        return;
    }
    flags = xombe_fec_decode(&d->fec, in, bits, &st);
    d->err_t = st.eps_t;
    d->err_1 = st.eps1;
    d->ops->start_bits(d, bits, flags);
    d->c.decision.e0 = st.eps0;
    d->c.decision.e1 = st.eps1;
}

int xombe_dec_half_bits(xombe_dec *d, const uint8_t bits[7], unsigned flags, int second_half, int n, int16_t *pcm)
{
    int N = xombe__dec_half_begin(d, second_half, n);
    if (N < 0)
        return -1;
    second_half = second_half != 0;
    if (!second_half)
        d->ops->start_bits(d, bits, flags);
    xombe__dec_synth(d, second_half, N, pcm);
    return N;
}

int xombe_dec_half(xombe_dec *d, const uint8_t in[9], unsigned flags, int second_half, int n, int16_t *pcm)
{
    int N = xombe__dec_half_begin(d, second_half, n);
    if (N < 0)
        return -1;
    second_half = second_half != 0;
    if (!second_half)
        start_frame(d, in, flags);
    xombe__dec_synth(d, second_half, N, pcm);
    return N;
}

int xombe_dec_frame_bits(xombe_dec *d, const uint8_t bits[XOMBE_BITS_BYTES], unsigned flags, int16_t pcm[160])
{
    int fresh = xombe_dec_half_bits(d, bits, flags, 0, 80, pcm) >= 0;
    if (!fresh)
        memset(pcm, 0, 80 * sizeof(*pcm));
    xombe_dec_half_bits(d, bits, flags, 1, 80, pcm + 80);
    return fresh;
}

int xombe_dec_frame(xombe_dec *d, const uint8_t in[XOMBE_FRAME_BYTES], unsigned flags, int16_t pcm[160])
{
    int fresh = xombe_dec_half(d, in, flags, 0, 80, pcm) >= 0;
    if (!fresh)
        memset(pcm, 0, 80 * sizeof(*pcm));
    xombe_dec_half(d, in, flags, 1, 80, pcm + 80);
    return fresh;
}

void xombe__dec_start_air(xombe_dec_core *d, xombe_fec *fec, const uint8_t air[9], unsigned flags, uint8_t bits[7])
{
    xombe_fec_stats st;
    memset(bits, 0, 7);
    if (flags & XOMBE_FLAG_LOST) {
        xombe__dec_start_frame(d, bits, XOMBE_FLAG_LOST);
        return;
    }
    flags = xombe_fec_decode(fec, air, bits, &st);
    d->eps_t = st.eps_t;
    d->eps_1 = st.eps1;
    xombe__dec_start_frame(d, bits, flags);
    d->decision.e0 = st.eps0;
    d->decision.e1 = st.eps1;
}
