#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "xombe.h"

/* public API contracts: handles, sizes, argument checks, call order, and delays, public headers only. every
   codec test runs for each codec built */

static int failures, checks;

static void check(int ok, const char *what)
{
    checks++;
    if (!ok) {
        failures++;
        fprintf(stderr, "FAIL %s\n", what);
    }
}

static const char *codec_name(enum xombe_codec c)
{
    return c == XOMBE_DSTAR ? "D-STAR" : "DMR";
}

static void checkc(int ok, enum xombe_codec c, const char *what)
{
    char msg[256];
    snprintf(msg, sizeof(msg), "%s: %s", codec_name(c), what);
    check(ok, msg);
}

static const enum xombe_codec codecs[] = {
    XOMBE_DMR,
#if XOMBE_HAVE_DSTAR
    XOMBE_DSTAR,
#endif
};
#define NCODECS ((int)(sizeof(codecs) / sizeof(codecs[0])))

static const uint8_t dmr_silence_bits[XOMBE_BITS_BYTES] = {0xF8, 0x01, 0xA9, 0x9F, 0x8C, 0xE0, 0x80};
static const uint8_t dstar_null_frame[XOMBE_FRAME_BYTES] = {0x9E, 0x8D, 0x32, 0x88, 0x26, 0x1A, 0x3F, 0x61, 0xE8};
static const uint8_t dstar_null_bits[XOMBE_BITS_BYTES] = {0xF8, 0x54, 0xAA, 0x9D, 0xE8, 0x0A, 0x00};
static const uint8_t fill[XOMBE_FRAME_BYTES] = {0xA5, 0xA5, 0xA5, 0xA5, 0xA5, 0xA5, 0xA5, 0xA5, 0xA5};

/* each codec's silence frame (DMR silence, D-STAR null frame) and its voice bits */
static uint8_t silence_frame[2][XOMBE_FRAME_BYTES];
static const uint8_t *const silence_bits[2] = {dmr_silence_bits, dstar_null_bits};

static void setup_frames(void)
{
    uint8_t f[XOMBE_FRAME_BYTES];
    xombe_fec_encode(XOMBE_DMR, dmr_silence_bits, silence_frame[XOMBE_DMR]);
    memcpy(silence_frame[XOMBE_DSTAR], dstar_null_frame, XOMBE_FRAME_BYTES);
    memcpy(f, fill, sizeof(f));
    check(xombe_fec_encode(XOMBE_DSTAR, dstar_null_bits, f) == (XOMBE_HAVE_DSTAR ? 0 : -1) &&
              !memcmp(f, XOMBE_HAVE_DSTAR ? dstar_null_frame : fill, sizeof(f)),
          "D-STAR null frame from its voice bits, refused without D-STAR");
}

/* pseudo-random frame n: every tenth the silence frame, every third a clean one from random voice bits */
static void stream_frame(enum xombe_codec c, int n, uint8_t fr[XOMBE_FRAME_BYTES])
{
    int i;
    for (i = 0; i < XOMBE_FRAME_BYTES; i++)
        fr[i] = (uint8_t)((n * 131 + i * 71 + (n * i) % 13) & 0xff);
    if (n % 10 == 3) {
        memcpy(fr, silence_frame[c], XOMBE_FRAME_BYTES);
    } else if (n % 3 == 0) {
        uint8_t bits[XOMBE_BITS_BYTES];
        for (i = 0; i < XOMBE_BITS_BYTES; i++)
            bits[i] = (uint8_t)(n * 37 + i * 101);
        bits[0] &= 0x7f;
        bits[6] &= c == XOMBE_DMR ? 0x80 : 0x00;
        xombe_fec_encode(c, bits, fr);
    }
}

static int popcount_bytes(const uint8_t *a, const uint8_t *b, int n)
{
    int i, k = 0;
    for (i = 0; i < n; i++) {
        unsigned v = (unsigned)(a[i] ^ b[i]);
        while (v) {
            k += v & 1;
            v >>= 1;
        }
    }
    return k;
}

static void test_fec(enum xombe_codec c)
{
    xombe_fec f;
    xombe_fec_stats st;
    uint8_t bits[XOMBE_BITS_BYTES], back[XOMBE_BITS_BYTES], fr[XOMBE_FRAME_BYTES];
    int n, i, round = 1, flipped = 1, rate = 1;
    checkc(xombe_fec_init(&f, c) == 0, c, "fec init");
    for (n = 0; n < 72; n++) {
        int diff;
        for (i = 0; i < XOMBE_BITS_BYTES; i++)
            bits[i] = (uint8_t)(n * 53 + i * 29 + n * i);
        bits[6] &= c == XOMBE_DMR ? 0x80 : 0x00;
        checkc(xombe_fec_encode(c, bits, fr) == 0, c, "fec encode");
        round &= xombe_fec_decode(&f, fr, back, n % 2 ? NULL : &st) == 0 && !memcmp(bits, back, sizeof(bits));
        if (n % 2 == 0)
            round &= st.eps0 == 0 && st.eps1 == 0 && st.eps_t == 0 && !st.bf && !st.hb;
        fr[n / 8] ^= (uint8_t)(0x80 >> (n % 8));
        xombe_fec_decode(&f, fr, back, &st);
        diff = popcount_bytes(bits, back, XOMBE_BITS_BYTES);
        flipped &= diff <= 1 && st.eps_t <= 1 && (diff == 0 || st.eps_t == 0) && st.eps_t == st.eps0 + st.eps1;
        rate &= st.eps_r > 0;
    }
    checkc(round, c, "fec round trip of clean frames, no corrections, NULL stats accepted");
    checkc(flipped, c, "any single flipped bit is corrected, or reaches the voice bits as one bit");
    checkc(rate, c, "the running error rate counts the corrections");
}

static void test_decoder(enum xombe_codec c)
{
    const uint8_t *silence = silence_frame[c];
    static xombe_dec_mem mem;
    int16_t pcm[XOMBE_HALF_MAX];
    xombe_dec *d, *h;
    void *heap;

    checkc(xombe_dec_init(NULL, c) == NULL, c, "decoder init rejects NULL");
    checkc(xombe_dec_size() <= XOMBE_DEC_SIZE, c, "decoder size within XOMBE_DEC_SIZE");
    d = xombe_dec_init(&mem, c);
    checkc(d == (xombe_dec *)&mem, c, "decoder handle is the caller's memory");
    {
        xombe_decision r;
        xombe_dec_decision_get(d, &r);
        checkc(r.mute_src == -1 && r.e0 == -1 && r.e1 == -1 && r.run == 0, c, "decision before the first frame");
    }

    checkc(xombe_dec_set_noise(d, (enum xombe_noise_mode)7, 1) == -1, c, "set_noise rejects an unknown mode");
    checkc(xombe_dec_set_noise(d, XOMBE_NOISE_PRODUCTION, 1) == 0, c, "set_noise accepts production");
    checkc(xombe_dec_set_noise(d, XOMBE_NOISE_CONFORMANCE, 0) == 0, c, "set_noise accepts conformance");
    {
        static xombe_dec_mem m[3];
        xombe_dec *def = xombe_dec_init(&m[0], c), *prod = xombe_dec_init(&m[1], c);
        xombe_dec *conf = xombe_dec_init(&m[2], c);
        int16_t a[160], b[160], cc[160];
        int f, same_prod = 1, same_conf = 1;
        xombe_dec_set_noise(prod, XOMBE_NOISE_PRODUCTION, 0);
        xombe_dec_set_noise(conf, XOMBE_NOISE_CONFORMANCE, 0);
        for (f = 0; f < 4; f++) {
            xombe_dec_frame(def, silence, 0, a);
            xombe_dec_frame(prod, silence, 0, b);
            xombe_dec_frame(conf, silence, 0, cc);
            same_prod &= !memcmp(a, b, sizeof(a));
            same_conf &= !memcmp(a, cc, sizeof(a));
        }
        checkc(same_prod && !same_conf, c, "default noise source is production, seed 0");
    }

    d = xombe_dec_init(&mem, c);
    checkc(xombe_dec_half(d, silence, 0, 1, 80, pcm) == -1, c, "second half before a first half is refused");
    checkc(xombe_dec_half(d, silence, 0, 0, 200, pcm) == XOMBE_HALF_MAX, c, "first half, length clamped up to 84");
    checkc(xombe_dec_half(d, silence, 0, 0, 80, pcm) == -1, c, "two first halves in a row are refused");
    checkc(xombe_dec_half(d, silence, 0, 1, 10, pcm) == XOMBE_HALF_MIN, c, "second half, length clamped to 76");
    checkc(xombe_dec_half(d, silence, 0, 1, 80, pcm) == -1, c, "two second halves in a row are refused");
    checkc(xombe_dec_half(d, silence, 0, 0, 80, pcm) == 80 && xombe_dec_half(d, NULL, 0, 1, 80, pcm) == 80, c,
           "second half reads no frame");
    checkc(xombe_dec_half_bits(d, NULL, 0, 1, 80, pcm) == -1, c, "bits second half before a first half is refused");
    checkc(xombe_dec_half_bits(d, silence_bits[c], 0, 0, 80, pcm) == 80 &&
               xombe_dec_half_bits(d, NULL, 0, 1, 80, pcm) == 80,
           c, "bits second half reads no bits");
    {
        static xombe_dec_mem lone_mem, twin_mem;
        xombe_dec *lone = xombe_dec_init(&lone_mem, c), *twin = xombe_dec_init(&twin_mem, c);
        int16_t full[160], want[160];
        int i, zeros = 1;
        xombe_dec_half(lone, silence, 0, 0, 80, pcm);
        xombe_dec_half(twin, silence, 0, 0, 80, want);
        xombe_dec_half(twin, NULL, 0, 1, 80, want + 80);
        for (i = 0; i < 160; i++)
            full[i] = 0x5555;
        checkc(xombe_dec_frame(lone, silence, 0, full) == 0, c, "frame after an unpaired first half returns 0");
        for (i = 0; i < 80; i++)
            zeros &= full[i] == 0;
        checkc(zeros && !memcmp(full + 80, want + 80, 80 * sizeof(*full)), c,
               "frame after an unpaired first half: zeros, then the pending second half");
        checkc(xombe_dec_frame(lone, silence, 0, full) == 1, c, "the next frame returns 1");
        xombe_dec_frame(twin, silence, 0, want);
        checkc(!memcmp(full, want, sizeof(full)), c, "the next frame decodes normally");
        xombe_dec_half_bits(lone, silence_bits[c], 0, 0, 80, pcm);
        checkc(xombe_dec_frame_bits(lone, silence_bits[c], 0, full) == 0, c,
               "bits frame after an unpaired first half returns 0");
    }
    {
        xombe_decision r;
        int16_t f[160];
        xombe_dec_frame(d, silence, 0, f);
        xombe_dec_decision_get(d, &r);
        checkc(r.cls == XOMBE_CLASS_SILENCE && r.conceal == XOMBE_CONCEAL_NONE && r.e0 == 0 && r.e1 == 0, c,
               "decision of a silence frame");
    }
    {
        /* a lost frame plays as a bad frame, without reading bits. D-STAR's bad voice frames never mute, so its
           run stays short of the mute */
        static xombe_dec_mem lm, bm;
        xombe_dec *l = xombe_dec_init(&lm, c), *b = xombe_dec_init(&bm, c);
        static const uint8_t zero[XOMBE_BITS_BYTES];
        int16_t pl[160], pb[160];
        int i, same = 1, erasure = 1, last_lost = c == XOMBE_DMR ? 6 : 3;
        for (i = 0; i < 8; i++) {
            xombe_decision r;
            unsigned lost = i >= 2 && i <= last_lost;
            xombe_dec_frame_bits(l, lost ? NULL : silence_bits[c], lost ? XOMBE_FLAG_LOST : 0, pl);
            xombe_dec_frame_bits(b, lost ? zero : silence_bits[c], lost ? XOMBE_FLAG_BF : 0, pb);
            same &= !memcmp(pl, pb, sizeof(pl));
            xombe_dec_decision_get(l, &r);
            if (lost)
                erasure &= r.cls == XOMBE_CLASS_ERASURE && r.conceal != XOMBE_CONCEAL_NONE && r.run == i - 1 &&
                           r.e0 == -1 && r.e1 == -1 && !r.bf && !r.hb;
        }
        checkc(same, c, "XOMBE_FLAG_LOST plays as XOMBE_FLAG_BF");
        checkc(erasure, c, "decision of a lost frame: erasure, concealed, run counts it");
    }

    heap = malloc(xombe_dec_size());
    h = xombe_dec_init(heap, c);
    if (h) {
        int16_t a[160], b[160];
        int i;
        d = xombe_dec_init(&mem, c);
        for (i = 0; i < 20; i++) {
            xombe_dec_frame(d, silence, 0, a);
            xombe_dec_frame(h, silence, 0, b);
        }
        checkc(!memcmp(a, b, sizeof(a)), c, "malloc'd and static decoders agree");
    }
    checkc(h != NULL, c, "decoder in malloc'd memory");
    free(heap);
}

/* frame calls equal the caller's own FEC into the *_bits calls, and two half calls */
static void test_decoder_bits(enum xombe_codec c)
{
    static xombe_dec_mem am, bm, hm;
    xombe_dec *a = xombe_dec_init(&am, c), *b = xombe_dec_init(&bm, c), *h = xombe_dec_init(&hm, c);
    xombe_fec fec;
    int n, same = 1, halves = 1, decided = 1, seen = 0;
    xombe_fec_init(&fec, c);
    for (n = 0; n < 300; n++) {
        uint8_t fr[XOMBE_FRAME_BYTES], bits[XOMBE_BITS_BYTES];
        int16_t pa[160], pb[160], ph[160];
        unsigned lost = n % 41 >= 7 && n % 41 < 11 ? XOMBE_FLAG_LOST : 0, flags = XOMBE_FLAG_LOST;
        xombe_fec_stats st;
        xombe_decision ra, rb;
        stream_frame(c, n, fr);
        same &= xombe_dec_frame(a, lost ? NULL : fr, lost, pa) == 1;
        if (!lost) {
            flags = xombe_fec_decode(&fec, fr, bits, &st);
            xombe_dec_set_errors(b, st.eps_t, st.eps1);
        }
        same &= xombe_dec_frame_bits(b, lost ? NULL : bits, flags, pb) == 1;
        halves &= xombe_dec_half(h, lost ? NULL : fr, lost, 0, 80, ph) == 80 &&
                  xombe_dec_half(h, NULL, 0, 1, 80, ph + 80) == 80 && !memcmp(pa, ph, sizeof(pa));
        same &= !memcmp(pa, pb, sizeof(pa));
        xombe_dec_decision_get(a, &ra);
        xombe_dec_decision_get(b, &rb);
        decided &= ra.cls == rb.cls && ra.conceal == rb.conceal && ra.run == rb.run && ra.mute_src == rb.mute_src &&
                   ra.bf == rb.bf && ra.hb == rb.hb && rb.e0 == -1 && rb.e1 == -1;
        decided &= lost ? ra.e0 == -1 && ra.e1 == -1 && !ra.bf && !ra.hb
                        : ra.e0 == st.eps0 && ra.e1 == st.eps1 && ra.bf == st.bf && ra.hb == st.hb;
        if (ra.conceal == XOMBE_CONCEAL_NONE)
            seen |= 1 << ra.cls;
        else
            seen |= 1 << (4 + ra.conceal);
    }
    checkc(same, c, "xombe_dec_frame equals xombe_fec_decode, xombe_dec_set_errors, and xombe_dec_frame_bits");
    checkc(halves, c, "xombe_dec_frame equals two xombe_dec_half calls");
    checkc(decided, c, "frame decisions: the bits one with the Golay corrections");
    checkc((seen & 0x65) == 0x65, c, "the stream plays voice, silence, repeats, and mutes");
}

/* tone frames: same ID and level rules as forced tones, refused calls write nothing, decoded as tones */
static void test_tone_frames(enum xombe_codec c)
{
    static const int bad[][2] = {{0, 90}, {4, 90}, {123, 90}, {127, 90}, {164, 90}, {-1, 90}, {255, 90},
                                 {5, -1}, {5, 128}, {129, 1000}};
    static const int good[][2] = {{5, 0}, {122, 127}, {128, 64}, {163, 1}};
    static xombe_dec_mem dm;
    xombe_dec *d = xombe_dec_init(&dm, c);
    uint8_t f[XOMBE_FRAME_BYTES];
    int16_t pcm[160];
    unsigned i;
    int refused = 1, accepted = 1, knox = 1;
    for (i = 0; i < sizeof(bad) / sizeof(bad[0]); i++) {
        memcpy(f, fill, sizeof(f));
        refused &= xombe_tone_frame(c, bad[i][0], bad[i][1], f) == -1 && !memcmp(f, fill, sizeof(f));
    }
    for (i = 0; i < sizeof(good) / sizeof(good[0]); i++) {
        xombe_decision r;
        accepted &= xombe_tone_frame(c, good[i][0], good[i][1], f) == 0;
        xombe_dec_frame(d, f, 0, pcm);
        xombe_dec_decision_get(d, &r);
        accepted &= r.cls == XOMBE_CLASS_TONE && r.conceal == XOMBE_CONCEAL_NONE && r.e0 == 0 && r.e1 == 0;
    }
    for (i = 144; i <= 159; i++) {
        int r;
        memcpy(f, fill, sizeof(f));
        r = xombe_tone_frame(c, (int)i, 90, f);
        knox &= c == XOMBE_DMR ? r == 0 : r == -1 && !memcmp(f, fill, sizeof(f));
    }
    checkc(refused, c, "xombe_tone_frame refuses IDs outside 5..122, 128..163 and levels outside 0..127");
    checkc(accepted, c, "xombe_tone_frame takes range ends, its frames decode as clean tones");
    checkc(knox, c, "Knox tone frames: DMR has them, D-STAR refuses them and writes nothing");
}

/* calls naming a codec the build lacks, or none at all */
static void test_codecs(void)
{
    static xombe_dec_mem dm;
    xombe_fec f, g;
    uint8_t fr[XOMBE_FRAME_BYTES];
    memset(&f, 0x5A, sizeof(f));
    g = f;
    check(xombe_dec_init(&dm, (enum xombe_codec)7) == NULL, "decoder init refuses an unknown codec");
    check(xombe_fec_init(&f, (enum xombe_codec)7) == -1 && !memcmp(&f, &g, sizeof(f)),
          "fec init refuses an unknown codec, f unchanged");
    memcpy(fr, fill, sizeof(fr));
    check(xombe_fec_encode((enum xombe_codec)7, dmr_silence_bits, fr) == -1 && !memcmp(fr, fill, sizeof(fr)),
          "fec encode refuses an unknown codec, writes nothing");
    check(xombe_tone_frame((enum xombe_codec)7, 128, 64, fr) == -1 && !memcmp(fr, fill, sizeof(fr)),
          "tone frame refuses an unknown codec, writes nothing");
    check((xombe_dec_init(&dm, XOMBE_DSTAR) != NULL) == XOMBE_HAVE_DSTAR, "D-STAR decoder only in a D-STAR build");
    check((xombe_fec_init(&f, XOMBE_DSTAR) == 0) == XOMBE_HAVE_DSTAR, "D-STAR FEC only in a D-STAR build");
    check((xombe_tone_frame(XOMBE_DSTAR, 128, 64, fr) == 0) == XOMBE_HAVE_DSTAR,
          "D-STAR tone frames only in a D-STAR build");
#if XOMBE_HAVE_DSTAR
    {
        xombe_fec_stats st;
        uint8_t bits[XOMBE_BITS_BYTES];
        xombe_fec_init(&f, XOMBE_DSTAR);
        check(xombe_fec_decode(&f, dstar_null_frame, bits, &st) == 0 && !memcmp(bits, dstar_null_bits, 7),
              "D-STAR null frame's voice bits");
    }
#endif
}

#if XOMBE_HAVE_ENCODER
/* harmonic signal that comes and goes */
static void test_signal(int n, int16_t pcm[160])
{
    int i;
    for (i = 0; i < 160; i++) {
        float t = (float)(n * 160 + i), x = 0.0f;
        int h;
        if ((n / 20) % 2 == 0)
            for (h = 1; h <= 8; h++)
                x += 1500.0f / (float)h * sinf(0.0393f * (float)h * t + 0.3f * (float)h);
        pcm[i] = (int16_t)(x + (float)((n * 160 + i) * 7919 % 61) - 30.0f);
    }
}

#if XOMBE_HAVE_DSTAR
#define ALL_OUTPUTS (XOMBE_ENC_DMR | XOMBE_ENC_DSTAR)
#else
#define ALL_OUTPUTS XOMBE_ENC_DMR
#endif

static void test_encoder(void)
{
    static xombe_enc_mem mem;
    xombe_enc_cfg cfg;
    int16_t pcm[160];
    uint8_t f[XOMBE_FRAME_BYTES];
    int i, got = 0;

    check(xombe_enc_size() <= XOMBE_ENC_SIZE, "encoder size within XOMBE_ENC_SIZE");
    check(xombe_enc_init(NULL, NULL) == NULL, "encoder init rejects NULL");
    check(xombe_enc_init(&mem, NULL) != NULL, "encoder init with default settings");

    xombe_enc_cfg_default(&cfg);
    cfg.ns_strength_permille = -1;
    check(xombe_enc_init(&mem, &cfg) == NULL, "ns_strength_permille below 0 is refused");
    cfg.ns_strength_permille = 1001;
    check(xombe_enc_init(&mem, &cfg) == NULL, "ns_strength_permille above 1000 is refused");
    cfg.ns_strength_permille = 0;
    check(xombe_enc_init(&mem, &cfg) != NULL, "ns_strength_permille 0 is accepted");
    cfg.ns_strength_permille = 1000;
    check(xombe_enc_init(&mem, &cfg) != NULL, "ns_strength_permille 1000 is accepted");
    xombe_enc_cfg_default(&cfg);
    cfg.voicing_scale_permille = 0;
    check(xombe_enc_init(&mem, &cfg) == NULL, "voicing_scale_permille 0 is refused");
    cfg.voicing_scale_permille = -350;
    check(xombe_enc_init(&mem, &cfg) == NULL, "voicing_scale_permille below 0 is refused");
    xombe_enc_cfg_default(&cfg);
    cfg.dmr_dtx = 5;
    cfg.dstar_dtx = 5;
    check(xombe_enc_init(&mem, &cfg) != NULL, "nonzero dtx means on");

    memset(pcm, 0, sizeof(pcm));
    xombe_enc_init(&mem, NULL);
    for (i = 0; i < XOMBE_ENC_DELAY + 2; i++)
        got += xombe_enc_frame((xombe_enc *)&mem, pcm, NULL, f, NULL);
    check(got == 2, "encoder output starts after XOMBE_ENC_DELAY frames");

    xombe_enc_cfg_default(&cfg);
    cfg.outputs = 0x80u;
    check(xombe_enc_init(&mem, &cfg) == NULL, "an unknown output is refused");
    cfg.outputs = 0;
    check(xombe_enc_init(&mem, &cfg) != NULL, "outputs 0 means DMR");
    check(xombe_enc_delay((xombe_enc *)&mem) == XOMBE_ENC_DELAY, "DMR-only delay");
    cfg.outputs = XOMBE_ENC_DSTAR;
    check((xombe_enc_init(&mem, &cfg) != NULL) == XOMBE_HAVE_DSTAR, "D-STAR output only in a D-STAR build");
    xombe_enc_cfg_default(&cfg);
    check(cfg.outputs == XOMBE_ENC_DMR && cfg.dmr_dtx && !cfg.dstar_dtx && cfg.dmr_quant == XOMBE_ENC_QUANT_CLOSED,
          "defaults: DMR output, DMR DTX on, D-STAR DTX off, closed loop");

    {
        /* pcm NULL is a frame of zeros, for both outputs and both forms */
        static xombe_enc_mem m2, m3, m4;
        xombe_enc *a, *b, *c2, *d2;
        uint8_t fa[2][XOMBE_FRAME_BYTES], fb[2][XOMBE_FRAME_BYTES], ba[2][XOMBE_BITS_BYTES];
        uint8_t bb[2][XOMBE_BITS_BYTES];
        int same = 1;
        xombe_enc_cfg_default(&cfg);
        cfg.outputs = ALL_OUTPUTS;
        a = xombe_enc_init(&mem, &cfg);
        b = xombe_enc_init(&m2, &cfg);
        c2 = xombe_enc_init(&m3, &cfg);
        d2 = xombe_enc_init(&m4, &cfg);
        for (i = 0; i < 8; i++) {
            int ra = xombe_enc_frame(a, NULL, NULL, fa[0], fa[1]), rb = xombe_enc_frame(b, pcm, NULL, fb[0], fb[1]);
            int rc = xombe_enc_frame_bits(c2, NULL, NULL, ba[0], ba[1]);
            int rd = xombe_enc_frame_bits(d2, pcm, NULL, bb[0], bb[1]);
            same &= ra == rb && rc == rd && !memcmp(fa, fb, sizeof(fa)) && !memcmp(ba, bb, sizeof(ba));
        }
        check(same, "pcm NULL encodes zeros");
    }
}

/* an output the encoder does not produce gets its codec's silence frame, in either form */
static void test_encoder_absent(void)
{
    static xombe_enc_mem mem;
    int16_t pcm[160];
    uint8_t f[XOMBE_FRAME_BYTES], b[XOMBE_BITS_BYTES];
    int n, ok = 1;
    xombe_enc_init(&mem, NULL);
    for (n = 0; n < 10; n++) {
        test_signal(n, pcm);
        xombe_enc_frame((xombe_enc *)&mem, pcm, NULL, NULL, f);
        xombe_enc_frame_bits((xombe_enc *)&mem, pcm, NULL, NULL, b);
        ok &= !memcmp(f, dstar_null_frame, sizeof(f)) && !memcmp(b, dstar_null_bits, sizeof(b));
    }
    check(ok, "a DMR-only encoder writes D-STAR null frames and their bits");
#if XOMBE_HAVE_DSTAR
    {
        xombe_enc_cfg cfg;
        xombe_enc_cfg_default(&cfg);
        cfg.outputs = XOMBE_ENC_DSTAR;
        xombe_enc_init(&mem, &cfg);
        for (n = 0, ok = 1; n < 10; n++) {
            test_signal(n, pcm);
            xombe_enc_frame((xombe_enc *)&mem, pcm, NULL, f, NULL);
            xombe_enc_frame_bits((xombe_enc *)&mem, pcm, NULL, b, NULL);
            ok &= !memcmp(f, silence_frame[XOMBE_DMR], sizeof(f)) && !memcmp(b, dmr_silence_bits, sizeof(b));
        }
        check(ok, "a D-STAR-only encoder writes DMR silence frames and their bits");
    }
#endif
}

/* frames are the voice bits through the FEC, for both codecs, through the delay fill, tones, and voice */
static void test_encoder_bits(void)
{
    static const xombe_enc_tone tones[] = {{128, 64}, {150, 90}, {40, 100}};
    static xombe_enc_mem ma, mb;
    xombe_enc_cfg cfg;
    xombe_enc *a, *b;
    int16_t pcm[160];
    int n, k, same = 1;
    xombe_enc_cfg_default(&cfg);
    cfg.outputs = ALL_OUTPUTS;
    a = xombe_enc_init(&ma, &cfg);
    b = xombe_enc_init(&mb, &cfg);
    for (n = 0; n < 100; n++) {
        uint8_t f[2][XOMBE_FRAME_BYTES], bits[2][XOMBE_BITS_BYTES], g[XOMBE_FRAME_BYTES];
        const xombe_enc_tone *t = n % 17 == 5 ? &tones[(n / 17) % 3] : NULL;
        test_signal(n, pcm);
        same &= xombe_enc_frame(a, pcm, t, f[0], f[1]) == xombe_enc_frame_bits(b, pcm, t, bits[0], bits[1]);
        for (k = 0; k < NCODECS; k++) {
            xombe_fec_encode(codecs[k], bits[codecs[k]], g);
            same &= !memcmp(g, f[codecs[k]], sizeof(g));
        }
    }
    check(same, "xombe_enc_frame equals xombe_fec_encode of xombe_enc_frame_bits, both codecs");
}

/* forced tone requests: ID and level ranges, refused calls write and consume nothing */
static void test_encoder_tone(void)
{
    static const xombe_enc_tone bad[] = {{0, 90}, {4, 90}, {123, 90}, {127, 90}, {164, 90}, {-1, 90},
                                         {255, 90}, {5, -1}, {5, 128}, {129, 1000}};
    static const xombe_enc_tone good[] = {{5, 0}, {122, 127}, {128, 64}, {163, 1}, {144, 90}, {159, 90}};
    static xombe_enc_mem ma, mb;
    xombe_enc_cfg cfg;
    xombe_enc *a, *b;
    int16_t pcm[160];
    uint8_t fa[2][XOMBE_FRAME_BYTES], fb[2][XOMBE_FRAME_BYTES], ba[2][XOMBE_BITS_BYTES];
    unsigned i, n;
    int refused = 1, untouched = 1, same = 1, accepted = 1;
    xombe_enc_cfg_default(&cfg);
    cfg.outputs = ALL_OUTPUTS;
    a = xombe_enc_init(&ma, &cfg);
    b = xombe_enc_init(&mb, &cfg);
    for (n = 0; n < 40; n++) {
        for (i = 0; i < 160; i++)
            pcm[i] = (int16_t)(3000.0f * sinf(0.05f * (float)(n * 160 + i)) + (float)((n * 160 + i) * 7919 % 61));
        if (n % 4 == 0) {
            for (i = 0; i < sizeof(bad) / sizeof(bad[0]); i++) {
                memset(fa, 0xA5, sizeof(fa));
                memset(ba, 0xA5, sizeof(ba));
                refused &= xombe_enc_frame(a, pcm, &bad[i], fa[0], fa[1]) == -1 &&
                           xombe_enc_frame_bits(a, pcm, &bad[i], ba[0], ba[1]) == -1;
                untouched &= !memcmp(fa[0], fill, 9) && !memcmp(fa[1], fill, 9) && !memcmp(ba[0], fill, 7) &&
                             !memcmp(ba[1], fill, 7);
            }
        }
        same &= xombe_enc_frame(a, pcm, NULL, fa[0], fa[1]) == xombe_enc_frame(b, pcm, NULL, fb[0], fb[1]) &&
                !memcmp(fa, fb, sizeof(fa));
    }
    check(refused, "xombe_enc_frame refuses IDs outside 5..122, 128..163 and levels outside 0..127");
    check(untouched, "a refused request writes nothing");
    check(same, "refused requests consume nothing");
    for (i = 0; i < sizeof(good) / sizeof(good[0]); i++)
        accepted &= xombe_enc_frame(a, pcm, &good[i], fa[0], fa[1]) == 1;
    check(accepted, "range ends and Knox IDs are accepted");
    check(xombe_enc_frame(a, pcm, &good[1], NULL, NULL) == 1 &&
              xombe_enc_frame_bits(a, pcm, &good[1], NULL, NULL) == 1,
          "NULL outputs are accepted");
}

/* FNV-1a 32 of the DMR voice bits of 100 test_signal frames, open loop, per kernel backend */
#if XOMBE_KERNELS_Q15
#define OPEN_LOOP_HASH 0xb7487b6du
#else
#define OPEN_LOOP_HASH 0x99dd5352u
#endif

static uint32_t fnv1a(uint32_t h, const uint8_t *p, size_t n)
{
    while (n--)
        h = (h ^ *p++) * 16777619u;
    return h;
}

/* quantizer setting: range, open loop unchanged, D-STAR output independent of it */
static void test_encoder_quant(void)
{
    static xombe_enc_mem m1, m2;
    xombe_enc_cfg cfg;
    xombe_enc *e1, *e2;
    int16_t pcm[160];
    uint8_t a1[7], a2[7];
    uint32_t h = 2166136261u;
    int n, differ = 0;

    xombe_enc_cfg_default(&cfg);
    cfg.dmr_quant = 2;
    check(xombe_enc_init(&m1, &cfg) == NULL, "an unknown quantizer is refused");
    cfg.dmr_quant = -1;
    check(xombe_enc_init(&m1, &cfg) == NULL, "a negative quantizer is refused");
    cfg.dmr_quant = XOMBE_ENC_QUANT_OPEN;
    e1 = xombe_enc_init(&m1, &cfg);
    e2 = xombe_enc_init(&m2, NULL);
    check(e1 && e2, "encoder init with each quantizer");
    if (!e1 || !e2)
        return;
    for (n = 0; n < 100; n++) {
        test_signal(n, pcm);
        if (xombe_enc_frame_bits(e1, pcm, NULL, a1, NULL) == 1)
            h = fnv1a(h, a1, 7);
        xombe_enc_frame_bits(e2, pcm, NULL, a2, NULL);
        differ += memcmp(a1, a2, 7) != 0;
    }
    if (h != OPEN_LOOP_HASH)
        fprintf(stderr, "open-loop hash 0x%08lx\n", (unsigned long)h);
    check(h == OPEN_LOOP_HASH, "open-loop output unchanged");
    check(differ > 0, "closed-loop output differs from open loop");
#if XOMBE_HAVE_DSTAR
    {
        uint8_t d1[XOMBE_FRAME_BYTES], d2[XOMBE_FRAME_BYTES], f1[XOMBE_FRAME_BYTES], f2[XOMBE_FRAME_BYTES];
        int same = 1;
        cfg.outputs = XOMBE_ENC_DMR | XOMBE_ENC_DSTAR;
        e1 = xombe_enc_init(&m1, &cfg);
        cfg.dmr_quant = XOMBE_ENC_QUANT_CLOSED;
        e2 = xombe_enc_init(&m2, &cfg);
        for (n = 0, differ = 0; n < 100; n++) {
            test_signal(n, pcm);
            xombe_enc_frame(e1, pcm, NULL, f1, d1);
            xombe_enc_frame(e2, pcm, NULL, f2, d2);
            same &= !memcmp(d1, d2, sizeof(d1));
            differ += memcmp(f1, f2, sizeof(f1)) != 0;
        }
        check(same && differ > 0, "D-STAR output of a dual encoder independent of the quantizer");
    }
#endif
}

#if XOMBE_HAVE_DSTAR
static void test_encoder_dstar(void)
{
    const uint8_t *silence = silence_frame[XOMBE_DMR];
    static xombe_enc_mem m1, m2, m3, m4;
    static uint8_t a1[100][XOMBE_FRAME_BYTES];
    xombe_enc_cfg cfg;
    xombe_enc *e1, *e2, *e3, *e4;
    int16_t pcm[160];
    uint8_t a2[XOMBE_FRAME_BYTES], a3[XOMBE_FRAME_BYTES], d2[XOMBE_FRAME_BYTES], d3[XOMBE_FRAME_BYTES];
    uint8_t d4[XOMBE_FRAME_BYTES];
    int n, same = 1, nulls = 0, dtx_nulls = 0, silent = 1, r1 = 0, r2 = 0, r3 = 0;
    int lag = XOMBE_ENC_DSTAR_DELAY - XOMBE_ENC_DELAY;

    xombe_enc_cfg_default(&cfg);
    e1 = xombe_enc_init(&m1, &cfg);
    cfg.outputs = XOMBE_ENC_DMR | XOMBE_ENC_DSTAR;
    e2 = xombe_enc_init(&m2, &cfg);
    cfg.outputs = XOMBE_ENC_DSTAR;
    e3 = xombe_enc_init(&m3, &cfg);
    cfg.dstar_dtx = 1;
    e4 = xombe_enc_init(&m4, &cfg);
    check(e1 && e2 && e3 && e4, "encoder init with each output set");
    if (!e1 || !e2 || !e3 || !e4)
        return;
    check(xombe_enc_delay(e1) == XOMBE_ENC_DELAY && xombe_enc_delay(e2) == XOMBE_ENC_DSTAR_DELAY &&
              xombe_enc_delay(e3) == XOMBE_ENC_DSTAR_DELAY && lag > 0,
          "the D-STAR output adds look-ahead to both outputs");
    for (n = 0; n < 100; n++) {
        test_signal(n, pcm);
        r1 += xombe_enc_frame(e1, pcm, NULL, a1[n], NULL);
        r2 += xombe_enc_frame(e2, pcm, NULL, a2, d2);
        r3 += xombe_enc_frame(e3, pcm, NULL, a3, d3);
        same &= !memcmp(n >= lag ? a1[n - lag] : silence, a2, XOMBE_FRAME_BYTES);
        silent &= !memcmp(a3, silence, XOMBE_FRAME_BYTES);
        if (n >= XOMBE_ENC_DSTAR_DELAY)
            nulls += !memcmp(d2, dstar_null_frame, sizeof(d2)) + !memcmp(d3, dstar_null_frame, sizeof(d3));
        if (xombe_enc_frame(e4, pcm, NULL, NULL, d4))
            dtx_nulls += !memcmp(d4, dstar_null_frame, sizeof(d4));
    }
    check(r1 == 100 - XOMBE_ENC_DELAY && r2 == 100 - XOMBE_ENC_DSTAR_DELAY && r3 == r2,
          "both outputs share XOMBE_ENC_DSTAR_DELAY");
    check(same, "with D-STAR, the DMR output is the DMR-only one, delayed");
    check(silent, "an encoder without DMR writes DMR silence frames");
    check(nulls == 0, "with D-STAR DTX off, every D-STAR frame is a voice frame");
    check(dtx_nulls > 0, "with D-STAR DTX on, inactive frames are null frames");
}
#endif
#endif

#if XOMBE_HAVE_DSTAR
/* transcoder from codec c: init, lost frames, clean output, two instances agree, and the decision of the input is
   the source decoder's */
static void test_xc(enum xombe_codec c)
{
    enum xombe_codec to = c == XOMBE_DMR ? XOMBE_DSTAR : XOMBE_DMR;
    static xombe_xc_mem x1, x2;
    static xombe_dec_mem dm;
    xombe_xc *x, *y;
    xombe_xc_decision d;
    xombe_dec *dd;
    xombe_fec fs;
    uint8_t fr[XOMBE_FRAME_BYTES], o1[XOMBE_FRAME_BYTES], o2[XOMBE_FRAME_BYTES], bits[XOMBE_BITS_BYTES];
    int16_t pcm[160];
    int n, agree = 1, clean = 1, decided = 1, lost_ok = 1, outs = 1;
    x = xombe_xc_init(&x1, c);
    if (c == XOMBE_DSTAR && !XOMBE_HAVE_ENCODER) {
        checkc(x == NULL, c, "transcoder from D-STAR needs the encoder");
        return;
    }
    y = xombe_xc_init(&x2, c);
    dd = xombe_dec_init(&dm, c);
    checkc(x == (xombe_xc *)&x1 && y != NULL, c, "transcoder handle is the caller's memory");
    checkc(xombe_xc_size() <= XOMBE_XC_SIZE, c, "transcoder size within XOMBE_XC_SIZE");
    checkc(xombe_xc_init(NULL, c) == NULL, c, "transcoder init rejects NULL");
    xombe_xc_decision_get(x, &d);
    checkc(d.in.e0 == -1 && d.in.e1 == -1 && d.in.mute_src == -1, c, "transcoder decision before the first frame");
    checkc(xombe_xc_frame(x, silence_frame[c], 0, o1) == 1 && xombe_xc_frame(y, silence_frame[c], 0, o2) == 1, c,
           "transcoder takes a silence frame");
    xombe_dec_frame(dd, silence_frame[c], 0, pcm);
    xombe_xc_decision_get(x, &d);
    checkc(d.in.cls == XOMBE_CLASS_SILENCE && d.out == (c == XOMBE_DMR ? XOMBE_XC_OUT_VOICE : XOMBE_XC_OUT_SILENCE),
           c, "a silence frame transcodes as a D-STAR all-unvoiced voice frame or a DMR silence frame");
    xombe_fec_init(&fs, to);
    for (n = 0; n < 300; n++) {
        unsigned flags = n % 37 == 5 ? XOMBE_FLAG_LOST : 0;
        xombe_fec_stats st;
        xombe_decision r;
        stream_frame(c, n, fr);
        agree &= xombe_xc_frame(x, flags ? NULL : fr, flags, o1) == 1;
        agree &= xombe_xc_frame(y, fr, flags, o2) == 1 && !memcmp(o1, o2, sizeof(o1));
        xombe_fec_decode(&fs, o1, bits, &st);
        clean &= st.eps0 == 0 && st.eps1 == 0;
        xombe_xc_decision_get(x, &d);
        if (flags)
            lost_ok &= d.in.cls == XOMBE_CLASS_ERASURE && d.in.conceal != XOMBE_CONCEAL_NONE && d.in.e0 == -1 &&
                       !d.in.bf && !d.in.hb;
        xombe_dec_frame(dd, fr, flags, pcm);
        xombe_dec_decision_get(dd, &r);
        decided &= d.in.cls == r.cls && d.in.conceal == r.conceal && d.in.run == r.run &&
                   d.in.mute_src == r.mute_src && d.in.e0 == r.e0 && d.in.e1 == r.e1 && d.in.bf == r.bf &&
                   d.in.hb == r.hb;
        if (c == XOMBE_DMR)
            outs &= d.out == XOMBE_XC_OUT_VOICE || d.out == XOMBE_XC_OUT_TONE || d.out == XOMBE_XC_OUT_KNOX ||
                    d.out == XOMBE_XC_OUT_RESEND;
        else
            outs &= d.out == XOMBE_XC_OUT_VOICE || d.out == XOMBE_XC_OUT_PULSE || d.out == XOMBE_XC_OUT_TONE ||
                    d.out == XOMBE_XC_OUT_SILENCE;
    }
    checkc(agree, c, "two transcoders agree, a lost frame's input is not read");
    checkc(clean, c, "transcoder output frames are clean codewords of the other codec");
    checkc(lost_ok, c, "decision of a lost frame: erasure, concealed, no error counts");
    checkc(decided, c, "the transcoder's input decision is the source decoder's");
    checkc(outs, c, "output classes are those of its direction");
    checkc(xombe_xc_init_shared(NULL, c, &x2) == NULL && xombe_xc_init_shared(&x2, c, NULL) == NULL, c,
           "shared init rejects NULL memory or workspace");
}

static void test_xc_codecs(void)
{
    static xombe_xc_mem m;
    check(xombe_xc_init(&m, (enum xombe_codec)7) == NULL, "transcoder init refuses an unknown codec");
#if XOMBE_HAVE_ENCODER
    {
        /* two transcoders on one workspace, interleaved, equal two with their own */
        static xombe_xc_mem c1, r1;
        static xombe_qws_mem ws;
        static xombe_xc_shared_mem c2, r2;
        xombe_xc *ca = xombe_xc_init(&c1, XOMBE_DMR), *ra = xombe_xc_init(&r1, XOMBE_DSTAR), *cb, *rb;
        uint8_t fr[XOMBE_FRAME_BYTES], o1[XOMBE_FRAME_BYTES], o2[XOMBE_FRAME_BYTES], q1[XOMBE_FRAME_BYTES];
        uint8_t q2[XOMBE_FRAME_BYTES];
        int n, same = 1;
        cb = xombe_xc_init_shared(&c2, XOMBE_DMR, &ws);
        rb = xombe_xc_init_shared(&r2, XOMBE_DSTAR, &ws);
        check(cb != NULL && rb != NULL, "shared init in XOMBE_XC_SHARED_SIZE");
        for (n = 0; cb && rb && n < 200; n++) {
            stream_frame(XOMBE_DMR, n, fr);
            xombe_xc_frame(ca, fr, 0, o1);
            xombe_xc_frame(cb, fr, 0, o2);
            xombe_xc_frame(ra, o1, 0, q1);
            xombe_xc_frame(rb, o2, 0, q2);
            same &= !memcmp(o1, o2, sizeof(o1)) && !memcmp(q1, q2, sizeof(q1));
        }
        check(same, "transcoders sharing a workspace equal transcoders with their own");
    }
#endif
}
#endif

int main(void)
{
    int k;
    check(xombe_version() == XOMBE_VERSION, "library version matches the header");
    setup_frames();
    test_codecs();
    for (k = 0; k < NCODECS; k++) {
        test_fec(codecs[k]);
        test_decoder(codecs[k]);
        test_decoder_bits(codecs[k]);
        test_tone_frames(codecs[k]);
#if XOMBE_HAVE_DSTAR
        test_xc(codecs[k]);
#endif
    }
#if XOMBE_HAVE_DSTAR
    test_xc_codecs();
#endif
#if XOMBE_HAVE_ENCODER
    test_encoder();
    test_encoder_absent();
    test_encoder_bits();
    test_encoder_tone();
    test_encoder_quant();
#if XOMBE_HAVE_DSTAR
    test_encoder_dstar();
#endif
#endif
    printf("API: %d/%d checks passed\n", checks - failures, checks);
    return failures ? 1 : 0;
}
