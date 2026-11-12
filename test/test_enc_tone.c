#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tone_fit.h"
#if XOMBE_HAVE_DSTAR
#include "../src/dstar/dstar_enc.h"
#endif

/* forced tones (xombe_enc_frame): position, decoded tone, unchanged analysis
   state, no leak into later frames, and the D-STAR mirror */

#define MAXF 4000
#define RUN 30 /* forced frames in a row for the decode measurements */
#define FIT_AT 12

static int failures, checks;

static void check(int ok, const char *what)
{
    checks++;
    if (!ok) {
        failures++;
        printf("FAIL: %s\n", what);
    }
}

static int16_t *load_raw(const char *dir, const char *name, int *frames)
{
    char path[4096];
    FILE *f;
    long size;
    int16_t *x;
    snprintf(path, sizeof(path), "%s/%s", dir, name);
    f = fopen(path, "rb");
    if (!f) {
        perror(path);
        exit(2);
    }
    fseek(f, 0, SEEK_END);
    size = ftell(f);
    fseek(f, 0, SEEK_SET);
    x = malloc((size_t)size + 2);
    *frames = (int)(fread(x, sizeof(int16_t), (size_t)size / 2, f) / 160);
    fclose(f);
    if (*frames > MAXF)
        *frames = MAXF;
    return x;
}

/* per input frame: requested tone (id 0: none) and the outputs that belong to it */
typedef struct {
    int frames;
    xombe_enc_tone req[MAXF];
    uint8_t a2[MAXF][7];
    uint8_t ds[MAXF][9];
} run;

static xombe_enc_mem mem_a, mem_b;

/* xombe_enc_frame with DMR voice bits and a D-STAR frame (ds NULL: none) */
static int enc_frame(xombe_enc *e, const int16_t *pcm, const xombe_enc_tone *t, uint8_t a2[7], uint8_t *ds)
{
    uint8_t dsb[7];
    int got = xombe_enc_frame_bits(e, pcm, t, a2, ds ? dsb : NULL);
    if (ds)
        xombe_fec_encode(XOMBE_DSTAR, dsb, ds);
    return got;
}

/* encodes r->frames of x with r->req, flushing so every input frame has its output */
static void encode(const xombe_enc_cfg *cfg, const int16_t *x, run *r)
{
    static const int16_t zeros[160];
    xombe_enc *e = xombe_enc_init(&mem_a, cfg);
    int dual = (cfg->outputs & ~XOMBE_ENC_DMR) != 0, delay = xombe_enc_delay(e), c;
    memset(r->a2, 0, sizeof(r->a2));
    memset(r->ds, 0, sizeof(r->ds));
    for (c = 0; c < r->frames + delay; c++) {
        const int16_t *pcm = c < r->frames ? x + 160 * c : zeros;
        const xombe_enc_tone *t = c < r->frames && r->req[c].id ? &r->req[c] : NULL;
        uint8_t a[7], d[9] = {0};
        int got;
        got = enc_frame(e, pcm, t, a, dual ? d : NULL);
        if (got == 1) {
            memcpy(r->a2[c - delay], a, 7);
            memcpy(r->ds[c - delay], d, 9);
        }
    }
}

static void tone_bytes(int id, int ad, uint8_t bytes[7])
{
    uint8_t d[49];
    xombe__tone_frame_bits(id, ad, d);
    xombe__pack49(d, bytes);
}

static int is_tone(const uint8_t bytes[7])
{
    uint8_t d[49];
    xombe__unpack49(bytes, d);
    return xombe__frame_class(d) == XOMBE_CLASS_TONE;
}

#if XOMBE_HAVE_DSTAR
static int ds_code(int id)
{
    int T;
    for (T = 0; T < 256; T++)
        if (xombe__ds_tone_id[T] == id)
            return T;
    return -1;
}

static int ds_is_tone(const uint8_t f[9])
{
    xombe_ds_frame fr;
    xombe__ds_fec_decode(f, &fr);
    return fr.cls == XOMBE_DS_TONE;
}
#endif

static xombe_enc_cfg cfg_of(unsigned outputs, int tones)
{
    xombe_enc_cfg c;
    xombe_enc_cfg_default(&c);
    c.outputs = outputs;
    c.tones = tones;
    return c;
}

static void label(char what[160], unsigned outputs, int tones, const char *s)
{
    const char *o = outputs == XOMBE_ENC_DMR ? "AMBE+2" : outputs == XOMBE_ENC_DSTAR ? "D-STAR" : "dual";
    snprintf(what, 160, "%s, tones %s: %s", o, tones ? "on" : "off", s);
}

/* forced frames land at their own input's output, in both outputs, and nothing else turns into a tone */
static void test_position(const int16_t *x, int frames, unsigned outputs, int tones)
{
    static run plain, forced;
    static const int at[][3] = {{40, 20, 90},   {41, 129, 70},  {77, 150, 100}, {120, 5, 0},
                                {121, 122, 127}, {122, 163, 3}, {200, 144, 60}, {201, 160, 110}};
    xombe_enc_cfg cfg = cfg_of(outputs, tones);
    int i, n = (int)(sizeof(at) / sizeof(at[0])), first = at[0][0], pos_ok = 1, before_ok = 1, leak = 0;
    int plain_tones = 0;
    char what[160];
    memset(&plain, 0, sizeof(plain));
    plain.frames = frames;
    encode(&cfg, x, &plain);
    forced = plain;
    for (i = 0; i < n; i++) {
        forced.req[at[i][0]].id = at[i][1];
        forced.req[at[i][0]].ad = at[i][2];
    }
    encode(&cfg, x, &forced);
    for (i = 0; i < frames; i++) {
        plain_tones += is_tone(plain.a2[i]);
        if (forced.req[i].id) {
            uint8_t want[7];
            tone_bytes(forced.req[i].id, forced.req[i].ad, want);
            if (outputs & XOMBE_ENC_DMR)
                pos_ok &= !memcmp(forced.a2[i], want, 7);
#if XOMBE_HAVE_DSTAR
            if (outputs & XOMBE_ENC_DSTAR) {
                int T = ds_code(forced.req[i].id);
                int ad = forced.req[i].ad;
                xombe_ds_frame fr;
                xombe__ds_fec_decode(forced.ds[i], &fr);
                if (T >= 0)
                    pos_ok &= fr.cls == XOMBE_DS_TONE && fr.tid == forced.req[i].id && fr.AD == (ad < 8 ? 8 : ad);
                else
                    pos_ok &= fr.cls == XOMBE_DS_VOICE;
            }
#endif
        } else {
            leak += is_tone(forced.a2[i]);
#if XOMBE_HAVE_DSTAR
            if (outputs & XOMBE_ENC_DSTAR)
                leak += ds_is_tone(forced.ds[i]);
#endif
        }
        if (i < first)
            before_ok &= !memcmp(plain.a2[i], forced.a2[i], 7) && !memcmp(plain.ds[i], forced.ds[i], 9);
    }
    label(what, outputs, tones, "speech input without detected tones");
    check(plain_tones == 0, what);
    label(what, outputs, tones, "each forced tone is its input frame's output");
    check(pos_ok, what);
    label(what, outputs, tones, "frames before the first request are unchanged");
    check(before_ok, what);
    label(what, outputs, tones, "no tone frame without a request");
    check(leak == 0, what);
}

/* a forced tone inside a detected tone run changes nothing but its own frame: tone frames leave the AMBE+2
   predictor at class tone and the D-STAR mirror alone, so the rest of the stream is bit-identical unless the
   request leaked into the detector's confirmation and hold */
static void test_in_tone_run(const int16_t *x, int frames, unsigned outputs)
{
    static run plain, forced;
    xombe_enc_cfg cfg = cfg_of(outputs, 1);
    int i, len = 0, placed = 0, same = 1, j;
    memset(&plain, 0, sizeof(plain));
    plain.frames = frames;
    encode(&cfg, x, &plain);
    forced = plain;
    for (i = 0; i < frames && placed < 6; i++) {
        int both = is_tone(plain.a2[i]);
#if XOMBE_HAVE_DSTAR
        if (outputs & XOMBE_ENC_DSTAR)
            both &= ds_is_tone(plain.ds[i]);
#endif
        len = both ? len + 1 : 0;
        if (len == 5) {
            /* last frame of a natural run's first five, also the frame after its end when the run is short */
            forced.req[i - 2].id = 130 + placed;
            forced.req[i - 2].ad = 80;
            placed++;
        }
    }
    check(placed >= 3, "tone input has detected tone runs");
    encode(&cfg, x, &forced);
    for (i = 0; i < frames; i++) {
        if (forced.req[i].id)
            continue;
        same &= !memcmp(plain.a2[i], forced.a2[i], 7);
        for (j = 0; j < 9; j++)
            same &= plain.ds[i][j] == forced.ds[i][j];
    }
    check(same, outputs & XOMBE_ENC_DSTAR ? "dual: a forced tone in a tone run leaves every other frame unchanged"
                                          : "AMBE+2: a forced tone in a tone run leaves every other frame unchanged");
}

#define SAME_FIELD(a, b, f) (!memcmp(&(a)->f, &(b)->f, sizeof((a)->f)))

/* VAD, hangover, noise and comfort-noise tracking, and the tone detector's state */
static int same_analysis(const struct xombe_enc *a, const struct xombe_enc *b)
{
    size_t lo = offsetof(struct xombe_enc, vad_nl), hi = offsetof(struct xombe_enc, tone_spec);
    int same = !memcmp((const char *)a + lo, (const char *)b + lo, hi - lo) && a->tone_prev_id == b->tone_prev_id &&
               a->tone_confirmed == b->tone_confirmed && a->tone_hold_used == b->tone_hold_used &&
               a->tone_last_id == b->tone_last_id && a->tone_last_ad == b->tone_last_ad && a->p_prev == b->p_prev;
#if XOMBE_HAVE_DSTAR
    if (a->cfg.outputs & XOMBE_ENC_DSTAR) {
        const xombe_ds_enc *da = XOMBE_DS_ENC(a), *db = XOMBE_DS_ENC(b);
        same &= SAME_FIELD(da, db, ns_mem) && SAME_FIELD(&da->tones, &db->tones, detected) &&
                SAME_FIELD(&da->tones, &db->tones, run) && SAME_FIELD(&da->tones, &db->tones, confirmed) &&
                SAME_FIELD(&da->tones, &db->tones, held) && da->tones.hold_left == db->tones.hold_left &&
                SAME_FIELD(&da->tones, &db->tones, level_db);
    }
#endif
    return same;
}

/* lockstep with an unforced encoder: up to and including the call that emits the forced frame, the analysis
   state is the same, a request changes only the output and the prediction memory */
static void test_state(const int16_t *x, unsigned outputs, int tones, int dtx, int start)
{
    char what[160];
    xombe_enc_cfg cfg = cfg_of(outputs, tones);
    xombe_enc *a, *b;
    xombe_enc_tone t = {150, 90};
    int c, same = 1, delay;
    cfg.dmr_dtx = dtx;
    a = xombe_enc_init(&mem_a, &cfg);
    b = xombe_enc_init(&mem_b, &cfg);
    delay = xombe_enc_delay(a);
    for (c = 0; c <= start + delay + 3; c++) {
        const xombe_enc_tone *r = c >= start && c < start + 4 ? &t : NULL;
        uint8_t a2[9], ds[9];
        xombe_enc_frame(a, x + 160 * c, r, a2, ds);
        xombe_enc_frame(b, x + 160 * c, NULL, a2, ds);
        same &= same_analysis(a, b);
    }
    label(what, outputs, tones, "a request leaves VAD, noise, comfort-noise, and tone detector state as without it");
    check(same, what);
}

static void fit(const double *x, int id, double *freq, double *amp, int *nf)
{
    float lo, hi;
    *nf = xombe__tone_freqs(id, &lo, &hi);
    freq[0] = lo;
    freq[1] = hi;
    search_freqs(x, freq, *nf, *nf == 1 ? 0.005 : 0.05);
    sinusoid_fit(x, FIT_LEN, freq, *nf, amp);
}

/* frequency and level errors of a fitted tone, the T5 limits */
static int tone_ok(int id, int ad, const double *freq, const double *amp, int nf, double *worst_df, double *worst_dl)
{
    float lo, hi;
    double nominal[2];
    int c, ok = 1;
    xombe__tone_freqs(id, &lo, &hi);
    nominal[0] = lo;
    nominal[1] = hi;
    for (c = 0; c < nf; c++) {
        double df = nf == 1 ? fabs(freq[c] - nominal[c]) : fabs(freq[c] - nominal[c]) / nominal[c];
        double dl = 20.0 * log10(amp[c] / 32768.0 + 1e-30) - expected_db(ad, freq[c]);
        ok &= nf == 1 ? df <= 0.1 : df <= (id >= 160 ? 0.04 : 0.01);
        ok &= fabs(dl) <= 0.5;
        if (nf == 2 && df > *worst_df)
            *worst_df = df;
        if (fabs(dl) > fabs(*worst_dl))
            *worst_dl = dl;
    }
    return ok;
}

/* RUN forced frames in speech, decoded through xombe_dec (and xombe_dsd), measured FIT_AT frames into the run */
static void test_decode(const int16_t *x, int frames, unsigned outputs)
{
    static run r;
    static double y[MAXF * 160];
    static const int ids[] = {5, 20, 64, 122, 128, 129, 143, 144, 150, 156, 159, 160, 161, 163};
    static const int ads[] = {40, 90, 115};
    const int start = 60;
    xombe_enc_cfg cfg = cfg_of(outputs, 1);
    double wdf = 0.0, wdl = 0.0;
    int ii, ai, i, ok = 1, cases = 0;
#if XOMBE_HAVE_DSTAR
    double wdf_ds = 0.0, wdl_ds = 0.0, knox_df = 0.0, knox_res = 0.0;
    int ok_ds = 1, knox_voice = 1;
#endif
    for (ii = 0; ii < (int)(sizeof(ids) / sizeof(ids[0])); ii++) {
        for (ai = 0; ai < 3; ai++) {
            int id = ids[ii], ad = ads[ai], nf;
            double freq[2], amp[2];
            xombe_dec dec;
            memset(&r, 0, sizeof(r));
            r.frames = start + RUN + 1 < frames ? start + RUN + 1 : frames;
            for (i = start; i < start + RUN; i++) {
                r.req[i].id = id;
                r.req[i].ad = ad;
            }
            encode(&cfg, x, &r);
            xombe_dec_set_noise(xombe_dec_init(&dec, XOMBE_DMR), XOMBE_NOISE_CONFORMANCE, 0);
            for (i = 0; i < start + RUN; i++) {
                int16_t pcm[160];
                int k;
                xombe_dec_frame_bits(&dec, r.a2[i], 0, pcm);
                for (k = 0; k < 160; k++)
                    y[i * 160 + k] = pcm[k];
            }
            fit(y + (start + FIT_AT) * 160, id, freq, amp, &nf);
            if (!tone_ok(id, ad, freq, amp, nf, &wdf, &wdl)) {
                ok = 0;
                printf("  AMBE+2 id %d AD %d: %.2f/%.2f Hz\n", id, ad, freq[0], freq[1]);
            }
            cases++;
#if XOMBE_HAVE_DSTAR
            if (outputs & XOMBE_ENC_DSTAR) {
                static xombe_dec_mem dm;
                xombe_dec *dsd = xombe_dec_init(&dm, XOMBE_DSTAR);
                int voice = 1;
                xombe_dec_set_noise(dsd, XOMBE_NOISE_CONFORMANCE, 0);
                for (i = 0; i < start + RUN; i++) {
                    int16_t pcm[160];
                    xombe_decision dc;
                    int k;
                    xombe_dec_frame(dsd, r.ds[i], 0, pcm);
                    xombe_dec_decision_get(dsd, &dc);
                    if (i >= start)
                        voice &= dc.cls == XOMBE_CLASS_VOICE;
                    for (k = 0; k < 160; k++)
                        y[i * 160 + k] = pcm[k];
                }
                if (ds_code(id) < 0) {
                    /* Knox: voice frames with the two tones on the nearest D-STAR harmonics */
                    const double *seg = y + (start + FIT_AT) * 160;
                    double e = 0.0, res, df = 0.0;
                    float lo, hi;
                    int c2;
                    xombe__tone_freqs(id, &lo, &hi);
                    fit(seg, id, freq, amp, &nf);
                    for (i = 0; i < FIT_LEN; i++)
                        e += seg[i] * seg[i];
                    res = sinusoid_fit(seg, FIT_LEN, freq, nf, amp) / e;
                    for (c2 = 0; c2 < 2; c2++) {
                        double d = fabs(freq[c2] - (c2 ? hi : lo)) / (c2 ? hi : lo);
                        df = d > df ? d : df;
                    }
                    knox_df = df > knox_df ? df : knox_df;
                    knox_res = res > knox_res ? res : knox_res;
                    if (!voice || df > 0.03 || res > 0.6) {
                        knox_voice = 0;
                        printf("  D-STAR Knox id %d AD %d: voice %d, frequency error %.2f%%, residual %.2f\n", id, ad,
                               voice, 100.0 * df, res);
                    }
                } else {
                    fit(y + (start + FIT_AT) * 160, id, freq, amp, &nf);
                    if (!tone_ok(id, ad, freq, amp, nf, &wdf_ds, &wdl_ds)) {
                        ok_ds = 0;
                        printf("  D-STAR id %d AD %d: %.2f/%.2f Hz, %.1f/%.1f\n", id, ad, freq[0], freq[1],
                               20.0 * log10(amp[0] / 32768.0), expected_db(ad, freq[0]));
                    }
                }
            }
#endif
        }
    }
    printf("forced tones (%s): %d cases, AMBE+2 worst dual frequency error %.2f%%, worst level error %+.3f dB\n",
           outputs & XOMBE_ENC_DSTAR ? "dual" : "AMBE+2", cases, 100.0 * wdf, wdl);
    check(ok, "forced tones decode at their frequency and level (T5 limits)");
#if XOMBE_HAVE_DSTAR
    if (outputs & XOMBE_ENC_DSTAR) {
        printf("forced tones, D-STAR: worst dual frequency error %.2f%%, worst level error %+.3f dB\n", 100.0 * wdf_ds,
               wdl_ds);
        check(ok_ds, "forced D-STAR tones decode at their frequency and level (T5 limits)");
        printf("forced Knox tones, D-STAR: worst frequency error %.2f%%, worst residual %.2f of the energy\n",
               100.0 * knox_df, knox_res);
        check(knox_voice, "forced Knox tones are D-STAR voice frames carrying the tone");
    }
#endif
}

/* the packers' frames are the forced encoder's, for every ID at levels 0..127 in steps of 7 and 127, and decode
   as tone frames of that ID and level */
static void test_packers(const int16_t *x, int frames, unsigned outputs)
{
    static run r;
    static int16_t loop[MAXF * 160];
    xombe_enc_cfg cfg = cfg_of(outputs, 1);
    int id, ad, i, n = 0, same = 1, decoded = 1;
    memset(&r, 0, sizeof(r));
    for (i = 0; i < MAXF; i++)
        memcpy(loop + 160 * i, x + 160 * (i % frames), sizeof(int16_t) * 160);
    for (id = 5; id <= 163; id++)
        for (ad = 0; ad <= 127 + 6; ad += 7)
            if (xombe__tone_id_accepted(id) && n < MAXF) {
                r.req[n].id = id;
                r.req[n].ad = ad > 127 ? 127 : ad;
                n++;
            }
    r.frames = n;
    encode(&cfg, loop, &r);
    for (i = 0; i < n; i++) {
        static xombe_dec dec;
        uint8_t a2[7], f[9], d[9], d49[49];
        int16_t pcm[160];
        xombe_decision dc;
        xombe_fec fec;
        int ok = xombe_tone_frame(XOMBE_DMR, r.req[i].id, r.req[i].ad, f) == 0, pid, pad;
        xombe_fec_init(&fec, XOMBE_DMR);
        xombe_fec_decode(&fec, f, a2, NULL);
        if (outputs & XOMBE_ENC_DMR)
            same &= ok && !memcmp(a2, r.a2[i], 7);
        xombe_dec_set_noise(xombe_dec_init(&dec, XOMBE_DMR), XOMBE_NOISE_CONFORMANCE, 0);
        xombe_dec_frame(&dec, f, 0, pcm);
        xombe_dec_decision_get(&dec, &dc);
        xombe__unpack49(a2, d49);
        decoded &= dc.cls == XOMBE_CLASS_TONE && xombe__tone_frame_parse(d49, &pid, &pad) && pid == r.req[i].id &&
                   pad == r.req[i].ad;
#if XOMBE_HAVE_DSTAR
        if (outputs & XOMBE_ENC_DSTAR) {
            int T = ds_code(r.req[i].id);
            ok = xombe_tone_frame(XOMBE_DSTAR, r.req[i].id, r.req[i].ad, d);
            same &= T >= 0 ? ok == 0 && !memcmp(d, r.ds[i], 9) : ok == -1;
            if (T >= 0) {
                static xombe_dec_mem dm;
                xombe_dec *dsd = xombe_dec_init(&dm, XOMBE_DSTAR);
                xombe_decision dd;
                xombe_ds_frame fr;
                xombe_dec_set_noise(dsd, XOMBE_NOISE_CONFORMANCE, 0);
                xombe_dec_frame(dsd, d, 0, pcm);
                xombe_dec_decision_get(dsd, &dd);
                xombe__ds_fec_decode(d, &fr);
                decoded &= dd.cls == XOMBE_CLASS_TONE && fr.tid == r.req[i].id &&
                           fr.AD == (r.req[i].ad < 8 ? 8 : r.req[i].ad);
            }
        }
#else
        (void)d;
#endif
    }
    check(n < MAXF, "every ID in the packer sweep");
    check(same, outputs & XOMBE_ENC_DSTAR ? "packers equal the forced dual encoder's frames"
                                          : "xombe_tone_frame equals the forced DMR encoder's frames");
    check(decoded, "packer frames decode as tone frames of their ID and level");
}

#if XOMBE_HAVE_DSTAR
/* the D-STAR encoder's mirror equals the decoder's, replayed from the emitted frames, after every frame, with
   forced tones in speech */
static void test_ds_mirror(const int16_t *x, int frames)
{
    static run r;
    xombe_enc_cfg cfg = cfg_of(XOMBE_ENC_DMR | XOMBE_ENC_DSTAR, 1);
    xombe_enc *e = xombe_enc_init(&mem_a, &cfg);
    xombe_ds_mirror dm;
    int c, bad = 0, n = 0, knox = 0, tone = 0, knox_bad = 0, delay = xombe_enc_delay(e);
    uint32_t s = 12345u;
    xombe__ds_mirror_init(&dm);
    memset(&r, 0, sizeof(r));
    for (c = 30; c < frames; c += 7 + (int)(s % 23)) {
        int k, len = 1 + (int)((s >> 8) % 3), id = (int)((s >> 12) % 154) + 5;
        if (id > 122)
            id += 5;
        for (k = 0; k < len && c + k < frames; k++) {
            r.req[c + k].id = id;
            r.req[c + k].ad = (int)((s >> 16) % 128);
        }
        knox += id >= 144 && id <= 159;
        tone += id < 144 || id > 159;
        s = s * 1664525u + 1013904223u;
    }
    for (c = 0; c < frames; c++) {
        static xombe_ds_enc before;
        uint8_t a2[7], ds[9], want[9];
        const xombe_enc_tone *t = c >= delay ? &r.req[c - delay] : NULL;
        int is_knox = t && t->id >= 144 && t->id <= 159;
        if (is_knox) {
            /* the encoder's own Knox frame at the requested level */
            int k = 0, b[9];
            before = *XOMBE_DS_ENC(e);
            while (before.tones.knox[k].id != t->id)
                k++;
            xombe__ds_knox_amps(&before, k, (float)t->ad, b);
            xombe__ds_encode_fields(b, want);
        }
        if (enc_frame(e, x + 160 * c, r.req[c].id ? &r.req[c] : NULL, a2, ds) == 1) {
            const xombe_ds_mirror *m = &XOMBE_DS_ENC(e)->mirror;
            xombe_ds_frame fr;
            xombe_ds_params p;
            xombe__ds_fec_decode(ds, &fr);
            xombe__ds_params_from_frame(&dm, &fr, &p);
            bad += m->L != dm.L || m->G != dm.G || memcmp(m->a, dm.a, sizeof(m->a));
            knox_bad += is_knox && memcmp(ds, want, 9);
            n++;
        }
    }
    printf("D-STAR mirror with %d forced tone and %d Knox requests: %d/%d frames equal\n", tone, knox, n - bad, n);
    check(knox > 0 && tone > 0, "mirror test forces both kinds");
    check(bad == 0, "the D-STAR mirror stays equal to the decoder's across forced tones");
    check(knox_bad == 0, "a forced Knox tone is the D-STAR encoder's Knox voice frame");
}
#endif

int main(int argc, char **argv)
{
    const char *vec = argc > 1 ? argv[1] : "vectors";
    int fs, ft;
    int16_t *speech = load_raw(vec, "e3/es1_orig.raw", &fs), *tones = load_raw(vec, "e7/tones.raw", &ft);
    test_position(speech, fs, XOMBE_ENC_DMR, 1);
    test_position(speech, fs, XOMBE_ENC_DMR, 0);
    test_in_tone_run(tones, ft, XOMBE_ENC_DMR);
    test_state(speech, XOMBE_ENC_DMR, 1, 1, 50);
    test_state(tones, XOMBE_ENC_DMR, 1, 1, 50);
    test_state(speech, XOMBE_ENC_DMR, 0, 0, 80);
    test_decode(speech, fs, XOMBE_ENC_DMR);
    test_packers(speech, fs, XOMBE_ENC_DMR);
#if XOMBE_HAVE_DSTAR
    test_position(speech, fs, XOMBE_ENC_DMR | XOMBE_ENC_DSTAR, 1);
    test_position(speech, fs, XOMBE_ENC_DMR | XOMBE_ENC_DSTAR, 0);
    test_position(speech, fs, XOMBE_ENC_DSTAR, 1);
    test_in_tone_run(tones, ft, XOMBE_ENC_DMR | XOMBE_ENC_DSTAR);
    test_state(speech, XOMBE_ENC_DMR | XOMBE_ENC_DSTAR, 1, 1, 50);
    test_state(tones, XOMBE_ENC_DMR | XOMBE_ENC_DSTAR, 1, 0, 50);
    test_decode(speech, fs, XOMBE_ENC_DMR | XOMBE_ENC_DSTAR);
    test_ds_mirror(speech, fs);
    test_packers(speech, fs, XOMBE_ENC_DMR | XOMBE_ENC_DSTAR);
#endif
    free(speech);
    free(tones);
    printf("forced tones: %d/%d checks passed\n", checks - failures, checks);
    return failures ? 1 : 0;
}
