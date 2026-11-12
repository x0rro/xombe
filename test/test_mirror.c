#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../src/internal.h"
#include "../src/tables/tables.h"

#define FS 8000

static long e1_frames, e1_bad, mirror_frames, mirror_bad, e7_fail;
static const xombe_enc_tone *forced; /* per input frame, id 0: none, NULL: no requests */
static long forced_frames;
static double mirror_max_err;
static uint32_t rng = 20261019u;

static uint32_t rnd(void)
{
    rng = rng * 1664525u + 1013904223u;
    return rng >> 8;
}

static double urand(void)
{
    return (rnd() & 0xFFFFFF) / 16777216.0;
}

static int16_t *load_raw(const char *dir, const char *name, long *n)
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
    *n = (long)fread(x, sizeof(int16_t), (size_t)size / 2, f);
    fclose(f);
    return x;
}

static int e1_frame_ok(const uint8_t d[49])
{
    int c = xombe__get_bits(d, 0, 6), b[9], i;
    if (c == 60 || c == 61)
        return 0;
    if (c == 63) {
        int id = xombe__get_bits(d, 12, 8);
        for (i = 1; i < 4; i++)
            if (xombe__get_bits(d, 12 + 8 * i, 8) != id)
                return 0;
        return xombe__tone_id_accepted(id) && !d[45] && !d[46] && !d[47] && !d[48];
    }
    xombe__get_fields(d, b);
    if (c == 62)
        return b[0] == XOMBE_SILENCE_B0 && b[1] == XOMBE_SILENCE_B1;
    return b[0] <= 119 && (b[1] < 18 || b[1] > 31);
}

static int same_pred(const xombe_pred *a, const xombe_pred *b)
{
    int l;
    if (a->cls != b->cls || a->L != b->L || a->f0 != b->f0 || a->gamma != b->gamma || memcmp(a->v2, b->v2, 8))
        return 0;
    for (l = 0; l <= XOMBE_LMAX; l++)
        if (a->a[l] != b->a[l])
            return 0;
    return 1;
}

/* E2(a): decoder's dequantized f0, L, and amplitudes equal the encoder's reconstruction, the dequantization of
   the emitted fields from the encoder's predictor before the frame. the encoder keeps that reconstruction as its
   predictor after a voice frame, so the two predictors must also agree */
static int mirror_ok(const xombe_pred *before, const xombe_enc *enc, const xombe_dec *dec, const uint8_t d[49],
                     int cls)
{
    const xombe_set *p2 = &dec->c.p2;
    xombe_params rec;
    int b[9], l;
    if (enc->pi.cls != cls || !same_pred(&enc->pi, &dec->c.pi))
        return 0;
    if (cls == XOMBE_CLASS_TONE)
        return 1;
    xombe__get_fields(d, b);
    xombe__dequant(b, before, cls == XOMBE_CLASS_SILENCE, &rec);
    if (p2->L != rec.L || p2->f0 != rec.f0)
        return 0;
    for (l = 1; l <= p2->L; l++) {
        double err = fabs((double)p2->a[l] - rec.a[l]);
        if (err > mirror_max_err)
            mirror_max_err = err;
        if (err > 1e-6)
            return 0;
    }
    return 1;
}

typedef struct {
    long frames;
    int16_t *pcm;
    int *cls;
    int *b0;
    int *voiced_bands;
    int state_finite; /* encoder and decoder state free of NaN and infinity after the run */
} encode_result;

static int all_finite(const float *v, int n)
{
    int i;
    for (i = 0; i < n; i++)
        if (!isfinite(v[i]))
            return 0;
    return 1;
}

static int pred_finite(const xombe_pred *p)
{
    return isfinite(p->f0) && isfinite(p->gamma) && all_finite(p->a, XOMBE_LMAX + 1);
}

static int state_finite(const xombe_enc *e, const xombe_dec *d)
{
    return pred_finite(&e->pi) && pred_finite(&d->c.pi) && isfinite(d->c.lambda) && isfinite(d->c.avg_f0) &&
           isfinite(d->hpf.y1) && isfinite(d->hpf.y2) && isfinite(e->hpf.y1) && isfinite(e->hpf.y2) &&
           isfinite(d->pitch_phase) && isfinite(d->epoch) && all_finite(d->tail, XOMBE_TAIL) &&
           all_finite(d->uv_overlap, XOMBE_TAIL) && isfinite(e->ns_floor) &&
           all_finite(e->ns_noise, XOMBE_NS_BINS) && all_finite(e->cn_power, XOMBE_SPEC_BINS) &&
           all_finite(e->vad_nl, 8) && isfinite(e->p_prev);
}

static void e7_check(int ok, const char *what)
{
    if (!ok) {
        e7_fail++;
        printf("E7 FAIL: %s\n", what);
    }
}

static void run_case(const char *name, const int16_t *x, long n, const xombe_enc_cfg *cfg, encode_result *res)
{
    static xombe_enc enc;
    static xombe_dec dec;
    long frames = n / 160, f, out = 0, bad = 0;
    xombe_enc_init(&enc, cfg);
    xombe_dec_set_noise(xombe_dec_init(&dec, XOMBE_DMR), XOMBE_NOISE_CONFORMANCE, 0);
    if (res) {
        res->frames = frames;
        res->pcm = malloc(sizeof(int16_t) * 160 * (size_t)(frames + 1));
        res->cls = malloc(sizeof(int) * (size_t)(frames + 1));
        res->b0 = malloc(sizeof(int) * (size_t)(frames + 1));
        res->voiced_bands = malloc(sizeof(int) * (size_t)(frames + 1));
    }
    for (f = 0; out < frames; f++) {
        uint8_t bits[7], d[49];
        int16_t pcm[160];
        xombe_pred before = enc.pi;
        const xombe_enc_tone *t = forced && f < frames && forced[f].id ? &forced[f] : NULL;
        int got = xombe_enc_frame_bits(&enc, f < frames ? x + 160 * f : NULL, t, bits, NULL);
        forced_frames += t != NULL;
        int cls, b[9];
        if (!got)
            continue;
        xombe__unpack49(bits, d);
        e1_frames++;
        if (!e1_frame_ok(d)) {
            bad++;
            e1_bad++;
        }
        xombe_dec_frame_bits(&dec, bits, 0, pcm);
        cls = xombe__frame_class(d);
        mirror_bad += !mirror_ok(&before, &enc, &dec, d, cls);
        mirror_frames++;
        if (res) {
            xombe__get_fields(d, b);
            memcpy(res->pcm + 160 * out, pcm, sizeof(pcm));
            res->cls[out] = cls;
            res->b0[out] = b[0];
            res->voiced_bands[out] = cls == XOMBE_CLASS_VOICE ? xombe__count_voiced(xombe__vuv[b[1]]) : 0;
        }
        out++;
    }
    if (bad)
        printf("  %s: %ld E1 violations\n", name, bad);
    if (res)
        res->state_finite = state_finite(&enc, &dec);
}

static void free_result(encode_result *r)
{
    free(r->pcm);
    free(r->cls);
    free(r->b0);
    free(r->voiced_bands);
}

static void test_files(const char *vec)
{
    static const char *const files[] = {
        "e3/es1_orig.raw", "e3/es2_orig.raw", "e3/es3_orig.raw", "e3/sn10_orig.raw", "e3/ttsf_orig.raw",
        "e3/ttsc_orig.raw", "e3/ttsl_orig.raw", "e5/low5.raw", "e5/hp5.raw", "e5/es1_lead0.5s.raw",
        "e7/glide.raw", "e7/mix.raw", "e7/saw150.raw", "e7/tones.raw"};
    xombe_enc_cfg cfg[4];
    unsigned i, v;
    for (v = 0; v < 4; v++)
        xombe_enc_cfg_default(&cfg[v]);
    cfg[1].dmr_dtx = 0;
    cfg[2].ns = 0;
    cfg[2].tones = 0;
    cfg[3].voicing_scale_permille = 600;
    for (i = 0; i < sizeof(files) / sizeof(files[0]); i++) {
        long n;
        int16_t *x = load_raw(vec, files[i], &n);
        for (v = 0; v < 4; v++)
            if (v == 0 || i % 3 == v - 1)
                run_case(files[i], x, n, &cfg[v], NULL);
        free(x);
    }
}

static uint32_t lcg(uint32_t *s)
{
    *s = *s * 1664525u + 1013904223u;
    return *s >> 8;
}

/* E2 across forced tones: bursts of one to three tone frames at random IDs and levels in speech and noise */
static void test_forced(const char *vec)
{
    static const char *const files[] = {"e3/es1_orig.raw", "e3/sn10_orig.raw", "e3/ttsl_orig.raw", "e5/low5.raw"};
    xombe_enc_cfg cfg;
    uint32_t r = 7u;
    unsigned i;
    for (i = 0; i < sizeof(files) / sizeof(files[0]); i++) {
        long n, f, frames;
        int16_t *x = load_raw(vec, files[i], &n);
        xombe_enc_tone *req;
        frames = n / 160;
        req = calloc((size_t)frames + 1, sizeof(*req));
        for (f = 20; f < frames; f += 5 + (long)(lcg(&r) % 30)) {
            int k, len = 1 + (int)(lcg(&r) % 3), id = 5 + (int)(lcg(&r) % 154);
            if (id > 122)
                id += 5;
            for (k = 0; k < len && f + k < frames; k++) {
                req[f + k].id = id;
                req[f + k].ad = (int)(lcg(&r) % 128);
            }
        }
        forced = req;
        xombe_enc_cfg_default(&cfg);
        cfg.tones = (int)(i & 1);
        cfg.dmr_dtx = i != 2;
        run_case(files[i], x, n, &cfg, NULL);
        forced = NULL;
        free(req);
        free(x);
    }
}

static void random_frame(int b[9], int silence)
{
    b[0] = silence ? XOMBE_SILENCE_B0 : (int)(rnd() % 120);
    b[1] = silence ? XOMBE_SILENCE_B1 : (int)(rnd() % 18);
    b[2] = (int)(rnd() % 32);
    b[3] = (int)(rnd() % 512);
    b[4] = (int)(rnd() % 128);
    b[5] = (int)(rnd() % 32);
    b[6] = (int)(rnd() % 16);
    b[7] = (int)(rnd() % 16);
    b[8] = (int)(rnd() % 8);
}

static int same_recon(const xombe_params *a, const xombe_params *b)
{
    int l;
    if (a->L != b->L || fabsf(a->gamma - b->gamma) > 1e-6f)
        return 0;
    for (l = 1; l <= a->L; l++)
        if (fabsf(a->a[l] - b->a[l]) > 1e-5f)
            return 0;
    return 1;
}

static int clamp_hit(const xombe_params *p)
{
    int l;
    if (p->gamma >= XOMBE_GAIN_MAX)
        return 1;
    for (l = 1; l <= p->L; l++)
        if (p->a[l] <= XOMBE_AMP_MIN || p->a[l] >= XOMBE_AMP_MAX)
            return 1;
    return 0;
}

static long idem_frames, idem_exact, idem_equivalent, idem_clamped, idem_bad;

static void test_idempotence(int count)
{
    xombe_pred pi;
    int i;
    xombe__pred_init(&pi);
    for (i = 0; i < count; i++) {
        int b[9], q[9], silence = (rnd() % 5) == 0, k, same = 1;
        xombe_params p, back;
        uint8_t v2[8];
        random_frame(b, silence);
        xombe__dequant(b, &pi, silence, &p);
        memcpy(q, b, sizeof(q));
        xombe__quantize(&pi, silence, p.L, p.a, q);
        for (k = 2; k < 9; k++)
            same &= q[k] == b[k];
        idem_frames++;
        if (same) {
            idem_exact++;
        } else {
            xombe__dequant(q, &pi, silence, &back);
            if (same_recon(&p, &back))
                idem_equivalent++;
            else if (clamp_hit(&p))
                idem_clamped++;
            else {
                idem_bad++;
                if (idem_bad <= 10)
                    printf("  idempotence: frame %d b=%d,%d,%d,%d,%d,%d,%d got %d,%d,%d,%d,%d,%d,%d\n", i, b[2], b[3],
                           b[4], b[5], b[6], b[7], b[8], q[2], q[3], q[4], q[5], q[6], q[7], q[8]);
            }
        }
        if (silence) {
            xombe__pred_store_other(&pi, XOMBE_CLASS_SILENCE);
        } else {
            xombe__v2_from_b1(b[1], v2);
            xombe__pred_store_voice(&pi, &p, v2);
        }
    }
}

static double level_dbfs(const int16_t *x, long n)
{
    double s = 0.0;
    long i;
    for (i = 0; i < n; i++)
        s += (double)x[i] * x[i];
    return s > 0.0 ? 10.0 * log10(s / n / (32768.0 * 32768.0)) : -200.0;
}

static int16_t clip16(double v)
{
    return (int16_t)(v > 32767.0 ? 32767 : v < -32768.0 ? -32768 : floor(v + 0.5));
}

static void am_sweep(int16_t *x, long n, long lead, double f_start, double f_end, double *f_at)
{
    double phase = 0.0;
    long i;
    for (i = 0; i < n; i++) {
        double t = (double)(i - lead) / FS, dur = (double)(n - lead) / FS, f, v = 0.0, env;
        int k;
        if (i < lead) {
            x[i] = 0;
            if (f_at)
                f_at[i] = 0.0;
            continue;
        }
        f = f_start * pow(f_end / f_start, t / dur);
        if (f_at)
            f_at[i] = f;
        phase += 2.0 * XOMBE_PI * f / FS;
        for (k = 1; k * f < 3900.0; k++)
            v += cos(k * phase);
        env = 1.0 + 0.5 * sin(2.0 * XOMBE_PI * 3.5 * t);
        x[i] = clip16(600.0 * env * v / sqrt((double)(3900.0 / f)));
    }
}

static void sweep_report(const char *name, double f_start, double f_end, int16_t *x, long n)
{
    double *f_at = malloc(sizeof(double) * (size_t)n);
    encode_result r;
    xombe_enc_cfg cfg;
    long f, in_range = 0, good = 0, octave_low = 0;
    am_sweep(x, n, FS / 2, f_start, f_end, f_at);
    xombe_enc_cfg_default(&cfg);
    run_case(name, x, n, &cfg, &r);
    for (f = 0; f < r.frames; f++) {
        double truth = f_at[160 * f + 120 < n ? 160 * f + 120 : n - 1], coded;
        if (truth < 66.0 || truth > 398.0 || r.cls[f] != XOMBE_CLASS_VOICE)
            continue;
        in_range++;
        coded = xombe__w0[r.b0[f]] * FS;
        if (fabs(coded - truth) <= 0.05 * truth)
            good++;
        else if (fabs(coded - truth / 2.0) <= 0.05 * truth / 2.0)
            octave_low++;
    }
    printf("E7 %s (informative): coded pitch within 5%% on %ld of %ld voice frames in 66..398 Hz (%ld an octave "
           "low)\n", name, good, in_range, octave_low);
    free_result(&r);
    free(f_at);
}

/* decoded over input power in dB on the pause frames of x from frame from on, frames whose clean source and two
   neighbours on each side are below 20 dB (16-bit units), which also absorbs the decoder's delay
   pcm[i] is aligned with x[i] */
static double pause_suppression(const int16_t *x, const int16_t *clean, const int16_t *pcm, long n, long from)
{
    long nfr = n / 160, i, j, k;
    double pin = 0.0, pout = 0.0;
    for (i = from; i < nfr - 2; i++) {
        int quiet = 1;
        for (j = i - 2; j <= i + 2 && quiet; j++) {
            double e = 0.0;
            for (k = 0; k < 160; k++)
                e += (double)clean[160 * j + k] * clean[160 * j + k];
            quiet = 10.0 * log10(e / 160.0 + 1e-3) < 20.0;
        }
        if (!quiet)
            continue;
        for (k = 160 * i; k < 160 * (i + 1); k++) {
            pin += (double)x[k] * x[k];
            pout += (double)pcm[k] * pcm[k];
        }
    }
    return 10.0 * log10((pout + 1e-9) / (pin + 1e-9));
}

/* VAD recovery (inputs as E5's vad_* from tools/e5.py --vad-gen)
   three triggers push the noise levels under the background of sn10 (es1 + white noise, 10 dB SNR, 0.5 s
   noise lead-in): 120 zero samples prepended, digital-silence gaps in two pauses (0 from 0.16 s, constant +1
   from 3.68 s, 8 frames each), and the noise raised 8 dB from 2.0 s
   from 7.5 s on, pause frames must be suppressed at most 3 dB less than in the unmodified clip, i.e. the noise
   levels have recovered (with the digital-silence guard, vad_gap never locks) */
static void test_vad_recovery(const char *vec, const char *inputs)
{
    long n, nc, i, m;
    int16_t *sn = load_raw(inputs, "sn10.d0.raw", &n), *clean = load_raw(vec, "e5/es1_lead0.5s.raw", &nc);
    int16_t *x = malloc(sizeof(int16_t) * (size_t)(n + 120)), *cl = calloc((size_t)(n + 120), sizeof(int16_t));
    const char *names[4] = {"sn10", "vad_on", "vad_gap", "vad_step"};
    const long from = 375;
    double ref = 0.0, g = pow(10.0, 8.0 / 20.0);
    int t;
    encode_result r;
    xombe_enc_cfg cfg;
    xombe_enc_cfg_default(&cfg);
    if (nc < n)
        n = nc;
    for (t = 0; t < 4; t++) {
        const int16_t *in = x;
        long off = t == 1 ? 120 : 0;
        double sup;
        memset(x, 0, sizeof(int16_t) * (size_t)(n + 120));
        for (i = 0; i < n; i++) {
            cl[off + i] = clean[i];
            x[off + i] = t == 3 && i >= 2 * FS ? clip16(clean[i] + g * (sn[i] - clean[i])) : sn[i];
        }
        if (t == 2) {
            for (i = 8 * 160; i < 16 * 160; i++)
                x[i] = 0;
            for (i = 184 * 160; i < 192 * 160; i++)
                x[i] = 1;
        }
        run_case(names[t], in, n + off, &cfg, &r);
        m = 160 * r.frames - off;
        sup = pause_suppression(x + off, cl + off, r.pcm + off, m < n ? m : n, from);
        if (t == 0) {
            ref = sup;
            printf("E7 VAD recovery: sn10 pause frames from 7.5 s decoded at %+.1f dB\n", sup);
        } else {
            printf("E7 VAD recovery: %s %+.1f dB (%+.1f vs sn10, limit +3)\n", names[t], sup, sup - ref);
            e7_check(sup - ref <= 3.0, "noise level not recovered after a VAD trigger");
        }
        free_result(&r);
    }
    free(sn);
    free(clean);
    free(x);
    free(cl);
}

/* E7 robustness, mandatory:
   digital silence decodes to -85 dBFS or less
   clipped full-scale speech, sweeps, and stationary inputs keep E1/E2 (checked on every frame by run_case)
   a stationary tone is absorbed (tone frames, or comfort noise after about 3 s without the tone detector)
   ten minutes of mixed abuse leave the state finite and the codec quiet again on digital silence
   informative: the sweeps' pitch accuracy */
static void test_e7(const char *vec, const char *inputs, double stress_seconds)
{
    const long lead = FS / 2, lead_frames = lead / 160;
    long n = 10 * FS, i, es1_n;
    int16_t *x = calloc((size_t)(30 * FS), sizeof(int16_t)), *es1 = load_raw(vec, "e3/es1_orig.raw", &es1_n);
    encode_result r;
    xombe_enc_cfg cfg;
    double lvl;
    long tone_frames = 0, other_frames = 0, absorbed_at = -1;
    xombe_enc_cfg_default(&cfg);

    run_case("digital silence", x, n, &cfg, &r);
    lvl = level_dbfs(r.pcm, 160 * r.frames);
    printf("E7 digital silence: decoded %.1f dBFS (limit -85)\n", lvl);
    e7_check(lvl <= -85.0, "digital silence above -85 dBFS");
    free_result(&r);

    for (i = 0; i < es1_n; i++)
        x[i] = clip16(16.0 * es1[i]);
    run_case("clipped full-scale speech", x, es1_n, &cfg, &r);
    printf("E7 clipped full-scale speech: decoded %.1f dBFS, %ld frames\n", level_dbfs(r.pcm, 160 * r.frames),
           r.frames);
    free_result(&r);

    sweep_report("rising AM comb sweep 50->500 Hz", 50.0, 500.0, x, 20 * FS);
    sweep_report("falling AM comb sweep 500->50 Hz", 500.0, 50.0, x, 20 * FS);

    /* tone detector on: every frame from two frames after the onset is a tone frame */
    memset(x, 0, sizeof(int16_t) * 30 * FS);
    for (i = lead; i < 8 * FS; i++)
        x[i] = clip16(8000.0 * sin(2.0 * XOMBE_PI * 1000.0 * i / FS));
    run_case("stationary 1 kHz tone", x, 8 * FS, &cfg, &r);
    for (i = lead_frames + 2; i < r.frames; i++) {
        tone_frames += r.cls[i] == XOMBE_CLASS_TONE;
        other_frames += r.cls[i] != XOMBE_CLASS_TONE;
    }
    printf("E7 stationary 1 kHz tone (7.5 s): %ld tone frames, %ld others after onset + 2 frames\n", tone_frames,
           other_frames);
    e7_check(other_frames == 0, "stationary tone not coded as tone frames throughout");
    free_result(&r);

    /* tone detector off: the steady comb becomes background (comfort noise) about 3 s after its onset, at the
       latest 4 s after it, and stays background */
    cfg.tones = 0;
    am_sweep(x, 8 * FS, lead, 150.0, 150.0, NULL);
    for (i = lead; i < 8 * FS; i++)
        x[i] = clip16(x[i] / (1.0 + 0.5 * sin(2.0 * XOMBE_PI * 3.5 * (double)(i - lead) / FS)));
    run_case("stationary 150 Hz comb", x, 8 * FS, &cfg, &r);
    for (i = r.frames - 1; i >= lead_frames && r.cls[i] == XOMBE_CLASS_SILENCE; i--)
        absorbed_at = i;
    printf("E7 stationary 150 Hz comb, tone detector off: background from %.2f s after onset\n",
           absorbed_at < 0 ? -1.0 : (absorbed_at - lead_frames) * 0.02);
    e7_check(absorbed_at >= 0 && absorbed_at - lead_frames <= 4 * 50, "steady comb not absorbed within 4 s");
    free_result(&r);
    free(x);

    test_vad_recovery(vec, inputs);

    if (stress_seconds > 0.0) {
        long total = (long)(stress_seconds * FS), pos = 0, bad_samples = 0, quiet_start = -1, quiet_len = 0;
        int16_t *s = malloc(sizeof(int16_t) * (size_t)total);
        int seg = 0;
        xombe_enc_cfg_default(&cfg);
        while (pos < total) {
            long len = (long)(2 + rnd() % 9) * FS, k;
            int kind = seg++ % 6;
            if (pos + len > total)
                len = total - pos;
            for (k = 0; k < len; k++) {
                double v;
                switch (kind) {
                case 0: v = es1[k % es1_n]; break;
                case 1: v = 3000.0 * (urand() - 0.5); break;
                case 2: v = 16.0 * es1[k % es1_n]; break;
                case 3: v = 9000.0 * sin(2.0 * XOMBE_PI * (300.0 + 2700.0 * k / len) * k / FS); break;
                case 4: v = 0.0; break;
                default: v = es1[k % es1_n] + 2000.0 * (urand() - 0.5); break;
                }
                s[pos + k] = clip16(v);
            }
            if (kind == 4 && len >= 2 * FS) {
                quiet_start = pos;
                quiet_len = len;
            }
            pos += len;
        }
        run_case("stress", s, total, &cfg, &r);
        for (i = 0; i < 160 * r.frames; i++)
            bad_samples += r.pcm[i] == 32767 || r.pcm[i] == -32768;
        printf("E7 stress: %.0f s, %ld frames, %ld saturated output samples (informative), state %s\n",
               stress_seconds, r.frames, bad_samples, r.state_finite ? "finite" : "NOT FINITE");
        e7_check(r.state_finite, "encoder or decoder state not finite after the stress run");
        if (quiet_start >= 0) {
            /* last digital-silence stretch, from 1 s into it to its end */
            long a = (quiet_start + FS) / 160, b = (quiet_start + quiet_len) / 160;
            lvl = level_dbfs(r.pcm + 160 * a, 160 * (b - a));
            printf("E7 stress: last digital silence (%.1f s) decodes to %.1f dBFS (limit -85)\n",
                   (b - a) * 0.02, lvl);
            e7_check(lvl <= -85.0, "digital silence after the stress run above -85 dBFS");
        }
        free_result(&r);
        free(s);
    }
    free(es1);
    printf("E7: %s\n", e7_fail ? "FAIL" : "PASS");
}

int main(int argc, char **argv)
{
    const char *vec = argc > 1 ? argv[1] : "vectors", *inputs = NULL;
    int i, e7 = 0;
    double stress = 0.0;
    for (i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "--e7"))
            e7 = 1;
        else if (!strcmp(argv[i], "--stress") && i + 1 < argc)
            stress = atof(argv[++i]);
        else if (!strcmp(argv[i], "--inputs") && i + 1 < argc)
            inputs = argv[++i];
    }
    if (e7 && !inputs) {
        fprintf(stderr, "test_mirror: --e7 needs --inputs DIR, the dithered E3 inputs (scripts/gen_inputs.py)\n");
        return 2;
    }
    test_files(vec);
    test_forced(vec);
    test_idempotence(5000);
    if (e7)
        test_e7(vec, inputs, stress);
    printf("E1: %ld/%ld frames valid\n", e1_frames - e1_bad, e1_frames);
    printf("E2 mirror: %ld/%ld frames equal (max |error| %.3g log2)\n", mirror_frames - mirror_bad, mirror_frames,
           mirror_max_err);
    printf("E2 forced tones: %ld requests among the frames above\n", forced_frames);
    printf("E2 idempotence: %ld frames, %ld identical, %ld equivalent (unused HOC field), %ld clamped, %ld mismatches\n",
           idem_frames, idem_exact, idem_equivalent, idem_clamped, idem_bad);
    printf("E1/E2: %s\n", e1_bad || mirror_bad || idem_bad ? "FAIL" : "PASS");
    return e1_bad || mirror_bad || idem_bad || e7_fail ? 1 : 0;
}
