/* AMBE+2 tone detector, reference vs merged (xombe__tone_detect over the generic search), on every segment position the
   encoder reads: decoded tone frames (T5), synthetic E6 cases, the tone vectors, and speech and other non-tone
   vectors. --full takes every AD in T5 (default every 6th), --list prints each disagreement, --time the slowest
   segment of each detector. exits 1 on a regression (see verdict) */
#include "../src/enc/tone_detect.c" /* its fit helpers are static */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "ref.h"

#define MAXFILES 64
#define T5_FRAMES 8
#define T5_FIRST 3 /* decoded frames before the first segment */

/* reference: AMBE+2 detector before the merge into xombe__tone_search, per-ID filter gain as tabulated then */
static float ref_hpf_db[XOMBE_TONE_IDS];

#define DUAL_REFINE_SPAN 0.025f
#define DUAL_TRIALS 24 /* 2 passes x 2 components x 3 steps x 2 signs */

/* the reference's float front end and correlation, as the float build has them */
static void ref_correlate(const float *x, int len, float hz, component *c)
{
    float w = 2.0f * XOMBE_PI_F * hz / XOMBE_FS, hc = cosf(0.5f * w), hs = sinf(0.5f * w);
    float cw = hc * hc - hs * hs, s = 0.0f, d = 0.0f, prev;
    int n = len - 1;
    if (cw >= 0.0f) {
        float l = -4.0f * hs * hs;
        for (; n > 0; n -= 2) {
            d += x[n] + l * s;
            s += d;
            d += x[n - 1] + l * s;
            s += d;
        }
        if (n == 0) {
            d += x[0] + l * s;
            s += d;
        }
        prev = s - d;
    } else {
        float l = 4.0f * hc * hc;
        for (; n > 0; n -= 2) {
            d = x[n] + l * s - d;
            s = d - s;
            d = x[n - 1] + l * s - d;
            s = d - s;
        }
        if (n == 0) {
            d = x[0] + l * s - d;
            s = d - s;
        }
        prev = d - s;
    }
    c->hz = hz;
    c->rc = s - cw * prev;
    c->rs = 2.0f * hs * hc * prev;
    c->h[0] = hc;
    c->h[1] = hs;
    c->m[0] = cosf((float)(2 * len - 1) * 0.5f * w);
    c->m[1] = sinf((float)(2 * len - 1) * 0.5f * w);
}

static void ref_fit_component(const float *x, int len, float hz, component *c)
{
    ref_correlate(x, len, hz, c);
    gram(c, c, len, c->g);
}

static float ref_prepare(xombe_enc *e, const int16_t *seg, float *x, float pw[XOMBE_TONE_FFT / 2 + 1])
{
    static const float min_amp = 32768.0f * 0.0031622777f;
    xombe_cpx *X = e->tone_spec;
    float mean, energy = 0.0f;
    long sum = 0;
    int n, k;
    for (n = 0; n < XOMBE_TONE_SEG; n++)
        sum += seg[n];
    mean = (float)sum / XOMBE_TONE_SEG;
    for (n = 0; n < XOMBE_TONE_SEG; n++) {
        x[n] = (float)seg[n] - mean;
        energy += x[n] * x[n];
    }
    if (energy / XOMBE_TONE_SEG < 0.5f * min_amp * min_amp)
        return 0.0f;
    memset(e->tone_spec, 0, sizeof(e->tone_spec));
    for (n = 0; n < XOMBE_TONE_SEG; n++) {
        float v = x[n] * xombe__tone_win[n];
        if (n & 1)
            X[n >> 1].im = v;
        else
            X[n >> 1].re = v;
    }
    xombe__rfft_packed(X, XOMBE_TONE_FFT, xombe__fft_tw, XOMBE_TW_N);
    for (k = 0; k <= XOMBE_TONE_FFT / 2; k++)
        pw[k] = X[k].re * X[k].re + X[k].im * X[k].im;
    return energy;
}

static float ref_band_total(const float *pw)
{
    float total = 0.0f;
    int k;
    for (k = (int)ceilf(GATE_LO_HZ / BIN_HZ); k <= (int)floorf(GATE_HI_HZ / BIN_HZ); k++)
        total += pw[k];
    return total;
}

static int ref_gate_pass(const float *pw, const float nominal[2], float total)
{
    float near = 0.0f;
    int k, r0lo, r0hi, r1lo, r1hi;
    gate_run(nominal[0], &r0lo, &r0hi);
    gate_run(nominal[1], &r1lo, &r1hi);
    for (k = 0; k <= XOMBE_TONE_FFT / 2; k++)
        if ((k >= r0lo && k <= r0hi) || (k >= r1lo && k <= r1hi))
            near += pw[k];
    return !(near < DUAL_GATE_SHARE * total);
}

/* returns 1 if seen[n] repeats an earlier pair, already rejected since none beats the current best */
static int visited(const float seen[][2], int n)
{
    int i;
    for (i = 0; i < n; i++)
        if (seen[i][0] == seen[n][0] && seen[i][1] == seen[n][1])
            return 1;
    return 0;
}

static int ref_detect(xombe_enc *e, const int16_t *seg, int *id, int *ad)
{
    float x[XOMBE_TONE_SEG], energy, total, best_dual = -1.0f, e1, amp1[2] = {0.0f, 0.0f};
    component single;
    float pw[XOMBE_TONE_FFT / 2 + 1];
    int k, peak = 2, single_id = -1, dual_id = -1, dual_ad = 0, cand;
    energy = ref_prepare(e, seg, x, pw);
    if (energy == 0.0f)
        return 0;
    for (k = 3; k < XOMBE_TONE_FFT / 2; k++)
        if (pw[k] > pw[peak])
            peak = k;
    {
        float a = 0.5f * logf(pw[peak - 1] + 1e-30f), b = 0.5f * logf(pw[peak] + 1e-30f);
        float c = 0.5f * logf(pw[peak + 1] + 1e-30f), den = a - 2.0f * b + c;
        float delta = den != 0.0f ? 0.5f * (a - c) / den : 0.0f;
        ref_fit_component(x, XOMBE_TONE_SEG, ((float)peak + delta) * BIN_HZ, &single);
    }
    e1 = fit(energy, &single, 1, XOMBE_TONE_SEG, amp1);
    cand = (int)floorf(single.hz / SINGLE_ID_HZ + 0.5f);
    if (e1 >= EXPLAINED_MIN && cand >= 5 && cand <= 122)
        single_id = cand;
    total = ref_band_total(pw);
    for (cand = 128; cand <= 163; cand++) {
        float lo, hi, nominal[2], amp[2], ex;
        component comp[2];
        float seen[DUAL_TRIALS + 1][2];
        int pass, c, s, nseen;
        if (xombe__tone_freqs(cand, &lo, &hi) != 2)
            continue;
        nominal[0] = lo;
        nominal[1] = hi;
        if (!ref_gate_pass(pw, nominal, total))
            continue;
        ref_fit_component(x, XOMBE_TONE_SEG, nominal[0], &comp[0]);
        ref_fit_component(x, XOMBE_TONE_SEG, nominal[1], &comp[1]);
        ex = fit(energy, comp, 2, XOMBE_TONE_SEG, amp);
        seen[0][0] = comp[0].hz;
        seen[0][1] = comp[1].hz;
        nseen = 1;
        for (pass = 0; pass < 2; pass++)
            for (c = 0; c < 2; c++)
                for (s = 0; s < 3; s++) {
                    static const float steps[3] = {0.01f, 0.005f, 0.0025f};
                    int sign;
                    for (sign = -1; sign <= 1; sign += 2) {
                        component keep = comp[c];
                        float hz = keep.hz + (float)sign * steps[s] * nominal[c], trial_amp[2] = {0.0f, 0.0f}, te;
                        if (fabsf(hz - nominal[c]) > DUAL_REFINE_SPAN * nominal[c] * (1.0f + 1e-6f))
                            continue;
                        seen[nseen][c] = hz;
                        seen[nseen][1 - c] = comp[1 - c].hz;
                        if (visited(seen, nseen))
                            continue;
                        nseen++;
                        ref_fit_component(x, XOMBE_TONE_SEG, hz, &comp[c]);
                        te = fit(energy, comp, 2, XOMBE_TONE_SEG, trial_amp);
                        if (te > ex) {
                            ex = te;
                            amp[0] = trial_amp[0];
                            amp[1] = trial_amp[1];
                        } else {
                            comp[c] = keep;
                        }
                    }
                }
        if (ex >= EXPLAINED_MIN && fabsf(20.0f * log10f((amp[0] + 1e-30f) / (amp[1] + 1e-30f))) <= DUAL_TWIST_DB &&
            ex > best_dual) {
            best_dual = ex;
            dual_id = cand;
            dual_ad = ad_from_level(0.5f * (peak_db(amp[0]) + peak_db(amp[1])) - ref_hpf_db[cand]);
        }
    }
    if (dual_id >= 0 && (single_id < 0 || best_dual > e1 + DUAL_PREFERENCE)) {
        *id = dual_id;
        *ad = dual_ad;
        return 1;
    }
    if (single_id >= 0) {
        *id = single_id;
        *ad = ad_from_level(peak_db(amp1[0]) - ref_hpf_db[single_id]);
        return 1;
    }
    return 0;
}

/* the detector under test */
static int new_detect(xombe_enc *e, const int16_t *seg, int *id, int *ad)
{
    return xombe__tone_detect(e, seg, id, ad);
}

typedef struct {
    const char *name;
    long segs, det[2], id_diff, ad_diff, max_dad, only[2], frames[2], truth, miss[2], wrong[2], ad_off[2];
} stats;

static struct xombe_enc enc_ref, enc_new;
static int list, timing;
static double slow_t[2];
static char slow_at[2][128];

static double now(void)
{
    return (double)clock() / CLOCKS_PER_SEC;
}

static void reset_decide(void)
{
    enc_ref.tone_prev_id = enc_new.tone_prev_id = -1;
    enc_ref.tone_confirmed = enc_new.tone_confirmed = 0;
    enc_ref.tone_hold_used = enc_new.tone_hold_used = 0;
    enc_ref.tone_last_id = enc_new.tone_last_id = 0;
    enc_ref.tone_last_ad = enc_new.tone_last_ad = 0;
}

/* both detectors on one segment, true_id > 0 names the tone it holds (0: no tone, -1: unknown) */
static void one(stats *s, const int16_t *seg, int true_id, int true_ad, const char *where, long pos)
{
    struct xombe_enc *es[2] = {&enc_ref, &enc_new};
    int det[2], id[2] = {-1, -1}, ad[2] = {0, 0}, k;
    for (k = 0; k < 2; k++) {
        double t0 = timing ? now() : 0.0;
        det[k] = (k ? new_detect : ref_detect)(es[k], seg, &id[k], &ad[k]);
        if (timing) {
            double t = now() - t0;
            if (t > slow_t[k]) {
                slow_t[k] = t;
                snprintf(slow_at[k], sizeof(slow_at[k]), "%s @%ld", where, pos);
            }
        }
        if (!det[k])
            id[k] = -1;
        s->det[k] += det[k];
        {
            int tid, tad;
            s->frames[k] += xombe__tone_decide(es[k], det[k], id[k], ad[k], &tid, &tad);
        }
        if (true_id > 0) {
            if (!det[k])
                s->miss[k]++;
            else if (id[k] != true_id)
                s->wrong[k]++;
            else if (true_ad >= 0 && ad[k] != true_ad)
                s->ad_off[k]++;
        }
    }
    s->segs++;
    s->truth += true_id > 0;
    if (det[0] != det[1]) {
        s->only[det[1]]++;
    } else if (det[0] && id[0] != id[1]) {
        s->id_diff++;
    } else if (det[0] && ad[0] != ad[1]) {
        int d = abs(ad[0] - ad[1]);
        s->ad_diff++;
        s->max_dad = d > s->max_dad ? d : s->max_dad;
    } else {
        return;
    }
    if (list)
        printf("  %s @%ld: ref %d/%d new %d/%d (true %d/%d)\n", where, pos, id[0], ad[0], id[1], ad[1], true_id,
               true_ad);
}

static void report(const stats *s)
{
    printf("%-8s %7ld segs  det %6ld %6ld  tone frames %6ld %6ld  ref-only %4ld new-only %4ld  ID diff %4ld  "
           "AD diff %4ld (max %ld)",
           s->name, s->segs, s->det[0], s->det[1], s->frames[0], s->frames[1], s->only[0], s->only[1], s->id_diff,
           s->ad_diff, s->max_dad);
    if (s->truth)
        printf("\n%-8s %7ld true  missed %4ld %4ld  wrong ID %4ld %4ld  AD off %4ld %4ld", "", s->truth, s->miss[0],
               s->miss[1], s->wrong[0], s->wrong[1], s->ad_off[0], s->ad_off[1]);
    printf("\n");
}

static int16_t *load(const char *path, long *n)
{
    FILE *f = fopen(path, "rb");
    long size;
    int16_t *x;
    if (!f) {
        perror(path);
        exit(2);
    }
    fseek(f, 0, SEEK_END);
    size = ftell(f);
    fseek(f, 0, SEEK_SET);
    *n = size / 2;
    x = calloc((size_t)*n + XOMBE_TONE_SEG, sizeof(int16_t)); /* zero tail, as the encoder flushes */
    if (fread(x, sizeof(int16_t), (size_t)*n, f) != (size_t)*n)
        exit(2);
    fclose(f);
    return x;
}

/* every segment at 160k of x[0..n), zeros past the end; tone true_id over [on, off) */
static void run_signal(stats *s, const int16_t *x, long n, const char *name, int true_id, int true_ad, long on,
                       long off)
{
    long p;
    reset_decide();
    for (p = 0; p < n; p += XOMBE_FRAME) {
        int t = true_id > 0 && p >= on && p + XOMBE_TONE_SEG <= off ? true_id : true_id > 0 ? -1 : true_id;
        one(s, x + p, t, true_ad, name, p);
    }
}

static void run_file(stats *s, const char *path, int tone)
{
    long n;
    int16_t *x = load(path, &n);
    run_signal(s, x, n, path, tone ? -1 : 0, 0, 0, 0);
    free(x);
}

static void run_t5(stats *s, int ad_step)
{
    static int16_t x[T5_FRAMES * XOMBE_FRAME];
    int id, ad, f;
    for (id = 5; id <= 163; id++) {
        float lo, hi;
        if (!xombe__tone_id_accepted(id) || !xombe__tone_freqs(id, &lo, &hi))
            continue;
        for (ad = 30; ad <= 120; ad += ad_step) {
            xombe_dec dec;
            uint8_t d[49], bytes[7];
            char where[32];
            xombe__tone_frame_bits(id, ad, d);
            xombe__pack49(d, bytes);
            xombe_dec_set_noise(xombe_dec_init(&dec, XOMBE_DMR), XOMBE_NOISE_CONFORMANCE, 0);
            for (f = 0; f < T5_FRAMES; f++)
                xombe_dec_frame_bits(&dec, bytes, 0, x + f * XOMBE_FRAME);
            snprintf(where, sizeof(where), "T5 id %d ad %d", id, ad);
            reset_decide();
            for (f = T5_FIRST; (f + 1) * XOMBE_FRAME + XOMBE_TONE_SEG - XOMBE_FRAME <= T5_FRAMES * XOMBE_FRAME; f++)
                one(s, x + f * XOMBE_FRAME, id, ad, where, f * XOMBE_FRAME);
        }
    }
}

/* tools/e6.py gen: DTMF at -16 dBFS, Knox and call progress at -20, single tones at 300, 1000, 2000, and 3000 Hz
   from -6 to -46 dBFS, each 0.4 s between 0.5 s and 0.3 s of silence, phases from a fixed LCG */
static void run_e6(stats *s)
{
    static int16_t x[8000 + 1600 + XOMBE_TONE_SEG];
    static const int singles[4] = {300, 1000, 2000, 3000};
    unsigned lcg = 20261019u;
    int c, n, len = 3200, on = 4000, total = on + len + 2400;
    for (c = 0; c < 36 + 20; c++) {
        double f[2], lev, ph[2];
        int id, nf, i;
        char where[32];
        float lo, hi;
        if (c < 36) {
            id = 128 + c;
            xombe__tone_freqs(id, &lo, &hi);
            f[0] = lo;
            f[1] = hi;
            nf = 2;
            lev = id < 144 ? -16.0 : -20.0;
        } else {
            f[0] = singles[(c - 36) / 5];
            lev = -6.0 - 10.0 * ((c - 36) % 5);
            id = (int)floor(f[0] / 31.25 + 0.5);
            nf = 1;
        }
        for (i = 0; i < nf; i++) {
            lcg = lcg * 1664525u + 1013904223u;
            ph[i] = 2.0 * XOMBE_PI * (lcg >> 8) / 16777216.0;
        }
        memset(x, 0, sizeof(x));
        for (n = 0; n < len; n++) {
            double v = 0.0;
            for (i = 0; i < nf; i++)
                v += 32768.0 * pow(10.0, lev / 20.0) * cos(2.0 * XOMBE_PI * f[i] / XOMBE_FS * n + ph[i]);
            v = floor(v + 0.5);
            x[on + n] = (int16_t)(v > 32767 ? 32767 : v < -32768 ? -32768 : v);
        }
        snprintf(where, sizeof(where), "E6 id %d %.0f dBFS", id, lev);
        run_signal(s, x, total, where, id, -1, on, on + len);
    }
}

int main(int argc, char **argv)
{
    static const char *tone_vec[] = {"e7/tones.raw", "dstar/tn.in.raw", "dstar/kx.in.raw", "dstar/tn.dec.raw",
                                     "dstar/kx.dec.raw"};
    static const char *speech_vec[] = {
        "e3/es1_orig.raw",       "e3/es2_orig.raw",          "e3/es3_orig.raw",          "e3/sn10_orig.raw",
        "e3/ttsc_orig.raw",      "e3/ttsf_orig.raw",         "e3/ttsl_orig.raw",
        "perceptual/es2_orig.raw", "perceptual/es3_orig.raw", "perceptual/es4_orig.raw",
        "perceptual/es1_nsdtx_off_ref.raw", "perceptual/es1_nsdtx_on_ref.raw", "perceptual/es2_ref.raw",
        "perceptual/es3_ref.raw", "perceptual/es4_ref.raw", "e5/es1_lead0.5s.raw", "e5/hp5.raw", "e5/low5.raw",
        "e7/glide.raw", "e7/mix.raw", "e7/saw150.raw", "dstar/sp.in.raw", "dstar/nz.in.raw", "dstar/sp.dec.raw",
        "dstar/sper.dec.raw", "dstar/nz.dec.raw"};
    const char *extra_tone[MAXFILES], *extra_speech[MAXFILES];
    stats t5, e6, tv, sp;
    char path[4096];
    int i, full = 0, ntone = 0, nspeech = 0, bad = 0;
    if (argc < 2) {
        fprintf(stderr, "usage: tone_cmp VEC [--full] [--list] [--time] [--tones FILE]... [--speech FILE]...\n");
        return 2;
    }
    for (i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "--full"))
            full = 1;
        else if (!strcmp(argv[i], "--list"))
            list = 1;
        else if (!strcmp(argv[i], "--time"))
            timing = 1;
        else if (!strcmp(argv[i], "--tones") && i + 1 < argc && ntone < MAXFILES)
            extra_tone[ntone++] = argv[++i];
        else if (!strcmp(argv[i], "--speech") && i + 1 < argc && nspeech < MAXFILES)
            extra_speech[nspeech++] = argv[++i];
        else
            return 2;
    }
    memset(&t5, 0, sizeof(t5));
    memset(&e6, 0, sizeof(e6));
    memset(&tv, 0, sizeof(tv));
    memset(&sp, 0, sizeof(sp));
    t5.name = "T5";
    e6.name = "E6";
    tv.name = "tones";
    sp.name = "speech";
    for (i = 0; i < XOMBE_TONE_IDS; i++) {
        float lo, hi;
        ref_hpf_db[i] = xombe__tone_freqs(i, &lo, &hi) ? (float)(0.5 * (ref_hpf_gain_db(lo) + ref_hpf_gain_db(hi)))
                                                         : 0.0f;
    }
    run_t5(&t5, full ? 1 : 6);
    run_e6(&e6);
    for (i = 0; i < (int)(sizeof(tone_vec) / sizeof(tone_vec[0])); i++) {
        snprintf(path, sizeof(path), "%s/%s", argv[1], tone_vec[i]);
        run_file(&tv, path, 1);
    }
    for (i = 0; i < ntone; i++)
        run_file(&tv, extra_tone[i], 1);
    for (i = 0; i < (int)(sizeof(speech_vec) / sizeof(speech_vec[0])); i++) {
        snprintf(path, sizeof(path), "%s/%s", argv[1], speech_vec[i]);
        run_file(&sp, path, 0);
    }
    for (i = 0; i < nspeech; i++)
        run_file(&sp, extra_speech[i], 0);
    printf("columns: reference, new\n");
    report(&t5);
    report(&e6);
    report(&tv);
    report(&sp);
    if (timing)
        printf("slowest: reference %.1f us (%s), new %.1f us (%s)\n", 1e6 * slow_t[0], slow_at[0], 1e6 * slow_t[1],
               slow_at[1]);
    /* verdict: no new speech detections, and no T5 or E6 tone the reference finds that the new one misses */
    bad |= sp.only[1] != 0 || sp.frames[1] > sp.frames[0];
    bad |= t5.miss[1] + t5.wrong[1] > t5.miss[0] + t5.wrong[0];
    bad |= e6.miss[1] + e6.wrong[1] > e6.miss[0] + e6.wrong[0];
    printf("tone_cmp: %s\n", bad ? "FAIL" : "PASS");
    return bad;
}
