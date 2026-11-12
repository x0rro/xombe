#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "xombe.h"
#include "../tools/io.h"

#define MAX_CALLS 20000
#define LEVEL_FLOOR_DB (-80.0)

enum { CAT_VOICE, CAT_SILENCE, CAT_TONE, CAT_REPEAT, CAT_MUTE, NCAT };
static const char *const cat_names[NCAT] = {"voice", "silence", "tone", "repeat", "mute"};
static const double lim_median[NCAT] = {0.5, 0.5, 0.25, 0.75, 0.25};
static const double lim_p90[NCAT] = {7.0, 4.0, 1.0, 7.0, 1.0};
static const double lim_share[NCAT] = {0.25, 0.12, 0.03, 0.25, 0.08};

typedef struct {
    enum xombe_class cls;
    enum xombe_conceal conceal;
    int run;
    long mute_src;
    int tone_id, tone_ad;
} frame_info;

static double deltas[NCAT][MAX_CALLS];
static int ndeltas[NCAT];
static double snr_mute[MAX_CALLS], snr_cn[MAX_CALLS];
static int nsnr_mute, nsnr_cn;
static int steady_fail, steady_mute_n, steady_tone_n;
static double steady_mute_max, steady_tone_max;
static long frames_total, decision_fail;

static enum xombe_class parse_class(const char *p)
{
    if (!strncmp(p, "\"voice\"", 7))
        return XOMBE_CLASS_VOICE;
    if (!strncmp(p, "\"erasure\"", 9))
        return XOMBE_CLASS_ERASURE;
    if (!strncmp(p, "\"silence\"", 9))
        return XOMBE_CLASS_SILENCE;
    return XOMBE_CLASS_TONE;
}

static enum xombe_conceal parse_conceal(const char *p)
{
    if (!strncmp(p, "\"none\"", 6))
        return XOMBE_CONCEAL_NONE;
    if (!strncmp(p, "\"repeat\"", 8))
        return XOMBE_CONCEAL_REPEAT;
    return XOMBE_CONCEAL_MUTE;
}

static double level_db(const int16_t *x, int n)
{
    double s = 0.0;
    int i;
    for (i = 0; i < n; i++)
        s += (double)x[i] * x[i];
    if (s <= 0.0)
        return -1e9;
    return 10.0 * log10(s / n) - 20.0 * log10(32768.0);
}

static int cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return x < y ? -1 : x > y;
}

static double percentile(double *v, int n, double q)
{
    double pos;
    int i;
    if (n == 0)
        return 0.0;
    qsort(v, (size_t)n, sizeof(double), cmp_double);
    pos = q * (n - 1);
    i = (int)floor(pos);
    if (i >= n - 1)
        return v[n - 1];
    return v[i] + (pos - i) * (v[i + 1] - v[i]);
}

static int expect_field(const stream_reader *r, const char *e, const char *key, const char **out)
{
    *out = e ? stream_key(e, key) : NULL;
    if (*out)
        return 1;
    fprintf(stderr, "%s:%ld: missing \"expect\".\"%s\"\n", r->path, r->lineno, key);
    return 0;
}

static int run_file(const char *path, int verbose)
{
    stream_reader r;
    stream_half h;
    int have_prev = 0, ev, ok = 1;
    long frame = -1;
    frame_info cur, prev;
    xombe_dec_mem mem;
    xombe_dec *dec;
    memset(&cur, 0, sizeof(cur));
    memset(&prev, 0, sizeof(prev));
    if (!stream_open(&r, path))
        return 1;
    dec = xombe_dec_init(&mem, XOMBE_DMR);
    xombe_dec_set_noise(dec, XOMBE_NOISE_CONFORMANCE, 0);
    while (ok && (ev = stream_next(&r, &h)) != STREAM_END) {
        if (ev == STREAM_ERROR) {
            ok = 0;
            break;
        }
        if (ev == STREAM_SET) {
            xombe_dec_set_errors(dec, r.eps_t, r.eps_1);
            continue;
        }
        if (ev == STREAM_INIT) {
            dec = xombe_dec_init(&mem, XOMBE_DMR);
            xombe_dec_set_noise(dec, XOMBE_NOISE_CONFORMANCE, 0);
            frame = -1;
            have_prev = 0;
            continue;
        }
        {
            const uint8_t *d = h.d;
            int16_t pcm[XOMBE_HALF_MAX], ref[XOMBE_HALF_MAX];
            int i, n, nref = 0;
            const char *q;
            double ld, lr, dl, num = 0.0, den = 0.0;
            int cat;
            n = xombe_dec_half_bits(dec, h.bits, h.flags, h.second_half, h.nsamples, pcm);
            if (n < 0) {
                fprintf(stderr, "%s:%ld: half-frame call out of order\n", path, r.lineno);
                ok = 0;
                break;
            }
            if (!h.second_half) {
                xombe_decision got;
                const char *e = stream_key(h.line, "expect"), *cls, *conceal;
                long run;
                frame++;
                frames_total++;
                if (have_prev || frame > 0)
                    prev = cur;
                have_prev = frame > 0;
                memset(&cur, 0, sizeof(cur));
                if (!expect_field(&r, e, "class", &cls) || !expect_field(&r, e, "conceal", &conceal) ||
                    !stream_int(&r, e, "run", &run)) {
                    ok = 0;
                    break;
                }
                cur.cls = parse_class(cls);
                cur.conceal = parse_conceal(conceal);
                cur.run = (int)run;
                cur.mute_src = -2;
                if ((q = stream_key(e, "mute_src")))
                    cur.mute_src = *q == '"' ? -1 : atol(q);
                for (i = 0; i < 8; i++)
                    cur.tone_id = (cur.tone_id << 1) | d[12 + i];
                for (i = 0; i < 6; i++)
                    cur.tone_ad = (cur.tone_ad << 1) | d[6 + i];
                cur.tone_ad = (cur.tone_ad << 1) | d[44];
                xombe_dec_decision_get(dec, &got);
                if (got.cls != cur.cls || got.conceal != cur.conceal || got.run != cur.run ||
                    (cur.conceal == XOMBE_CONCEAL_MUTE && got.mute_src != cur.mute_src)) {
                    decision_fail++;
                    if (verbose && decision_fail <= 20)
                        fprintf(stderr, "%s frame %ld: expect %d/%d/%d/%ld got %d/%d/%d/%lld\n", path, frame,
                                cur.cls, cur.conceal, cur.run, cur.mute_src, got.cls, got.conceal, got.run,
                                (long long)got.mute_src);
                }
            }
            q = stream_key(h.line, "pcm");
            if (q && (q = strchr(q, '['))) {
                q++;
                while (nref < XOMBE_HALF_MAX && *q && *q != ']') {
                    char *endp;
                    ref[nref++] = (int16_t)strtol(q, &endp, 10);
                    if (endp == q)
                        break;
                    q = endp;
                    while (*q == ',' || *q == ' ')
                        q++;
                }
            }
            if (nref != n) {
                fprintf(stderr, "%s frame %ld: length %d, reference %d\n", path, frame, n, nref);
                decision_fail++;
                continue;
            }
            ld = level_db(pcm, n);
            lr = level_db(ref, n);
            dl = (ld > LEVEL_FLOOR_DB ? ld : LEVEL_FLOOR_DB) - (lr > LEVEL_FLOOR_DB ? lr : LEVEL_FLOOR_DB);
            if (cur.conceal == XOMBE_CONCEAL_MUTE)
                cat = CAT_MUTE;
            else if (cur.conceal == XOMBE_CONCEAL_REPEAT)
                cat = CAT_REPEAT;
            else
                cat = cur.cls == XOMBE_CLASS_SILENCE ? CAT_SILENCE : cur.cls == XOMBE_CLASS_TONE ? CAT_TONE
                                                                                                  : CAT_VOICE;
            deltas[cat][ndeltas[cat]++] = dl;
            for (i = 0; i < n; i++) {
                num += (double)ref[i] * ref[i];
                den += ((double)ref[i] - pcm[i]) * ((double)ref[i] - pcm[i]);
            }
            if (have_prev) {
                int steady_mute = cur.conceal == XOMBE_CONCEAL_MUTE && prev.conceal == XOMBE_CONCEAL_MUTE &&
                                  prev.mute_src == cur.mute_src;
                int steady_tone = cur.conceal == XOMBE_CONCEAL_NONE && cur.cls == XOMBE_CLASS_TONE &&
                                  prev.conceal == XOMBE_CONCEAL_NONE && prev.cls == XOMBE_CLASS_TONE &&
                                  prev.tone_id == cur.tone_id && prev.tone_ad == cur.tone_ad &&
                                  ((cur.tone_id >= 5 && cur.tone_id <= 122) ||
                                   (cur.tone_id >= 128 && cur.tone_id <= 143));
                int steady_cn = cur.conceal == XOMBE_CONCEAL_NONE && cur.cls == XOMBE_CLASS_SILENCE &&
                                (prev.conceal == XOMBE_CONCEAL_MUTE ||
                                 (prev.conceal == XOMBE_CONCEAL_NONE && prev.cls == XOMBE_CLASS_SILENCE));
                if (steady_mute) {
                    steady_mute_n++;
                    if (fabs(dl) > steady_mute_max)
                        steady_mute_max = fabs(dl);
                    if (fabs(dl) > 2.0)
                        steady_fail++;
                    if (lr >= LEVEL_FLOOR_DB && nsnr_mute < MAX_CALLS)
                        snr_mute[nsnr_mute++] = 10.0 * log10(num / (den > 0 ? den : 1e-9));
                }
                if (steady_tone) {
                    steady_tone_n++;
                    if (fabs(dl) > steady_tone_max)
                        steady_tone_max = fabs(dl);
                    if (fabs(dl) > 1.0)
                        steady_fail++;
                }
                if (steady_cn && nsnr_cn < MAX_CALLS && num > 0)
                    snr_cn[nsnr_cn++] = 10.0 * log10(num / (den > 0 ? den : 1e-9));
            }
        }
    }
    if (!stream_close(&r))
        ok = 0;
    return ok ? 0 : 1;
}

int main(int argc, char **argv)
{
    int i, metrics = 0, fail = 0;
    for (i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--metrics")) {
            metrics = 1;
            continue;
        }
        if (run_file(argv[i], 1))
            return 2;
    }
    printf("T4(a) decisions: %ld/%ld frames match\n", frames_total - decision_fail, frames_total);
    if (decision_fail)
        fail = 1;
    if (metrics) {
        int c;
        double m;
        printf("T4(b) level, pooled  |median|   P90 |dL|  share>3dB   (limits)\n");
        for (c = 0; c < NCAT; c++) {
            double absd[MAX_CALLS], med, p90, share;
            int k, over = 0, ok;
            for (k = 0; k < ndeltas[c]; k++) {
                absd[k] = fabs(deltas[c][k]);
                over += absd[k] > 3.0;
            }
            med = fabs(percentile(deltas[c], ndeltas[c], 0.5));
            p90 = percentile(absd, ndeltas[c], 0.9);
            share = ndeltas[c] ? (double)over / ndeltas[c] : 0.0;
            ok = med <= lim_median[c] && p90 <= lim_p90[c] && share <= lim_share[c];
            printf("  %-8s n=%5d  %6.3f  %8.3f  %8.1f%%   (%.2f / %.1f / %.0f%%) %s\n", cat_names[c], ndeltas[c],
                   med, p90, 100.0 * share, lim_median[c], lim_p90[c], 100.0 * lim_share[c], ok ? "ok" : "FAIL");
            if (!ok)
                fail = 1;
        }
        printf("T4(c) steady state: mute calls %d (max |dL| %.2f dB, limit 2.0), tone calls %d (max %.2f dB, limit 1.0)"
               ", %d violations %s\n",
               steady_mute_n, steady_mute_max, steady_tone_n, steady_tone_max, steady_fail,
               steady_fail ? "FAIL" : "ok");
        if (steady_fail)
            fail = 1;
        m = percentile(snr_mute, nsnr_mute, 0.5);
        printf("T4(d) noise sync: steady mute median SNR %.1f dB over %d calls (>= 15) %s\n", m, nsnr_mute,
               m >= 15.0 ? "ok" : "FAIL");
        if (m < 15.0)
            fail = 1;
        m = percentile(snr_cn, nsnr_cn, 0.5);
        printf("T4(d) noise sync: steady comfort noise median SNR %.1f dB over %d calls (>= 12) %s\n", m, nsnr_cn,
               m >= 12.0 ? "ok" : "FAIL");
        if (m < 12.0)
            fail = 1;
    }
    printf("T4%s: %s\n", metrics ? "" : "(a)", fail ? "FAIL" : "PASS");
    return fail;
}
