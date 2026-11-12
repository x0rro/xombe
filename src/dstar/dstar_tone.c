#include <math.h>
#include <string.h>
#include "dstar_enc.h"
#include "../tables/dstar_enc_tables.h"
#include "../common/tone_table.h"
#include "../common/prof.h"

/* D-STAR encoder tone pipeline
     detect   least-squares detector (xombe__tone_search) on XOMBE_DS_TONE_LEN samples from 160t + TONE_OFF, over
              single tones and the dual tones with a D-STAR code (DTMF, call progress, and 163), then if none a
              second pass over the Knox pairs alone (T = -ID), so a Knox candidate never changes a D-STAR decision
     run rule a tone frame lies in a run of at least MINRUN detections of the same T (MINRUN_SINGLE for single
              tones, whose runs also continue across one 31.25 Hz step), after which HOLD frames without a
              detection repeat the last tone
     in noise with no tone on the frame or either neighbour after the run rule, a DTMF detection on the
              XOMBE_DS_NZ_LEN-sample window from 160n + NZ_OFF at threshold NZ_EMIN is a tone frame when its
              residual is within NZ_RES_DB of the noise floor
   detection level is the AMBE+2 level 127 + L/0.711 + TONE_ADOFF, L the detector's level in dB */

#define TONE_EMIN 0.99f /* minimum explained energy fraction */
#define TONE_OFF 20
#define TONE_ADOFF 0.6f
#define MINRUN 2
#define MINRUN_SINGLE 3
#define HOLD 1
#define SINGLE_T_MIN 37 /* D-STAR T of the lowest single tone */
#define NZ_OFF (-80)
#define NZ_EMIN 0.85f
#define NZ_RES_DB 10.0f
#define KNOX_REST 5.5f /* encoder's other Knox harmonics this far below the tone level, log2 */
#define KNOX_GAIN 0.9f /* target offset, log2 */
#define RING XOMBE_DS_TONE_RING

/* encoded frame is 4 behind the newest, the run rule decides the frame after it from two detections ahead */
XOMBE_STATIC_ASSERT(XOMBE_ENC_DELAY + XOMBE_ENC_LOOKAHEAD == 4, tone_pipeline_matches_the_delay);
XOMBE_STATIC_ASSERT(-NZ_OFF <= XOMBE_ENC_RAW_LEAD, tone_in_noise_window_in_raw_lead);
XOMBE_STATIC_ASSERT(XOMBE_DS_TONE_LEN <= XOMBE_TONE_SEG && XOMBE_DS_NZ_LEN <= XOMBE_TONE_SEG, ds_tone_segments_fit);

/* DTMF 128..143, call progress 160..162, and 163, high component first */
static int dual_pairs(xombe_tone_pair p[20])
{
    int id, n = 0;
    for (id = 128; id <= 162; id++) {
        float lo, hi;
        if (id > 143 && id < 160)
            continue;
        xombe__tone_freqs(id, &lo, &hi);
        p[n].id = (short)id;
        p[n].hz[0] = hi;
        p[n].hz[1] = lo;
        n++;
    }
    /* D-STAR T 19 (AMBE+2 163) is detected as 340 + 550 Hz, not the AMBE+2 table's 350 + 490 Hz */
    p[n].id = 163;
    p[n].hz[0] = 550.0f;
    p[n].hz[1] = 340.0f;
    return n + 1;
}

/* tone-in-noise candidates, the 16 DTMF pairs in D-STAR T order */
static void dtmf_pairs(xombe_tone_pair p[16])
{
    int T;
    for (T = 0; T < 16; T++) {
        float lo, hi;
        p[T].id = xombe__ds_tone_id[T];
        xombe__tone_freqs(p[T].id, &lo, &hi);
        p[T].hz[0] = hi;
        p[T].hz[1] = lo;
    }
}

static void set_tone(xombe_ds_tone *d, int T, float level)
{
    d->on = 1;
    d->T = T;
    d->level = 127.0f + level / XOMBE_TONE_DB_PER_AD + TONE_ADOFF;
}

/* run continues on the same T, or single tones one step apart */
static int same_tone(int a, int b)
{
    return a == b || (a >= SINGLE_T_MIN && b >= SINGLE_T_MIN && (a - b <= 1 && b - a <= 1));
}

void xombe__ds_tone_init(xombe_ds_tones *s)
{
    int id;
    memset(s, 0, sizeof(*s));
    /* Knox IDs that the AMBE+2 tone table gives a harmonic pair for */
    for (id = XOMBE_TONE_KNOX_FIRST; id < XOMBE_TONE_KNOX_FIRST + XOMBE_TONE_KNOX_COUNT; id++) {
        float f0;
        int l1, l2;
        xombe__tone_harmonics(id, &f0, &l1, &l2);
        if (l1 <= 0 || l2 <= 0)
            continue;
        s->knox[s->nknox].id = (short)id;
        s->knox[s->nknox].hz[0] = f0 * (float)l2;
        s->knox[s->nknox].hz[1] = f0 * (float)l1;
        s->knox_f0[s->nknox] = f0 / XOMBE_FS;
        s->knox_l[s->nknox][0] = (uint8_t)l1;
        s->knox_l[s->nknox][1] = (uint8_t)l2;
        s->nknox++;
    }
}

static void detect(xombe_enc *e, xombe_ds_tones *s, int m, const int16_t *seg)
{
    xombe_tone_pair pairs[20];
    xombe_tone_search ts;
    xombe_tone_seg sg;
    xombe_tone_hit h;
    xombe_ds_tone *d = &s->detected[m % RING];
    const xombe_ds_tone *p = &s->detected[(m + RING - 1) % RING];
    int run = 0;
    ts.len = XOMBE_DS_TONE_LEN;
    ts.win = xombe__ds_tone_win;
    ts.emin = TONE_EMIN;
    ts.singles = 1;
    ts.fine_peak = 1;
    ts.pairs = pairs;
    ts.npairs = dual_pairs(pairs);
    d->on = 0;
    xombe__tone_prepare(e, seg, &ts, &sg);
    if (xombe__tone_search(&sg, &ts, &h)) {
        set_tone(d, xombe__ds_tone_code(h.id), h.level);
    } else if (s->nknox) {
        ts.singles = 0;
        ts.pairs = s->knox;
        ts.npairs = s->nknox;
        if (xombe__tone_search(&sg, &ts, &h))
            set_tone(d, -h.id, h.level);
    }
    if (d->on) {
        run = 1;
        if (m > 0 && p->on && same_tone(d->T, p->T)) {
            run = s->run[(m + RING - 1) % RING] + 1;
            run = run > MINRUN_SINGLE ? MINRUN_SINGLE : run;
        }
    }
    s->run[m % RING] = run;
}

/* run rule on frame m, detections of m + 1 and m + 2 already in the ring */
static void decide(xombe_ds_tones *s, int m)
{
    const xombe_ds_tone *d = &s->detected[m % RING];
    xombe_ds_tone *o = &s->confirmed[m % RING];
    o->on = 0;
    if (d->on) {
        int need = d->T >= SINGLE_T_MIN ? MINRUN_SINGLE : MINRUN, len = s->run[m % RING], k;
        for (k = 1; k <= MINRUN_SINGLE - 1 && len < need; k++) {
            const xombe_ds_tone *a = &s->detected[(m + k - 1) % RING], *b = &s->detected[(m + k) % RING];
            if (!b->on || !same_tone(b->T, a->T))
                break;
            len++;
        }
        if (len >= need) {
            *o = *d;
            s->held = *d;
            s->hold_left = HOLD;
        } else {
            s->hold_left = 0;
        }
    } else if (s->hold_left > 0) {
        *o = s->held;
        s->hold_left--;
    }
}

/* frame n's final tone from the run rule or the tone-in-noise pass, raw is frame n's input */
static void finish(xombe_enc *e, xombe_ds_tones *s, int n, const int16_t *raw)
{
    int64_t sum = 0;
    float floor_db;
    int i;
    for (i = 0; i < XOMBE_FRAME; i++)
        sum += (int32_t)raw[i] * raw[i];
    s->level_db[s->level_pos] = 10.0f * log10f((float)sum / XOMBE_FRAME + 1e-9f);
    s->level_pos = (s->level_pos + 1) % (XOMBE_DS_NZ_FLOOR + 1);
    if (s->n_levels <= XOMBE_DS_NZ_FLOOR)
        s->n_levels++;
    s->ndet = s->detected[n % RING];
    s->ntrk = s->confirmed[n % RING];
    s->tone = s->ntrk;
    if (s->tone.on || (n > 0 && s->confirmed[(n + RING - 1) % RING].on) || s->confirmed[(n + 1) % RING].on)
        return;
    {
        xombe_tone_pair pairs[16];
        xombe_tone_search ts;
        xombe_tone_seg sg;
        xombe_tone_hit h;
        dtmf_pairs(pairs);
        ts.len = XOMBE_DS_NZ_LEN;
        ts.win = xombe__ds_nz_win;
        ts.emin = NZ_EMIN;
        ts.singles = 0;
        ts.fine_peak = 0;
        ts.pairs = pairs;
        ts.npairs = 16;
        xombe__tone_prepare(e, raw + NZ_OFF, &ts, &sg);
        if (!xombe__tone_search(&sg, &ts, &h))
            return;
        floor_db = s->level_db[0];
        for (i = 1; i < s->n_levels; i++)
            floor_db = s->level_db[i] < floor_db ? s->level_db[i] : floor_db;
        if (10.0f * log10f((h.res > 1e-9f ? h.res : 1e-9f) / XOMBE_DS_NZ_LEN) > floor_db + NZ_RES_DB)
            return;
        set_tone(&s->tone, xombe__ds_tone_code(h.id), h.level);
    }
}

void xombe__ds_tone_push(xombe_enc *e)
{
    xombe_ds_tones *s = &XOMBE_DS_ENC(e)->tones;
    int t = s->n_frames++;
    if (t >= 1)
        XOMBE_PF(11, detect(e, s, t - 1, xombe__fe_raw(e, 1) + TONE_OFF));
    if (t >= 3)
        decide(s, t - 3);
    if (t >= 4)
        XOMBE_PF(12, finish(e, s, t - 4, xombe__fe_raw(e, 4)));
}

void xombe__ds_enc_force_tone(xombe_enc *e, int id, int ad)
{
    xombe_ds_enc_frame *fr = &XOMBE_DS_ENC(e)->fr;
    int T = xombe__ds_tone_code(id);
    fr->voice = 0;
    fr->tone.on = 1;
    if (T >= 0) {
        fr->cls = XOMBE_DS_CLS_TONE;
        fr->tone.T = T;
        fr->tone.level = XOMBE_DS_AD_LEVEL(ad);
    } else {
        fr->cls = XOMBE_DS_CLS_KNOX;
        fr->tone.T = -id;
        fr->tone.level = (float)ad;
    }
}

void xombe__ds_knox_amps(xombe_ds_enc *s, int k, float AD, int b[9])
{
    const xombe_ds_tones *t = &s->tones;
    float f0 = t->knox_f0[k], lm = XOMBE_TONE_AD_LOG10 * (AD - 127.0f) * XOMBE_LOG2_10, at[XOMBE_LMAX + 1], f0d;
    int L = xombe__L_from_f0(f0), b0 = xombe__ds_knox_b0[t->knox[k].id - XOMBE_TONE_KNOX_FIRST], Ld, l, j;
    /* target peak, log2 min(16384*10^(0.03555*(AD - 127)), 2^14) */
    lm = (lm > 0.0f ? 0.0f : lm) + 14.0f;
    if (t->knox_l[k][0] > L && t->knox_l[k][1] > L)
        lm -= KNOX_REST;
    f0d = xombe__ds_f0(b0);
    Ld = xombe__ds_L(b0);
    for (l = 1; l <= Ld; l++)
        at[l] = lm - KNOX_REST;
    /* each tone on its nearest D-STAR harmonic */
    for (j = 0; j < 2; j++) {
        int h = (int)floorf((float)t->knox_l[k][j] * f0 / f0d + 0.5f);
        if (h >= 1 && h <= Ld)
            at[h] = lm;
    }
    for (l = 1; l <= Ld; l++)
        at[l] += KNOX_GAIN;
    xombe__ds_enc_amps(s, b0, XOMBE_DS_B1_VOICED, at, b);
}
