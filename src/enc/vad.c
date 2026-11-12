#include <math.h>
#include <string.h>
#include "../internal.h"

#define VAD_BAND_OFFSET 1e-3f
#define VAD_FLOOR_REL 1e-3f
#define VAD_FLOOR_ABS 1.0f
#define ABSORB_SPREAD XOMBE_DB(4.0f)
/* digital silence (any constant code): FLAT_RUN equal consecutive input samples in a window */
#define FLAT_RUN 16
/* history: log2 powers in Q10, their sum in linear Q20. smoothed powers lie in [VAD_BAND_OFFSET, 2^30] (int16
   input's total power), so logs fit int16 and 300 linear values int64, and a value leaving the history removes
   exactly what it added */
#define HIST_LOG_ONE 1024
#define HIST_LIN_ONE 1048576.0f
#define VAD_LOW XOMBE_HZ(100.0f)
#define VAD_HIGH XOMBE_HZ(4000.0f)

static void band_powers(const float *P, float E[8])
{
    int k, b;
    for (b = 0; b < 8; b++)
        E[b] = VAD_BAND_OFFSET;
    for (k = 0; k < XOMBE_SPEC_BINS; k++) {
        float f = (float)k / XOMBE_FFT;
        if (f < VAD_LOW || f >= VAD_HIGH)
            continue;
        b = (int)(f * XOMBE_BANDS_PER_UNIT);
        E[b] += P[k];
    }
}

static int16_t hist_log(float x)
{
    return (int16_t)floorf(log2f(x) * HIST_LOG_ONE + 0.5f);
}

static int64_t hist_lin(int16_t q)
{
    return (int64_t)(exp2f((float)q / HIST_LOG_ONE) * HIST_LIN_ONE + 0.5f);
}

static float noise_floor(const float nl[8])
{
    float m = 0.0f;
    int b;
    for (b = 0; b < 8; b++)
        if (nl[b] > m)
            m = nl[b];
    m *= VAD_FLOOR_REL;
    return m > VAD_FLOOR_ABS ? m : VAD_FLOOR_ABS;
}

static float vad_metric(const float E[8], const float nl[8], float fl)
{
    float best = 0.0f, second = 0.0f;
    int b;
    for (b = 0; b < 8; b++) {
        float snr = log2f((E[b] > fl ? E[b] : fl) / (nl[b] > fl ? nl[b] : fl));
        if (snr < 0.0f)
            snr = 0.0f;
        if (snr > best) {
            second = best;
            best = snr;
        } else if (snr > second) {
            second = snr;
        }
    }
    return 0.5f * (best + second);
}

/* minimum-statistics release (Martin, 2001), returns 1 if any band moved. a band above its noise level in every
   half of the 3 s history is locked (its background reads as speech, so the noise update can't raise it), and its
   level rises to the minimum smoothed power over the history. clean speech dips to the background in pauses, so
   the level overshoots only when a pauseless sound fills the whole history */
static int vad_release(xombe_enc *e)
{
    int b, i, moved = 0;
    for (b = 0; b < 8; b++) {
        int lo = INT16_MAX;
        if (e->vad_above[b] < XOMBE_VAD_HIST)
            continue;
        for (i = 0; i < XOMBE_VAD_HIST; i++)
            if (e->vad_hist[i][b] < lo)
                lo = e->vad_hist[i][b];
        if ((float)lo > log2f(e->vad_nl[b]) * HIST_LOG_ONE) {
            e->vad_nl[b] = exp2f((float)lo / HIST_LOG_ONE);
            moved = 1;
        }
    }
    return moved;
}

/* returns 1 if x[0..n) holds FLAT_RUN equal consecutive samples */
int xombe__vad_flat(const int16_t *x, int n)
{
    int i, run = 1;
    for (i = 1; i < n; i++) {
        run = x[i] == x[i - 1] ? run + 1 : 1;
        if (run >= FLAT_RUN)
            return 1;
    }
    return 0;
}

/* flat: window holds digital silence (xombe__vad_flat), which can't lower the noise levels */
void xombe__vad_half(xombe_enc *e, const float *P, int flat, int *speech, float *v)
{
    float E[8], fl;
    int16_t *row;
    int b, i, redo = 0;
    band_powers(P, E);
    if (!e->vad_ready) {
        memcpy(e->vad_nl, E, sizeof(E));
        memcpy(e->vad_smooth, E, sizeof(E));
        e->vad_ready = 1;
    } else {
        for (b = 0; b < 8; b++)
            e->vad_smooth[b] = 0.8f * e->vad_smooth[b] + 0.2f * E[b];
    }
    fl = noise_floor(e->vad_nl);
    *v = vad_metric(E, e->vad_nl, fl);
    *speech = *v >= XOMBE_VAD_THRESHOLD;
    for (b = 0; b < 8; b++) {
        if ((E[b] > fl ? E[b] : fl) > (e->vad_nl[b] > fl ? e->vad_nl[b] : fl)) {
            if (e->vad_above[b] < XOMBE_VAD_HIST)
                e->vad_above[b]++;
        } else {
            e->vad_above[b] = 0;
        }
    }
    if (!*speech) {
        for (b = 0; b < 8; b++) {
            float nl;
            if (E[b] < e->vad_nl[b] / 4.0f)
                nl = 0.5f * e->vad_nl[b] + 0.5f * E[b];
            else
                nl = 0.9f * e->vad_nl[b] + 0.1f * E[b];
            if (!flat || nl > e->vad_nl[b])
                e->vad_nl[b] = nl;
        }
    }
    /* stationarity test floors stored smoothed powers at each half's own fl */
    row = e->vad_hist[e->vad_hist_pos];
    for (b = 0; b < 8; b++) {
        if (e->vad_hist_n == XOMBE_VAD_HIST)
            e->vad_hist_sum[b] -= hist_lin(row[b]);
        row[b] = hist_log(e->vad_smooth[b]);
        e->vad_hist_sum[b] += hist_lin(row[b]);
    }
    e->vad_hist_fl[e->vad_hist_pos] = hist_log(fl);
    e->vad_hist_pos = (e->vad_hist_pos + 1) % XOMBE_VAD_HIST;
    if (e->vad_hist_n < XOMBE_VAD_HIST)
        e->vad_hist_n++;
    if (e->vad_hist_n == XOMBE_VAD_HIST) {
        int steady = 1;
        for (b = 0; b < 8 && steady; b++) {
            int lo = INT16_MAX, hi = INT16_MIN;
            for (i = 0; i < XOMBE_VAD_HIST; i++) {
                int x = e->vad_hist[i][b] > e->vad_hist_fl[i] ? e->vad_hist[i][b] : e->vad_hist_fl[i];
                if (x < lo)
                    lo = x;
                if (x > hi)
                    hi = x;
            }
            steady = hi - lo <= ABSORB_SPREAD * HIST_LOG_ONE;
        }
        if (steady) {
            for (b = 0; b < 8; b++) {
                float mean = (float)e->vad_hist_sum[b] / (XOMBE_VAD_HIST * HIST_LIN_ONE);
                if (mean > e->vad_nl[b])
                    e->vad_nl[b] = mean;
            }
            redo = 1;
        }
        redo |= vad_release(e);
    }
    /* redecide after absorption or release, without a second noise update */
    if (redo) {
        fl = noise_floor(e->vad_nl);
        *v = vad_metric(E, e->vad_nl, fl);
        *speech = *v >= XOMBE_VAD_THRESHOLD;
    }
}
