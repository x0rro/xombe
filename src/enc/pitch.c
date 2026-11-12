#include <math.h>
#include <string.h>
#include "kernels.h"

/* period grid from 21 samples in half-sample steps, index i is period (GRID_BASE + i)/2 */
#define GRID_BASE 42
#define DIGITAL_SILENCE_E0 1e-2f
/* tracking range 0.8..1.2 of a period, as TRACK_*_NUM/TRACK_DEN */
#define TRACK_LO_NUM 4
#define TRACK_HI_NUM 6
#define TRACK_DEN 5
#define BACKWARD_ACCEPT 0.48f
#define HALVING_PERIOD 100.0f
#define HALVING_MARGIN 0.1f
#define INITIAL_PERIOD 100.0f
#define REFINE_TOP XOMBE_HZ(2000.0f) /* refinement fits harmonics below this */
/* refinement tries period +/-1 sample in 1/REFINE_STEPS steps, none below REFINE_MIN_PERIOD */
#define REFINE_STEPS 8
#define REFINE_MIN_PERIOD 20.0f
/* submultiple P/n replaces the forward estimate P when its cumulative error c has c <= LOOSE and c/c(P) <=
   LOOSE_RATIO, c <= TIGHT and c/c(P) <= TIGHT_RATIO, or c <= ANY */
#define SUBMULT_LOOSE 0.85f
#define SUBMULT_LOOSE_RATIO 1.7f
#define SUBMULT_TIGHT 0.4f
#define SUBMULT_TIGHT_RATIO 3.5f
#define SUBMULT_ANY 0.05f

void xombe__pitch_init(xombe_enc *e)
{
    e->p_prev = INITIAL_PERIOD;
    e->e_prev1 = e->e_prev2 = 0.0f;
}

static float grid_period(int i)
{
    return 0.5f * (float)(GRID_BASE + i);
}

static int grid_index(float P)
{
    return (int)(2.0f * P + 0.5f) - GRID_BASE;
}

/* b > 0 */
static int floor_div(int a, int b)
{
    return a >= 0 ? a / b : -((-a + b - 1) / b);
}

static float error_at(const xombe_enc *e, int slot, float P)
{
    return e->psilent[slot] ? 1.0f : xombe__kern_error_at(e, slot, P);
}

void xombe__pitch_errors(xombe_enc *e, const float *centre, int slot)
{
    int i;
    xombe__kern_pitch_acf(e, centre, slot);
    e->psilent[slot] = e->pe0[slot] < DIGITAL_SILENCE_E0;
    for (i = 0; i < XOMBE_PITCH_GRID; i++)
        e->perr[slot][i] = error_at(e, slot, grid_period(i));
}

/* in exact arithmetic, since 0.8*P and 1.2*P land on the grid for some P and those ends are included */
static void grid_range(int i, int *lo, int *hi)
{
    *lo = -floor_div(TRACK_DEN * GRID_BASE - TRACK_LO_NUM * (GRID_BASE + i), TRACK_DEN);
    *hi = floor_div(TRACK_HI_NUM * (GRID_BASE + i) - TRACK_DEN * GRID_BASE, TRACK_DEN);
    if (*lo < 0)
        *lo = 0;
    if (*hi > XOMBE_PITCH_GRID - 1)
        *hi = XOMBE_PITCH_GRID - 1;
}

static int argmin_range(const float *v, int lo, int hi)
{
    int i, best = lo;
    for (i = lo + 1; i <= hi; i++)
        if (v[i] < v[best])
            best = i;
    return best;
}

XOMBE_STATIC_ASSERT(XOMBE_PITCH_GRID <= 256, grid_index_fits_uchar);

/* argmin_range(v, lo[i], hi[i]) for every i, lo[i] when empty, first on ties. linear time with a monotonic deque,
   both ends being nondecreasing in i */
static void sliding_argmin(const float *v, const unsigned char *lo, const unsigned char *hi, unsigned char *out)
{
    unsigned char dq[XOMBE_PITCH_GRID];
    int i, head = 0, tail = 0, next = 0;
    for (i = 0; i < XOMBE_PITCH_GRID; i++) {
        if (lo[i] > hi[i]) {
            out[i] = lo[i];
            continue;
        }
        if (next < lo[i])
            next = lo[i];
        for (; next <= hi[i]; next++) {
            while (tail > head && v[dq[tail - 1]] > v[next])
                tail--;
            dq[tail++] = (unsigned char)next;
        }
        while (dq[head] < lo[i])
            head++;
        out[i] = dq[head];
    }
}

float xombe__pitch_track(xombe_enc *e, int cur, int next1, int next2)
{
    const float *E0 = e->perr[cur], *E1 = e->perr[next1], *E2 = e->perr[next2];
    float min2[XOMBE_PITCH_GRID], ce_f[XOMBE_PITCH_GRID], ce_b, period;
    unsigned char rlo[XOMBE_PITCH_GRID], rhi[XOMBE_PITCH_GRID], am[XOMBE_PITCH_GRID];
    int i, lo, hi, pb, pf, n;
    if (e->psilent[cur])
        return e->p_prev;
    grid_range(grid_index(e->p_prev), &lo, &hi);
    pb = lo <= hi ? argmin_range(E0, lo, hi) : 0;
    ce_b = E0[pb] + e->e_prev1 + e->e_prev2;
    for (i = 0; i < XOMBE_PITCH_GRID; i++) {
        grid_range(i, &lo, &hi);
        rlo[i] = (unsigned char)lo;
        rhi[i] = (unsigned char)hi;
    }
    sliding_argmin(E2, rlo, rhi, am);
    for (i = 0; i < XOMBE_PITCH_GRID; i++)
        min2[i] = E2[am[i]];
    sliding_argmin(E1, rlo, rhi, am);
    for (i = 0; i < XOMBE_PITCH_GRID; i++)
        ce_f[i] = E0[i] + E1[am[i]] + min2[am[i]];
    pf = argmin_range(ce_f, 0, XOMBE_PITCH_GRID - 1);
    for (n = (GRID_BASE + pf) / GRID_BASE; n >= 2; n--) {
        /* grid index nearest P/n, halves up: floor((GRID_BASE + pf)/n - GRID_BASE + 1/2) */
        int idx = floor_div(2 * (GRID_BASE + pf - n * GRID_BASE) + n, 2 * n);
        float c, ratio;
        if (idx < 0)
            idx = 0;
        c = ce_f[idx];
        ratio = c / (ce_f[pf] > 1e-9f ? ce_f[pf] : 1e-9f);
        if ((c <= SUBMULT_LOOSE && ratio <= SUBMULT_LOOSE_RATIO) || (c <= SUBMULT_TIGHT && ratio <= SUBMULT_TIGHT_RATIO) ||
            c <= SUBMULT_ANY) {
            pf = idx;
            break;
        }
    }
    i = (ce_b <= BACKWARD_ACCEPT || ce_b <= ce_f[pf]) ? pb : pf;
    period = grid_period(i);
    e->e_prev2 = e->e_prev1;
    e->e_prev1 = E0[i];
    e->p_prev = period;
    if (period > HALVING_PERIOD && error_at(e, cur, 0.5f * period) <= E0[i] + HALVING_MARGIN)
        period /= 2.0f;
    return period;
}

int xombe__pitch_window_len(float period)
{
    int n = 2 * (int)floorf((2.5f * period - 1.0f) / 2.0f + 0.5f) + 1;
    if (n < XOMBE_BASE_WINDOW)
        n = XOMBE_BASE_WINDOW;
    if (n > XOMBE_MAX_WINDOW)
        n = XOMBE_MAX_WINDOW;
    return n;
}

float xombe__pitch_refine(xombe_enc *e, const xombe_spec *S, float period)
{
    float best = 0.0f, band_r[8], band_e[8], best_p = period;
    int d, found = 0, kmax = (int)(REFINE_TOP * XOMBE_FFT);
    xombe_lf_spec sp;
    xombe__kern_lf_prepare(e, S, 1, kmax, &sp);
    for (d = -REFINE_STEPS; d <= REFINE_STEPS; d++) {
        float p = period + (float)d / REFINE_STEPS, r;
        if (p < REFINE_MIN_PERIOD)
            continue;
        r = xombe__kern_line_fit(e, &sp, p, (int)ceilf(0.5f * XOMBE_FFT / p), kmax, band_r, band_e);
        if (!found || r < best) {
            found = 1;
            best = r;
            best_p = p;
        }
    }
    return best_p;
}
