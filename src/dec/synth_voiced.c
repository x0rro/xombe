#include <string.h>
#include "dkernels.h"

#define WT_SIZE XOMBE_WT_SIZE
/* epochs per half-frame: tracks stay below 0.0625 cycles/sample (500 Hz), so at most 6 in 84 samples */
#define EPOCH_MAX 8
/* octave jumps below OCTAVE_MAX, the epoch track follows the octave instead of jumping
   doubled: previous f0 between OCTAVE_UP_LO and OCTAVE_UP_HI of the current
   halved: current f0 between OCTAVE_DOWN_LO and OCTAVE_DOWN_HI of the previous */
#define OCTAVE_UP_LO 0.4f
#define OCTAVE_UP_HI 0.6f
#define OCTAVE_DOWN_LO 0.41667f
#define OCTAVE_DOWN_HI 0.625f
#define OCTAVE_MAX XOMBE_HZ(250.0f)

/* epochs of a pitch track over a half-frame of N samples: phase (in periods) starts at th0, frequency ramps
   linearly from fs to fs + delta, so the phase at n is th0 + n*(fs + n*delta/(2N)), the running sum of the
   midpoint frequencies fs + (n + 1/2)*delta/N, evaluated in closed form so float rounding does not accumulate
   each integer the phase reaches is an epoch. within its sample the phase is th + u*inc + h*u*(u - 1), inc the
   sample's increment, h = delta/(2N), the epoch placed by linear interpolation u = (k - th)/inc, then one
   Newton step on the curvature term. returns the epoch count (floor of the phase at N), that phase in *th_end */
static int epoch_track(float th0, float fs, float delta, int N, float t[EPOCH_MAX], float *th_end)
{
    float h = 0.5f * delta / (float)N, th = th0;
    int n, k = 1, count = 0;
    for (n = 0; n < N; n++) {
        float m = (float)(n + 1), next = th0 + m * (fs + h * m), inc = next - th;
        while ((float)k <= next && count < EPOCH_MAX) {
            float u = ((float)k - th) / inc;
            u -= h * u * (u - 1.0f) / (inc + h * (2.0f * u - 1.0f));
            t[count++] = (float)n + u;
            k++;
        }
        th = next;
    }
    *th_end = th;
    return count;
}

void xombe__synth_voiced(xombe_dec *d, const xombe_set *cur, const float *M, const uint8_t *vh, int cur_voiced,
                         int N, float *out)
{
    float wcur[WT_SIZE], old[XOMBE_HALF_MAX + XOMBE_TAIL], fresh[XOMBE_HALF_MAX + XOMBE_TAIL];
    int pv = d->prev_voiced, cv = cur_voiced, any = pv || cv;
    float fp = d->prev_f0, fc = cur->f0, fs = 0.0f, fe = 0.0f;
    float pc = d->pitch_phase, ec_start = d->epoch, ec = d->epoch;
    float delta, phi = pc, phi_new = 0.0f, inv_N = 1.0f / (float)N, t_old[EPOCH_MAX], t_new[EPOCH_MAX];
    float per_p = 1.0f / fp, per_c = 1.0f / fc;
    int up = 0, down = 0, K = 0, K_new = 0, k, n;

    if (pv && cv) {
        fs = fp;
        if (OCTAVE_UP_LO * fc < fp && fp < OCTAVE_UP_HI * fc && fp < OCTAVE_MAX) {
            up = 1;
            fe = fc / 2.0f;
        } else if (OCTAVE_DOWN_LO * fp < fc && fc < OCTAVE_DOWN_HI * fp && fc < OCTAVE_MAX) {
            down = 1;
            fe = 2.0f * fc;
        } else {
            fe = fc;
        }
    } else if (pv) {
        fs = fe = fp;
    } else if (cv) {
        fs = fe = fc;
    }
    delta = fe - fs;
    if (any)
        K = epoch_track(pc, fs, delta, N, t_old, &phi);

    memset(old, 0, sizeof(old));
    memcpy(old, d->tail, sizeof(d->tail));
    if (pv)
        for (k = 0; k < K; k++)
            xombe__dk_add_grain(old, N, t_old[k], fp, per_p, d->wprev);
    for (n = 0; n < N; n++)
        out[n] += (1.0f - (float)n * inv_N) * old[n];

    memset(fresh, 0, sizeof(fresh));
    if (cv) {
        float pp = pc, dp = delta, fsp = fs;
        xombe__dk_wavetable(cur, M, vh, wcur);
        /* across an octave jump the new track runs at twice (up) or half (down) the old one's rate */
        if (up) {
            pp = 2.0f * pc;
            if (pp >= 1.0f) {
                pp -= 1.0f;
                ec = -(pp / 2.0f) * per_p;
            }
            dp = 2.0f * delta;
            fsp = 2.0f * fs;
        } else if (down) {
            pp = pc / 2.0f;
            dp = delta / 2.0f;
            fsp = fs / 2.0f;
        }
        K_new = epoch_track(pp, fsp, dp, N, t_new, &phi_new);
        xombe__dk_add_grain(fresh, N, ec, fc, per_c, wcur);
        for (k = 0; k < K_new; k++)
            xombe__dk_add_grain(fresh, N, t_new[k], fc, per_c, wcur);
    }
    for (n = 0; n < N; n++)
        out[n] += ((float)n * inv_N) * fresh[n];
    memcpy(d->tail, fresh + N, sizeof(d->tail));

    if (!any) {
        d->pitch_phase = 0.0f;
        d->epoch = 0.0f;
    } else {
        d->pitch_phase = cv ? phi_new - K_new : phi - K;
        if (cv && K_new > 0)
            d->epoch = t_new[K_new - 1] - N;
        else if (pv && K > 0)
            d->epoch = t_old[K - 1] - N;
        else
            d->epoch = ec_start - N;
    }
    if (cv)
        memcpy(d->wprev, wcur, sizeof(wcur));
    d->prev_f0 = cur->f0;
    d->prev_voiced = cv;
}
