#include <math.h>
#include <string.h>
#include "../internal.h"

/* 100 Hz suppression bins, NS_BINS_PER_UNIT per cycle/sample */
#define NS_BINS_PER_UNIT 80.0f
/* voiced harmonics take half the suppression below NS_VOICED_FULL, a share rising with frequency up to
   NS_VOICED_RAMP */
#define NS_VOICED_FULL XOMBE_HZ(400.0f)
#define NS_VOICED_RAMP XOMBE_HZ(800.0f)
#define NS_GAIN_SLOPE 1.35f
#define NS_FAST 0.25f
#define NS_SLOW 0.1f
#define NS_FLOOR_UP 1.0444f
#define NS_FLOOR_DOWN 0.9591f
#define NS_MIN_LOG2 (-16.0f)

void xombe__ns_bins(const float *P, float Xm[XOMBE_NS_BINS])
{
    float acc[XOMBE_NS_BINS];
    int k;
    memset(acc, 0, sizeof(acc));
    for (k = 0; k < XOMBE_SPEC_BINS; k++) {
        int b = (int)((float)k / XOMBE_FFT * NS_BINS_PER_UNIT);
        if (b < XOMBE_NS_BINS)
            acc[b] += P[k];
    }
    for (k = 0; k < XOMBE_NS_BINS; k++)
        Xm[k] = sqrtf(acc[k]);
}

float xombe__ns_shape(float f)
{
    static const float fx[6] = {0.0f, XOMBE_HZ(500.0f), XOMBE_HZ(1400.0f), XOMBE_HZ(3000.0f), XOMBE_HZ(3900.0f),
                                XOMBE_HZ(4000.0f)};
    static const float fy[6] = {0.25f, 0.25f, 1.0f, 1.0f, 0.5f, 0.5f};
    int i;
    if (f <= fx[0])
        return fy[0];
    for (i = 1; i < 6; i++)
        if (f <= fx[i])
            return fy[i - 1] + (fy[i] - fy[i - 1]) * (f - fx[i - 1]) / (fx[i] - fx[i - 1]);
    return fy[5];
}

void xombe__ns_update(xombe_enc *e, const float Xm[XOMBE_NS_BINS], int speech, int confident_run,
                      int no_band_voiced)
{
    float s = e->ns_strength;
    int k;
    if (!e->ns_ready) {
        memcpy(e->ns_noise, Xm, sizeof(e->ns_noise));
        e->ns_floor = 0.5f * s;
        e->ns_ready = 1;
    } else if (confident_run >= XOMBE_CONFIDENT_RUN) {
        for (k = 0; k < XOMBE_NS_BINS; k++)
            e->ns_noise[k] = (1.0f - NS_FAST) * e->ns_noise[k] + NS_FAST * Xm[k];
    } else if (no_band_voiced) {
        for (k = 0; k < XOMBE_NS_BINS; k++)
            if (e->ns_noise[k] > Xm[k] / 3.0f)
                e->ns_noise[k] = (1.0f - NS_SLOW) * e->ns_noise[k] + NS_SLOW * Xm[k];
    }
    e->ns_floor *= speech ? NS_FLOOR_UP : NS_FLOOR_DOWN;
    if (e->ns_floor < 0.5f * s)
        e->ns_floor = 0.5f * s;
    if (e->ns_floor > s)
        e->ns_floor = s;
}

void xombe__ns_bin_voicing(const uint8_t vb[8], int ambe2_rule, uint8_t vbin[XOMBE_NS_BINS])
{
    int k;
    for (k = 0; k < XOMBE_NS_BINS; k++) {
        float fc = ((float)k + 0.5f) / NS_BINS_PER_UNIT;
        int j = (int)(fc * XOMBE_BANDS_PER_UNIT);
        vbin[k] = (uint8_t)(ambe2_rule ? xombe__voiced_at(vb, fc) : vb[j > 7 ? 7 : j]);
    }
}

void xombe__ns_gain(const xombe_enc *e, xombe_ns_mem *m, const float Xm[XOMBE_NS_BINS], int speech, float v,
                    const uint8_t vbin[XOMBE_NS_BINS], float G[XOMBE_NS_BINS])
{
    float g[XOMBE_NS_BINS], gf[XOMBE_NS_BINS], mu;
    int k;
    for (k = 0; k < XOMBE_NS_BINS; k++) {
        float fc = ((float)k + 0.5f) / NS_BINS_PER_UNIT;
        float xp = vbin[k] ? Xm[k] : Xm[k] * (1.0f + xombe__ns_shape(fc)) * 0.5f;
        float gk = 0.0f;
        if (xp > e->ns_noise[k]) {
            gk = NS_GAIN_SLOPE * (1.0f - e->ns_noise[k] / xp);
            if (gk > 1.0f)
                gk = 1.0f;
        }
        g[k] = gk;
    }
    for (k = 0; k < XOMBE_NS_BINS; k++) {
        float lo = g[k > 0 ? k - 1 : 0], hi = g[k < XOMBE_NS_BINS - 1 ? k + 1 : XOMBE_NS_BINS - 1];
        gf[k] = 0.25f * lo + 0.5f * g[k] + 0.25f * hi;
    }
    mu = speech ? 0.0f : (XOMBE_VAD_THRESHOLD - v) / XOMBE_VAD_THRESHOLD;
    if (mu < 0.0f)
        mu = 0.0f;
    if (mu > 1.0f)
        mu = 1.0f;
    for (k = 0; k < XOMBE_NS_BINS; k++) {
        float fc = ((float)k + 0.5f) / NS_BINS_PER_UNIT;
        float gs = m->have_prev ? 0.5f * (gf[k] + m->prev_gain[k]) : gf[k];
        float fk = vbin[k] ? e->ns_floor : e->ns_floor * xombe__ns_shape(fc);
        float mx = gs > fk ? gs : fk;
        if (mx > 1.0f)
            mx = 1.0f;
        G[k] = (1.0f - mu) * mx + mu * fk;
    }
    memcpy(m->prev_gain, gf, sizeof(gf));
    m->have_prev = 1;
}

static float log2_floored(float x)
{
    float d = x > 0.0f ? log2f(x) : NS_MIN_LOG2;
    return d < NS_MIN_LOG2 ? NS_MIN_LOG2 : d;
}

void xombe__ns_harmonic_offsets(const float G[XOMBE_NS_BINS], int L, float f0, const uint8_t *vh, float *da)
{
    int l;
    for (l = 1; l <= L; l++) {
        float f = (float)l * f0, pos = f * NS_BINS_PER_UNIT - 0.5f, gl, d;
        int i;
        if (pos <= 0.0f) {
            gl = G[0];
        } else if (pos >= XOMBE_NS_BINS - 1) {
            gl = G[XOMBE_NS_BINS - 1];
        } else {
            i = (int)pos;
            gl = G[i] + (G[i + 1] - G[i]) * (pos - (float)i);
        }
        d = log2_floored(gl);
        if (vh[l]) {
            if (f < NS_VOICED_FULL)
                d *= 0.5f;
            else if (f < NS_VOICED_RAMP)
                d *= f / NS_VOICED_RAMP;
        }
        da[l] = d;
    }
}

void xombe__ns_apply_voice(const float G[XOMBE_NS_BINS], float *a, int L, float f0, const uint8_t v2[8])
{
    uint8_t vh[XOMBE_LMAX + 1];
    float da[XOMBE_LMAX + 1];
    int l;
    xombe__harmonic_voicing(v2, f0, L, vh);
    xombe__ns_harmonic_offsets(G, L, f0, vh, da);
    for (l = 1; l <= L; l++) {
        a[l] += da[l];
        if (a[l] < XOMBE_AMP_MIN)
            a[l] = XOMBE_AMP_MIN;
    }
}

void xombe__ns_apply_silence(const xombe_enc *e, float *a, int L, float f0)
{
    int l;
    for (l = 1; l <= L; l++) {
        float f = (float)l * f0;
        a[l] += log2_floored(e->ns_floor * xombe__ns_shape(f));
        if (a[l] < XOMBE_AMP_MIN)
            a[l] = XOMBE_AMP_MIN;
    }
}
