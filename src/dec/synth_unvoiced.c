#include <string.h>
#include "dkernels.h"
#include "../tables/dsp_tables.h"

static int16_t noise_draw(xombe_dec *d)
{
    if (d->noise_mode == XOMBE_NOISE_PRODUCTION) {
        uint32_t x = d->prng;
        x ^= x << 13;
        x ^= x >> 17;
        x ^= x << 5;
        d->prng = x;
        return (int16_t)((int)(x >> 16) - 32768);
    }
    return xombe__lcg_draw(&d->lcg);
}

void xombe__noise_init(xombe_dec *d)
{
    int i;
    d->lcg = 0;
    for (i = 0; i < XOMBE_TAIL; i++)
        d->noise_hist[i] = noise_draw(d);
    memset(d->uv_overlap, 0, sizeof(d->uv_overlap));
}

static void noise_reset(xombe_dec *d)
{
    memset(d->noise_hist, 0, sizeof(d->noise_hist));
}

static void noise_segment(xombe_dec *d, int N, int16_t *seg)
{
    int16_t buf[XOMBE_TAIL + XOMBE_HALF_MAX];
    int i;
    memcpy(buf, d->noise_hist, sizeof(d->noise_hist));
    for (i = 0; i < N; i++)
        buf[XOMBE_TAIL + i] = noise_draw(d);
    memcpy(seg, buf, sizeof(buf[0]) * 2 * N);
    memcpy(d->noise_hist, buf + N, sizeof(d->noise_hist));
}

static float fade(int n, float rs)
{
    return xombe__sqrt_n[n] * rs;
}

void xombe__synth_unvoiced(xombe_dec *d, const xombe_set *cur, const float *M, const uint8_t *vh, int N,
                           float *out)
{
    int16_t seg[2 * XOMBE_HALF_MAX];
    float y[XOMBE_UV_FFT], rs = xombe__rsqrt_half[N - XOMBE_HALF_MIN];
    int n;
    if (cur->cls == XOMBE_CLASS_TONE) {
        noise_reset(d);
        for (n = 0; n < N; n++)
            out[n] += d->uv_overlap[n] * fade(N - n, rs);
        memset(d->uv_overlap, 0, sizeof(d->uv_overlap));
        return;
    }
    noise_segment(d, N, seg);
    xombe__dk_unvoiced(seg, N, cur, M, vh, d->c.avg_f0, y);
    for (n = 0; n < N; n++)
        out[n] += d->uv_overlap[n] * fade(N - n, rs) + y[n] * fade(n, rs);
    memset(d->uv_overlap, 0, sizeof(d->uv_overlap));
    memcpy(d->uv_overlap, y + N, sizeof(float) * N);
}
