#include "../internal.h"

static void cfft(xombe_cpx *x, int m, const xombe_cpx *tw, int stride, int inverse)
{
    int i, j, len;
    for (i = 1, j = 0; i < m; i++) {
        int bit = m >> 1;
        for (; j & bit; bit >>= 1)
            j ^= bit;
        j |= bit;
        if (i < j) {
            xombe_cpx t = x[i];
            x[i] = x[j];
            x[j] = t;
        }
    }
    for (len = 2; len <= m; len <<= 1) {
        int half = len >> 1, step = stride * (m / len);
        for (i = 0; i < m; i += len) {
            for (j = 0; j < half; j++) {
                xombe_cpx w = tw[j * step];
                xombe_cpx a = x[i + j], b = x[i + j + half], t;
                if (inverse)
                    w.im = -w.im;
                t.re = b.re * w.re - b.im * w.im;
                t.im = b.re * w.im + b.im * w.re;
                x[i + j].re = a.re + t.re;
                x[i + j].im = a.im + t.im;
                x[i + j + half].re = a.re - t.re;
                x[i + j + half].im = a.im - t.im;
            }
        }
    }
}

void xombe__rfft(const float *in, xombe_cpx *out, int n, const xombe_cpx *tw, int tw_n_real)
{
    int k;
    for (k = 0; k < n / 2; k++) {
        out[k].re = in[2 * k];
        out[k].im = in[2 * k + 1];
    }
    xombe__rfft_packed(out, n, tw, tw_n_real);
}

void xombe__rfft_packed(xombe_cpx *out, int n, const xombe_cpx *tw, int tw_n_real)
{
    int m = n / 2, k, s = tw_n_real / n;
    cfft(out, m, tw, 2 * s, 0);
    out[m] = out[0];
    for (k = 0; k <= m / 2; k++) {
        xombe_cpx z1 = out[k], z2 = out[m - k], e, o, w, r1, r2;
        e.re = 0.5f * (z1.re + z2.re);
        e.im = 0.5f * (z1.im - z2.im);
        o.re = 0.5f * (z1.im + z2.im);
        o.im = -0.5f * (z1.re - z2.re);
        if (k < m / 2 || k == 0) {
            w = tw[k * s];
        } else {
            w.re = 0.0f;
            w.im = -1.0f;
        }
        r1.re = e.re + w.re * o.re - w.im * o.im;
        r1.im = e.im + w.re * o.im + w.im * o.re;
        r2.re = e.re - w.re * o.re + w.im * o.im;
        r2.im = -(e.im - (w.re * o.im + w.im * o.re));
        out[k] = r1;
        out[m - k] = r2;
    }
}

void xombe__irfft(const xombe_cpx *in, float *out, int n, const xombe_cpx *tw, int tw_n_real)
{
    xombe_cpx z[XOMBE_IRFFT_MAX / 2];
    int m = n / 2, k, s = tw_n_real / n;
    float scale = 1.0f / (float)m;
    for (k = 0; k < m; k++) {
        xombe_cpx a = in[k], b = in[m - k], fe, fo, d, w;
        fe.re = 0.5f * (a.re + b.re);
        fe.im = 0.5f * (a.im - b.im);
        d.re = 0.5f * (a.re - b.re);
        d.im = 0.5f * (a.im + b.im);
        w = tw[k * s];
        w.im = -w.im;
        fo.re = d.re * w.re - d.im * w.im;
        fo.im = d.re * w.im + d.im * w.re;
        z[k].re = fe.re - fo.im;
        z[k].im = fe.im + fo.re;
    }
    cfft(z, m, tw, 2 * s, 1);
    for (k = 0; k < m; k++) {
        out[2 * k] = z[k].re * scale;
        out[2 * k + 1] = z[k].im * scale;
    }
}
