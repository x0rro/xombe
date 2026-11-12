#include <math.h>
#include "../internal.h"

/* sin(pi*r/m) for m > 0, r reduced exactly to |r| <= m/2 first so the argument carries a single rounding
   however large r is (the reduced |r| must be exact in float: m < 2^24) */
float xombe__sin_ratio(long r, long m)
{
    r %= 2 * m;
    if (r > m)
        r -= 2 * m;
    else if (r < -m)
        r += 2 * m;
    if (2 * r > m)
        r = m - r;
    else if (2 * r < -m)
        r = -m - r;
    return sinf(XOMBE_PI_F * ((float)r / (float)m));
}

/* cos(2 pi i/(n - 1)) as sin(pi (4i + n - 1)/(2(n - 1))) */
void xombe__hamming(float *w, int n)
{
    int i;
    for (i = 0; i < n; i++)
        w[i] = 0.54f - 0.46f * xombe__sin_ratio(4L * i + (n - 1), 2L * (n - 1));
}

/* over i = 0..n-1 the cosine sums to 1 and its square to (n + 1)/2 (n > 3) */
float xombe__hamming_sum2(int n)
{
    return 0.2916f * (float)n + 0.1058f * (float)(n + 1) - 0.4968f;
}

/* sin(n*x)/sin(x) at x = pi*p/m, |p| < m */
static float dirichlet(int n, long p, long m)
{
    return p == 0 ? (float)n : xombe__sin_ratio((long)n * p, m) / xombe__sin_ratio(p, m);
}

/* xombe__sin_ratio(t, g) for g a power of two, without divisions */
static float sin_ratio_pow2(long t, long g)
{
    t &= 2 * g - 1;
    if (t > g)
        t -= 2 * g;
    if (2 * t > g)
        t = g - t;
    else if (2 * t < -g)
        t = -g - t;
    return sinf(XOMBE_PI_F * ((float)t * (1.0f / (float)g)));
}

/* transform of the centred xombe__hamming(n) at q/g cycles/sample (0 <= q < g/2): Dirichlet kernels at
   x = pi*q/g and x +/- pi/(n - 1), over the common denominator g*(n - 1) */
float xombe__hamming_dtft(int n, long q, long g)
{
    long p = q * (n - 1), m = g * (n - 1);
    /* for g a power of two the centre kernel's ratios reduce exactly to multiples of 1/g: sin_ratio(n*p, m)
       reduces n*p to (n - 1) times n*q reduced mod 2g, and the float quotient r/m equals the exact product
       t*(1/g) */
    if (q > 0 && (g & (g - 1)) == 0)
        return 0.54f * (sin_ratio_pow2((long)n * q, g) / sin_ratio_pow2(q, g)) +
               0.23f * (dirichlet(n, p - g, m) + dirichlet(n, p + g, m));
    return 0.54f * dirichlet(n, p, m) + 0.23f * (dirichlet(n, p - g, m) + dirichlet(n, p + g, m));
}
