#include <math.h>
#include <string.h>
#include "cquant.h"
#include "../tables/tables.h"
#include "../tables/dsp_tables.h"
#if XOMBE_KERNELS_Q15
#include "../dsp/q15ops.h"
#if XOMBE_HAVE_DSTAR
#include "../dstar/dstar_internal.h"
#endif
#endif

#define CQ_INF 3.0e38f
#define CQ_AMIN (-15.0f)
#define CQ_SQRT2 1.41421356f
/* relative slack of the pair search's row bound against float rounding (float eps is 6e-8) */
#define CQ_BOUND_SLACK 1e-5f
#define CQ_BANDS 32   /* bound bands of the pair search's row order */
#define CQ_ROW_COST 4 /* cap units of a row, a pair is 1 (Cortex-M33: ~280 cycles a row, ~85 a pair) */

#if XOMBE_KERNELS_Q15
#define CQ_BIAS 0x40000000u /* 2^30: |v| <= 2^29 shifted to unsigned, so v >> k floors without a signed shift */
#endif

#ifdef XOMBE_CQ_STATS
long xombe__cq_stat[4]; /* pair searches, rows evaluated, pairs evaluated, columns scanned */
#define CQ_STAT(k) (xombe__cq_stat[k]++)
#else
#define CQ_STAT(k) ((void)0)
#endif

static int usable(const unsigned char *mask, int i)
{
    return !mask || mask[i];
}

static float clip(float v, float lo, float hi)
{
    return v < lo ? lo : v > hi ? hi : v;
}

void xombe__cq_amp_weights(const float *t, int L, float gamma, float floor_, float *w)
{
    float mx = t[1];
    int l;
    for (l = 2; l <= L; l++)
        if (t[l] > mx)
            mx = t[l];
    w[0] = 0.0f;
    for (l = 1; l <= L; l++) {
        float v = exp2f(gamma * (t[l] - mx));
        w[l] = v > floor_ ? v : floor_;
    }
}

#if XOMBE_KERNELS_Q15
/* int16 pair (lo, hi) of values in range */
static uint32_t pair16(int32_t lo, int32_t hi)
{
    return ((uint32_t)lo & 0xffffu) | (uint32_t)hi << 16;
}

/* 2^e, -126 <= e <= 127 */
static float pow2i(int e)
{
    uint32_t b = (uint32_t)(e + 127) << 23;
    float f;
    memcpy(&f, &b, sizeof(f));
    return f;
}

/* E with |v| < 2^E (2^(E-1) <= |v| for normal v) */
static int expo(float v)
{
    uint32_t b;
    memcpy(&b, &v, sizeof(b));
    return (int)(b >> 23 & 255) - 126;
}

/* (int)(v rounded to nearest), |v| < 2^22 so v +/- 1/2 is exact and the error is at most 1/2 */
static int32_t round_i(float v)
{
    return (int32_t)(v + (v < 0.0f ? -0.5f : 0.5f));
}

/* rounds c[0..n) to ci = c*2^s, s the largest at most min(smax, 13 - E(max |c|)) with
   sum (|c| 2^s + 1/2)*mx*2^t <= lim, mx the image columns' maxima (units of 2^-t). returns s, and in *err a bound
   on |sum c*m - 2^-(s+t) sum ci*mi| over the rows (each mi within 1/2 of m*2^t): sum |c| 2^-(t+1) + 2^-(s+1) sum mx
   s = -100: c too large */
static int quant_coef(const float *c, int n, const float *mx, int t, int smax, float lim, int32_t *ci, float *err)
{
    float cm = 0.0f, sc = 0.0f, sm = 0.0f, u, f;
    int k, s;
    for (k = 0; k < n; k++) {
        float a = fabsf(c[k]);
        cm = a > cm ? a : cm;
        sc += a;
        sm += mx[k];
    }
    s = 13 - expo(cm);
    s = s > smax ? smax : s;
    s = s > 100 ? 100 : s < -100 ? -100 : s;
    u = pow2i(t);
    for (;; s--) {
        float tot = 0.0f;
        f = pow2i(s);
        for (k = 0; k < n; k++)
            tot += (fabsf(c[k]) * f + 0.5f) * mx[k] * u;
        if (tot <= lim || s == -100)
            break;
    }
    for (k = 0; k < n; k++)
        ci[k] = round_i(c[k] * f);
    *err = (sc * pow2i(-t - 1) + pow2i(-s - 1) * sm) * 1.0001f;
    return s;
}

/* sorts v[0..n-1] ascending with its tags, ties in input order (the stable order): buckets by value, the
   bucket index being monotonic in v, then insertion within them */
static void sort_proj(float *v, short *tag, int n)
{
    float tv[128], lo, hi, sc;
    short tt[128];
    unsigned char bk[128];
    int cnt[129], i, j;
    if (n < 2)
        return;
    lo = hi = v[0];
    for (i = 1; i < n; i++) {
        lo = v[i] < lo ? v[i] : lo;
        hi = v[i] > hi ? v[i] : hi;
    }
    if (!(hi > lo))
        return;
    sc = (float)(n - 1) / (hi - lo);
    memset(cnt, 0, sizeof(int) * (n + 1));
    for (i = 0; i < n; i++) {
        int b = (int)((v[i] - lo) * sc);
        b = b < n - 1 ? b : n - 1;
        bk[i] = (unsigned char)b;
        cnt[b + 1]++;
    }
    for (i = 0; i < n; i++)
        cnt[i + 1] += cnt[i];
    for (i = 0; i < n; i++) {
        int p = cnt[bk[i]]++;
        tv[p] = v[i];
        tt[p] = tag[i];
    }
    for (i = 1; i < n; i++) {
        float x = tv[i];
        short g = tt[i];
        for (j = i; j > 0 && tv[j - 1] > x; j--) {
            tv[j] = tv[j - 1];
            tt[j] = tt[j - 1];
        }
        tv[j] = x;
        tt[j] = g;
    }
    memcpy(v, tv, sizeof(float) * n);
    memcpy(tag, tt, sizeof(short) * n);
}
#else
/* sorts v[0..n-1] ascending with its tags (bottom-up merge sort, stable) */
static void sort_proj(float *v, short *tag, int n)
{
    float tv[128];
    short tt[128];
    int w, i;
    for (w = 1; w < n; w *= 2) {
        for (i = 0; i < n; i += 2 * w) {
            int a = i, am = i + w < n ? i + w : n, b = am, bm = i + 2 * w < n ? i + 2 * w : n, o = i;
            while (a < am || b < bm)
                if (b >= bm || (a < am && v[a] <= v[b])) {
                    tv[o] = v[a];
                    tt[o++] = tag[a++];
                } else {
                    tv[o] = v[b];
                    tt[o++] = tag[b++];
                }
        }
        memcpy(v, tv, sizeof(float) * n);
        memcpy(tag, tt, sizeof(short) * n);
    }
}

#endif

/* pair-search bounds (prba_pairs): Lc, Mz = Lc^-1 Q43, the principal direction of the usable p(j), and the rows
   sorted along it. returns 0 if Q44 is not safely positive definite, the search is then exhaustive */
static int prepare_bounds(xombe_cq_ws *ws)
{
    const xombe_amp_tables *ts = ws->ts;
    float tr = 0.0f, C[4][4], d[4] = {0.5f, 0.5f, 0.5f, 0.5f};
    int a, b, k, j, it, cnt = 0;
#if XOMBE_KERNELS_Q15
    const xombe_cq_q15_tab *qt = ws->qt;
    if (!qt)
        return 0;
#endif
    for (a = 0; a < 4; a++)
        tr += ws->Q[3 + a][3 + a];
    for (a = 0; a < 4; a++)
        for (b = 0; b < 4; b++) {
            float s = b <= a ? ws->Q[3 + a][3 + b] : 0.0f;
            if (b > a) {
                ws->Lc[a][b] = 0.0f;
                continue;
            }
            for (k = 0; k < b; k++)
                s -= ws->Lc[a][k] * ws->Lc[b][k];
            if (a == b) {
                if (!(s > 1e-6f * tr))
                    return 0;
                ws->Lc[a][a] = sqrtf(s);
            } else
                ws->Lc[a][b] = s / ws->Lc[b][b];
        }
    for (b = 0; b < 3; b++)
        for (a = 0; a < 4; a++) {
            float s = ws->Q[3 + a][b];
            for (k = 0; k < a; k++)
                s -= ws->Lc[a][k] * ws->Mz[k][b];
            ws->Mz[a][b] = s / ws->Lc[a][a];
        }
    /* principal axis of p = Lc^T x4 over the usable rows, covariance Lc^T cov(x4) Lc, p itself not needed
       (dir*p = (Lc dir)*x4, |p|^2 = e4q) */
    {
#if XOMBE_KERNELS_Q15
        const float *sx = qt->sx;
        float cx[4][4], T[4][4], u[4], emax = 0.0f;
        for (j = 0; j < ts->n_prba58; j++)
            if (usable(ts->ok_prba58, j) && ws->e4q[j] > emax)
                emax = ws->e4q[j];
        cnt = qt->cnt;
#define sxx qt->sxx
#else
        float sx[4] = {0, 0, 0, 0}, sxx[4][4], cx[4][4], T[4][4], u[4], emax = 0.0f;
        for (a = 0; a < 4; a++)
            for (b = 0; b < 4; b++)
                sxx[a][b] = 0.0f;
        for (j = 0; j < ts->n_prba58; j++) {
            const float *x = ts->prba58[j];
            if (!usable(ts->ok_prba58, j))
                continue;
            cnt++;
            for (a = 0; a < 4; a++) {
                sx[a] += x[a];
                for (b = 0; b <= a; b++)
                    sxx[a][b] += x[a] * x[b];
            }
            if (ws->e4q[j] > emax)
                emax = ws->e4q[j];
        }
#endif
        if (!cnt)
            return 0;
        ws->pmax = 1.001f * sqrtf(emax);
        for (a = 0; a < 4; a++)
            for (b = 0; b <= a; b++)
                cx[a][b] = cx[b][a] = sxx[a][b] - sx[a] * sx[b] / (float)cnt;
#undef sxx
        /* T = cov(x4) Lc, C = Lc^T T */
        for (a = 0; a < 4; a++)
            for (b = 0; b < 4; b++) {
                float s = 0.0f;
                for (k = b; k < 4; k++)
                    s += cx[a][k] * ws->Lc[k][b];
                T[a][b] = s;
            }
        for (a = 0; a < 4; a++)
            for (b = 0; b < 4; b++) {
                float s = 0.0f;
                for (k = a; k < 4; k++)
                    s += ws->Lc[k][a] * T[k][b];
                C[a][b] = s;
            }
        /* power iteration, any unit direction gives a valid bound, the principal one a tight one */
        for (it = 0; it < 16; it++) {
            float nd[4], nn = 0.0f;
            for (a = 0; a < 4; a++) {
                nd[a] = C[a][0] * d[0] + C[a][1] * d[1] + C[a][2] * d[2] + C[a][3] * d[3];
                nn += nd[a] * nd[a];
            }
            if (!(nn > 0.0f))
                break;
            nn = 1.0f / sqrtf(nn);
            for (a = 0; a < 4; a++)
                d[a] = nd[a] * nn;
        }
        {
            float nn = sqrtf(d[0] * d[0] + d[1] * d[1] + d[2] * d[2] + d[3] * d[3]);
            for (a = 0; a < 4; a++)
                ws->dir[a] = d[a] / nn;
        }
        /* u = Lc dir, Lc lower triangular */
        for (a = 0; a < 4; a++) {
            float s = 0.0f;
            for (k = 0; k <= a; k++)
                s += ws->Lc[a][k] * ws->dir[k];
            u[a] = s;
        }
        ws->n4 = 0;
        for (j = 0; j < ts->n_prba58; j++)
            if (usable(ts->ok_prba58, j)) {
                const float *x = ts->prba58[j];
                ws->ps[ws->n4] = u[0] * x[0] + u[1] * x[1] + u[2] * x[2] + u[3] * x[3];
                ws->pord[ws->n4++] = (short)j;
            }
    }
    sort_proj(ws->ps, ws->pord, ws->n4);
    /* per prba24 row: dir*(Mz x3), so dir*z(i) = dir*zq - dz(i) */
    {
        float md[3];
        int i;
        for (b = 0; b < 3; b++)
            md[b] = ws->Mz[0][b] * ws->dir[0] + ws->Mz[1][b] * ws->dir[1] + ws->Mz[2][b] * ws->dir[2] +
                    ws->Mz[3][b] * ws->dir[3];
#if XOMBE_KERNELS_Q15
        int32_t mi[3];
        uint32_t p0, p1;
        float e;
        int sm = quant_coef(md, 3, qt->x24, qt->tx24, 100, 1073741824.0f, mi, &e);
        if (sm == -100)
            return 0;
        p0 = pair16(mi[0], mi[1]);
        p1 = pair16(mi[2], 0);
        for (i = 0; i < ts->n_prba24; i++)
            ws->dz[i] = q15_smlad(qt->r24[i][3], p0, q15_smuad(qt->r24[i][4], p1));
        /* |dz| <= 2^30, so its float conversion is within 2^6 units */
        ws->dz_s = pow2i(-(sm + qt->tx24));
        ws->dz_e = e + 64.0f * ws->dz_s;
#else
        for (i = 0; i < ts->n_prba24; i++) {
            const float *x = ts->prba24[i];
            ws->dz[i] = md[0] * x[0] + md[1] * x[1] + md[2] * x[2];
        }
#endif
    }
    /* per prba24 row: Schur quadratic x3^T(Q33 - Mz^T Mz)x3 (unusable rows: CQ_INF), and a rounding magnitude
       bound mw >= x3^T(|Q33| + 2|Mz^T Mz|)x3 for every row */
    {
        float Sc[3][3], sa = 0.0f, n3 = 0.0f;
        int i;
        for (a = 0; a < 3; a++)
            for (b = 0; b < 3; b++) {
                float mm = ws->Mz[0][a] * ws->Mz[0][b] + ws->Mz[1][a] * ws->Mz[1][b] + ws->Mz[2][a] * ws->Mz[2][b] +
                           ws->Mz[3][a] * ws->Mz[3][b];
                Sc[a][b] = ws->Q[a][b] - mm;
                sa += fabsf(ws->Q[a][b]) + 2.0f * fabsf(mm);
            }
#if XOMBE_KERNELS_Q15
        {
            const float c[6] = {Sc[0][0], Sc[1][1], Sc[2][2], Sc[0][1], Sc[0][2], Sc[1][2]};
            int32_t ci[6];
            uint32_t p0, p1, p2;
            int sc = quant_coef(c, 6, qt->m24, qt->tm24, 100, 536870912.0f, ci, &ws->sq_d);
            if (sc == -100)
                return 0;
            p0 = pair16(ci[0], ci[1]);
            p1 = pair16(ci[2], ci[3]);
            p2 = pair16(ci[4], ci[5]);
            for (i = 0; i < ts->n_prba24; i++) {
                const uint32_t *r = qt->r24[i];
                ws->sq[i] = q15_smlad(r[0], p0, q15_smlad(r[1], p1, q15_smuad(r[2], p2)));
            }
            ws->sq_x = sc + qt->tm24;
            n3 = qt->n3;
        }
#else
        for (i = 0; i < ts->n_prba24; i++) {
            const float *x = ts->prba24[i];
            float xx = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
            if (!usable(ts->ok_prba24, i)) {
                ws->sq[i] = CQ_INF;
                continue;
            }
            ws->sq[i] = Sc[0][0] * x[0] * x[0] + Sc[1][1] * x[1] * x[1] + Sc[2][2] * x[2] * x[2] +
                        2.0f * (Sc[0][1] * x[0] * x[1] + Sc[0][2] * x[0] * x[2] + Sc[1][2] * x[1] * x[2]);
            if (xx > n3)
                n3 = xx;
        }
#endif
        ws->n3max = sqrtf(n3);
        ws->mw = sa * n3; /* |x|^T S|x| <= (sum S) |x|^2 for S >= 0 */
    }
    /* column bound, roles swapped: Q33 = L3 L3^T, M3 = L3^-1 Q34,
       S'(j) = x4^T(Q44 - M3^T M3)x4 - 2 x4*(q4 - M3^T zq') - |zq'|^2, zq' = L3^-1 q3 */
    ws->col_ok = 0;
    {
        float tr3 = ws->Q[0][0] + ws->Q[1][1] + ws->Q[2][2], Sc[4][4], sa = 0.0f, n4 = 0.0f;
        for (a = 0; a < 3; a++)
            for (b = 0; b <= a; b++) {
                float s = ws->Q[a][b];
                for (k = 0; k < b; k++)
                    s -= ws->L3[a][k] * ws->L3[b][k];
                if (a == b) {
                    if (!(s > 1e-6f * tr3))
                        return 1;
                    ws->L3[a][a] = sqrtf(s);
                } else
                    ws->L3[a][b] = s / ws->L3[b][b];
            }
        for (b = 0; b < 4; b++)
            for (a = 0; a < 3; a++) {
                float s = ws->Q[a][3 + b];
                for (k = 0; k < a; k++)
                    s -= ws->L3[a][k] * ws->M3[k][b];
                ws->M3[a][b] = s / ws->L3[a][a];
            }
        for (a = 0; a < 4; a++)
            for (b = 0; b < 4; b++) {
                float mm = ws->M3[0][a] * ws->M3[0][b] + ws->M3[1][a] * ws->M3[1][b] + ws->M3[2][a] * ws->M3[2][b];
                Sc[a][b] = ws->Q[3 + a][3 + b] - mm;
                sa += fabsf(ws->Q[3 + a][3 + b]) + 2.0f * fabsf(mm);
            }
#if XOMBE_KERNELS_Q15
        {
            const float c[10] = {Sc[0][0], Sc[1][1], Sc[2][2], Sc[3][3], Sc[0][1],
                                 Sc[0][2], Sc[0][3], Sc[1][2], Sc[1][3], Sc[2][3]};
            int32_t ci[10];
            uint32_t p[5];
            int sc = quant_coef(c, 10, qt->m58, qt->tm58, 100, 536870912.0f, ci, &ws->sq4_d);
            if (sc == -100)
                return 1;
            for (k = 0; k < 5; k++)
                p[k] = pair16(ci[2 * k], ci[2 * k + 1]);
            for (j = 0; j < ts->n_prba58; j++) {
                const uint32_t *r = qt->r58[j];
                ws->sq4[j] = q15_smlad(r[0], p[0], q15_smlad(r[1], p[1], q15_smlad(r[2], p[2], q15_smlad(r[3], p[3],
                             q15_smuad(r[4], p[4])))));
            }
            ws->sq4_x = sc + qt->tm58;
            n4 = qt->n4;
        }
#else
        for (j = 0; j < ts->n_prba58; j++) {
            const float *x = ts->prba58[j];
            float s = 0.0f, xx = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3];
            for (a = 0; a < 4; a++)
                for (b = 0; b < 4; b++)
                    s += x[a] * Sc[a][b] * x[b];
            ws->sq4[j] = s;
            if (usable(ts->ok_prba58, j) && xx > n4)
                n4 = xx;
        }
#endif
        ws->n4max = sqrtf(n4);
        ws->mw4 = sa * n4;
        ws->col_ok = 1;
    }
    return 1;
}

#if XOMBE_KERNELS_Q15
static const xombe_cq_q15_tab *q15_tab(const xombe_amp_tables *ts)
{
    if (ts == &xombe__amp_ambe2)
        return &xombe__cq_q15_ambe2;
#if XOMBE_HAVE_DSTAR
    if (ts == &xombe__amp_dstar)
        return &xombe__cq_q15_dstar;
#endif
    return 0;
}

/* the map's block form: within block k, A[l] = c0[k] + b1[l] c1[k], so its weighted moments need only the
   block sums of w, w b1, and w b1^2 */
void xombe__cq_prepare(xombe_cq_ws *ws, const xombe_amp_tables *ts, int L, const float *w)
{
    float W0[4], W1[4], W2[4];
    int k, l, g, h, i, m, j, pos = 1, off = 0, same = ws->ts == ts && ws->L == L;
    ws->ts = ts;
    ws->qt = q15_tab(ts);
    ws->L = L;
    ws->work = 0;
    ws->sw = 0.0f;
    ws->w[0] = 0.0f;
    for (l = 1; l <= L; l++) {
        ws->w[l] = w ? w[l] : 1.0f;
        ws->sw += ws->w[l];
    }
    for (k = 0; k < 4 && !same; k++) {
        int J = ts->blk[L][k], M = J < 6 ? J : 6;
        const float *kern = xombe__dct_blk + xombe__dct_blk_off[J];
        ws->pos[k] = pos;
        ws->J[k] = J;
        ws->M[k] = M;
        /* block means C0, C1 of (R_2k, R_2k+1), R = S8 G with S8[i][m] = alpha_m cos(pi m (i + 1/2)/8) */
        for (g = 0; g < 7; g++) {
            float s0 = 2.0f * xombe__dct8[g + 1][2 * k], s1 = 2.0f * xombe__dct8[g + 1][2 * k + 1];
            ws->c0[k][g] = 0.5f * (s0 + s1);
            ws->c1[k][g] = (0.5f / CQ_SQRT2) * (s0 - s1);
        }
        for (i = 0; i < J; i++) {
            ws->b1[pos + i] = 2.0f * kern[J + i];
            for (m = 2; m < 6; m++)
                ws->H[pos + i][m - 2] = m < M ? 2.0f * kern[m * J + i] : 0.0f;
        }
        pos += J;
    }
    for (k = 0; k < 4; k++) {
        float s0 = 0.0f, s1 = 0.0f, s2 = 0.0f;
        for (l = ws->pos[k]; l < ws->pos[k] + ws->J[k]; l++) {
            float wb = ws->w[l] * ws->b1[l];
            s0 += ws->w[l];
            s1 += wb;
            s2 += wb * ws->b1[l];
        }
        W0[k] = s0;
        W1[k] = s1;
        W2[k] = s2;
    }
    for (g = 0; g < 7; g++) {
        float s = 0.0f;
        for (k = 0; k < 4; k++)
            s += ws->c0[k][g] * W0[k] + ws->c1[k][g] * W1[k];
        ws->mu[g] = s / ws->sw;
    }
    for (g = 0; g < 7; g++)
        for (h = g; h < 7; h++) {
            float s = 0.0f;
            for (k = 0; k < 4; k++) {
                float dg = ws->c0[k][g] - ws->mu[g], dh = ws->c0[k][h] - ws->mu[h];
                s += dg * dh * W0[k] + (dg * ws->c1[k][h] + ws->c1[k][g] * dh) * W1[k] +
                     ws->c1[k][g] * ws->c1[k][h] * W2[k];
            }
            ws->Q[g][h] = ws->Q[h][g] = s;
        }
    {
        const float q00 = ws->Q[3][3], q11 = ws->Q[4][4], q22 = ws->Q[5][5], q33 = ws->Q[6][6];
        const float q01 = 2.0f * ws->Q[3][4], q02 = 2.0f * ws->Q[3][5], q03 = 2.0f * ws->Q[3][6];
        const float q12 = 2.0f * ws->Q[4][5], q13 = 2.0f * ws->Q[4][6], q23 = 2.0f * ws->Q[5][6];
        for (j = 0; j < ts->n_prba58; j++) {
            const float *x = ts->prba58[j];
            ws->e4q[j] = x[0] * (q00 * x[0] + q01 * x[1] + q02 * x[2] + q03 * x[3]) +
                         x[1] * (q11 * x[1] + q12 * x[2] + q13 * x[3]) + x[2] * (q22 * x[2] + q23 * x[3]) +
                         q33 * x[3] * x[3];
        }
    }
    ws->chol_ok = prepare_bounds(ws);
    /* HOC rows' w-only terms, G = H^T W H and H^T w per block (H is zero past M - 2 columns) */
    for (k = 0; k < 4; k++) {
        float G[4][4], h1[4];
        int n, pk = ws->pos[k];
        ws->hoc_off[k] = off;
        off += ts->n_hoc[k];
        if (ws->M[k] <= 2)
            continue;
        for (g = 0; g < 4; g++) {
            float s = 0.0f;
            for (i = 0; i < ws->J[k]; i++)
                s += ws->w[pk + i] * ws->H[pk + i][g];
            h1[g] = s;
            for (h = g; h < 4; h++) {
                s = 0.0f;
                for (i = 0; i < ws->J[k]; i++)
                    s += ws->w[pk + i] * ws->H[pk + i][g] * ws->H[pk + i][h];
                G[g][h] = s;
            }
        }
        for (n = 0; n < ts->n_hoc[k]; n++) {
            const float *c = ts->hoc[k] + 4 * n;
            ws->hq[ws->hoc_off[k] + n] = G[0][0] * c[0] * c[0] + G[1][1] * c[1] * c[1] + G[2][2] * c[2] * c[2] +
                                         G[3][3] * c[3] * c[3] +
                                         2.0f * (G[0][1] * c[0] * c[1] + G[0][2] * c[0] * c[2] + G[0][3] * c[0] * c[3] +
                                                 G[1][2] * c[1] * c[2] + G[1][3] * c[1] * c[3] + G[2][3] * c[2] * c[3]);
            ws->hl[ws->hoc_off[k] + n] = h1[0] * c[0] + h1[1] * c[1] + h1[2] * c[2] + h1[3] * c[3];
        }
    }
}
#else
void xombe__cq_prepare(xombe_cq_ws *ws, const xombe_amp_tables *ts, int L, const float *w)
{
    float mu[7];
    int k, l, g, h, i, m, j, pos = 1, off = 0, same = ws->ts == ts && ws->L == L;
    ws->ts = ts;
    ws->L = L;
    ws->work = 0;
    ws->sw = 0.0f;
    ws->w[0] = 0.0f;
    for (l = 1; l <= L; l++) {
        ws->w[l] = w ? w[l] : 1.0f;
        ws->sw += ws->w[l];
    }
    for (k = 0; k < 4 && !same; k++) {
        int J = ts->blk[L][k], M = J < 6 ? J : 6;
        const float *kern = xombe__dct_blk + xombe__dct_blk_off[J];
        float c0[7], c1[7];
        ws->pos[k] = pos;
        ws->J[k] = J;
        ws->M[k] = M;
        for (g = 0; g < 7; g++) {
            float s0 = 2.0f * xombe__dct8[g + 1][2 * k], s1 = 2.0f * xombe__dct8[g + 1][2 * k + 1];
            c0[g] = 0.5f * (s0 + s1);
            c1[g] = (0.5f / CQ_SQRT2) * (s0 - s1);
        }
        for (i = 0; i < J; i++) {
            float b1 = 2.0f * kern[J + i];
            for (g = 0; g < 7; g++)
                ws->A[pos + i][g] = c0[g] + b1 * c1[g];
            for (m = 2; m < 6; m++)
                ws->H[pos + i][m - 2] = m < M ? 2.0f * kern[m * J + i] : 0.0f;
        }
        pos += J;
    }
    for (g = 0; g < 7; g++) {
        float s = 0.0f;
        for (l = 1; l <= L; l++)
            s += ws->w[l] * ws->A[l][g];
        mu[g] = s / ws->sw;
    }
    for (l = 1; l <= L; l++)
        for (g = 0; g < 7; g++)
            ws->Wm[l][g] = ws->w[l] * (ws->A[l][g] - mu[g]);
    for (g = 0; g < 7; g++)
        for (h = g; h < 7; h++) {
            float s = 0.0f;
            for (l = 1; l <= L; l++)
                s += (ws->A[l][g] - mu[g]) * ws->Wm[l][h];
            ws->Q[g][h] = ws->Q[h][g] = s;
        }
    for (j = 0; j < ts->n_prba58; j++) {
        const float *x = ts->prba58[j];
        float e = 0.0f;
        for (g = 0; g < 4; g++) {
            float s = 0.0f;
            for (h = 0; h < 4; h++)
                s += ws->Q[3 + g][3 + h] * x[h];
            e += x[g] * s;
        }
        ws->e4q[j] = e;
    }
    ws->chol_ok = prepare_bounds(ws);
    /* HOC candidates in the harmonic domain, block-local */
    for (k = 0; k < 4 && !same; k++) {
        int J = ws->J[k], M = ws->M[k], n;
        ws->hv_off[k] = off;
        if (M <= 2)
            continue;
        for (n = 0; n < ts->n_hoc[k]; n++) {
            const float *hc = ts->hoc[k] + 4 * n;
            for (i = 0; i < J; i++) {
                float s = 0.0f;
                for (m = 2; m < M; m++)
                    s += ws->H[ws->pos[k] + i][m - 2] * hc[m - 2];
                ws->hv[off + n * J + i] = s;
            }
        }
        off += ts->n_hoc[k] * J;
    }
}

#endif

/* (b3, b4) pair search
   E(i, j) = e3(i) + e4(j) + v(i)*x4(j), e3 = x3^T Q33 x3 - 2 x3*q3, e4 = x4^T Q44 x4 - 2 x4*q4, v = 2 Q43 x3,
   evaluated in that order. the exhaustive search over usable pairs is pruned by bounds that never exclude a
   pair it would keep: with Q44 = Lc Lc^T, p = Lc^T x4, and z(i) = Lc^-1(q4 - Q43 x3),
     E(i, j) = S(i) + |p(j) - z(i)|^2,  S(i) = e3(i) - |z(i)|^2,
   so S(i) bounds row i and (dir*p(j) - dir*z(i))^2 bounds pair (i, j) for a unit direction dir (the rows'
   principal axis, sorted along it). each bound is lowered by a rounding slack, which also covers dir*z(i)
   taken as dir*zq - dir*(Mz x3) */

typedef struct {
    float e;
    int i, j;
} cq_pair;

static void pair_insert(cq_pair *top, int n, int n58, float e, int i, int j)
{
    int key = i * n58 + j, p;
    if (!(e < top[n - 1].e || (e == top[n - 1].e && key < top[n - 1].i * n58 + top[n - 1].j)))
        return;
    for (p = n - 1; p > 0; p--) {
        const cq_pair *q = &top[p - 1];
        if (!(e < q->e || (e == q->e && key < q->i * n58 + q->j)))
            break;
        top[p] = top[p - 1];
    }
    top[p].e = e;
    top[p].i = i;
    top[p].j = j;
}

static float row_e3(const xombe_cq_ws *ws, const float *x3, const float q3[3])
{
    float e = 0.0f, d = 0.0f;
    int a, b;
    for (a = 0; a < 3; a++) {
        float s = 0.0f;
        for (b = 0; b < 3; b++)
            s += ws->Q[a][b] * x3[b];
        e += x3[a] * s;
        d += x3[a] * q3[a];
    }
    return e - 2.0f * d;
}

static void row_v(const xombe_cq_ws *ws, const float *x3, float v[4])
{
    int b;
    for (b = 0; b < 4; b++)
        v[b] = 2.0f * (x3[0] * ws->Q[0][3 + b] + x3[1] * ws->Q[1][3 + b] + x3[2] * ws->Q[2][3 + b]);
}

static float pair_e(const xombe_cq_ws *ws, float e3, const float v[4], int j)
{
    const float *x = ws->ts->prba58[j];
    return (e3 + ws->e4[j]) + (v[0] * x[0] + v[1] * x[1] + v[2] * x[2] + v[3] * x[3]);
}

static void eval_row(const xombe_cq_ws *ws, int i, float e3, const float v[4], cq_pair *top, int n)
{
    const xombe_amp_tables *ts = ws->ts;
    int j, n58 = ts->n_prba58;
    for (j = 0; j < n58; j++)
        if (usable(ts->ok_prba58, j)) {
            CQ_STAT(2);
            pair_insert(top, n, n58, pair_e(ws, e3, v, j), i, j);
        }
}

/* row i with bound sb = S(i) - slack: columns outward from dir*z(i) = zdq - dz(i) along the sorted projections,
   each side until the projection bound exceeds the n-th best, e3 and v computed at the row's first pair
   returns the number of columns skipped by their own bound */
static int eval_row_proj(const xombe_cq_ws *ws, int i, const float q3[3], float sb, float zdq, const float *ps,
                         const short *pord, int n4, cq_pair *top, int n, int *left_)
{
    const xombe_amp_tables *ts = ws->ts;
    const float *x3 = ts->prba24[i];
#if XOMBE_KERNELS_Q15
    /* dz within dz_e: the true centre is in zd +/- ze, distances shrink by ze on either side */
    float e3 = 0.0f, v[4], zd = zdq - (float)ws->dz[i] * ws->dz_s, ze = ws->dz_e;
#else
    float e3 = 0.0f, v[4], zd = zdq - ws->dz[i];
#endif
    int lo = 0, hi = n4, t, n58 = ts->n_prba58, skip = 0, ready = 0, left = *left_ - CQ_ROW_COST;
    CQ_STAT(1);
    while (lo < hi) {
        int mid = (lo + hi) >> 1;
        if (ps[mid] < zd)
            lo = mid + 1;
        else
            hi = mid;
    }
    for (t = lo; t < n4; t++) {
#if XOMBE_KERNELS_Q15
        float d = ps[t] - zd - ze > 0.0f ? ps[t] - zd - ze : 0.0f;
#else
        float d = ps[t] - zd;
#endif
        int j = pord[t];
        float e;
        CQ_STAT(3);
        if (sb + d * d > top[n - 1].e)
            break;
        if (ws->colb[j] > top[n - 1].e) {
            skip++;
            continue;
        }
        if (!ready) {
            e3 = row_e3(ws, x3, q3);
            row_v(ws, x3, v);
            ready = 1;
        }
        if (--left < 0)
            break;
        CQ_STAT(2);
        e = pair_e(ws, e3, v, j);
        if (e <= top[n - 1].e)
            pair_insert(top, n, n58, e, i, j);
    }
    for (t = lo - 1; t >= 0 && left >= 0; t--) {
#if XOMBE_KERNELS_Q15
        float d = zd - ze - ps[t] > 0.0f ? zd - ze - ps[t] : 0.0f;
#else
        float d = zd - ps[t];
#endif
        int j = pord[t];
        float e;
        CQ_STAT(3);
        if (sb + d * d > top[n - 1].e)
            break;
        if (ws->colb[j] > top[n - 1].e) {
            skip++;
            continue;
        }
        if (!ready) {
            e3 = row_e3(ws, x3, q3);
            row_v(ws, x3, v);
            ready = 1;
        }
        if (--left < 0)
            break;
        CQ_STAT(2);
        e = pair_e(ws, e3, v, j);
        if (e <= top[n - 1].e)
            pair_insert(top, n, n58, e, i, j);
    }
    *left_ = left;
    return skip;
}

#if XOMBE_KERNELS_Q15
static float row_sb(const xombe_cq_ws *ws, int i)
{
    return (float)ws->sb[i] * ws->sb_s - ws->sb_k;
}

static void pairs_all(const xombe_cq_ws *ws, const float q[7], int n, cq_pair *top)
{
    const xombe_amp_tables *ts = ws->ts;
    int i;
    for (i = 0; i < ts->n_prba24; i++) {
        float v[4];
        if (!usable(ts->ok_prba24, i))
            continue;
        CQ_STAT(1);
        row_v(ws, ts->prba24[i], v);
        eval_row(ws, i, row_e3(ws, ts->prba24[i], q), v, top, n);
    }
}

/* n best pairs by weighted mean-removed error of A3 x3 + A4 x4 against y, ascending, ties by row-major index
   the bounds' linear terms in int16 at the largest scale that keeps every sum within 2^30, the row bound's
   quadratic part (scale 2^-sq_x) shifted down to it */
static void prba_pairs(xombe_cq_ws *ws, const float *y, int n, cq_pair *top)
{
    const xombe_amp_tables *ts = ws->ts;
    const xombe_cq_q15_tab *qt = ws->qt;
    float q[7], zq[4], my = 0.0f, e4abs = 0.0f, pm = ws->pmax, zdq, sb0;
    int l, g, i, j, m, k, i0 = -1, left = 1 << 30, budget = ws->cap > 0 ? ws->cap : 1 << 30;
    short cand[512], ord[512]; /* band << 9 | row, then rows by band */
    CQ_STAT(0);
    for (l = 1; l <= ws->L; l++)
        my += ws->w[l] * y[l];
    my /= ws->sw;
    for (g = 0; g < 7; g++)
        q[g] = 0.0f;
    for (k = 0; k < 4; k++) {
        float u0 = 0.0f, u1 = 0.0f;
        for (l = ws->pos[k]; l < ws->pos[k] + ws->J[k]; l++) {
            float wu = ws->w[l] * (y[l] - my);
            u0 += wu;
            u1 += wu * ws->b1[l];
        }
        for (g = 0; g < 7; g++)
            q[g] += (ws->c0[k][g] - ws->mu[g]) * u0 + ws->c1[k][g] * u1;
    }
    for (j = 0; j < ts->n_prba58; j++) {
        const float *x = ts->prba58[j];
        ws->e4[j] = ws->e4q[j] - 2.0f * (x[0] * q[3] + x[1] * q[4] + x[2] * q[5] + x[3] * q[6]);
        if (usable(ts->ok_prba58, j) && fabsf(ws->e4[j]) > e4abs)
            e4abs = fabsf(ws->e4[j]);
    }
    for (m = 0; m < n; m++) {
        top[m].e = CQ_INF;
        top[m].i = top[m].j = 1 << 20;
    }
    if (ws->brute || !ws->chol_ok) {
        pairs_all(ws, q, n, top);
        return;
    }
    /* zq = Lc^-1 q4 */
    for (m = 0; m < 4; m++) {
        float s = q[3 + m];
        for (g = 0; g < m; g++)
            s -= ws->Lc[m][g] * zq[g];
        zq[m] = s / ws->Lc[m][m];
    }
    /* S(i) = x3^T(Q33 - Mz^T Mz)x3 - 2 x3*h - |zq|^2, h = q3 - Mz^T zq, the slack covers the rounding of S and of
       E (|e3| <= x3^T|Q33|x3 + 2|x3||q3|, |v*x4| <= |e4| + |p|^2 + 2|z||p| <= |e4| + 2|p|^2 + |z|^2,
       |z|^2 <= 2|zq|^2 + 2|Mz x3|^2), the integer parts' errors are added to it */
    {
        float h[3], zzq = 0.0f, nq = 0.0f, nh = 0.0f, K, lh[3], el;
        int32_t hi[4];
        uint32_t h01, h2;
        int sl, x;
        for (m = 0; m < 4; m++)
            zzq += zq[m] * zq[m];
        for (g = 0; g < 3; g++) {
            h[g] = q[g] - (ws->Mz[0][g] * zq[0] + ws->Mz[1][g] * zq[1] + ws->Mz[2][g] * zq[2] + ws->Mz[3][g] * zq[3]);
            nq += q[g] * q[g];
            nh += h[g] * h[g];
            lh[g] = -2.0f * h[g];
        }
        K = 2.0f * ws->n3max * (sqrtf(nq) + sqrtf(nh)) + 3.0f * zzq + 3.0f * e4abs + 3.0f * pm * pm;
        /* colb[j] = S'(j) - slack, or -CQ_INF without the bound. |sq4| <= 2^29 and |lin| <= 2^30 convert to
           float within 2^5 and 2^6 units */
        if (ws->col_ok) {
            float zq3[3], h4[4], zz3 = 0.0f, n4q = 0.0f, nh4 = 0.0f, K4, l4[4], e4l, s4, sx4;
            uint32_t p0, p1;
            int s;
            for (m = 0; m < 3; m++) {
                float t = q[m];
                for (g = 0; g < m; g++)
                    t -= ws->L3[m][g] * zq3[g];
                zq3[m] = t / ws->L3[m][m];
                zz3 += zq3[m] * zq3[m];
            }
            for (g = 0; g < 4; g++) {
                h4[g] = q[3 + g] - (ws->M3[0][g] * zq3[0] + ws->M3[1][g] * zq3[1] + ws->M3[2][g] * zq3[2]);
                n4q += q[3 + g] * q[3 + g];
                nh4 += h4[g] * h4[g];
                l4[g] = -2.0f * h4[g];
            }
            s = quant_coef(l4, 4, qt->x58, qt->tx58, 100, 1073741824.0f, hi, &e4l);
            if (s == -100) {
                pairs_all(ws, q, n, top);
                return;
            }
            s4 = pow2i(-(s + qt->tx58));
            sx4 = pow2i(-ws->sq4_x);
            K4 = zz3 + CQ_BOUND_SLACK * (ws->mw4 + ws->mw + K + 2.0f * ws->n4max * (sqrtf(n4q) + sqrtf(nh4)) +
                                         3.0f * zz3) +
                 ws->sq4_d + e4l + 32.0f * sx4 + 64.0f * s4;
            p0 = pair16(hi[0], hi[1]);
            p1 = pair16(hi[2], hi[3]);
            for (j = 0; j < ts->n_prba58; j++) {
                const uint32_t *r = qt->r58[j];
                ws->colb[j] = (float)ws->sq4[j] * sx4 + (float)q15_smlad(r[5], p0, q15_smuad(r[6], p1)) * s4 - K4;
            }
        } else
            for (j = 0; j < ts->n_prba58; j++)
                ws->colb[j] = -CQ_INF;
        K = zzq + CQ_BOUND_SLACK * (ws->mw + K);
        /* sb = (sq >> k) + lin at 2^-x, x = sq_x - k: within sq_d + el, 1 unit for the floor of the shift, and
           2^7 units for the float conversion of |sb| < 2^31 */
        sl = quant_coef(lh, 3, qt->x24, qt->tx24, ws->sq_x - qt->tx24, 1073741824.0f, hi, &el);
        x = sl + qt->tx24;
        k = ws->sq_x - x;
        if (sl == -100 || k > 30) {
            pairs_all(ws, q, n, top);
            return;
        }
        ws->sb_s = pow2i(-x);
        ws->sb_k = K + ws->sq_d + el + (k ? ws->sb_s : 0.0f) + 128.0f * ws->sb_s;
        h01 = pair16(hi[0], hi[1]);
        h2 = pair16(hi[2], 0);
        {
            const uint32_t(*r24)[5] = qt->r24;
            const unsigned char *ok = ts->ok_prba24;
            int32_t mn = INT32_MAX, bk = (int32_t)(CQ_BIAS >> k), *sb = ws->sb;
            const int32_t *sq = ws->sq;
            if (k)
                for (i = 0; i < ts->n_prba24; i++)
                    sb[i] = q15_smlad(r24[i][3], h01,
                                      q15_smlad(r24[i][4], h2, (int32_t)(((uint32_t)sq[i] + CQ_BIAS) >> k) - bk));
            else
                for (i = 0; i < ts->n_prba24; i++)
                    sb[i] = q15_smlad(r24[i][3], h01, q15_smlad(r24[i][4], h2, sq[i]));
            for (i = 0; i < ts->n_prba24; i++) {
                if (ok && !ok[i])
                    sb[i] = INT32_MAX;
                else if (sb[i] < mn) {
                    mn = sb[i];
                    i0 = i;
                }
            }
        }
    }
    if (i0 < 0)
        return;
    sb0 = row_sb(ws, i0);
    zdq = ws->dir[0] * zq[0] + ws->dir[1] * zq[1] + ws->dir[2] * zq[2] + ws->dir[3] * zq[3];
    eval_row_proj(ws, i0, q, sb0, zdq, ws->ps, ws->pord, ws->n4, top, n, &left);
    left = budget - ((1 << 30) - left); /* the first row always runs to completion */
    /* columns with bound within the n-th best so far, in projection order (the n-th best only drops), refiltered
       when the rows have skipped more columns by their bound than the list holds. other rows within the n-th
       best go by ascending bound in CQ_BANDS bands (index order within a band), so a cap cuts the least
       promising. without a cap the order changes only how fast the n-th best drops, not the n best
       the integer pre-filter sb <= T keeps a superset of the rows with bound <= t0 (T is past (t0 + sb_k)/sb_s by
       more than its float rounding, and INT32_MAX - 1 at most so it drops the unusable rows) */
    {
        float t0 = top[n - 1].e, scale = t0 > sb0 ? (float)CQ_BANDS / (t0 - sb0) : 0.0f, X = (t0 + ws->sb_k) / ws->sb_s;
        float Xm = X + (fabsf(X) + (fabsf(t0) + fabsf(ws->sb_k)) / ws->sb_s) * 1e-6f + 512.0f;
        int32_t T = Xm >= 2147483000.0f ? INT32_MAX - 1 : Xm <= -2147483000.0f ? INT32_MIN : (int32_t)Xm;
        int nf = 0, nh = 0, t, skip = 0, cnt[CQ_BANDS + 1] = {0};
        for (t = 0; t < ws->n4; t++)
            if (ws->colb[ws->pord[t]] <= top[n - 1].e) {
                ws->fps[nf] = ws->ps[t];
                ws->ford[nf++] = ws->pord[t];
            }
        for (i = 0; i < ts->n_prba24; i++)
            if (ws->sb[i] <= T && i != i0) {
                float sb = row_sb(ws, i);
                int b;
                if (!(sb <= t0))
                    continue;
                b = (int)((sb - sb0) * scale);
                b = b < CQ_BANDS ? b : CQ_BANDS - 1;
                cand[nh++] = (short)(b << 9 | i);
                cnt[b + 1]++;
            }
        for (k = 0; k < CQ_BANDS; k++)
            cnt[k + 1] += cnt[k];
        for (t = 0; t < nh; t++)
            ord[cnt[cand[t] >> 9]++] = (short)(cand[t] & 511);
        for (t = 0; t < nh && left > 0; t++) {
            float sb;
            i = ord[t];
            sb = row_sb(ws, i);
            if (sb > top[n - 1].e)
                continue;
            skip += eval_row_proj(ws, i, q, sb, zdq, ws->fps, ws->ford, nf, top, n, &left);
            if (skip > nf) {
                int f = 0;
                for (k = 0; k < nf; k++)
                    if (ws->colb[ws->ford[k]] <= top[n - 1].e) {
                        ws->fps[f] = ws->fps[k];
                        ws->ford[f++] = ws->ford[k];
                    }
                nf = f;
                skip = 0;
            }
        }
    }
    if (budget - left > ws->work)
        ws->work = budget - left;
}
#else
static void prba_pairs(xombe_cq_ws *ws, const float *y, int n, cq_pair *top)
{
    const xombe_amp_tables *ts = ws->ts;
    float u[XOMBE_LMAX + 1], q[7], zq[4], my = 0.0f, e4abs = 0.0f, pm = ws->pmax, zdq, sb0 = CQ_INF;
    int l, g, i, j, m, L = ws->L, i0 = -1, left = 1 << 30, budget = ws->cap > 0 ? ws->cap : 1 << 30;
    short cand[512], ord[512];
    CQ_STAT(0);
    for (l = 1; l <= L; l++)
        my += ws->w[l] * y[l];
    my /= ws->sw;
    for (l = 1; l <= L; l++)
        u[l] = y[l] - my;
    for (g = 0; g < 7; g++) {
        float s = 0.0f;
        for (l = 1; l <= L; l++)
            s += ws->Wm[l][g] * u[l];
        q[g] = s;
    }
    for (j = 0; j < ts->n_prba58; j++) {
        const float *x = ts->prba58[j];
        ws->e4[j] = ws->e4q[j] - 2.0f * (x[0] * q[3] + x[1] * q[4] + x[2] * q[5] + x[3] * q[6]);
        if (usable(ts->ok_prba58, j) && fabsf(ws->e4[j]) > e4abs)
            e4abs = fabsf(ws->e4[j]);
    }
    for (m = 0; m < n; m++) {
        top[m].e = CQ_INF;
        top[m].i = top[m].j = 1 << 20;
    }
    if (ws->brute || !ws->chol_ok) {
        for (i = 0; i < ts->n_prba24; i++) {
            float v[4];
            if (!usable(ts->ok_prba24, i))
                continue;
            CQ_STAT(1);
            row_v(ws, ts->prba24[i], v);
            eval_row(ws, i, row_e3(ws, ts->prba24[i], q), v, top, n);
        }
        return;
    }
    for (m = 0; m < 4; m++) {
        float s = q[3 + m];
        for (g = 0; g < m; g++)
            s -= ws->Lc[m][g] * zq[g];
        zq[m] = s / ws->Lc[m][m];
    }
    {
        float h[3], zzq = 0.0f, nq = 0.0f, nh = 0.0f, K;
        for (m = 0; m < 4; m++)
            zzq += zq[m] * zq[m];
        for (g = 0; g < 3; g++) {
            h[g] = q[g] - (ws->Mz[0][g] * zq[0] + ws->Mz[1][g] * zq[1] + ws->Mz[2][g] * zq[2] + ws->Mz[3][g] * zq[3]);
            nq += q[g] * q[g];
            nh += h[g] * h[g];
        }
        K = 2.0f * ws->n3max * (sqrtf(nq) + sqrtf(nh)) + 3.0f * zzq + 3.0f * e4abs + 3.0f * pm * pm;
        if (ws->col_ok) {
            float zq3[3], h4[4], zz3 = 0.0f, n4q = 0.0f, nh4 = 0.0f, K4;
            for (m = 0; m < 3; m++) {
                float t = q[m];
                for (g = 0; g < m; g++)
                    t -= ws->L3[m][g] * zq3[g];
                zq3[m] = t / ws->L3[m][m];
                zz3 += zq3[m] * zq3[m];
            }
            for (g = 0; g < 4; g++) {
                h4[g] = q[3 + g] - (ws->M3[0][g] * zq3[0] + ws->M3[1][g] * zq3[1] + ws->M3[2][g] * zq3[2]);
                n4q += q[3 + g] * q[3 + g];
                nh4 += h4[g] * h4[g];
            }
            K4 = zz3 + CQ_BOUND_SLACK * (ws->mw4 + ws->mw + K + 2.0f * ws->n4max * (sqrtf(n4q) + sqrtf(nh4)) +
                                         3.0f * zz3);
            for (j = 0; j < ts->n_prba58; j++) {
                const float *x = ts->prba58[j];
                ws->colb[j] = ws->sq4[j] - 2.0f * (x[0] * h4[0] + x[1] * h4[1] + x[2] * h4[2] + x[3] * h4[3]) - K4;
            }
        } else
            for (j = 0; j < ts->n_prba58; j++)
                ws->colb[j] = -CQ_INF;
        K = zzq + CQ_BOUND_SLACK * (ws->mw + K);
        for (i = 0; i < ts->n_prba24; i++) {
            const float *x3 = ts->prba24[i];
            float sb = ws->sq[i] - 2.0f * (x3[0] * h[0] + x3[1] * h[1] + x3[2] * h[2]) - K;
            ws->sb[i] = sb;
            if (sb < sb0) {
                sb0 = sb;
                i0 = i;
            }
        }
    }
    if (i0 < 0)
        return;
    zdq = ws->dir[0] * zq[0] + ws->dir[1] * zq[1] + ws->dir[2] * zq[2] + ws->dir[3] * zq[3];
    eval_row_proj(ws, i0, q, sb0, zdq, ws->ps, ws->pord, ws->n4, top, n, &left);
    left = budget - ((1 << 30) - left);
    {
        float t0 = top[n - 1].e, scale = t0 > sb0 ? (float)CQ_BANDS / (t0 - sb0) : 0.0f;
        int nf = 0, nh = 0, t, k, skip = 0, cnt[CQ_BANDS + 1] = {0};
        for (t = 0; t < ws->n4; t++)
            if (ws->colb[ws->pord[t]] <= top[n - 1].e) {
                ws->fps[nf] = ws->ps[t];
                ws->ford[nf++] = ws->pord[t];
            }
        for (i = 0; i < ts->n_prba24; i++)
            if (i != i0 && ws->sb[i] <= t0) {
                int b = (int)((ws->sb[i] - sb0) * scale);
                b = b < CQ_BANDS ? b : CQ_BANDS - 1;
                cand[nh++] = (short)(b << 9 | i);
                cnt[b + 1]++;
            }
        for (k = 0; k < CQ_BANDS; k++)
            cnt[k + 1] += cnt[k];
        for (t = 0; t < nh; t++)
            ord[cnt[cand[t] >> 9]++] = (short)(cand[t] & 511);
        for (t = 0; t < nh && left > 0; t++) {
            i = ord[t];
            if (ws->sb[i] > top[n - 1].e)
                continue;
            skip += eval_row_proj(ws, i, q, ws->sb[i], zdq, ws->fps, ws->ford, nf, top, n, &left);
            if (skip > nf) {
                int f = 0;
                for (k = 0; k < nf; k++)
                    if (ws->colb[ws->ford[k]] <= top[n - 1].e) {
                        ws->fps[f] = ws->fps[k];
                        ws->ford[f++] = ws->ford[k];
                    }
                nf = f;
                skip = 0;
            }
        }
    }
    if (budget - left > ws->work)
        ws->work = budget - left;
}

#endif

#if XOMBE_KERNELS_Q15
/* per block, the HOC row minimizing the weighted mean-removed error against y with the other blocks fixed,
   up to passes sweeps until nothing changes. a block whose others are unchanged since its last search is
   skipped. with y centred, row c of block k costs
     e = O2 + Y2 - 2 c*b + c^T G c - (O1 + Y1 - c*H^T w)^2 / sw
   Y1, Y2, b = sum w y, sum w y^2, H^T W y over block k, O1, O2 the same sums of the other blocks' residuals */
static void hoc_search(const xombe_cq_ws *ws, const float *y, int hoc[4], int passes)
{
    const xombe_amp_tables *ts = ws->ts;
    float Y1[4], Y2[4], B[4][4], P1[4], P2[4], cy = 0.0f, isw = 1.0f / ws->sw;
    int p, k, kk, i, n, l, step = 0, last[4] = {-1, -1, -1, -1}, chg[4] = {-1, -1, -1, -1};
    for (l = 1; l <= ws->L; l++)
        cy += ws->w[l] * y[l];
    cy *= isw;
    for (k = 0; k < 4; k++) {
        float s1 = 0.0f, s2 = 0.0f, b[4] = {0.0f, 0.0f, 0.0f, 0.0f};
        int pos = ws->pos[k], m;
        for (i = 0; i < ws->J[k]; i++) {
            float d = y[pos + i] - cy, wd = ws->w[pos + i] * d;
            s1 += wd;
            s2 += wd * d;
            for (m = 0; m < 4; m++)
                b[m] += wd * ws->H[pos + i][m];
        }
        Y1[k] = P1[k] = s1;
        Y2[k] = P2[k] = s2;
        memcpy(B[k], b, sizeof(b));
        if (ws->M[k] > 2) {
            const float *c = ts->hoc[k] + 4 * hoc[k];
            P1[k] = s1 - ws->hl[ws->hoc_off[k] + hoc[k]];
            P2[k] = s2 - 2.0f * (c[0] * b[0] + c[1] * b[1] + c[2] * b[2] + c[3] * b[3]) +
                    ws->hq[ws->hoc_off[k] + hoc[k]];
        }
    }
    for (p = 0; p < passes; p++) {
        int changed = 0;
        for (k = 0; k < 4; k++) {
            const float *b = B[k], *hq = ws->hq + ws->hoc_off[k], *hl = ws->hl + ws->hoc_off[k];
            float a1 = Y1[k], best = CQ_INF;
            int bi = 0;
            if (ws->M[k] <= 2)
                continue;
            if (last[k] >= 0) {
                for (kk = 0; kk < 4 && (kk == k || chg[kk] < last[k]); kk++)
                    ;
                if (kk == 4)
                    continue;
            }
            last[k] = step++;
            for (kk = 0; kk < 4; kk++)
                if (kk != k)
                    a1 += P1[kk];
            /* O2 + Y2 is the same for every row */
            for (n = 0; n < ts->n_hoc[k]; n++) {
                const float *c = ts->hoc[k] + 4 * n;
                float r1, e;
                if (!usable(ts->ok_hoc[k], n))
                    continue;
                r1 = a1 - hl[n];
                e = hq[n] - 2.0f * (c[0] * b[0] + c[1] * b[1] + c[2] * b[2] + c[3] * b[3]) - r1 * r1 * isw;
                if (e < best) {
                    best = e;
                    bi = n;
                }
            }
            if (bi != hoc[k]) {
                const float *c = ts->hoc[k] + 4 * bi;
                hoc[k] = bi;
                chg[k] = last[k];
                changed = 1;
                P1[k] = Y1[k] - hl[bi];
                P2[k] = Y2[k] - 2.0f * (c[0] * b[0] + c[1] * b[1] + c[2] * b[2] + c[3] * b[3]) + hq[bi];
            }
        }
        if (!changed)
            break;
    }
}

/* block k's HOC row n in the harmonic domain, v[0..J) */
static void hoc_vals(const xombe_cq_ws *ws, int k, int n, float *v)
{
    const float *c = ws->ts->hoc[k] + 4 * n;
    int i, pos = ws->pos[k];
    for (i = 0; i < ws->J[k]; i++) {
        const float *H = ws->H[pos + i];
        v[i] = H[0] * c[0] + H[1] * c[1] + H[2] * c[2] + H[3] * c[3];
    }
}
#else
static const float *hv_row(const xombe_cq_ws *ws, int k, int n)
{
    return ws->hv + ws->hv_off[k] + n * ws->J[k];
}

static void hoc_search(const xombe_cq_ws *ws, const float *y, int hoc[4], int passes)
{
    const xombe_amp_tables *ts = ws->ts;
    float rest[XOMBE_LMAX + 1];
    int p, k, kk, i, n, L = ws->L, step = 0, last[4] = {-1, -1, -1, -1}, chg[4] = {-1, -1, -1, -1};
    for (p = 0; p < passes; p++) {
        int changed = 0;
        for (k = 0; k < 4; k++) {
            float c = 0.0f, so1 = 0.0f, so2 = 0.0f, best = CQ_INF;
            int bi = 0, pos = ws->pos[k], J = ws->J[k], l;
            if (ws->M[k] <= 2)
                continue;
            if (last[k] >= 0) {
                for (kk = 0; kk < 4 && (kk == k || chg[kk] < last[k]); kk++)
                    ;
                if (kk == 4)
                    continue;
            }
            last[k] = step++;
            for (kk = 0; kk < 4; kk++) {
                const float *h = kk != k && ws->M[kk] > 2 ? hv_row(ws, kk, hoc[kk]) : 0;
                for (i = 0; i < ws->J[kk]; i++) {
                    l = ws->pos[kk] + i;
                    rest[l] = h ? y[l] - h[i] : y[l];
                }
            }
            /* error is shift invariant, centre rest on its weighted mean before the one-pass sums */
            for (l = 1; l <= L; l++)
                c += ws->w[l] * rest[l];
            c /= ws->sw;
            for (l = 1; l <= L; l++) {
                rest[l] -= c;
                if (l < pos || l >= pos + J) {
                    so1 -= ws->w[l] * rest[l];
                    so2 += ws->w[l] * rest[l] * rest[l];
                }
            }
            for (n = 0; n < ts->n_hoc[k]; n++) {
                const float *v = hv_row(ws, k, n);
                float s1 = so1, s2 = so2, e;
                if (!usable(ts->ok_hoc[k], n))
                    continue;
                for (i = 0; i < J; i++) {
                    float d = v[i] - rest[pos + i], wd = ws->w[pos + i] * d;
                    s1 += wd;
                    s2 += wd * d;
                }
                e = s2 - s1 * s1 / ws->sw;
                if (e < best) {
                    best = e;
                    bi = n;
                }
            }
            if (bi != hoc[k]) {
                hoc[k] = bi;
                chg[k] = last[k];
                changed = 1;
            }
        }
        if (!changed)
            break;
    }
}

#endif

/* G2..G8 (g[2..8]) and block DCT coefficients C[k][0..M-1] of r */
static void forward(const xombe_cq_ws *ws, const float *r, float g[9], float C[4][6])
{
    float R[8];
    int k, m, i;
    for (k = 0; k < 4; k++) {
        int J = ws->J[k], M = ws->M[k], pos = ws->pos[k];
        const float *kern = xombe__dct_blk + xombe__dct_blk_off[J];
        for (m = 0; m < M; m++) {
            float s = 0.0f;
            for (i = 0; i < J; i++)
                s += r[pos + i] * kern[m * J + i];
            C[k][m] = s / (float)J;
        }
        R[2 * k] = C[k][0] + CQ_SQRT2 * C[k][1];
        R[2 * k + 1] = C[k][0] - CQ_SQRT2 * C[k][1];
    }
    for (m = 0; m < 8; m++) {
        float s = 0.0f;
        for (i = 0; i < 8; i++)
            s += xombe__dct8[m][i] * R[i];
        g[m + 1] = s / 8.0f;
    }
}

/* open-loop HOC choice */
static void hoc_open(const xombe_cq_ws *ws, float C[4][6], int hoc[4])
{
    int k;
    for (k = 0; k < 4; k++)
        hoc[k] = ws->M[k] <= 2 ? 0
                               : xombe__nearest_row(ws->ts->hoc[k], ws->ts->n_hoc[k], 4, ws->M[k] - 2, &C[k][2],
                                                    ws->ts->ok_hoc[k]);
}

/* k usable gain indices nearest g, by distance then index */
static int nearest_gains(const xombe_amp_tables *ts, float g, float Gprev, int k, int *out)
{
    float d[XOMBE_CQ_MAX_GAIN];
    int i, p, n = 0;
    for (i = 0; i < ts->n_dg; i++) {
        float e;
        if (!usable(ts->ok_dg, i))
            continue;
        e = fabsf(xombe__gain_ts(ts, i, Gprev) - g);
        if (n == k && !(e < d[k - 1]))
            continue;
        p = n < k ? n++ : k - 1;
        for (; p > 0 && e < d[p - 1]; p--) {
            d[p] = d[p - 1];
            out[p] = out[p - 1];
        }
        d[p] = e;
        out[p] = i;
    }
    return n;
}

#if XOMBE_KERNELS_Q15
/* t[l] = A3 x3 + A4 x4, plus block k's HOC row hk[k] unless hk is NULL */
static void shape(const xombe_cq_ws *ws, int i3, int i4, const int *hk, float *t)
{
    const float *x3 = ws->ts->prba24[i3], *x4 = ws->ts->prba58[i4];
    int k, i;
    for (k = 0; k < 4; k++) {
        const float *c0 = ws->c0[k], *c1 = ws->c1[k];
        float a0 = (c0[0] * x3[0] + c0[1] * x3[1] + c0[2] * x3[2]) +
                   (c0[3] * x4[0] + c0[4] * x4[1] + c0[5] * x4[2] + c0[6] * x4[3]),
              a1 = (c1[0] * x3[0] + c1[1] * x3[1] + c1[2] * x3[2]) +
                   (c1[3] * x4[0] + c1[4] * x4[1] + c1[5] * x4[2] + c1[6] * x4[3]);
        float v[XOMBE_CQ_J_MAX];
        int pos = ws->pos[k], hv = hk && ws->M[k] > 2;
        if (hv)
            hoc_vals(ws, k, hk[k], v);
        for (i = 0; i < ws->J[k]; i++)
            t[pos + i] = a0 + ws->b1[pos + i] * a1 + (hv ? v[i] : 0.0f);
    }
}

static void sub_prba(const xombe_cq_ws *ws, const float *r, int i3, int i4, float *y)
{
    int l;
    shape(ws, i3, i4, 0, y);
    for (l = 1; l <= ws->L; l++)
        y[l] = r[l] - y[l];
}

static void sub_hoc(const xombe_cq_ws *ws, const float *r, const int hoc[4], float *y)
{
    float v[XOMBE_CQ_J_MAX];
    int k, i;
    for (k = 0; k < 4; k++) {
        int pos = ws->pos[k], hv = ws->M[k] > 2;
        if (hv)
            hoc_vals(ws, k, hoc[k], v);
        for (i = 0; i < ws->J[k]; i++)
            y[pos + i] = hv ? r[pos + i] - v[i] : r[pos + i];
    }
}
#else
static void sub_prba(const xombe_cq_ws *ws, const float *r, int i3, int i4, float *y)
{
    const float *x3 = ws->ts->prba24[i3], *x4 = ws->ts->prba58[i4];
    int l;
    for (l = 1; l <= ws->L; l++) {
        const float *A = ws->A[l];
        y[l] = r[l] - (A[0] * x3[0] + A[1] * x3[1] + A[2] * x3[2]) -
               (A[3] * x4[0] + A[4] * x4[1] + A[5] * x4[2] + A[6] * x4[3]);
    }
}

static void sub_hoc(const xombe_cq_ws *ws, const float *r, const int hoc[4], float *y)
{
    int k, i;
    for (k = 0; k < 4; k++) {
        const float *h = ws->M[k] > 2 ? hv_row(ws, k, hoc[k]) : 0;
        for (i = 0; i < ws->J[k]; i++)
            y[ws->pos[k] + i] = h ? r[ws->pos[k] + i] - h[i] : r[ws->pos[k] + i];
    }
}

#endif

void xombe__cq_closed(xombe_cq_ws *ws, const float *r, const float *target, const float *base, float Gprev,
                      int n_pairs, int n_gain, int iters, xombe_cq_result *res)
{
    const xombe_amp_tables *ts = ws->ts;
    float g[9], C[4][6], y[XOMBE_LMAX + 1], best = CQ_INF;
    cq_pair top[XOMBE_CQ_MAX_PAIRS];
    int hoc[4], shoc[4], it, p, l, k, L = ws->L, pi = -1, pj = -1, searched = 0;
    forward(ws, r, g, C);
    hoc_open(ws, C, hoc);
    /* a pair search depends only on HOC (y = r - HOC's part). every search takes n_pairs, the best pair is the
       first whatever n is, so the last search is reused, HOC unchanged, by the next iteration (which then
       stops) and by the final step */
    for (it = 0; it < iters; it++) {
        if (!searched || memcmp(hoc, shoc, sizeof(hoc))) {
            sub_hoc(ws, r, hoc, y);
            prba_pairs(ws, y, n_pairs, top);
            memcpy(shoc, hoc, sizeof(hoc));
            searched = 1;
        }
        if (top[0].i == pi && top[0].j == pj)
            break;
        pi = top[0].i;
        pj = top[0].j;
        sub_prba(ws, r, pi, pj, y);
        hoc_search(ws, y, hoc, 2);
    }
    if (!searched || memcmp(hoc, shoc, sizeof(hoc))) {
        sub_hoc(ws, r, hoc, y);
        prba_pairs(ws, y, n_pairs, top);
    }
    res->err = CQ_INF;
    for (p = 0; p < n_pairs; p++) {
        float tv[XOMBE_LMAX + 1], s[XOMBE_LMAX + 1], mt = 0.0f, gc = 0.0f;
        int hk[4], gi[XOMBE_CQ_MAX_GAIN], ng, q;
        if (top[p].e >= CQ_INF)
            break;
        sub_prba(ws, r, top[p].i, top[p].j, y);
        memcpy(hk, hoc, sizeof(hk));
        hoc_search(ws, y, hk, 2);
        /* tv = A3 x3 + A4 x4 + sum H h */
#if XOMBE_KERNELS_Q15
        shape(ws, top[p].i, top[p].j, hk, tv);
#else
        {
            const float *x3 = ts->prba24[top[p].i], *x4 = ts->prba58[top[p].j];
            for (k = 0; k < 4; k++) {
                const float *h = ws->M[k] > 2 ? hv_row(ws, k, hk[k]) : 0;
                int i;
                for (i = 0; i < ws->J[k]; i++) {
                    const float *A;
                    l = ws->pos[k] + i;
                    A = ws->A[l];
                    tv[l] = (A[0] * x3[0] + A[1] * x3[1] + A[2] * x3[2]) +
                            (A[3] * x4[0] + A[4] * x4[1] + A[5] * x4[2] + A[6] * x4[3]) + (h ? h[i] : 0.0f);
                }
            }
        }
#endif
        for (l = 1; l <= L; l++)
            mt += tv[l];
        mt /= (float)L;
        for (l = 1; l <= L; l++) {
            s[l] = tv[l] - mt + base[l];
            gc += ws->w[l] * (r[l] - (tv[l] - mt));
        }
        gc = gc / ws->sw + xombe__half_log2[L];
        ng = nearest_gains(ts, gc, Gprev, n_gain, gi);
        for (q = 0; q < ng; q++) {
            float G = xombe__gain_ts(ts, gi[q], Gprev), e = 0.0f;
            for (l = 1; l <= L; l++) {
                float d = clip(s[l] + G, CQ_AMIN, XOMBE_CQ_AMAX) - target[l];
                e += ws->w[l] * d * d;
            }
            if (e < best) {
                best = e;
                res->err = e;
                res->b[0] = res->b[1] = 0;
                res->b[2] = gi[q];
                res->b[3] = top[p].i;
                res->b[4] = top[p].j;
                for (k = 0; k < 4; k++)
                    res->b[5 + k] = hk[k];
                res->u[0] = 0.0f;
                for (l = 1; l <= L; l++)
                    res->u[l] = s[l] + G;
            }
        }
    }
}

void xombe__cq_quantize(xombe_cq_ws *ws, const xombe_amp_tables *ts, const float *target, int L, const float *pt,
                        float Gprev, const float *w, int closed, int n_pairs, int n_gain, int iters, int clip_passes,
                        xombe_cq_result *res)
{
    float tg[XOMBE_LMAX + 1], pred[XOMBE_LMAX + 1], base[XOMBE_LMAX + 1], r[XOMBE_LMAX + 1];
    float teff[XOMBE_LMAX + 1], mp = 0.0f, best = CQ_INF;
    int l, cp;
    for (l = 1; l <= L; l++) {
        tg[l] = clip(target[l], CQ_AMIN, XOMBE_CQ_AMAX);
        mp += pt[l];
    }
    mp /= (float)L;
    tg[0] = pred[0] = base[0] = r[0] = teff[0] = 0.0f;
    for (l = 1; l <= L; l++) {
        pred[l] = ts->rho * (pt[l] - mp);
        base[l] = pred[l] - xombe__half_log2[L];
        r[l] = tg[l] - pred[l];
    }
    xombe__cq_prepare(ws, ts, L, w);
    if (!closed) {
        float g[9], C[4][6], mt = 0.0f;
        int gi;
        for (l = 1; l <= L; l++)
            mt += tg[l];
        nearest_gains(ts, mt / (float)L + xombe__half_log2[L], Gprev, 1, &gi);
        forward(ws, r, g, C);
        res->b[0] = res->b[1] = 0;
        res->b[2] = gi;
        res->b[3] = xombe__nearest_row(&ts->prba24[0][0], ts->n_prba24, 3, 3, &g[2], ts->ok_prba24);
        res->b[4] = xombe__nearest_row(&ts->prba58[0][0], ts->n_prba58, 4, 4, &g[5], ts->ok_prba58);
        hoc_open(ws, C, res->b + 5);
        return;
    }
    memcpy(teff, tg, sizeof(float) * (L + 1));
    for (cp = 0; cp <= clip_passes; cp++) {
        xombe_cq_result cur;
        int any = 0;
        if (cp > 0)
            for (l = 1; l <= L; l++)
                r[l] = teff[l] - pred[l];
        xombe__cq_closed(ws, r, tg, base, Gprev, n_pairs, n_gain, iters, &cur);
        if (cur.err < best) {
            best = cur.err;
            *res = cur;
        }
        /* clip-aware refinement: a harmonic whose target sits on a clip bound and whose unclipped
           reconstruction is already beyond it costs nothing, stop pulling it back */
        for (l = 1; l <= L; l++) {
            int hi = tg[l] >= XOMBE_CQ_AMAX - 1e-6f && cur.u[l] > tg[l];
            int lo = tg[l] <= CQ_AMIN + 1e-6f && cur.u[l] < tg[l];
            teff[l] = hi || lo ? cur.u[l] : tg[l];
            any |= hi || lo;
        }
        if (!any || best == 0.0f)
            break;
    }
}

#if XOMBE_HAVE_ENCODER
int xombe__cq_search_ambe2(xombe_cq_ws *ws, const xombe_pred *pi, const float *target, int L, const float *w,
                           int n_pairs, int n_gain, xombe_cq_result *res, int *bclosed, int *bopen)
{
    float tg[XOMBE_LMAX + 1], pt[XOMBE_LMAX + 1], base[XOMBE_LMAX + 1], r[XOMBE_LMAX + 1];
    float t[XOMBE_LMAX + 1], a[XOMBE_LMAX + 1], mp = 0.0f, eo = 0.0f;
    int bo[9], l, k;
    tg[0] = 0.0f;
    for (l = 1; l <= L; l++)
        tg[l] = clip(target[l], XOMBE_AMP_MIN, XOMBE_AMP_MAX);
    xombe__quantize(pi, 0, L, tg, bo);
    xombe__predict(pi, L, pt);
    for (l = 1; l <= L; l++)
        mp += pt[l];
    mp /= (float)L;
    r[0] = base[0] = 0.0f;
    for (l = 1; l <= L; l++) {
        float pred = XOMBE_RHO * (pt[l] - mp);
        base[l] = pred - xombe__half_log2[L];
        r[l] = tg[l] - pred;
    }
    xombe__cq_prepare(ws, &xombe__amp_ambe2, L, w);
    xombe__cq_closed(ws, r, tg, base, pi->gamma, n_pairs, n_gain, XOMBE_CQ_ITERS, res);
    if (bclosed)
        memcpy(bclosed, res->b, sizeof(res->b));
    if (bopen)
        memcpy(bopen, bo, sizeof(bo));
    for (k = 2; k < 9 && bo[k] == res->b[k]; k++)
        ;
    if (k == 9)
        return 0; /* same fields count as the closed choice */
    /* open choice's error against the decoder's reconstruction (clamp 15 - 2^-11) */
    xombe__spectral_shape(bo, L, t);
    xombe__assemble(t, pt, L, xombe__gain(bo[2], pi->gamma, 0), a);
    for (l = 1; l <= L; l++) {
        float d = a[l] - tg[l];
        eo += ws->w[l] * d * d;
    }
    if (res->err <= eo)
        return 0;
    for (k = 2; k < 9; k++)
        res->b[k] = bo[k];
    res->err = eo;
    return 1;
}
#endif
