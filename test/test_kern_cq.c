/* closed-loop quantizer (src/common/cquant.c): the pruned pair search against the exhaustive one on random
   shaped targets for both table sets, and with the q15 kernels the codebook images against the codebooks */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "../src/common/cquant.h"
#include "../src/dstar/dstar_internal.h"
#include "../src/tables/dsp_tables.h"

#define CASES 300

static unsigned rng = 12345u;

static float urand(void) /* [0, 1) */
{
    rng = rng * 1664525u + 1013904223u;
    return (float)(rng >> 8) / 16777216.0f;
}

static int irand(int n)
{
    return (int)(urand() * (float)n) % n;
}

static xombe_cq_ws ws;

/* random fields' shape plus noise, quantized with and without pruning, uncapped */
static int brute_cases(const xombe_amp_tables *ts, const char *name)
{
    int c, bad = 0;
    for (c = 0; c < CASES; c++) {
        float t[XOMBE_LMAX + 1], target[XOMBE_LMAX + 1], base[XOMBE_LMAX + 1], w[XOMBE_LMAX + 1], G;
        float noise = 0.05f + 0.6f * urand();
        int b[9] = {0}, L = 9 + irand(48), l, k, pass;
        xombe_cq_result res[2];
        b[3] = irand(ts->n_prba24);
        b[4] = irand(ts->n_prba58);
        for (k = 0; k < 4; k++)
            b[5 + k] = irand(ts->n_hoc[k]);
        xombe__spectral_shape_ts(ts, b, L, t);
        G = 2.0f * urand() - 1.0f;
        target[0] = base[0] = 0.0f;
        for (l = 1; l <= L; l++) {
            target[l] = t[l] + 3.0f * G + noise * (2.0f * urand() - 1.0f);
            base[l] = -xombe__half_log2[L];
        }
        xombe__cq_amp_weights(target, L, XOMBE_CQ_QWEIGHT, XOMBE_CQ_WEIGHT_FLOOR, w);
        for (pass = 0; pass < 2; pass++) {
            memset(&ws, 0, sizeof(ws));
            ws.brute = pass;
            xombe__cq_prepare(&ws, ts, L, c & 1 ? w : NULL);
            xombe__cq_closed(&ws, target, target, base, 0.0f, XOMBE_CQ_N_PAIRS, XOMBE_CQ_N_GAIN, XOMBE_CQ_ITERS,
                             &res[pass]);
        }
        for (k = 2; k < 9 && res[0].b[k] == res[1].b[k]; k++)
            ;
        if (k < 9 && bad++ < 5)
            fprintf(stderr, "%s case %d (L %d): pruned b%d %d, exhaustive %d\n", name, c, L, k, res[0].b[k],
                    res[1].b[k]);
    }
    printf("kern_cq: %s %d/%d cases pruned = exhaustive\n", name, CASES - bad, CASES);
    return bad;
}

#if XOMBE_KERNELS_Q15
/* every image entry within 1/2 unit of its exact value and within its column's maximum */
static int image_check(const xombe_amp_tables *ts, const xombe_cq_q15_tab *qt, const char *name)
{
    static const int p24[6][2] = {{0, 0}, {1, 1}, {2, 2}, {0, 1}, {0, 2}, {1, 2}};
    static const int p58[10][2] = {{0, 0}, {1, 1}, {2, 2}, {3, 3}, {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
    int i, k, bad = 0;
    for (i = 0; i < ts->n_prba24; i++) {
        const float *x = ts->prba24[i];
        for (k = 0; k < 9; k++) {
            uint32_t v = qt->r24[i][k / 2];
            int q = (int16_t)(k & 1 ? v >> 16 : v & 0xffff), t = k < 6 ? qt->tm24 : qt->tx24;
            double e = k < 6 ? (p24[k][0] == p24[k][1] ? 1.0 : 2.0) * x[p24[k][0]] * x[p24[k][1]] : x[k - 6];
            double m = k < 6 ? qt->m24[k] : qt->x24[k - 6];
            bad += fabs(q - ldexp(e, t)) > 0.5 || fabs(ldexp(q, -t)) > m;
        }
    }
    for (i = 0; i < ts->n_prba58; i++) {
        const float *x = ts->prba58[i];
        for (k = 0; k < 14; k++) {
            uint32_t v = qt->r58[i][k / 2];
            int q = (int16_t)(k & 1 ? v >> 16 : v & 0xffff), t = k < 10 ? qt->tm58 : qt->tx58;
            double e = k < 10 ? (p58[k][0] == p58[k][1] ? 1.0 : 2.0) * x[p58[k][0]] * x[p58[k][1]] : x[k - 10];
            double m = k < 10 ? qt->m58[k] : qt->x58[k - 10];
            bad += fabs(q - ldexp(e, t)) > 0.5 || fabs(ldexp(q, -t)) > m;
        }
    }
    printf("kern_cq: %s images %s\n", name, bad ? "WRONG, rerun scripts/gen_cq_tables.py" : "match the codebooks");
    return bad;
}
#endif

int main(int argc, char **argv)
{
    int bad = 0;
    (void)argc;
    (void)argv;
#if XOMBE_KERNELS_Q15
    bad += image_check(&xombe__amp_ambe2, &xombe__cq_q15_ambe2, "AMBE+2");
    bad += image_check(&xombe__amp_dstar, &xombe__cq_q15_dstar, "D-STAR");
#endif
    bad += brute_cases(&xombe__amp_ambe2, "AMBE+2");
    bad += brute_cases(&xombe__amp_dstar, "D-STAR");
    return bad != 0;
}
