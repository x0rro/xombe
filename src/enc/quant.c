#include "../internal.h"
#include "../tables/tables.h"
#include "../tables/dsp_tables.h"

void xombe__quantize(const xombe_pred *pi, int silence, int L, const float *target, int b[9])
{
    float pt[XOMBE_LMAX + 1];
    float r[XOMBE_LMAX + 1], mean_t = 0.0f, mean_p = 0.0f, R[8], G[9], C[4][7];
    int l, k, m, i, pos = 1;
    xombe__predict(pi, L, pt);
    for (l = 1; l <= L; l++) {
        mean_t += target[l];
        mean_p += pt[l];
    }
    mean_t /= (float)L;
    mean_p /= (float)L;
    {
        float dgamma = mean_t + xombe__half_log2[L] - xombe__amp_ambe2.grho * pi->gamma;
        if (silence)
            dgamma += XOMBE_SILENCE_GAIN_OFFSET;
        b[2] = xombe__nearest_row(xombe__dg, 32, 1, 1, &dgamma, NULL);
    }
    for (l = 1; l <= L; l++)
        r[l] = target[l] - XOMBE_RHO * (pt[l] - mean_p);
    for (k = 0; k < 4; k++) {
        int J = xombe__blk[L][k], top = J < 6 ? J : 6;
        const float *kern = xombe__dct_blk + xombe__dct_blk_off[J];
        for (m = 1; m <= 6; m++) {
            float s = 0.0f;
            if (m <= top)
                for (i = 0; i < J; i++)
                    s += r[pos + i] * kern[(m - 1) * J + i];
            C[k][m] = m <= top ? s / (float)J : 0.0f;
        }
        if (J < 2)
            C[k][2] = 0.0f;
        pos += J;
        R[2 * k] = C[k][1] + 1.41421356f * C[k][2];
        R[2 * k + 1] = C[k][1] - 1.41421356f * C[k][2];
    }
    for (m = 1; m <= 8; m++) {
        float s = 0.0f;
        for (i = 0; i < 8; i++)
            s += R[i] * xombe__dct8[m - 1][i];
        G[m] = s / 8.0f;
    }
    b[3] = xombe__nearest_row(&xombe__prba24[0][0], 512, 3, 3, &G[2], NULL);
    b[4] = xombe__nearest_row(&xombe__prba58[0][0], 128, 4, 4, &G[5], NULL);
    for (k = 0; k < 4; k++) {
        int J = xombe__blk[L][k], used = (J < 6 ? J : 6) - 2, rows;
        const float *tab = xombe__hoc_table(k, &rows);
        b[5 + k] = used > 0 ? xombe__nearest_row(tab, rows, 4, used, &C[k][3], NULL) : 0;
    }
}
