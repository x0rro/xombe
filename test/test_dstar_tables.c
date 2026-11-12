#include <math.h>
#include <stdio.h>
#include <string.h>
#include "../src/dstar/dstar_internal.h"
#include "dstar_ref.h"

/* generated D-STAR tables (scripts/gen_dstar_tables.py) against the model's tables_dstar record (vectors/dstar):
   every float within a float rounding of the double, the masks, VUV0, and LMPRBL exact, the pitch grid against
   its closed form, and the table sets' shapes */

static int failures, checks;

static void check(int ok, const char *what, long i)
{
    checks++;
    if (!ok) {
        failures++;
        fprintf(stderr, "FAIL %s (item %ld)\n", what, i);
    }
}

static int close_to(float t, double v)
{
    return fabs((double)t - v) <= 6e-8 * (fabs(v) > 1.0 ? fabs(v) : 1.0);
}

static void floats(const char *what, long row, const float *t, int n, const double *v, int nv)
{
    int i;
    check(n == nv, what, row);
    for (i = 0; i < n && i < nv; i++)
        check(close_to(t[i], v[i]), what, row * 100 + i);
}

static void bytes(const char *what, long row, const unsigned char *t, int n, const double *v, int nv)
{
    int i;
    check(n == nv, what, row);
    for (i = 0; i < n && i < nv; i++)
        check((double)t[i] == v[i], what, row * 1000 + i);
}

int main(int argc, char **argv)
{
    static const char *const hoc_name[4] = {"HOC1", "HOC2", "HOC3", "HOC4"};
    static const char *const ok_hoc_name[4] = {"ok_HOC1", "ok_HOC2", "ok_HOC3", "ok_HOC4"};
    const xombe_amp_tables *ts = &xombe__amp_dstar;
    double p0 = 0, p1 = 0;
    int seen[32], st, k, c_tab, c_row, c_v;
    char path[1024];
    dsref r;
    if (argc != 2) {
        fprintf(stderr, "usage: test_dstar_tables VECTORS/dstar\n");
        return 2;
    }
    snprintf(path, sizeof(path), "%s/tables_dstar.txt.gz", argv[1]);
    if (!dsref_open(&r, path))
        return 1;
    c_tab = dsref_need(&r, "tab", DSREF_STR);
    c_row = dsref_need(&r, "row", DSREF_INT);
    c_v = dsref_need(&r, "v", DSREF_FVEC);
    memset(seen, 0, sizeof(seen));
    while (c_v >= 0 && (st = dsref_next(&r)) > 0) {
        const char *t = dsref_s(&r, c_tab);
        long row = dsref_i(&r, c_row);
        const double *v = dsref_fv(&r, c_v);
        int n = dsref_len(&r, c_v);
        if (!strcmp(t, "p0"))
            p0 = v[0], seen[0]++;
        else if (!strcmp(t, "p1"))
            p1 = v[0], seen[1]++;
        else if (!strcmp(t, "rho"))
            check(close_to(XOMBE_DS_RHO, v[0]) && ts->rho == XOMBE_DS_RHO, "rho", 0), seen[2]++;
        else if (!strcmp(t, "grho"))
            check(close_to(XOMBE_DS_GRHO, v[0]) && ts->grho == XOMBE_DS_GRHO, "grho", 0), seen[3]++;
        else if (!strcmp(t, "DG") && row < 64)
            floats("DG", row, &xombe__ds_dg[row], 1, v, n), seen[4]++;
        else if (!strcmp(t, "PRBA24") && row < 512)
            floats("PRBA24", row, xombe__ds_prba24[row], 3, v, n), seen[5]++;
        else if (!strcmp(t, "PRBA58") && row < 128)
            floats("PRBA58", row, xombe__ds_prba58[row], 4, v, n), seen[6]++;
        else if (!strcmp(t, "VUV0") && row < 16)
            bytes("VUV0", row, xombe__ds_vuv[row], 8, v, n), seen[7]++;
        else if (!strcmp(t, "LMPRBL") && row < 57)
            bytes("LMPRBL", row, ts->blk[row], 4, v, n), seen[8]++;
        else if (!strcmp(t, "ok_DG"))
            bytes("ok_DG", row, ts->ok_dg, ts->n_dg, v, n), seen[9]++;
        else if (!strcmp(t, "ok_PRBA24"))
            bytes("ok_PRBA24", row, ts->ok_prba24, ts->n_prba24, v, n), seen[10]++;
        else if (!strcmp(t, "ok_PRBA58"))
            bytes("ok_PRBA58", row, ts->ok_prba58, ts->n_prba58, v, n), seen[11]++;
        else {
            for (k = 0; k < 4; k++) {
                if (!strcmp(t, hoc_name[k]) && row < ts->n_hoc[k]) {
                    floats(hoc_name[k], row, ts->hoc[k] + 4 * row, 4, v, n), seen[12 + k]++;
                    break;
                }
                if (!strcmp(t, ok_hoc_name[k])) {
                    bytes(ok_hoc_name[k], row, ts->ok_hoc[k], ts->n_hoc[k], v, n), seen[16 + k]++;
                    break;
                }
            }
            check(k < 4, t, row);
        }
    }
    if (!dsref_close(&r))
        failures++;
    {
        static const int want[20] = {1, 1, 1, 1, 64, 512, 128, 16, 57, 1, 1, 1, 16, 16, 16, 8, 1, 1, 1, 1};
        for (k = 0; k < 20; k++)
            check(seen[k] == want[k], "row count", k);
    }
    check(close_to(XOMBE_DS_P0, p0) && close_to(XOMBE_DS_P1, p1), "p0 p1", 0);
    for (k = 0; k < XOMBE_DS_B0_COUNT; k++) {
        double f0 = pow(2.0, p0 + p1 * (k + 0.5));
        int L = (int)floor(0.4627 / f0);
        L = L < 9 ? 9 : L > XOMBE_LMAX ? XOMBE_LMAX : L;
        check(close_to(xombe__ds_w0[k], f0), "W0", k);
        check(xombe__ds_lt[k] == L, "LTAB", k);
    }
    check(ts->n_dg == 64 && ts->n_prba24 == 512 && ts->n_prba58 == 128, "table set sizes", 0);
    check(xombe__ds_tone_id[0] == 129 && xombe__ds_tone_id[16] == 160 && xombe__ds_tone_id[37] == 5 &&
          xombe__ds_tone_id[154] == 122 && xombe__ds_tone_id[20] == XOMBE_DS_TONE_NONE &&
          xombe__ds_tone_id[155] == XOMBE_DS_TONE_NONE, "tone map", 0);
    {   /* the silence frame decodes cleanly to b0 124 and re-encodes to itself */
        static const uint8_t nf[9] = XOMBE_DS_NULL_FRAME;
        xombe_ds_frame f;
        uint8_t o[9];
        xombe__ds_fec_decode(nf, &f);
        check(f.ea == 0 && f.eb == 0 && f.b[0] == XOMBE_DS_B0_SILENCE && f.cls == XOMBE_DS_NULL, "null frame", 0);
        check(xombe__ds_encode_fields(f.b, o) == 0 && !memcmp(o, nf, 9), "null frame encode", 0);
    }
    {   /* out-of-range fields are clamped to their widths and still written */
        int lo[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1}, hi[9] = {999, 99, 99, 999, 999, 99, 99, 99, 99};
        int zero[9] = {0}, top[9] = {127, 15, 63, 511, 127, 15, 15, 15, 7};
        uint8_t o[9], want[9];
        xombe__ds_encode_fields(zero, want);
        check(xombe__ds_encode_fields(lo, o) == -1 && !memcmp(o, want, 9), "fields clamped up to 0", 0);
        xombe__ds_encode_fields(top, want);
        check(xombe__ds_encode_fields(hi, o) == -1 && !memcmp(o, want, 9), "fields clamped down to their widths", 0);
    }
    printf("%d checks, %d failures\n", checks, failures);
    return failures ? 1 : 0;
}
