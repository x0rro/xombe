#define _POSIX_C_SOURCE 200112L
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "../src/state.h"
#include "../src/dstar/dstar_enc.h"
#include "dstar_ref.h"
#if defined(XOMBE_CALLGRIND)
#include <valgrind/callgrind.h>
#endif

/* D-STAR encoder back end (src/dstar/dstar_enc.c) on the model's recorded inputs (enc stage):
     dstar_dump_enc_be C.enc.txt.gz OUT          n cls b0t mn b1 b out
       b0t  xombe__ds_b0_from_f0(1/Pr) of the recorded refined period
       mn   xombe__ds_m_noise(N, b0t) of the recorded window and pitch
       b1   xombe__ds_b1_search of the recorded a and w
       b    the recorded targets at, b0 and b1 quantized against the mirror run over the stream (xombe__ds_enc_amps;
            the mirror commits these fields, Knox frames commit their recorded fields, tone and NULL frames skip)
       out  the 9 bytes of b
     dstar_dump_enc_be --spec C.enc.txt.gz OUT   n cls mv w a b1
       the voicing on the C front end's spectrum: C.in.raw through the input high-pass, spectrum at window N centred
       on 160n + 120, the metric and weights at the recorded b0t's f0, a with the recorded second-half VAD
       decision (the VAD gate), b1 from those
   non-voice frames are written with the reference's values, so rows stay matched
   --bench REPS: (default mode) time xombe__ds_enc_amps over the stream REPS times, best to stderr, with
   -DXOMBE_CALLGRIND the voice frames run between CALLGRIND_START/STOP_INSTRUMENTATION */

#define MAXF 20000
#define HPF_GAIN 0.9427f /* must match INPUT_HPF_GAIN in enc/encoder.c */
#define PAD 400

typedef struct {
    int n, voice, knox, N, b0t, b0, b1, L, vad2;
    char cls[8];
    float Pr, a[8], w[8], mn[8], at[XOMBE_LMAX + 1];
    int b[9];
    uint8_t out[9];
} rec;

static rec R[MAXF];
static int nrec;

static int load(const char *path)
{
    dsref r;
    int st;
    if (!dsref_open(&r, path))
        return 0;
    int c_n = dsref_need(&r, "n", DSREF_INT), c_cls = dsref_need(&r, "cls", DSREF_STR),
        c_N = dsref_need(&r, "N", DSREF_INT), c_Pr = dsref_need(&r, "Pr", DSREF_FLOAT),
        c_b0t = dsref_need(&r, "b0t", DSREF_INT), c_b0 = dsref_need(&r, "b0", DSREF_INT),
        c_L = dsref_need(&r, "L", DSREF_INT), c_w = dsref_need(&r, "w", DSREF_FVEC),
        c_mn = dsref_need(&r, "mn", DSREF_FVEC), c_a = dsref_need(&r, "a", DSREF_FVEC),
        c_b1 = dsref_need(&r, "b1", DSREF_INT), c_at = dsref_need(&r, "at", DSREF_FVEC),
        c_b = dsref_need(&r, "b", DSREF_IVEC), c_out = dsref_need(&r, "out", DSREF_HEX),
        c_vad2 = dsref_need(&r, "vad2", DSREF_INT);
    while ((st = dsref_next(&r)) > 0 && nrec < MAXF) {
        rec *q = &R[nrec++];
        int j;
        memset(q, 0, sizeof(*q));
        q->n = (int)dsref_i(&r, c_n);
        snprintf(q->cls, sizeof(q->cls), "%s", dsref_s(&r, c_cls));
        q->voice = !strcmp(q->cls, "voice");
        q->knox = !strcmp(q->cls, "knox");
        q->N = (int)dsref_i(&r, c_N);
        q->Pr = (float)dsref_f(&r, c_Pr);
        q->b0t = (int)dsref_i(&r, c_b0t);
        q->b0 = (int)dsref_i(&r, c_b0);
        q->b1 = (int)dsref_i(&r, c_b1);
        q->L = (int)dsref_i(&r, c_L);
        q->vad2 = (int)dsref_i(&r, c_vad2);
        if (q->voice) {
            if (dsref_len(&r, c_a) != 8 || dsref_len(&r, c_w) != 8 || dsref_len(&r, c_mn) != 8 ||
                dsref_len(&r, c_at) != q->L) {
                fprintf(stderr, "%s: frame %d: bad vector lengths\n", path, q->n);
                return 0;
            }
            for (j = 0; j < 8; j++) {
                q->a[j] = (float)dsref_fv(&r, c_a)[j];
                q->w[j] = (float)dsref_fv(&r, c_w)[j];
                q->mn[j] = (float)dsref_fv(&r, c_mn)[j];
            }
            for (j = 0; j < q->L; j++)
                q->at[j + 1] = (float)dsref_fv(&r, c_at)[j];
        }
        for (j = 0; j < dsref_len(&r, c_b) && j < 9; j++)
            q->b[j] = (int)dsref_iv(&r, c_b)[j];
        if (dsref_len(&r, c_out) == 9)
            memcpy(q->out, dsref_bytes(&r, c_out), 9);
    }
    return dsref_close(&r) && st == 0;
}

static xombe_ds_enc ds;

static void run_amps(int b[][9])
{
    int i;
    xombe__ds_enc_reset(&ds);
#if defined(XOMBE_CALLGRIND)
    CALLGRIND_START_INSTRUMENTATION;
#endif
    for (i = 0; i < nrec; i++) {
        rec *q = &R[i];
        if (q->voice)
            xombe__ds_enc_amps(&ds, q->b0, q->b1, q->at, b[i]);
        else if (q->knox)
            xombe__ds_mirror_commit(&ds.mirror, q->b, NULL);
    }
#if defined(XOMBE_CALLGRIND)
    CALLGRIND_STOP_INSTRUMENTATION;
#endif
}

static int dump_exact(const char *out, int bench)
{
    static int b[MAXF][9];
    dsref_writer w;
    int i, j;
    run_amps(b);
    if (!dsref_wopen(&w, out, "enc", "n:i cls:s b0t:i mn:f* b1:i b:i* out:x", 7))
        return 1;
    for (i = 0; i < nrec; i++) {
        rec *q = &R[i];
        dsref_wi(&w, q->n);
        dsref_ws(&w, q->cls);
        if (q->voice) {
            float mn[8];
            uint8_t o[9];
            xombe__ds_m_noise(q->N, q->b0t, mn);
            dsref_wi(&w, xombe__ds_b0_from_f0(1.0f / q->Pr));
            dsref_wfv_float(&w, mn, 8);
            dsref_wi(&w, xombe__ds_b1_search(q->a, q->w));
            dsref_wiv_int(&w, b[i], 9);
            xombe__ds_encode_fields(b[i], o);
            dsref_wx(&w, o, 9);
        } else {
            long bl[9];
            int nb = q->b[0] >= 124 && !q->knox ? 1 : 9;
            for (j = 0; j < nb; j++)
                bl[j] = q->b[j];
            dsref_wi(&w, q->b0t);
            dsref_wfv_float(&w, q->mn, 0);
            dsref_wi(&w, q->b1);
            dsref_wiv(&w, bl, nb);
            dsref_wx(&w, q->out, 9);
        }
        dsref_wend(&w);
    }
    if (bench > 0) {
        double best = 1e30;
        int k, nv = 0;
        for (i = 0; i < nrec; i++)
            nv += R[i].voice;
        for (k = 0; k < bench; k++) {
            struct timespec t0, t1;
            double dt;
            clock_gettime(CLOCK_MONOTONIC, &t0);
            run_amps(b);
            clock_gettime(CLOCK_MONOTONIC, &t1);
            dt = (double)(t1.tv_sec - t0.tv_sec) + 1e-9 * (double)(t1.tv_nsec - t0.tv_nsec);
            if (dt < best)
                best = dt;
        }
        fprintf(stderr, "amps: %d voice frames, %.1f us per frame\n", nv, 1e6 * best / (nv ? nv : 1));
    }
    return dsref_wclose(&w) ? 0 : 1;
}

static int dump_spec(const char *encpath, const char *out)
{
    static unsigned char mem[XOMBE_ENC_SIZE];
    static xombe_spec S;
    char pcmpath[1024];
    size_t n = strlen(encpath);
    long ns, i;
    int16_t *x;
    float *xh;
    xombe_hpf hpf;
    xombe_enc *e = xombe_enc_init(mem, NULL);
    dsref_writer w;
    if (n < 11 || strcmp(encpath + n - 11, ".enc.txt.gz") || n > sizeof(pcmpath) - 8) {
        fprintf(stderr, "%s: expected C.enc.txt.gz\n", encpath);
        return 1;
    }
    memcpy(pcmpath, encpath, n - 11);
    strcpy(pcmpath + n - 11, ".in.raw");
    if (!(x = dsref_read_pcm(pcmpath, &ns)))
        return 1;
    xh = calloc((size_t)(ns + 2 * PAD + 160), sizeof(float));
    xombe__hpf_init(&hpf, HPF_GAIN);
    for (i = 0; i < ns; i++)
        xh[PAD + i] = xombe__hpf_step(&hpf, (float)x[i]);
    for (; i < ns + PAD + 160; i++)
        xh[PAD + i] = xombe__hpf_step(&hpf, 0.0f);
    if (!dsref_wopen(&w, out, "enc", "n:i cls:s mv:f* w:f* a:f* b1:i", 7))
        return 1;
    for (i = 0; i < nrec; i++) {
        rec *q = &R[i];
        dsref_wi(&w, q->n);
        dsref_ws(&w, q->cls);
        if (q->voice) {
            float m[8], wv[8], mn[8], a[8], f0 = xombe__ds_f0(q->b0t);
            xombe__enc_set_window(e, q->N);
            xombe__kern_spectrum(xh + PAD + 160 * q->n + 120, e->win, q->N, &S);
            xombe__voicing_metric(e, &S, f0, m, wv);
            xombe__ds_m_noise(q->N, q->b0t, mn);
            xombe__ds_soft(m, mn, XOMBE_DS_S_VN, q->vad2, a);
            dsref_wfv_float(&w, m, 8);
            dsref_wfv_float(&w, wv, 8);
            dsref_wfv_float(&w, a, 8);
            dsref_wi(&w, xombe__ds_b1_search(a, wv));
        } else {
            dsref_wfv_float(&w, NULL, 0);
            dsref_wfv_float(&w, NULL, 0);
            dsref_wfv_float(&w, NULL, 0);
            dsref_wi(&w, q->b1);
        }
        dsref_wend(&w);
    }
    free(x);
    free(xh);
    return dsref_wclose(&w) ? 0 : 1;
}

int main(int argc, char **argv)
{
    int spec = 0, bench = 0, i = 1;
    for (; i < argc && !strncmp(argv[i], "--", 2); i++) {
        if (!strcmp(argv[i], "--spec"))
            spec = 1;
        else if (!strcmp(argv[i], "--bench") && i + 1 < argc)
            bench = atoi(argv[++i]);
        else
            break;
    }
    if (argc - i != 2) {
        fprintf(stderr, "usage: dstar_dump_enc_be [--spec] [--bench REPS] C.enc.txt.gz OUT\n");
        return 2;
    }
    if (!load(argv[i]))
        return 1;
    return spec ? dump_spec(argv[i], argv[i + 1]) : dump_exact(argv[i + 1], bench);
}
