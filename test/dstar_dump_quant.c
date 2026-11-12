#define _POSIX_C_SOURCE 200112L
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "../src/common/cquant.h"
#include "../src/dstar/dstar_internal.h"
#include "dstar_ref.h"
#if defined(XOMBE_CALLGRIND)
#include <valgrind/callgrind.h>
#endif

/* closed-loop quantizer against the reference quantizer cases:
     dstar_dump_quant [options] closed quant_closed.txt.gz OUT   (tset parent pass b err)
     dstar_dump_quant [options] dstar  quant_dstar.txt.gz OUT    (id b ahat err)
     dstar_dump_quant [options] ambe2  quant_ambe2.txt.gz OUT    (id bopen bclosed pick b ahat err)
   extra columns: errpy, the C error of the reference's fields (closed, dstar), for mismatch analysis, and w
   (dstar, ambe2), recomputed with xombe__cq_amp_weights from the target and checked against the record
   options:
     --brute       no pruning in the pair search
     --bench REPS  time the quantizer over every call after the dump, best of REPS, to stderr
     --top N       with --bench, the N slowest calls' ids
     --ids a,b,... only these calls, with -DXOMBE_CALLGRIND a callgrind dump after each
     --cap N       dump the search capped at N, fail if a pair search exceeds max(N + 4, 132) units or the
                   summed error grows by more than --cap-tol REL (default 1e-3) against the exact search */

typedef struct {
    int kind; /* 0 closed, 1 dstar, 2 ambe2 */
    int tset; /* closed: 0 dstar, 1 ambe2 */
    long id, pass;
    int L, n_pairs, n_gain, iters, clip_passes, has_w, mode_closed;
    float Gprev;
    float r[XOMBE_LMAX + 1], target[XOMBE_LMAX + 1], base[XOMBE_LMAX + 1], w[XOMBE_LMAX + 1];
    int Lp;
    float ap[XOMBE_LMAX + 1];
    int bref[9];
} qcall;

typedef struct {
    xombe_cq_result res;
    int pick, bopen[9], bclosed[9];
    float ahat[XOMBE_LMAX + 1], err, errpy;
} qout;

static xombe_cq_ws ws;
static int brute, qcap;

static void fill(float *dst, const double *src, int n)
{
    int i;
    dst[0] = 0.0f;
    for (i = 0; i < n; i++)
        dst[i + 1] = (float)src[i];
}

static const xombe_amp_tables *set_of(const qcall *c)
{
    return c->kind == 1 || (c->kind == 0 && c->tset == 0) ? &xombe__amp_dstar : &xombe__amp_ambe2;
}

/* decoder's reconstruction of fields b against the memory */
static void recon(const qcall *c, const int b[9], float *a)
{
    const xombe_amp_tables *ts = set_of(c);
    float t[XOMBE_LMAX + 1], pt[XOMBE_LMAX + 1];
    xombe__predict_from(c->Lp, c->ap, c->L, pt);
    xombe__spectral_shape_ts(ts, b, c->L, t);
    xombe__assemble_ts(ts, t, pt, c->L, xombe__gain_ts(ts, b[2], c->Gprev), a);
}

static float werr(const qcall *c, const float *a, float amax)
{
    float e = 0.0f;
    int l;
    for (l = 1; l <= c->L; l++) {
        float t = c->target[l] < -15.0f ? -15.0f : c->target[l] > amax ? amax : c->target[l];
        float d = a[l] - t;
        e += (c->has_w ? c->w[l] : 1.0f) * d * d;
    }
    return e;
}

/* search's clipped error of fields b on a prepared workspace */
static float closed_err(const qcall *c, const int b[9])
{
    const xombe_amp_tables *ts = set_of(c);
    float t[XOMBE_LMAX + 1], mt = 0.0f, G = xombe__gain_ts(ts, b[2], c->Gprev), e = 0.0f;
    int l;
    xombe__spectral_shape_ts(ts, b, c->L, t);
    for (l = 1; l <= c->L; l++)
        mt += t[l];
    mt /= (float)c->L;
    for (l = 1; l <= c->L; l++) {
        float a = t[l] - mt + c->base[l] + G, d;
        a = a < -15.0f ? -15.0f : a > XOMBE_CQ_AMAX ? XOMBE_CQ_AMAX : a;
        d = a - c->target[l];
        e += (c->has_w ? c->w[l] : 1.0f) * d * d;
    }
    return e;
}

static void run(const qcall *c, qout *o)
{
    const float *w = c->has_w ? c->w : NULL;
    ws.brute = brute;
    ws.cap = qcap;
    if (c->kind == 0) {
        xombe__cq_prepare(&ws, set_of(c), c->L, w);
        xombe__cq_closed(&ws, c->r, c->target, c->base, c->Gprev, c->n_pairs, c->n_gain, c->iters, &o->res);
    } else if (c->kind == 1) {
        float pt[XOMBE_LMAX + 1];
        xombe__predict_from(c->Lp, c->ap, c->L, pt);
        xombe__cq_quantize(&ws, &xombe__amp_dstar, c->target, c->L, pt, c->Gprev, w, c->mode_closed, c->n_pairs,
                           c->n_gain, c->iters, c->clip_passes, &o->res);
    } else {
        xombe_pred pi;
        xombe__pred_init(&pi);
        pi.L = c->Lp;
        pi.gamma = c->Gprev;
        memcpy(pi.a, c->ap, sizeof(pi.a));
        o->pick = xombe__cq_search_ambe2(&ws, &pi, c->target, c->L, w, c->n_pairs, c->n_gain, &o->res, o->bclosed,
                                         o->bopen);
    }
}

static int cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return x < y ? -1 : x > y;
}

static double now(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + 1e-9 * t.tv_nsec;
}

static int in_ids(const char *ids, long id)
{
    const char *p = ids;
    while (p && *p) {
        char *e;
        long v = strtol(p, &e, 10);
        if (e == p)
            break;
        if (v == id)
            return 1;
        p = *e == ',' ? e + 1 : e;
    }
    return 0;
}

int main(int argc, char **argv)
{
    static const char *cols[3] = {
        "tset:s parent:i pass:i b:i* err:f errpy:f",
        "id:i w:f* b:i* ahat:f* err:f errpy:f",
        "id:i w:f* bopen:i* bclosed:i* pick:s b:i* ahat:f* err:f"};
    const char *ids = NULL;
    int kind, reps = 0, ntop = 0, ai = 1, st, ok = 1, worst = 0, nchg = 0;
    long n = 0, cap = 0, i;
    double e_exact = 0.0, e_cap = 0.0, cap_tol = 1e-3;
    qcall *calls = NULL;
    dsref r;
    dsref_writer wr;
    while (ai < argc && argv[ai][0] == '-') {
        if (!strcmp(argv[ai], "--brute"))
            brute = 1;
        else if (!strcmp(argv[ai], "--bench") && ai + 1 < argc)
            reps = atoi(argv[++ai]);
        else if (!strcmp(argv[ai], "--top") && ai + 1 < argc)
            ntop = atoi(argv[++ai]);
        else if (!strcmp(argv[ai], "--ids") && ai + 1 < argc)
            ids = argv[++ai];
        else if (!strcmp(argv[ai], "--cap") && ai + 1 < argc)
            qcap = atoi(argv[++ai]);
        else if (!strcmp(argv[ai], "--cap-tol") && ai + 1 < argc)
            cap_tol = atof(argv[++ai]);
        else
            break;
        ai++;
    }
    if (argc - ai != 3) {
        fprintf(stderr, "usage: dstar_dump_quant [--brute] [--bench REPS [--top N]] [--ids a,b] [--cap N [--cap-tol REL]] "
                        "closed|dstar|ambe2 IN OUT\n");
        return 2;
    }
    kind = !strcmp(argv[ai], "closed") ? 0 : !strcmp(argv[ai], "dstar") ? 1 : !strcmp(argv[ai], "ambe2") ? 2 : -1;
    if (kind < 0 || !dsref_open(&r, argv[ai + 1]))
        return 2;
    {
        int c_L = dsref_need(&r, "L", DSREF_INT), c_np = dsref_need(&r, "n_pairs", DSREF_INT);
        int c_ng = dsref_need(&r, "n_gain", DSREF_INT), c_t = dsref_need(&r, "target", DSREF_FVEC);
        int c_w = dsref_need(&r, "w", DSREF_FVEC), c_b = dsref_need(&r, "b", DSREF_IVEC);
        int c_id = dsref_need(&r, kind == 0 ? "parent" : "id", DSREF_INT);
        int c_ts = -1, c_pass = -1, c_it = -1, c_G = -1, c_r = -1, c_base = -1, c_cp = -1, c_mode = -1;
        int c_Lp = -1, c_ap = -1;
        if (kind == 0) {
            c_ts = dsref_need(&r, "tset", DSREF_STR);
            c_pass = dsref_need(&r, "pass", DSREF_INT);
            c_it = dsref_need(&r, "iters", DSREF_INT);
            c_G = dsref_need(&r, "Gprev", DSREF_FLOAT);
            c_r = dsref_need(&r, "r", DSREF_FVEC);
            c_base = dsref_need(&r, "base", DSREF_FVEC);
        } else if (kind == 1) {
            c_it = dsref_need(&r, "iters", DSREF_INT);
            c_cp = dsref_need(&r, "clip_passes", DSREF_INT);
            c_mode = dsref_need(&r, "mode", DSREF_STR);
            c_G = dsref_need(&r, "stG", DSREF_FLOAT);
            c_Lp = dsref_need(&r, "stL", DSREF_INT);
            c_ap = dsref_need(&r, "sta", DSREF_FVEC);
        } else {
            c_mode = dsref_need(&r, "mode", DSREF_STR);
            c_G = dsref_need(&r, "piG", DSREF_FLOAT);
            c_Lp = dsref_need(&r, "piL", DSREF_INT);
            c_ap = dsref_need(&r, "pia", DSREF_FVEC);
        }
        while ((st = dsref_next(&r)) > 0) {
            qcall *c;
            int k, L;
            if (ids && !in_ids(ids, dsref_i(&r, c_id)))
                continue;
            if (n == cap) {
                cap = cap ? 2 * cap : 1024;
                calls = realloc(calls, cap * sizeof(*calls));
                if (!calls)
                    return 2;
            }
            c = &calls[n++];
            memset(c, 0, sizeof(*c));
            c->kind = kind;
            c->id = dsref_i(&r, c_id);
            c->L = L = (int)dsref_i(&r, c_L);
            c->n_pairs = (int)dsref_i(&r, c_np);
            c->n_gain = (int)dsref_i(&r, c_ng);
            c->iters = c_it >= 0 ? (int)dsref_i(&r, c_it) : 2;
            c->clip_passes = c_cp >= 0 ? (int)dsref_i(&r, c_cp) : 0;
            c->mode_closed = c_mode < 0 || !strcmp(dsref_s(&r, c_mode), "closed");
            c->Gprev = (float)dsref_f(&r, c_G);
            fill(c->target, dsref_fv(&r, c_t), L);
            c->has_w = dsref_len(&r, c_w) > 0;
            if (c->has_w)
                fill(c->w, dsref_fv(&r, c_w), L);
            for (k = 0; k < 7; k++)
                c->bref[2 + k] = (int)dsref_iv(&r, c_b)[k];
            if (kind == 0) {
                c->tset = strcmp(dsref_s(&r, c_ts), "dstar") ? 1 : 0;
                c->pass = dsref_i(&r, c_pass);
                fill(c->r, dsref_fv(&r, c_r), L);
                fill(c->base, dsref_fv(&r, c_base), L);
            } else {
                int la = dsref_len(&r, c_ap), l;
                c->Lp = (int)dsref_i(&r, c_Lp);
                fill(c->ap, dsref_fv(&r, c_ap), la);
                for (l = la + 1; l <= XOMBE_LMAX; l++)
                    c->ap[l] = c->ap[la];
            }
            if (c->L < 1 || c->L > XOMBE_LMAX || c->n_pairs > XOMBE_CQ_MAX_PAIRS || c->n_gain > XOMBE_CQ_MAX_GAIN) {
                fprintf(stderr, "dstar_dump_quant: bad call %ld\n", c->id);
                return 2;
            }
        }
    }
    if (!dsref_close(&r) || st < 0)
        return 2;
    if (!dsref_wopen(&wr, argv[ai + 2], argv[ai], cols[kind], 9))
        return 2;
    for (i = 0; i < n; i++) {
        const qcall *c = &calls[i];
        qout o;
        float a[XOMBE_LMAX + 1];
#if defined(XOMBE_CALLGRIND)
        CALLGRIND_TOGGLE_COLLECT;
#endif
        run(c, &o);
#if defined(XOMBE_CALLGRIND)
        CALLGRIND_TOGGLE_COLLECT;
        if (ids)
            CALLGRIND_DUMP_STATS;
#endif
        if (qcap > 0) {
            qout x;
            int k, d = 0, cq;
            if (ws.work > worst)
                worst = ws.work;
            cq = qcap;
            qcap = 0;
            run(c, &x);
            qcap = cq;
            for (k = 2; k < 9; k++)
                d |= x.res.b[k] != o.res.b[k];
            nchg += d;
            e_exact += x.res.err;
            e_cap += o.res.err;
        }
        if (kind == 0) {
            dsref_ws(&wr, c->tset ? "ambe2" : "dstar");
            dsref_wi(&wr, c->id);
            dsref_wi(&wr, c->pass);
            dsref_wiv_int(&wr, o.res.b + 2, 7);
            dsref_wf(&wr, o.res.err);
            dsref_wf(&wr, closed_err(c, c->bref));
        } else {
            recon(c, o.res.b, a);
            dsref_wi(&wr, c->id);
            if (c->has_w) { /* amp_weights(target, 0.5) against the recorded weights */
                float wc[XOMBE_LMAX + 1];
                xombe__cq_amp_weights(c->target, c->L, XOMBE_CQ_QWEIGHT, XOMBE_CQ_WEIGHT_FLOOR, wc);
                dsref_wfv_float(&wr, wc + 1, c->L);
            } else
                dsref_wfv_float(&wr, NULL, 0);
            if (kind == 2) {
                dsref_wiv_int(&wr, o.bopen + 2, 7);
                dsref_wiv_int(&wr, o.bclosed + 2, 7);
                dsref_ws(&wr, o.pick ? "open" : "closed");
            }
            dsref_wiv_int(&wr, o.res.b + 2, 7);
            dsref_wfv_float(&wr, a + 1, c->L);
            dsref_wf(&wr, werr(c, a, kind == 1 ? XOMBE_DS_AMP_MAX : XOMBE_AMP_MAX));
            if (kind == 1) {
                recon(c, c->bref, a);
                dsref_wf(&wr, werr(c, a, XOMBE_DS_AMP_MAX));
            }
        }
        dsref_wend(&wr);
    }
    ok = dsref_wclose(&wr);
    if (qcap > 0) {
        int bound = qcap + 4 > 132 ? qcap + 4 : 132;
        fprintf(stderr, "cap %d: %d of %ld calls changed, error %+.4f%%, most units in a search %d (bound %d)\n", qcap,
                nchg, n, 100.0 * (e_cap / e_exact - 1.0), worst, bound);
        if (worst > bound || e_cap > (1.0 + cap_tol) * e_exact)
            ok = 0;
    }
    if (reps > 0 && n > 0) {
        double best = 1e30, *ct = malloc(n * sizeof(double));
        int rep;
        for (i = 0; i < n; i++)
            ct[i] = 1e30;
        for (rep = 0; rep < reps; rep++) {
            double t0 = now(), t1;
            for (i = 0; i < n; i++) {
                qout o;
                double s = ntop ? now() : 0.0;
                run(&calls[i], &o);
                if (ntop) {
                    double d = now() - s;
                    if (d < ct[i])
                        ct[i] = d;
                }
            }
            t1 = now() - t0;
            if (t1 < best)
                best = t1;
        }
        fprintf(stderr, "bench %s: %ld calls, best of %d: %.3f us/call\n", argv[ai], n, reps, 1e6 * best / n);
        if (ntop) {
            /* per-call distribution (each call's best time) */
            double *sorted = malloc(n * sizeof(double));
            static const double pc[5] = {0.5, 0.9, 0.99, 0.999, 1.0};
            int p;
            memcpy(sorted, ct, n * sizeof(double));
            qsort(sorted, n, sizeof(double), cmp_double);
            fprintf(stderr, "  per call (us):");
            for (p = 0; p < 5; p++) {
                long k = (long)(pc[p] * (n - 1) + 0.5);
                fprintf(stderr, " p%g %.2f", 100.0 * pc[p], 1e6 * sorted[k]);
            }
            fprintf(stderr, "\n");
            free(sorted);
        }
        for (rep = 0; rep < ntop; rep++) {
            long m = -1;
            for (i = 0; i < n; i++)
                if (ct[i] >= 0 && (m < 0 || ct[i] > ct[m]))
                    m = i;
            if (m < 0)
                break;
            fprintf(stderr, "  slow: id %ld pass %ld L %d  %.2f us\n", calls[m].id, calls[m].pass, calls[m].L,
                    1e6 * ct[m]);
            ct[m] = -1.0;
        }
        free(ct);
    }
    free(calls);
#ifdef XOMBE_CQ_STATS
    {
        extern long xombe__cq_stat[4];
        double n = (double)xombe__cq_stat[0];
        fprintf(stderr, "stats: %ld pair searches, %.1f rows, %.1f pairs, and %.1f columns per search\n",
                xombe__cq_stat[0], xombe__cq_stat[1] / n, xombe__cq_stat[2] / n, xombe__cq_stat[3] / n);
    }
#endif
    return ok ? 0 : 1;
}
