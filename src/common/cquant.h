#ifndef XOMBE_CQUANT_H
#define XOMBE_CQUANT_H

/* weighted closed-loop spectral-amplitude quantizer (b2..b8) in float32, on either table set
   model (params.h): a = clip((Tv - mean Tv) + base + G), base = rho*(Pt - mean Pt) - 1/2 log2 L,
   G = min(dg[b2] + grho*Gprev, gmax), Tv = A3*prba24[b3] + A4*prba58[b4] + sum_k H_k*hoc_k[b5+k], matrices per L
   search, all errors weighted by w and mean-removed:
     1. HOC open loop (nearest in the block-DCT domain), then up to iters rounds of the best (b3, b4) pair over
        all 512 x 128 given HOC, then HOC per block given the pair (two passes)
     2. the best n_pairs pairs given that HOC, each with its own HOC search and the n_gain gains nearest the
        weighted continuous gain, scored by the exact clipped reconstruction
   entries an ok mask marks unusable are skipped. the search clips to [-15, 15] for both sets, not to amax
   with ws->cap > 0 a pair search stops after cap units of work (a pair 1, a row 4), keeping the best pairs so
   far. its lowest-bound row always completes, the others run in ascending bound order, so a search takes at
   most max(cap + 4, 132) units. cap 0 is exact
   arrays are indexed 1..L (0 unused), fields come back in b[2..8] */

#include "../internal.h"
#if XOMBE_KERNELS_Q15
#include "../tables/cq_q15_tables.h"
#endif

#define XOMBE_CQ_MAX_PAIRS 8
#define XOMBE_CQ_MAX_GAIN 8
#define XOMBE_CQ_HOC_ROWS 72 /* sum of the HOC row counts, AMBE+2's 32 + 16 + 16 + 8 */
#define XOMBE_CQ_J_MAX 17    /* longest spectral block */
#define XOMBE_CQ_AMAX 15.0f  /* search's upper clip, both sets */
#define XOMBE_CQ_CAP 1250    /* encoder and transcoder ws.cap */
/* encoder and transcoder settings: weights gamma and floor, then n_pairs, n_gain, iters, and clip_passes */
#define XOMBE_CQ_QWEIGHT 0.5f
#define XOMBE_CQ_WEIGHT_FLOOR 0.05f
#define XOMBE_CQ_N_PAIRS 4
#define XOMBE_CQ_N_GAIN 3
#define XOMBE_CQ_ITERS 2
#define XOMBE_CQ_CLIP_PASSES 2

#if XOMBE_KERNELS_Q15
/* q15 kernels: the same search with the pair-search bounds in int16 SIMD against the codebook images
   (cq_q15_tables.h) and the HOC errors as quadratic forms, every bound lowered by its proven rounding error, so
   pruning still never excludes a pair the exhaustive search keeps. pair errors and the final scoring stay float.
   about 12 KB, keep it off the stack on small targets. zero it before first use, xombe__cq_prepare keeps the
   maps of L across calls with the same set and L */
typedef struct {
    const xombe_amp_tables *ts;
    const xombe_cq_q15_tab *qt;            /* ts's images, NULL: exhaustive pair search */
    int L;
    int pos[4], J[4], M[4];                /* block k: harmonics pos..pos+J-1, M = min(J, 6) */
    float w[XOMBE_LMAX + 1], sw;           /* weights, sw = sum w */
    float b1[XOMBE_LMAX + 1];              /* Tv's map from (G2..G8) at harmonic l of block k: c0[k] + b1[l] c1[k] */
    float c0[4][7], c1[4][7];
    float mu[7];                           /* weighted column means of that map */
    float Q[7][7];                         /* (A - mean)^T W (A - mean) */
    float H[XOMBE_LMAX + 1][4];            /* HOC synthesis: 2cos(pi m (i + 1/2)/J), m = 2..M-1, block-local i */
    float hq[XOMBE_CQ_HOC_ROWS];           /* every HOC row c: c^T H^T W H c, over its block */
    float hl[XOMBE_CQ_HOC_ROWS];           /* and c*H^T w */
    int hoc_off[4];                        /* block k's first row in hq, hl */
    float e4q[128];                        /* prba58 rows' quadratic terms x4^T Q44 x4 */
    float e4[128];                         /* e4 of the current pair search */
    /* pair-search bounds (prba_pairs), integers v stand for v * 2^-exponent */
    int chol_ok;                           /* 0: Q44 not safely positive definite, or no images: exhaustive */
    float Lc[4][4];                        /* Q44 = Lc Lc^T */
    float Mz[4][3];                        /* Lc^-1 Q43 */
    float dir[4], pmax;                    /* principal direction of the rows p = Lc^T x4, max |p| */
    float ps[128];                         /* usable rows' projections dir*p, ascending */
    short pord[128];                       /* their prba58 rows */
    int n4;                                /* their count */
    int col_ok;                            /* the column bound (roles swapped) is available */
    float L3[3][3], M3[3][4];              /* Q33 = L3 L3^T, M3 = L3^-1 Q34 */
    int32_t sq4[128];                      /* per prba58 row: x4^T(Q44 - M3^T M3)x4 at 2^-sq4_x, |.| <= 2^29 */
    int sq4_x;
    float sq4_d, mw4, n4max;               /* its error bound, rounding magnitude, max |x4| */
    float colb[128];                       /* per prba58 row: bound S' - slack in the current search */
    float fps[128];                        /* columns left after the first row, in projection order */
    short ford[128];
    int32_t sq[512];                       /* per prba24 row: x3^T(Q33 - Mz^T Mz)x3 at 2^-sq_x, |.| <= 2^29 */
    int sq_x;
    float sq_d, mw, n3max;                 /* its error bound, rounding magnitude, max |x3| */
    int32_t sb[512];                       /* per prba24 row: bound S - slack = sb * sb_s - sb_k, INT32_MAX unusable */
    float sb_s, sb_k;
    int32_t dz[512];                       /* per prba24 row: dir*(Mz x3) = dz * dz_s within dz_e */
    float dz_s, dz_e;
    int brute;                             /* 1: no pruning in the pair search (testing) */
    int cap;                               /* pair-search cap in units, 0: none (set after zeroing) */
    int work;                              /* most units a pair search took since xombe__cq_prepare */
} xombe_cq_ws;
#else
/* one search's workspace: per-L maps, w-only terms, and HOC candidates, filled by xombe__cq_prepare
   about 19 KB, keep it off the stack on small targets. zero it before first use, xombe__cq_prepare keeps the
   maps of L across calls with the same set and L */
typedef struct {
    const xombe_amp_tables *ts;
    int L;
    int pos[4], J[4], M[4];
    float w[XOMBE_LMAX + 1], sw;
    float A[XOMBE_LMAX + 1][7];            /* Tv's linear map from (G2..G8) = (prba24, prba58) */
    float Wm[XOMBE_LMAX + 1][7];           /* w*(A - weighted column mean) */
    float Q[7][7];
    float H[XOMBE_LMAX + 1][4];
    float hv[XOMBE_CQ_HOC_ROWS * XOMBE_CQ_J_MAX]; /* every HOC row's block values (H_k @ hoc_k[n]) */
    int hv_off[4];
    float e4q[128];
    float e4[128];
    int chol_ok;
    float Lc[4][4];
    float Mz[4][3];
    float dir[4], pmax;
    float ps[128];
    short pord[128];
    int n4;
    int col_ok;
    float L3[3][3], M3[3][4];
    float sq4[128], mw4, n4max;
    float colb[128];
    float fps[128];
    short ford[128];
    float sq[512], mw, n3max;
    float sb[512];
    float dz[512];
    int brute;
    int cap;
    int work;
} xombe_cq_ws;
#endif

typedef struct {
    int b[9];                   /* b[2..8] */
    float err;                  /* sum w (a - clip(target))^2, a the search's clipped reconstruction */
    float u[XOMBE_LMAX + 1];    /* winner's unclipped reconstruction */
} xombe_cq_result;

/* w[l] = max(2^(gamma (t[l] - max t)), floor), l = 1..L */
void xombe__cq_amp_weights(const float *target, int L, float gamma, float floor_, float *w);

/* fills the workspace for set ts, L harmonics, and weights w[1..L] (NULL: unweighted) */
void xombe__cq_prepare(xombe_cq_ws *ws, const xombe_amp_tables *ts, int L, const float *w);

/* closed-loop search on a prepared workspace: residual r = target - rho(Pt - mean Pt), clipped target,
   base = rho(Pt - mean Pt) - 1/2 log2 L, previous gain Gprev
   n_pairs <= XOMBE_CQ_MAX_PAIRS, n_gain <= XOMBE_CQ_MAX_GAIN */
void xombe__cq_closed(xombe_cq_ws *ws, const float *r, const float *target, const float *base, float Gprev,
                      int n_pairs, int n_gain, int iters, xombe_cq_result *res);

/* quantizes target[1..L] (clipped here to [-15, XOMBE_CQ_AMAX]) on set ts against prediction pt[1..L]
   (xombe__predict_from of the memory) and its gain Gprev. closed 0 is open loop, res->err and res->u unset */
void xombe__cq_quantize(xombe_cq_ws *ws, const xombe_amp_tables *ts, const float *target, int L, const float *pt,
                        float Gprev, const float *w, int closed, int n_pairs, int n_gain, int iters, int clip_passes,
                        xombe_cq_result *res);

#if XOMBE_HAVE_ENCODER
/* closed-loop search on the AMBE+2 set against predictor pi vs the encoder's open-loop choice (xombe__quantize),
   closed wins unless open's weighted error against the decoder's reconstruction is smaller, ties go closed
   target clipped to [-15, 15 - 2^-11]. returns 1 if the open-loop choice won, res->b holds the pick, res->err
   its error (the closed search's own for a closed pick), bclosed/bopen (b[2..8]) may be NULL */
int xombe__cq_search_ambe2(xombe_cq_ws *ws, const xombe_pred *pi, const float *target, int L, const float *w,
                           int n_pairs, int n_gain, xombe_cq_result *res, int *bclosed, int *bopen);
#endif

#endif
