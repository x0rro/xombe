#ifndef XOMBE_XC_STATE_H
#define XOMBE_XC_STATE_H

/* transcoder state behind xombe.h's xombe_xc handle (xc.c), the DMR -> D-STAR half (xcode.c) and the
   D-STAR -> DMR half (xrev.c) */

#include "dstar_internal.h"
#include "dsd_state.h"
#include "../common/cquant.h"

struct xombe_xc_trace;
struct xombe_xr_trace;

/* DMR -> D-STAR */
typedef struct {
    xombe_dec_core d;         /* DMR frame logic and sets */
    xombe_fec fec;
    xombe_ds_mirror mirror;   /* D-STAR decoder's prediction memory */
    uint8_t last[XOMBE_DS_FRAME_BYTES];
    int have_last_b, last_b[9]; /* fields of last, if it was a voice frame */
    int out;                  /* enum xombe_xc_out of the last frame */
    struct xombe_xc_trace *trace; /* trace.h */
    xombe_cq_ws *ws;          /* quantizer workspace, own or shared, nothing in it lives between frames */
} xombe_xc_fwd;

enum { XOMBE_XR_IN_VOICE, XOMBE_XR_IN_TONE, XOMBE_XR_IN_NULL, XOMBE_XR_IN_REPEAT, XOMBE_XR_IN_MUTE };

/* D-STAR -> DMR */
typedef struct {
    xombe_dec_core d;         /* D-STAR frame logic and sets (with s) */
    xombe_dsd_logic s;
    xombe_fec fec;
    xombe_pred pi;            /* DMR decoder's prediction memory */
    int prev_b1;              /* previous D-STAR voice frame's b1, -1 after any other frame but a repeat */
    int have_last, last_db[9]; /* fields of the last good D-STAR voice frame */
    int have_tone, tone_id, tone_ad;
    int out;
    struct xombe_xr_trace *trace;
    xombe_cq_ws *ws;
} xombe_xc_rev;

struct xombe_xc {
    int from; /* enum xombe_codec */
    int (*frame)(struct xombe_xc *x, const uint8_t in[9], unsigned flags, uint8_t out[9]);
    union {
        xombe_xc_fwd f;
        xombe_xc_rev r;
    } u;
};

/* xombe_xc_init's layout: the state, then its own workspace */
struct xombe_xc_own {
    struct xombe_xc x;
    xombe_cq_ws ws;
};

/* each half's init (ws zeroed here), frame, and output class */
void xombe__xc_fwd_init(xombe_xc_fwd *x, xombe_cq_ws *ws);
int xombe__xc_fwd_frame(xombe_xc_fwd *x, const uint8_t in[9], unsigned flags, uint8_t out[9]);
#if XOMBE_HAVE_ENCODER
void xombe__xc_rev_init(xombe_xc_rev *x, xombe_cq_ws *ws);
int xombe__xc_rev_frame(xombe_xc_rev *x, const uint8_t in[9], unsigned flags, uint8_t out[9]);
#endif

#endif
