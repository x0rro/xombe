#ifndef XOMBE_DSTAR_TRACE_H
#define XOMBE_DSTAR_TRACE_H

/* transcoder per-frame records for tests, filled when the state's trace pointer is set */

#include "../state.h"

/* DMR -> D-STAR half of xombe_xc (xcode.c), dumped by test/dstar_dump_xc.c */
typedef struct xombe_xc_trace {
    int cls;              /* AMBE+2 class, XOMBE_CLASS_ERASURE for a concealed frame */
    int conceal;          /* enum xombe_conceal */
    int cls_out;          /* class of the set transcoded (the played set when concealed) */
    int has_rb, rb[9];    /* AMBE+2 fields, frames decoded normally */
    int has_par;          /* voice and silence sets: f0, L, a, V2 */
    float f0;             /* cycles/sample */
    int L;
    float a[XOMBE_LMAX + 1];
    uint8_t V2[8], V1[8];
    int has_v1;           /* frames decoded normally: V1 and G */
    float G;
    int pulse, tid, AD;   /* tid, AD -1 when none */
    int has_soft;
    float soft[4];        /* per-kHz voiced fraction */
    int tL, tb0;          /* quantizer target tgt[1..tL] and its pitch, tb0 -1 when not quantized */
    float tgt[XOMBE_LMAX + 1];
    int has_b, b[9];      /* D-STAR fields of the output (of the resent frame for a resend) */
} xombe_xc_trace;

/* D-STAR -> DMR half of xombe_xc (xrev.c), dumped by test/dstar_dump_xrev.c */
typedef struct xombe_xr_trace {
    int in_cls;           /* XOMBE_XR_IN_* (xc_state.h) */
    int scls;             /* class of the D-STAR decoder's played set */
    int has_db, db[9];    /* D-STAR fields, frames not lost */
    int has_par;          /* voice and silence sets: f0, L, a */
    float f0;
    int L;
    float a[XOMBE_LMAX + 1];
    int tid, AD;          /* -1 when none */
    int tL;               /* quantizer target tgt[1..tL], 0 when not quantized */
    float tgt[XOMBE_LMAX + 1];
    int has_b, b[9];      /* AMBE+2 fields of voice and pulse frames */
    uint8_t d49[49];
} xombe_xr_trace;

#endif
