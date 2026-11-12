#include <string.h>
#include "../internal.h"
#include "xc_state.h"

/* xombe_xc: the DMR -> D-STAR half (xcode.c) or the D-STAR -> DMR half (xrev.c) by the input codec, each reached
   only through its own init so a program using one leaves the other out */

XOMBE_STATIC_ASSERT(sizeof(struct xombe_xc_own) <= XOMBE_XC_SIZE, xc_state_fits_XOMBE_XC_SIZE);
XOMBE_STATIC_ASSERT(sizeof(struct xombe_xc) <= XOMBE_XC_SHARED_SIZE, xc_state_fits_XOMBE_XC_SHARED_SIZE);
XOMBE_STATIC_ASSERT(sizeof(xombe_cq_ws) <= XOMBE_QWS_SIZE, cq_ws_fits_XOMBE_QWS_SIZE);

size_t xombe_xc_size(void)
{
    return sizeof(struct xombe_xc_own);
}

static int fwd_frame(xombe_xc *x, const uint8_t in[9], unsigned flags, uint8_t out[9])
{
    return xombe__xc_fwd_frame(&x->u.f, in, flags, out);
}

xombe_xc *xombe__xc_init_shared_dmr(void *mem, void *ws)
{
    xombe_xc *x = mem;
    if (!x || !ws)
        return NULL;
    xombe__xc_fwd_init(&x->u.f, ws);
    x->from = XOMBE_DMR;
    x->frame = fwd_frame;
    return x;
}

#if XOMBE_HAVE_ENCODER
static int rev_frame(xombe_xc *x, const uint8_t in[9], unsigned flags, uint8_t out[9])
{
    return xombe__xc_rev_frame(&x->u.r, in, flags, out);
}
#endif

xombe_xc *xombe__xc_init_shared_dstar(void *mem, void *ws)
{
#if XOMBE_HAVE_ENCODER
    xombe_xc *x = mem;
    if (!x || !ws)
        return NULL;
    xombe__xc_rev_init(&x->u.r, ws);
    x->from = XOMBE_DSTAR;
    x->frame = rev_frame;
    return x;
#else
    (void)mem;
    (void)ws;
    return NULL;
#endif
}

xombe_xc *xombe__xc_init_dmr(void *mem)
{
    return mem ? xombe__xc_init_shared_dmr(mem, &((struct xombe_xc_own *)mem)->ws) : NULL;
}

xombe_xc *xombe__xc_init_dstar(void *mem)
{
    return mem ? xombe__xc_init_shared_dstar(mem, &((struct xombe_xc_own *)mem)->ws) : NULL;
}

int xombe_xc_frame(xombe_xc *x, const uint8_t in[XOMBE_FRAME_BYTES], unsigned flags, uint8_t out[XOMBE_FRAME_BYTES])
{
    return x->frame(x, in, flags, out);
}

void xombe_xc_decision_get(const xombe_xc *x, xombe_xc_decision *out)
{
    if (x->from == XOMBE_DSTAR) {
        out->in = x->u.r.d.decision;
        out->out = (enum xombe_xc_out)x->u.r.out;
    } else {
        out->in = x->u.f.d.decision;
        out->out = (enum xombe_xc_out)x->u.f.out;
    }
}
