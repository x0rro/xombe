#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../src/internal.h"
#include "../src/dstar/xc_state.h"
#include "../src/dstar/trace.h"
#include "dstar_ref.h"

/* the D-STAR -> rate-33 transcoder's records (stage 'xrev'), exact quantizer search (ws->cap 0) unless --capped
   keeps the transcoder's XOMBE_CQ_CAP */

static const char *cls_name(int c)
{
    return c == XOMBE_CLASS_VOICE ? "voice" : c == XOMBE_CLASS_SILENCE ? "silence" : c == XOMBE_CLASS_TONE ? "tone" : "?";
}

int main(int argc, char **argv)
{
    static const char *const in_cls[] = {"voice", "tone", "null", "repeat", "mute"};
    static const char *const ocls[] = {"voice", "pulse", "tone", "knox", "silence", "resend"}; /* enum xombe_xc_out */
    static xombe_xc_mem mem;
    long n, nf;
    uint8_t *fr;
    dsref_writer w;
    xombe_xr_trace tr;
    xombe_xc *x;
    int ok, capped = argc == 4 && !strcmp(argv[1], "--capped");
    if (argc != 3 + capped) {
        fprintf(stderr, "usage: dstar_dump_xrev [--capped] IN.dstar OUT\n");
        return 2;
    }
    argv += capped;
    if (!(fr = dsref_read_frames(argv[1], 9, &nf)))
        return 1;
    ok = dsref_wopen(&w, argv[2], "xrev",
                     "n:i in:x ea:i eb:i cls:s scls:s db:i* f0:f L:i a:f* ocls:s tid:i AD:i tgt:f* b:i* d49:b out:x",
                     7);
    if (!ok) {
        free(fr);
        return 1;
    }
    x = xombe_xc_init(&mem, XOMBE_DSTAR);
    x->u.r.trace = &tr;
    if (!capped)
        x->u.r.ws->cap = 0;
    for (n = 0; n < nf && ok; n++) {
        uint8_t out[9];
        xombe_xc_decision dc;
        xombe_xc_frame(x, fr + 9 * n, 0, out);
        xombe_xc_decision_get(x, &dc);
        dsref_wi(&w, n);
        dsref_wx(&w, fr + 9 * n, 9);
        dsref_wi(&w, dc.in.e0);
        dsref_wi(&w, dc.in.e1);
        dsref_ws(&w, in_cls[tr.in_cls]);
        dsref_ws(&w, cls_name(tr.scls));
        dsref_wiv_int(&w, tr.db, tr.has_db ? 9 : 0);
        dsref_wf(&w, tr.has_par ? tr.f0 : NAN);
        dsref_wi(&w, tr.has_par ? tr.L : 0);
        dsref_wfv_float(&w, tr.a + 1, tr.has_par ? tr.L : 0);
        dsref_ws(&w, ocls[dc.out]);
        dsref_wi(&w, tr.tid);
        dsref_wi(&w, tr.AD);
        dsref_wfv_float(&w, tr.tgt + 1, tr.tL);
        dsref_wiv_int(&w, tr.b, tr.has_b ? 9 : 0);
        dsref_wb(&w, tr.d49, 49);
        dsref_wx(&w, out, 9);
        dsref_wend(&w);
    }
    free(fr);
    ok &= dsref_wclose(&w);
    return ok ? 0 : 1;
}
