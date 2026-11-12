#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../src/internal.h"
#include "../src/dstar/xc_state.h"
#include "../src/dstar/trace.h"
#include "dstar_ref.h"

/* the rate-33 -> D-STAR transcoder's records (stage 'xc') */

static const char *cls_name(int c)
{
    return c == XOMBE_CLASS_VOICE ? "voice" : c == XOMBE_CLASS_SILENCE ? "silence" : c == XOMBE_CLASS_TONE ? "tone" : "?";
}

static void bits8(dsref_writer *w, const uint8_t *v, int on)
{
    dsref_wb(w, v, on ? 8 : 0);
}

int main(int argc, char **argv)
{
    long n, nf;
    uint8_t *fr;
    dsref_writer w;
    xombe_xc_mem mem;
    xombe_xc_trace tr;
    xombe_xc *x;
    int ok;
    if (argc != 3) {
        fprintf(stderr, "usage: dstar_dump_xc IN.r33 OUT\n");
        return 2;
    }
    if (!(fr = dsref_read_frames(argv[1], 9, &nf)))
        return 1;
    ok = dsref_wopen(&w, argv[2], "xc",
                     "n:i in:x e0:i e1:i cls:s cls_out:s rb:i* f0:f L:i a:f* V2:b V1:b G:f pulse:i tid:i AD:i "
                     "soft:f* tgt:f* tb0:i b:i* out:x", 7);
    if (!ok) {
        free(fr);
        return 1;
    }
    x = xombe_xc_init(&mem, XOMBE_DMR);
    x->u.f.trace = &tr;
    for (n = 0; n < nf && ok; n++) {
        uint8_t out[9];
        xombe_xc_decision dc;
        xombe_xc_frame(x, fr + 9 * n, 0, out);
        xombe_xc_decision_get(x, &dc);
        dsref_wi(&w, n);
        dsref_wx(&w, fr + 9 * n, 9);
        dsref_wi(&w, dc.in.e0);
        dsref_wi(&w, dc.in.e1);
        dsref_ws(&w, tr.conceal == XOMBE_CONCEAL_REPEAT ? "repeat"
                     : tr.conceal == XOMBE_CONCEAL_MUTE ? "mute" : cls_name(tr.cls));
        dsref_ws(&w, cls_name(tr.cls_out));
        dsref_wiv_int(&w, tr.rb, tr.has_rb ? 9 : 0);
        dsref_wf(&w, tr.has_par ? tr.f0 : NAN);
        dsref_wi(&w, tr.has_par ? tr.L : 0);
        dsref_wfv_float(&w, tr.a + 1, tr.has_par ? tr.L : 0);
        bits8(&w, tr.V2, tr.has_par);
        bits8(&w, tr.V1, tr.has_v1);
        dsref_wf(&w, tr.has_v1 ? tr.G : NAN);
        dsref_wi(&w, tr.pulse);
        dsref_wi(&w, tr.tid);
        dsref_wi(&w, tr.AD);
        dsref_wfv_float(&w, tr.soft, tr.has_soft ? 4 : 0);
        dsref_wfv_float(&w, tr.tgt + 1, tr.tb0 >= 0 ? tr.tL : 0);
        dsref_wi(&w, tr.tb0);
        dsref_wiv_int(&w, tr.b, tr.has_b ? 9 : 0);
        dsref_wx(&w, out, 9);
        dsref_wend(&w);
    }
    free(fr);
    ok &= dsref_wclose(&w);
    return ok ? 0 : 1;
}
