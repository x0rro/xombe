#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../src/dstar/dstar_internal.h"
#include "dstar_ref.h"

/* D-STAR dequantization against the mirror (stage 'deq') */

int main(int argc, char **argv)
{
    long n, nf;
    uint8_t *fr;
    dsref_writer w;
    xombe_ds_mirror m;
    int k, ok;
    if (argc != 3) {
        fprintf(stderr, "usage: dstar_dump_deq IN.dstar OUT\n");
        return 2;
    }
    if (!(fr = dsref_read_frames(argv[1], 9, &nf)))
        return 1;
    if (!dsref_wopen(&w, argv[2], "deq",
                     "n:i voice:i b0:i b1:i b2:i b3:i b4:i b5:i b6:i b7:i b8:i f0:f L:i vl:b Tv:f* G:f stL:i "
                     "stG:f sta:f*", 0)) {
        free(fr);
        return 1;
    }
    xombe__ds_mirror_init(&m);
    for (n = 0; n < nf; n++) {
        xombe_ds_frame f;
        xombe_ds_params p;
        int voice;
        xombe__ds_fec_decode(fr + 9 * n, &f);
        voice = xombe__ds_params_from_frame(&m, &f, &p);
        dsref_wi(&w, n);
        dsref_wi(&w, voice);
        for (k = 0; k < 9; k++)
            dsref_wi(&w, f.b[k]);
        if (voice) {
            dsref_wf(&w, p.f0);
            dsref_wi(&w, p.L);
            dsref_wb(&w, p.vl + 1, p.L);
            dsref_wfv_float(&w, p.Tv + 1, p.L);
            dsref_wf(&w, p.G);
        } else {
            dsref_wf(&w, NAN);
            dsref_wi(&w, 0);
            dsref_wb(&w, NULL, 0);
            dsref_wfv_float(&w, NULL, 0);
            dsref_wf(&w, NAN);
        }
        dsref_wi(&w, m.L);
        dsref_wf(&w, m.G);
        dsref_wfv_float(&w, m.a + 1, m.L);
        dsref_wend(&w);
    }
    free(fr);
    ok = dsref_wclose(&w);
    return ok ? 0 : 1;
}
