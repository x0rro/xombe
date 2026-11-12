#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../src/state.h"
#include "../src/dstar/dstar_enc.h"
#include "dstar_ref.h"

/* whole encoder with the D-STAR output on PCM, against the model's records (enc stage)
   columns n cls det trk tone b0 b1 out:
       cls   voice, tone, knox, or silence (a NULL frame)
       det, trk, tone   tone pipeline on the frame (T and level, empty when none): detection, after the
             run rule, and the final tone after the tone-in-noise pass
       b0, b1  the frame's first two fields (126 and -1 on tone frames, 124 and -1 on NULL frames)
   --dual: both outputs enabled (the D-STAR output does not depend on it) */

static void wtone(dsref_writer *w, const xombe_ds_tone *t)
{
    float v[2];
    v[0] = (float)t->T;
    v[1] = t->level;
    dsref_wfv_float(w, v, t->on ? 2 : 0);
}

int main(int argc, char **argv)
{
    static const char *cls_name[4] = {"voice", "silence", "tone", "knox"};
    static xombe_enc_mem mem;
    xombe_enc_cfg cfg;
    xombe_enc *e;
    dsref_writer w;
    int16_t pcm[XOMBE_FRAME];
    uint8_t d[9];
    long frames = 0, n = 0;
    int dual = argc > 1 && !strcmp(argv[1], "--dual");
    FILE *in;
    if (argc != 3 + dual) {
        fprintf(stderr, "usage: dstar_dump_enc [--dual] IN.raw OUT\n");
        return 2;
    }
    in = fopen(argv[1 + dual], "rb");
    if (!in) {
        perror(argv[1 + dual]);
        return 1;
    }
    xombe_enc_cfg_default(&cfg);
    cfg.outputs = XOMBE_ENC_DSTAR | (dual ? XOMBE_ENC_DMR : 0);
    e = xombe_enc_init(&mem, &cfg);
    if (!e || !dsref_wopen(&w, argv[2 + dual], "enc", "n:i cls:s det:f* trk:f* tone:f* b0:i b1:i out:x", 7))
        return 1;
    for (;;) {
        const xombe_ds_enc_frame *fr = &XOMBE_DS_ENC(e)->fr;
        size_t k = fread(pcm, sizeof(pcm[0]), XOMBE_FRAME, in);
        int got;
        if (k > 0) { /* a partial last frame is zero-padded, as the model counts frames */
            memset(pcm + k, 0, (XOMBE_FRAME - k) * sizeof(pcm[0]));
            frames++;
            got = xombe_enc_frame(e, pcm, NULL, NULL, d);
        } else if (n < frames)
            got = xombe_enc_frame(e, NULL, NULL, NULL, d);
        else
            break;
        if (!got)
            continue;
        dsref_wi(&w, n++);
        dsref_ws(&w, cls_name[fr->cls]);
        wtone(&w, &fr->det);
        wtone(&w, &fr->trk);
        wtone(&w, &fr->tone);
        dsref_wi(&w, fr->cls == XOMBE_DS_CLS_TONE ? 126 : fr->cls == XOMBE_DS_CLS_NULL ? 124 : fr->b0);
        dsref_wi(&w, fr->cls == XOMBE_DS_CLS_VOICE ? fr->b1 : fr->cls == XOMBE_DS_CLS_KNOX ? 15 : -1);
        dsref_wx(&w, d, 9);
        dsref_wend(&w);
    }
    fclose(in);
    return dsref_wclose(&w) ? 0 : 1;
}
