#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../src/internal.h"
#include "../src/dstar/dsd_state.h"
#include "dstar_ref.h"

/* D-STAR decoder records: --dec (default) frame level (stage 'dec'), --half half-frame sets (stage 'half'),
   --pcm the PCM */

static const char *const LOG[] = {"voice", "tone", "null", "repeat", "mute"};

static const char *cls_name(int c)
{
    return c == XOMBE_CLASS_VOICE ? "voice" : c == XOMBE_CLASS_SILENCE ? "silence" : c == XOMBE_CLASS_TONE ? "tone" : "?";
}

static void write_half(dsref_writer *w, long n, int h, const xombe_set *s, const int16_t *pcm, int N)
{
    char V[9], vl[XOMBE_LMAX + 2];
    uint8_t vh[XOMBE_LMAX + 1];
    float a[XOMBE_LMAX + 1], M[XOMBE_LMAX + 1];
    double e2 = 0.0;
    int l, tone = s->cls == XOMBE_CLASS_TONE;
    for (l = 0; l < 8; l++)
        V[l] = tone || s->vb[l] ? 'V' : 'U';
    V[8] = 0;
    memset(vh, 0, sizeof(vh));
    if (!tone)
        xombe__ds_harm_voicing(s->vb, s->f0, s->L, vh);
    for (l = 1; l <= s->L; l++) {
        if (tone) {
            int on = s->e[l] > XOMBE_AMP_OFF;
            M[l] = on ? exp2f(s->e[l]) : 0.0f;
            a[l] = on ? s->e[l] : (float)log2(1e-30);
            vh[l] = (uint8_t)on;
        } else {
            M[l] = exp2f(s->e[l]);
            a[l] = s->a[l];
        }
        vl[l - 1] = vh[l] ? 'V' : 'U';
    }
    vl[s->L] = 0;
    for (l = 0; l < N; l++)
        e2 += (double)pcm[l] * pcm[l];
    dsref_wi(w, n);
    dsref_wi(w, h);
    dsref_ws(w, cls_name(s->cls));
    dsref_wi(w, s->repeat);
    dsref_wf(w, s->f0);
    dsref_wi(w, s->L);
    dsref_ws(w, V);
    dsref_ws(w, vl);
    dsref_wfv_float(w, a + 1, s->L);
    dsref_wfv_float(w, M + 1, s->L);
    dsref_wf(w, sqrt(e2 / N));
    dsref_wend(w);
}

int main(int argc, char **argv)
{
    long n, nf;
    uint8_t *fr;
    dsref_writer w;
    xombe_dec_mem mem;
    xombe_dec *s;
    FILE *pcm_out = NULL;
    int mode = 0, ok = 1, a = 1;
    if (argc == 4 && !strcmp(argv[1], "--half"))
        mode = 1, a = 2;
    else if (argc == 4 && !strcmp(argv[1], "--pcm"))
        mode = 2, a = 2;
    else if (argc == 4 && !strcmp(argv[1], "--dec"))
        a = 2;
    if (argc != a + 2) {
        fprintf(stderr, "usage: dstar_dump_dec [--dec | --half | --pcm] IN.dstar OUT\n");
        return 2;
    }
    if (!(fr = dsref_read_frames(argv[a], 9, &nf)))
        return 1;
    if (mode == 2)
        ok = (pcm_out = fopen(argv[a + 1], "wb")) != NULL;
    else if (mode == 1)
        ok = dsref_wopen(&w, argv[a + 1], "half", "n:i h:i cls:s rep:i f0:f L:i V:s vl:s a:f* M:f* rms:f", 0);
    else
        ok = dsref_wopen(&w, argv[a + 1], "dec",
                         "n:i log:s ea:i eb:i bf:i hb:i eT:i e1:i Picls:s PiL:i Pif0:f PiG:f lam:f fbar:f Pia:f*", 0);
    if (!ok) {
        free(fr);
        return 1;
    }
    s = xombe_dec_init(&mem, XOMBE_DSTAR);
    xombe_dec_set_noise(s, XOMBE_NOISE_CONFORMANCE, 0);
    for (n = 0; n < nf && ok; n++) {
        int16_t pcm[160];
        float fbar = s->c.avg_f0 * XOMBE_FS; /* the model records it before the frame's halves */
        xombe_dec_half(s, fr + 9 * n, 0, 0, 80, pcm);
        if (mode == 0) {
            const xombe_pred *pi = &s->c.pi;
            dsref_wi(&w, n);
            dsref_ws(&w, LOG[s->ds.log]);
            dsref_wi(&w, s->c.decision.e0);
            dsref_wi(&w, s->c.decision.e1);
            dsref_wi(&w, s->c.decision.bf);
            dsref_wi(&w, s->c.decision.hb);
            dsref_wi(&w, s->ds.eT);
            dsref_wi(&w, s->ds.e1);
            dsref_ws(&w, cls_name(pi->cls));
            dsref_wi(&w, pi->L);
            dsref_wf(&w, pi->f0);
            dsref_wf(&w, pi->gamma);
            dsref_wf(&w, s->c.lambda);
            dsref_wf(&w, fbar);
            dsref_wfv_float(&w, pi->a + 1, pi->L);
            dsref_wend(&w);
        } else if (mode == 1) {
            write_half(&w, n, 1, &s->c.p1, pcm, 80);
        }
        xombe_dec_half(s, fr + 9 * n, 0, 1, 80, pcm + 80);
        if (mode == 1)
            write_half(&w, n, 2, &s->c.p2, pcm + 80, 80);
        if (mode == 2)
            ok = fwrite(pcm, 2, 160, pcm_out) == 160;
    }
    free(fr);
    if (mode == 2)
        ok &= fclose(pcm_out) == 0;
    else
        ok &= dsref_wclose(&w);
    return ok ? 0 : 1;
}
