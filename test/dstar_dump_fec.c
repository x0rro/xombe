#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../src/dstar/dstar_internal.h"
#include "dstar_ref.h"

/* D-STAR FEC and fields (stage 'fec'): every frame of a stream with the error rule 'ds', --enc re-encodes the
   fields and tone cases of a fec_enc file, --ids prints "n class toneid" per frame, toneid -1 if none (conformance) */

static const char *class_name(const xombe_ds_frame *f)
{
    switch (f->cls) {
    case XOMBE_DS_VOICE: return "voice";
    case XOMBE_DS_TONE: return "tone";
    default: return f->b[0] == 124 ? "null124" : f->b[0] == 125 ? "null125" : "null127";
    }
}

static int dump_stream(const char *in, const char *out)
{
    long n, nf;
    uint8_t *fr = dsref_read_frames(in, 9, &nf), enc[9];
    dsref_writer w;
    xombe_ds_err er;
    int k, ok;
    if (!fr)
        return 1;
    if (!dsref_wopen(&w, out, "fec",
                     "n:i in:x A:i B:i C:i pn:i da:i ea:i db:i eb:i d48:b b0:i b1:i b2:i b3:i b4:i b5:i b6:i b7:i "
                     "b8:i cls:s T:i AD:i bf:i hb:i eT:i e1:i enc:x", 0)) {
        free(fr);
        return 1;
    }
    xombe__ds_err_init(&er);
    for (n = 0; n < nf; n++) {
        xombe_ds_frame f;
        xombe_ds_flags fl;
        xombe__ds_fec_decode(fr + 9 * n, &f);
        xombe__ds_err_update(&er, f.ea, f.eb, &fl);
        xombe__ds_encode_fields(f.b, enc);
        dsref_wi(&w, n);
        dsref_wx(&w, fr + 9 * n, 9);
        dsref_wi(&w, (long)f.A);
        dsref_wi(&w, (long)f.B);
        dsref_wi(&w, (long)f.C);
        dsref_wi(&w, (long)f.pn);
        dsref_wi(&w, f.da);
        dsref_wi(&w, f.ea);
        dsref_wi(&w, f.db);
        dsref_wi(&w, f.eb);
        dsref_wb(&w, f.d, 48);
        for (k = 0; k < 9; k++)
            dsref_wi(&w, f.b[k]);
        dsref_ws(&w, class_name(&f));
        dsref_wi(&w, f.T);
        dsref_wi(&w, f.AD);
        dsref_wi(&w, fl.bf);
        dsref_wi(&w, fl.hb);
        dsref_wi(&w, fl.eT);
        dsref_wi(&w, fl.e1);
        dsref_wx(&w, enc, 9);
        dsref_wend(&w);
    }
    free(fr);
    ok = dsref_wclose(&w);
    return ok ? 0 : 1;
}

static int dump_enc(const char *in, const char *out)
{
    dsref r;
    dsref_writer w;
    int c_kind, c_b, c_T, c_level, st, k, ok;
    if (!dsref_open(&r, in))
        return 1;
    c_kind = dsref_need(&r, "kind", DSREF_STR);
    c_b = dsref_need(&r, "b", DSREF_IVEC);
    c_T = dsref_need(&r, "T", DSREF_INT);
    c_level = dsref_need(&r, "level", DSREF_FLOAT);
    if (c_kind < 0 || c_b < 0 || c_T < 0 || c_level < 0 ||
        !dsref_wopen(&w, out, "fec_enc", "kind:s b:i* T:i level:f d48:b out:x", 0)) {
        dsref_close(&r);
        return 1;
    }
    while ((st = dsref_next(&r)) > 0) {
        uint8_t d[48], o[9];
        int b[9];
        if (!strcmp(dsref_s(&r, c_kind), "tone")) {
            int T = (int)dsref_i(&r, c_T);
            xombe__ds_tone_bits(T, xombe__ds_tone_q((float)dsref_f(&r, c_level)), d);
            xombe__ds_tone_frame(T, xombe__ds_tone_q((float)dsref_f(&r, c_level)), o);
            xombe__ds_get_fields(d, b);
        } else {
            if (dsref_len(&r, c_b) != 9) {
                fprintf(stderr, "%s: b needs 9 fields\n", in);
                st = -1;
                break;
            }
            for (k = 0; k < 9; k++)
                b[k] = (int)dsref_iv(&r, c_b)[k];
            xombe__ds_put_fields(b, d);
            if (xombe__ds_encode_fields(b, o) != 0) {
                fprintf(stderr, "%s: field out of range\n", in);
                st = -1;
                break;
            }
        }
        dsref_ws(&w, dsref_s(&r, c_kind));
        dsref_wiv_int(&w, b, 9);
        dsref_wi(&w, dsref_i(&r, c_T));
        dsref_wf(&w, dsref_f(&r, c_level));
        dsref_wb(&w, d, 48);
        dsref_wx(&w, o, 9);
        dsref_wend(&w);
    }
    ok = dsref_close(&r) && st == 0;
    ok &= dsref_wclose(&w);
    return ok ? 0 : 1;
}

static int dump_ids(const char *in)
{
    long n, nf;
    uint8_t *fr = dsref_read_frames(in, 9, &nf);
    if (!fr)
        return 1;
    for (n = 0; n < nf; n++) {
        xombe_ds_frame f;
        xombe__ds_fec_decode(fr + 9 * n, &f);
        printf("%ld %s %d\n", n, class_name(&f), f.tid);
    }
    free(fr);
    return 0;
}

int main(int argc, char **argv)
{
    if (argc == 3 && !strcmp(argv[1], "--ids"))
        return dump_ids(argv[2]);
    if (argc == 4 && !strcmp(argv[1], "--enc"))
        return dump_enc(argv[2], argv[3]);
    if (argc == 3)
        return dump_stream(argv[1], argv[2]);
    fprintf(stderr, "usage: dstar_dump_fec IN.dstar OUT | --enc fec_enc.txt.gz OUT | --ids IN.dstar\n");
    return 2;
}
