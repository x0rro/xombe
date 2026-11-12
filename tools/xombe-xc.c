#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "io.h"
#include "xombe.h"

/* DMR rate-33 frames -> D-STAR frames, or back with --reverse (xombe_xc_frame)
   --amb: the DMR side is .amb, through the FEC (exact on clean frames) */

#if XOMBE_HAVE_ENCODER
#define NO_REVERSE ""
#else
#define NO_REVERSE " (not in this build)"
#endif

static void usage(FILE *f, int status)
{
    fprintf(f,
            "usage: xombe-xc [options] IN.r33 OUT.dstar\n"
            "       xombe-xc --reverse [options] IN.dstar OUT.r33\n"
            "modes:\n"
            "  (default)    DMR rate 33 -> D-STAR\n"
            "  --reverse    D-STAR -> DMR rate 33" NO_REVERSE "\n"
            "  r33: DMR rate-33 frames, 9 bytes (xombe-enc --air, xombe-dec --air)\n"
            "  dstar: D-STAR frames, 9 bytes (xombe-enc --dstar, xombe-dec --dstar). '-' is stdin or stdout\n"
            "DMR format:\n"
            "  --amb        .amb in place of r33, either direction (xombe-enc, xombe-dec)\n"
            "options:\n"
            "  --lost LIST  frames lost, e.g. 3,10-12 (XOMBE_FLAG_LOST)\n"
            "  --decisions  per-frame decisions to stderr: 'n class conceal run mute e0 e1 bf hb out', as\n"
            "               xombe-dec's of the input, then the output kind\n"
            "  --help       this text\n"
            "  no delay: one output frame per input frame. a truncated final frame is an error\n");
    exit(status);
}

int main(int argc, char **argv)
{
    static const char *const outk[] = {"voice", "pulse", "tone", "knox", "silence", "resend"};
    static xombe_xc_mem mem;
    const char *paths[2];
    xombe_xc *x = NULL;
    io_frames lost = {NULL, 0};
    FILE *in, *out;
    long n = 0;
    xombe_fec fec;
    int i, np = 0, show = 0, rev = 0, amb = 0, ok = 1, got;
    io_prog = "xombe-xc";
    for (i = 1; i < argc; i++) {
        const char *a = argv[i];
        if (!strcmp(a, "--help"))
            usage(stdout, 0);
        else if (!strcmp(a, "--decisions"))
            show = 1;
        else if (!strcmp(a, "--amb"))
            amb = 1;
        else if (!strcmp(a, "--reverse")) {
            if (!XOMBE_HAVE_ENCODER)
                io_fail(": not in this build (XOMBE_ENCODER=OFF)", a);
            rev = 1;
        } else if (!strcmp(a, "--lost") && i + 1 < argc) {
            if (!io_frames_parse(argv[++i], &lost))
                io_fail(": not a frame list (e.g. 3,10-12)", argv[i]);
        } else if (np < 2 && (a[0] != '-' || !a[1]))
            paths[np++] = a;
        else
            usage(stderr, 2);
    }
    if (np != 2)
        usage(stderr, 2);
    if (!(in = io_open(paths[0], "rb")))
        return 1;
    if (amb && !rev) {
        char hdr[4];
        if (fread(hdr, 1, 4, in) != 4 || memcmp(hdr, IO_AMB_HEADER, 4)) {
            fprintf(stderr, "%s: not an .amb file\n", paths[0]);
            if (in != stdin)
                fclose(in);
            return 1;
        }
    }
    if (!(out = io_open(paths[1], "wb"))) {
        if (in != stdin)
            fclose(in);
        return 1;
    }
    if (amb && rev && !io_write(out, IO_AMB_HEADER, 4, paths[1]))
        ok = 0;
    xombe_fec_init(&fec, XOMBE_DMR);
    x = xombe_xc_init(&mem, rev ? XOMBE_DSTAR : XOMBE_DMR);
    while (ok) {
        uint8_t fr[XOMBE_FRAME_BYTES], ds[XOMBE_FRAME_BYTES], bits[XOMBE_BITS_BYTES], rec[IO_AMB_RECORD];
        unsigned flags = io_frames_has(&lost, n) ? XOMBE_FLAG_LOST : 0;
        if (amb && !rev) {
            if ((got = io_read_record(in, rec, IO_AMB_RECORD, paths[0])) == 1) {
                io_amb_to_bytes(rec, bits);
                xombe_fec_encode(XOMBE_DMR, bits, fr);
            }
        } else
            got = io_read_record(in, fr, XOMBE_FRAME_BYTES, paths[0]);
        if (got != 1) {
            ok = got == 0;
            break;
        }
        xombe_xc_frame(x, fr, flags, ds);
        if (show) {
            xombe_xc_decision d;
            xombe_xc_decision_get(x, &d);
            io_print_decision(n, &d.in, outk[d.out]);
        }
        if (amb && rev) {
            xombe_fec_stats st;
            xombe_fec_decode(&fec, ds, bits, &st);
            io_amb_from_bytes(bits, rec);
            ok = io_write(out, rec, IO_AMB_RECORD, paths[1]);
        } else
            ok = io_write(out, ds, XOMBE_FRAME_BYTES, paths[1]);
        if (!ok)
            break;
        n++;
    }
    if (ferror(in)) {
        fprintf(stderr, "%s: read error\n", paths[0]);
        ok = 0;
    }
    if (in != stdin)
        fclose(in);
    ok &= io_close(out, paths[1]);
    free(lost.r);
    return ok ? 0 : 1;
}
