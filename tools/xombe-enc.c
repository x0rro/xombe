#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "io.h"
#include "xombe.h"

#if XOMBE_HAVE_DSTAR
#define NO_DSTAR ""
#else
#define NO_DSTAR " (not in this build)"
#endif

static void usage(FILE *f, int status)
{
    fprintf(f,
            "usage: xombe-enc [options] IN.raw OUT\n"
            "  IN.raw: 8 kHz s16le mono, a partial final frame is dropped. '-' is stdin or stdout\n"
            "modes:\n"
            "  (default)            AMBE+2: OUT is .amb, or .bits with --bits, or air frames with --air\n"
            "  --dstar              D-STAR only: OUT is 9-byte D-STAR frames (.dstar)" NO_DSTAR "\n"
            "  --dual DSTAR_OUT     both: AMBE+2 to OUT, D-STAR frames to DSTAR_OUT" NO_DSTAR "\n"
            "AMBE+2 output format:\n"
            "  --bits               49 bytes per frame, one bit per byte\n"
            "  --air                72-bit air frames, 9 bytes, FEC applied\n"
            "options:\n"
            "  --ns on|off          noise suppression (default on)\n"
            "  --ns-strength S      suppression strength, 0..1 (default 0.5)\n"
            "  --dtx on|off         AMBE+2 silence frames for background (default on)\n"
            "  --tones on|off       tone detection (default on)\n"
            "  --vscale S           voicing decision scale, 0.001..1000000 (default 0.35)\n"
            "  --dstar-dtx on|off   D-STAR null frames for background (default off)" NO_DSTAR "\n"
            "  --quant closed|open  AMBE+2 amplitude quantizer, open is faster (default closed)\n"
            "  --help               this text\n"
            "  --ns-strength and --vscale round S to thousandths. delay: 3 frames, 4 with D-STAR, flushed at the end\n");
    exit(status);
}

static int on_off(const char *s)
{
    if (!strcmp(s, "on"))
        return 1;
    if (!strcmp(s, "off"))
        return 0;
    io_fail(": expected on or off", s);
    return 0;
}

static int quantizer(const char *s)
{
    if (!strcmp(s, "closed"))
        return XOMBE_ENC_QUANT_CLOSED;
    if (!strcmp(s, "open"))
        return XOMBE_ENC_QUANT_OPEN;
    io_fail(": expected closed or open", s);
    return 0;
}

static float number(const char *s, float lo, float hi)
{
    float v;
    if (!io_parse_float(s, lo, hi, &v)) {
        fprintf(stderr, "%s: '%s' is not a number in %.7g..%.7g\n", io_prog, s, lo, hi);
        exit(2);
    }
    return v;
}

/* the library takes the scales in thousandths */
static int permille(const char *s, float lo, float hi)
{
    return (int)(number(s, lo, hi) * 1000.0f + 0.5f);
}

static int write_frame(FILE *out, enum io_format fmt, const uint8_t bits[XOMBE_BITS_BYTES], const char *path)
{
    uint8_t buf[IO_BITS_RECORD];
    if (fmt == IO_AMB) {
        io_amb_from_bytes(bits, buf);
        return io_write(out, buf, IO_AMB_RECORD, path);
    }
    if (fmt == IO_AIR) {
        xombe_fec_encode(XOMBE_DMR, bits, buf);
        return io_write(out, buf, IO_AIR_RECORD, path);
    }
    io_unpack_bits(bits, buf);
    return io_write(out, buf, IO_BITS_RECORD, path);
}

static int write_dstar(FILE *out, const uint8_t bits[XOMBE_BITS_BYTES], const char *path)
{
    uint8_t fr[XOMBE_FRAME_BYTES];
    xombe_fec_encode(XOMBE_DSTAR, bits, fr);
    return io_write(out, fr, sizeof(fr), path);
}

enum mode { MODE_AMBE2, MODE_DSTAR, MODE_DUAL };

/* one frame through the encoder, pcm NULL flushes, returns 1 if it wrote a frame. the DMR output goes to out in
   fmt, the D-STAR output to out (--dstar) or dout (--dual) */
static int encode_frame(xombe_enc *enc, enum mode mode, const int16_t *pcm, FILE *out, enum io_format fmt,
                        const char *path, FILE *dout, const char *dpath, int *ok)
{
    uint8_t dmr[XOMBE_BITS_BYTES], dstar[XOMBE_BITS_BYTES];
    if (xombe_enc_frame_bits(enc, pcm, NULL, mode != MODE_DSTAR ? dmr : NULL, mode != MODE_AMBE2 ? dstar : NULL) != 1)
        return 0;
    if (mode == MODE_DSTAR)
        *ok = write_dstar(out, dstar, path);
    else
        *ok = write_frame(out, fmt, dmr, path) && (mode != MODE_DUAL || write_dstar(dout, dstar, dpath));
    return 1;
}

int main(int argc, char **argv)
{
    xombe_enc_cfg cfg;
    static xombe_enc_mem mem;
    xombe_enc *enc;
    const char *paths[2], *dpath = NULL;
    enum io_format fmt = IO_AMB;
    enum mode mode = MODE_AMBE2;
    int np = 0, i, ok = 1, have_fmt = 0;
    long frames = 0, written = 0;
    FILE *in, *out, *dout = NULL;
    io_prog = "xombe-enc";
    xombe_enc_cfg_default(&cfg);
    for (i = 1; i < argc; i++) {
        const char *a = argv[i];
        int more = i + 1 < argc;
        if (!strcmp(a, "--help"))
            usage(stdout, 0);
        else if (!strcmp(a, "--air")) {
            fmt = IO_AIR;
            have_fmt = 1;
        } else if (!strcmp(a, "--bits")) {
            fmt = IO_BITS;
            have_fmt = 1;
        } else if (!strcmp(a, "--ns") && more)
            cfg.ns = on_off(argv[++i]);
        else if (!strcmp(a, "--ns-strength") && more)
            cfg.ns_strength_permille = permille(argv[++i], 0.0f, 1.0f);
        else if (!strcmp(a, "--dtx") && more)
            cfg.dmr_dtx = on_off(argv[++i]);
        else if (!strcmp(a, "--tones") && more)
            cfg.tones = on_off(argv[++i]);
        else if (!strcmp(a, "--quant") && more)
            cfg.dmr_quant = quantizer(argv[++i]);
        else if (!strcmp(a, "--vscale") && more)
            cfg.voicing_scale_permille = permille(argv[++i], 0.001f, 1e6f);
        else if (!strcmp(a, "--dstar") || (!strcmp(a, "--dual") && more) || (!strcmp(a, "--dstar-dtx") && more)) {
            if (!XOMBE_HAVE_DSTAR)
                io_fail(": not in this build (XOMBE_DSTAR=OFF)", a);
            if (!strcmp(a, "--dstar"))
                mode = MODE_DSTAR;
            else if (!strcmp(a, "--dual")) {
                mode = MODE_DUAL;
                dpath = argv[++i];
            } else
                cfg.dstar_dtx = on_off(argv[++i]);
        } else if (a[0] == '-' && a[1])
            usage(stderr, 2);
        else if (np < 2)
            paths[np++] = a;
        else
            usage(stderr, 2);
    }
    if (np != 2)
        usage(stderr, 2);
    if (mode == MODE_DSTAR && have_fmt)
        io_fail("--bits and --air select the AMBE+2 output's format, --dstar has none", NULL);
    cfg.outputs = mode == MODE_AMBE2 ? XOMBE_ENC_DMR
                  : mode == MODE_DSTAR ? XOMBE_ENC_DSTAR
                                       : XOMBE_ENC_DMR | XOMBE_ENC_DSTAR;
    enc = xombe_enc_init(&mem, &cfg);
    if (!enc)
        io_fail("invalid encoder settings", NULL);
    in = io_open(paths[0], "rb");
    if (!in)
        return 1;
    out = io_open(paths[1], "wb");
    if (!out)
        return 1;
    if (dpath && !(dout = io_open(dpath, "wb")))
        return 1;
    if (fmt == IO_AMB && mode != MODE_DSTAR)
        ok = io_write(out, IO_AMB_HEADER, 4, paths[1]);
    while (ok) {
        int16_t pcm[160];
        size_t got = io_read_pcm(in, pcm, 160);
        if (got != 160) {
            /* as the reference encoder's frame count */
            if (got)
                fprintf(stderr, "%s: partial final frame of %zu samples dropped\n", paths[0], got);
            break;
        }
        frames++;
        written += encode_frame(enc, mode, pcm, out, fmt, paths[1], dout, dpath, &ok);
    }
    if (ferror(in)) {
        fprintf(stderr, "%s: read error\n", paths[0]);
        ok = 0;
    }
    while (ok && written < frames)
        written += encode_frame(enc, mode, NULL, out, fmt, paths[1], dout, dpath, &ok);
    if (in != stdin)
        fclose(in);
    ok &= io_close(out, paths[1]);
    if (dout)
        ok &= io_close(dout, dpath);
    return ok ? 0 : 1;
}
