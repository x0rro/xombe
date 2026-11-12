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
            "usage: xombe-dec [options] IN OUT.raw\n"
            "  OUT.raw: 8 kHz s16le mono, 160 samples per frame. '-' is stdin or stdout\n"
            "input formats (modes):\n"
            "  (default)        AMBE+2 .amb file\n"
            "  --bits           AMBE+2, 49 bytes per frame, one bit per byte\n"
            "  --air            AMBE+2 72-bit air frames, 9 bytes, FEC applied\n"
            "  --dstar          D-STAR (AMBE 2400) frames, 9 bytes as received, FEC applied" NO_DSTAR "\n"
            "  --stream         AMBE+2 JSONL stream IN.jsonl[.gz], one half call per line\n"
            "options:\n"
            "  --lost LIST      frames lost, e.g. 3,10-12 (XOMBE_FLAG_LOST)\n"
            "  --conformance-noise\n"
            "                   noise source that reproduces the test vectors (default: production source)\n"
            "  --decisions      per-frame decisions to stderr: 'n class conceal run mute e0 e1 bf hb', mute: the\n"
            "                   muted set's source frame, init, or -. e0 e1: Golay corrections, -1 if lost or\n"
            "                   post-FEC (.amb, --bits)\n"
            "  --help           this text\n"
            "  a truncated final frame is an error\n");
    exit(status);
}

static void print_decision(const xombe_dec *d, long frame)
{
    xombe_decision r;
    xombe_dec_decision_get(d, &r);
    io_print_decision(frame, &r, NULL);
}

static int read_frame(FILE *in, const char *path, enum io_format fmt, uint8_t bits[XOMBE_BITS_BYTES], uint8_t air[IO_AIR_RECORD])
{
    uint8_t buf[IO_BITS_RECORD];
    size_t want = fmt == IO_AMB ? IO_AMB_RECORD : fmt == IO_BITS ? IO_BITS_RECORD : IO_AIR_RECORD;
    int i, got = io_read_record(in, fmt == IO_AIR ? air : buf, want, path);
    if (got != 1)
        return got;
    if (fmt == IO_AMB) {
        io_amb_to_bytes(buf, bits);
    } else if (fmt == IO_BITS) {
        for (i = 0; i < IO_BITS_RECORD; i++)
            buf[i] &= 1;
        io_pack_bits(buf, bits);
    }
    return 1;
}

static int run_stream(const char *path, const char *out_path, int conformance, int show)
{
    stream_reader r;
    stream_half h;
    FILE *out;
    long frame = -1;
    xombe_dec_mem mem;
    xombe_dec *dec;
    int ev, ok = 1;
    if (!stream_open(&r, path))
        return 1;
    out = io_open(out_path, "wb");
    if (!out) {
        stream_close(&r);
        return 1;
    }
    dec = xombe_dec_init(&mem, XOMBE_DMR);
    if (conformance)
        xombe_dec_set_noise(dec, XOMBE_NOISE_CONFORMANCE, 0);
    while (ok && (ev = stream_next(&r, &h)) != STREAM_END) {
        if (ev == STREAM_ERROR) {
            ok = 0;
        } else if (ev == STREAM_SET) {
            xombe_dec_set_errors(dec, r.eps_t, r.eps_1);
        } else if (ev == STREAM_INIT) {
            dec = xombe_dec_init(&mem, XOMBE_DMR);
            if (conformance)
                xombe_dec_set_noise(dec, XOMBE_NOISE_CONFORMANCE, 0);
            frame = -1;
        } else {
            int16_t pcm[XOMBE_HALF_MAX];
            int ns = xombe_dec_half_bits(dec, h.bits, h.flags, h.second_half, h.nsamples, pcm);
            if (ns < 0) {
                fprintf(stderr, "%s:%ld: half-frame call out of order\n", path, r.lineno);
                ok = 0;
                break;
            }
            ok = io_write_pcm(out, pcm, (size_t)ns, out_path);
            if (!h.second_half) {
                frame++;
                if (show)
                    print_decision(dec, frame);
            }
        }
    }
    ok &= stream_close(&r);
    ok &= io_close(out, out_path);
    return ok ? 0 : 1;
}

int main(int argc, char **argv)
{
    enum io_format fmt = IO_AMB;
    int conformance = 0, show = 0, stream = 0, dstar = 0, np = 0, nfmt = 0, i, ok = 1, got;
    io_frames lost = {NULL, 0};
    const char *paths[2];
    FILE *in, *out;
    xombe_dec_mem mem;
    xombe_dec *dec;
    long frame = 0;
    io_prog = "xombe-dec";
    for (i = 1; i < argc; i++) {
        const char *a = argv[i];
        if (!strcmp(a, "--help"))
            usage(stdout, 0);
        else if (!strcmp(a, "--bits")) {
            fmt = IO_BITS;
            nfmt++;
        } else if (!strcmp(a, "--air")) {
            fmt = IO_AIR;
            nfmt++;
        } else if (!strcmp(a, "--dstar")) {
            if (!XOMBE_HAVE_DSTAR)
                io_fail(": not in this build (XOMBE_DSTAR=OFF)", a);
            dstar = 1;
            fmt = IO_AIR;
            nfmt++;
        } else if (!strcmp(a, "--stream")) {
            stream = 1;
            nfmt++;
        } else if (!strcmp(a, "--lost") && i + 1 < argc) {
            if (!io_frames_parse(argv[++i], &lost))
                io_fail(": not a frame list (e.g. 3,10-12)", argv[i]);
        } else if (!strcmp(a, "--conformance-noise"))
            conformance = 1;
        else if (!strcmp(a, "--decisions"))
            show = 1;
        else if (a[0] == '-' && a[1])
            usage(stderr, 2);
        else if (np < 2)
            paths[np++] = a;
        else
            usage(stderr, 2);
    }
    if (np != 2)
        usage(stderr, 2);
    if (nfmt > 1)
        io_fail("--bits, --air, --dstar and --stream are exclusive", NULL);
    if (stream) {
        if (lost.n)
            io_fail("--lost does not apply to --stream", NULL);
        return run_stream(paths[0], paths[1], conformance, show);
    }
    in = io_open(paths[0], "rb");
    if (!in)
        return 1;
    if (fmt == IO_AMB) {
        char hdr[4];
        if (fread(hdr, 1, 4, in) != 4 || memcmp(hdr, IO_AMB_HEADER, 4)) {
            fprintf(stderr, "%s: not an .amb file\n", paths[0]);
            return 1;
        }
    }
    out = io_open(paths[1], "wb");
    if (!out)
        return 1;
    dec = xombe_dec_init(&mem, dstar ? XOMBE_DSTAR : XOMBE_DMR);
    if (conformance)
        xombe_dec_set_noise(dec, XOMBE_NOISE_CONFORMANCE, 0);
    for (;;) {
        uint8_t bits[XOMBE_BITS_BYTES], air[IO_AIR_RECORD];
        int16_t pcm[160];
        unsigned flags = io_frames_has(&lost, frame) ? XOMBE_FLAG_LOST : 0;
        if ((got = read_frame(in, paths[0], fmt, bits, air)) != 1) {
            ok = got == 0;
            break;
        }
        if (fmt == IO_AIR)
            xombe_dec_frame(dec, air, flags, pcm);
        else
            xombe_dec_frame_bits(dec, bits, flags, pcm);
        if (show)
            print_decision(dec, frame);
        if (!io_write_pcm(out, pcm, 160, paths[1])) {
            ok = 0;
            break;
        }
        frame++;
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
