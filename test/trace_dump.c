/* dumps the codec's stage taps (src/state.h) for one input, for test/trace_cmp.py to compare two builds or
   against another implementation, decoder stages in two 80-sample halves per frame, encoder stages in the default
   configuration
   OUT.trace: little-endian records, u32 frame, u8 half, u8 stage (enum xombe_stage), u16 n, f32 f0, n x f32
   encoder records carry the frame they belong to (the call that returned its bits, less XOMBE_ENC_DELAY)
   and half 0 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../src/internal.h"
#include "../tools/io.h"

typedef struct {
    FILE *f;
    const char *path;
    uint32_t frame;
    uint8_t half;
    int failed;
} dump;

static void put_u32(uint8_t *p, uint32_t v)
{
    p[0] = (uint8_t)v;
    p[1] = (uint8_t)(v >> 8);
    p[2] = (uint8_t)(v >> 16);
    p[3] = (uint8_t)(v >> 24);
}

static void put_f32(uint8_t *p, float x)
{
    uint32_t v;
    memcpy(&v, &x, sizeof(v));
    put_u32(p, v);
}

static void tap(void *ctx, enum xombe_stage stage, float f0, const float *v, int n)
{
    dump *t = ctx;
    uint8_t head[12], buf[4 * (XOMBE_HALF_MAX + XOMBE_LMAX)];
    int i;
    if (t->failed || n < 0 || n > (int)(sizeof(buf) / 4))
        return;
    put_u32(head, t->frame);
    head[4] = t->half;
    head[5] = (uint8_t)stage;
    head[6] = (uint8_t)n;
    head[7] = (uint8_t)(n >> 8);
    put_f32(head + 8, f0);
    for (i = 0; i < n; i++)
        put_f32(buf + 4 * i, v[i]);
    if (!io_write(t->f, head, sizeof(head), t->path) || !io_write(t->f, buf, 4 * (size_t)n, t->path))
        t->failed = 1;
}

static int run_dec(FILE *in, const char *in_path, dump *t)
{
    static xombe_dec_mem mem;
    xombe_dec *d = xombe_dec_init(&mem, XOMBE_DMR);
    uint8_t rec[IO_AMB_RECORD], bits[7];
    char head[sizeof(IO_AMB_HEADER) - 1];
    int16_t pcm[XOMBE_HALF_MAX];
    size_t got;
    xombe_dec_set_noise(d, XOMBE_NOISE_CONFORMANCE, 0);
    if (fread(head, 1, sizeof(head), in) != sizeof(head) || memcmp(head, IO_AMB_HEADER, sizeof(head)) != 0) {
        fprintf(stderr, "%s: not an .amb file\n", in_path);
        return 0;
    }
    xombe__dec_set_tap(d, tap, t);
    while ((got = fread(rec, 1, sizeof(rec), in)) == sizeof(rec)) {
        io_amb_to_bytes(rec, bits);
        t->half = 0;
        xombe_dec_half_bits(d, bits, 0, 0, 80, pcm);
        t->half = 1;
        xombe_dec_half_bits(d, bits, 0, 1, 80, pcm);
        t->frame++;
    }
    if (got) {
        fprintf(stderr, "%s: truncated frame\n", in_path);
        return 0;
    }
    return 1;
}

#if XOMBE_HAVE_ENCODER
static int run_enc(FILE *in, dump *t)
{
    static xombe_enc_mem mem;
    xombe_enc *e = xombe_enc_init(&mem, NULL);
    int16_t pcm[XOMBE_FRAME];
    uint8_t bits[7];
    uint32_t call = 0;
    int i;
    xombe__enc_set_tap(e, tap, t);
    t->half = 0;
    for (;; call++) {
        t->frame = call - XOMBE_ENC_DELAY;
        if (io_read_pcm(in, pcm, XOMBE_FRAME) != XOMBE_FRAME)
            break;
        xombe_enc_frame_bits(e, pcm, NULL, bits, NULL);
    }
    for (i = 0; i < XOMBE_ENC_DELAY; i++, call++) {
        t->frame = call - XOMBE_ENC_DELAY;
        xombe_enc_frame_bits(e, NULL, NULL, bits, NULL);
    }
    return 1;
}
#endif

int main(int argc, char **argv)
{
    dump t;
    FILE *in;
    int ok;
    if (argc != 4 || (strcmp(argv[1], "dec") != 0 && strcmp(argv[1], "enc") != 0)) {
        fprintf(stderr, "usage: trace_dump dec IN.amb OUT.trace | trace_dump enc IN.raw OUT.trace\n");
        return 2;
    }
#if !XOMBE_HAVE_ENCODER
    if (strcmp(argv[1], "enc") == 0) {
        fprintf(stderr, "trace_dump: built without the encoder\n");
        return 2;
    }
#endif
    memset(&t, 0, sizeof(t));
    if (!(in = io_open(argv[2], "rb")))
        return 1;
    t.path = argv[3];
    if (!(t.f = io_open(argv[3], "wb")))
        return 1;
#if XOMBE_HAVE_ENCODER
    ok = strcmp(argv[1], "dec") == 0 ? run_dec(in, argv[2], &t) : run_enc(in, &t);
#else
    ok = run_dec(in, argv[2], &t);
#endif
    fclose(in);
    ok &= io_close(t.f, argv[3]) && !t.failed;
    return ok ? 0 : 1;
}
