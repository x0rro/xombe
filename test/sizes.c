#include <stdio.h>
#include <string.h>
#include "xombe.h"

/* prints this build's row of the size table (scripts/size_table.sh) and checks every *_size() <= its *_SIZE.
   usage: sizes [--header] */

static int failures;

static void check(size_t got, size_t max, const char *what)
{
    if (got > max) {
        failures++;
        fprintf(stderr, "FAIL %s: %zu > %zu\n", what, got, max);
    }
}

#if XOMBE_HAVE_DSTAR
/* state and workspace bytes xombe_xc_init_shared writes for one input codec (it zeroes both in full): the last
   byte off its fill, over two fills */
static void shared_bytes(enum xombe_codec from, size_t *state, size_t *ws)
{
    static unsigned char mem[2][XOMBE_XC_SIZE]; /* larger than any shared state */
    static unsigned char qws[2][XOMBE_QWS_SIZE];
    static const unsigned char pat[2] = {0xA5, 0x5A};
    size_t s[2], w[2];
    int i;
    for (i = 0; i < 2; i++) {
        memset(mem[i], pat[i], sizeof(mem[i]));
        memset(qws[i], pat[i], sizeof(qws[i]));
        xombe_xc_init_shared(mem[i], from, qws[i]);
    }
    for (i = 0; i < 2; i++) {
        size_t k = sizeof(mem[i]);
        while (k > 0 && mem[i][k - 1] == pat[i])
            k--;
        s[i] = k;
        k = sizeof(qws[i]);
        while (k > 0 && qws[i][k - 1] == pat[i])
            k--;
        w[i] = k;
    }
    *state = s[0] > s[1] ? s[0] : s[1];
    *ws = w[0] > w[1] ? w[0] : w[1];
}
#endif

int main(int argc, char **argv)
{
    char row[1024];
    size_t n = 0;
    if (argc > 1 && !strcmp(argv[1], "--header")) {
        printf("| kernels | D-STAR | enc | enc delay DMR / D-STAR / both | dec | xc | xc shared DMR / D-STAR | QWS |\n");
        printf("|---|---|---|---|---|---|---|---|\n");
        return 0;
    }
    n += (size_t)snprintf(row + n, sizeof(row) - n, "| %s | %s |", XOMBE_KERNELS_Q15 ? "q15" : "float",
                          XOMBE_HAVE_DSTAR ? "on" : "off");
#if XOMBE_HAVE_ENCODER
    {
        static xombe_enc_mem mem;
        static const unsigned outs[3] = {XOMBE_ENC_DMR, XOMBE_ENC_DSTAR, XOMBE_ENC_DMR | XOMBE_ENC_DSTAR};
        xombe_enc_cfg cfg;
        int i;
        check(xombe_enc_size(), XOMBE_ENC_SIZE, "xombe_enc_size");
        n += (size_t)snprintf(row + n, sizeof(row) - n, " %zu |", xombe_enc_size());
        for (i = 0; i < 3; i++) {
            xombe_enc *e;
            xombe_enc_cfg_default(&cfg);
            cfg.outputs = outs[i];
            e = xombe_enc_init(&mem, &cfg);
            if (e)
                n += (size_t)snprintf(row + n, sizeof(row) - n, "%s%d", i ? " / " : " ", xombe_enc_delay(e));
            else
                n += (size_t)snprintf(row + n, sizeof(row) - n, "%s-", i ? " / " : " ");
            if (!e && !(outs[i] & XOMBE_ENC_DSTAR)) {
                failures++;
                fprintf(stderr, "FAIL encoder init, outputs %u\n", outs[i]);
            }
        }
        n += (size_t)snprintf(row + n, sizeof(row) - n, " |");
    }
#else
    n += (size_t)snprintf(row + n, sizeof(row) - n, " - | - |");
#endif
    check(xombe_dec_size(), XOMBE_DEC_SIZE, "xombe_dec_size");
    n += (size_t)snprintf(row + n, sizeof(row) - n, " %zu |", xombe_dec_size());
#if XOMBE_HAVE_DSTAR
    {
        size_t st, ws, qws;
        check(xombe_xc_size(), XOMBE_XC_SIZE, "xombe_xc_size");
        shared_bytes(XOMBE_DMR, &st, &ws);
        check(st, XOMBE_XC_SHARED_SIZE, "xc shared state from DMR");
        check(ws, XOMBE_QWS_SIZE, "xc workspace from DMR");
        qws = ws;
        n += (size_t)snprintf(row + n, sizeof(row) - n, " %zu | %zu", xombe_xc_size(), st);
#if XOMBE_HAVE_ENCODER
        shared_bytes(XOMBE_DSTAR, &st, &ws);
        check(st, XOMBE_XC_SHARED_SIZE, "xc shared state from D-STAR");
        check(ws, XOMBE_QWS_SIZE, "xc workspace from D-STAR");
        qws = ws > qws ? ws : qws;
        n += (size_t)snprintf(row + n, sizeof(row) - n, " / %zu |", st);
#else
        n += (size_t)snprintf(row + n, sizeof(row) - n, " / - |");
#endif
        n += (size_t)snprintf(row + n, sizeof(row) - n, " %zu |", qws);
    }
#else
    n += (size_t)snprintf(row + n, sizeof(row) - n, " - | - | - |");
#endif
    puts(row);
    if (failures)
        fprintf(stderr, "sizes: %d failures\n", failures);
    return failures != 0;
}
