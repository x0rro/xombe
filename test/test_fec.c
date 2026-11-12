#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../src/internal.h"

static int failures, checks;

static void check(int ok, const char *what, int index)
{
    checks++;
    if (!ok) {
        failures++;
        if (failures <= 20)
            fprintf(stderr, "FAIL %s (item %d)\n", what, index);
    }
}

static int hex_bytes(const char *s, uint8_t *out, int n)
{
    int i;
    for (i = 0; i < n; i++) {
        unsigned v;
        if (sscanf(s + 2 * i, "%2x", &v) != 1)
            return 0;
        out[i] = (uint8_t)v;
    }
    return 1;
}

static FILE *open_vec(const char *dir, const char *name)
{
    char path[4096];
    FILE *f;
    snprintf(path, sizeof(path), "%s/%s", dir, name);
    f = fopen(path, "r");
    if (!f) {
        perror(path);
        exit(2);
    }
    return f;
}

static void test_encode_vectors(const char *dir)
{
    char line[512], h49[64], h72[64];
    FILE *f = open_vec(dir, "fec_encode_vectors.txt");
    int n = 0;
    while (fgets(line, sizeof(line), f)) {
        uint8_t b[7], air[9], mine[9], back[7];
        xombe_fec fec;
        xombe_fec_stats st;
        unsigned flags;
        if (line[0] == '#' || sscanf(line, "%63s %63s", h49, h72) != 2)
            continue;
        hex_bytes(h49, b, 7);
        hex_bytes(h72, air, 9);
        xombe_fec_encode(XOMBE_DMR, b, mine);
        check(!memcmp(mine, air, 9), "encode vector", n);
        xombe_fec_init(&fec, XOMBE_DMR);
        flags = xombe_fec_decode(&fec, air, back, &st);
        check(!memcmp(back, b, 7) && st.eps0 == 0 && st.eps1 == 0 && flags == 0, "decode vector", n);
        n++;
    }
    fclose(f);
    check(n == 80, "80 encode vectors", n);
    printf("fec_encode_vectors: %d pairs\n", n);
}

static void test_air_map(const char *dir)
{
    char line[512], cw[16], label[32];
    FILE *f = open_vec(dir, "air72_map.txt");
    static const int offset[4] = {0, 24, 47, 58};
    static const int data_base[4] = {0, 12, 24, 35};
    static const int data_len[4] = {12, 12, 11, 14};
    uint8_t w[72], a[72];
    int n = 0, i;
    for (i = 0; i < 72; i++)
        w[i] = (uint8_t)i;
    xombe__interleave(w, a);
    while (fgets(line, sizeof(line), f)) {
        int air_bit, k, bit;
        char *p;
        if (line[0] == '#')
            continue;
        if (sscanf(line, "%d %15s", &air_bit, cw) != 2)
            continue;
        p = strchr(line, '[');
        if (!p || sscanf(p + 1, "%d] %31s", &bit, label) != 2)
            continue;
        k = cw[1] - '0';
        check(a[air_bit] == offset[k] + bit, "air map position", air_bit);
        if (bit < data_len[k]) {
            char want[16];
            snprintf(want, sizeof(want), "d%d", data_base[k] + bit);
            check(!strcmp(label, want), "air map data index", air_bit);
        } else {
            check(!strcmp(label, "parity"), "air map parity", air_bit);
        }
        n++;
    }
    fclose(f);
    check(n == 72, "72 air map entries", n);
    printf("air72_map: %d entries\n", n);
}

static uint32_t rng_state = 12345u;
static uint32_t rnd(void)
{
    rng_state = rng_state * 1103515245u + 12345u;
    return rng_state >> 8;
}

static void flip_codeword_bit(uint8_t air[9], int w_index)
{
    uint8_t w[72], a[72];
    int i;
    for (i = 0; i < 72; i++)
        w[i] = (uint8_t)(i == w_index);
    xombe__interleave(w, a);
    for (i = 0; i < 72; i++)
        if (a[i])
            air[i >> 3] ^= (uint8_t)(0x80 >> (i & 7));
}

static void test_injection(void)
{
    int trial, e0, e1, n = 0;
    for (trial = 0; trial < 400; trial++) {
        for (e0 = 0; e0 <= 3; e0++) {
            for (e1 = 0; e1 <= 3; e1++) {
                uint8_t b[7], air[9], back[7];
                int used[24] = {0}, used1[23] = {0}, i;
                xombe_fec fec;
                xombe_fec_stats st;
                for (i = 0; i < 7; i++)
                    b[i] = (uint8_t)rnd();
                b[6] &= 0x80;
                xombe_fec_encode(XOMBE_DMR, b, air);
                for (i = 0; i < e0;) {
                    int pos = (int)(rnd() % 24);
                    if (used[pos])
                        continue;
                    used[pos] = 1;
                    flip_codeword_bit(air, pos);
                    i++;
                }
                for (i = 0; i < e1;) {
                    int pos = (int)(rnd() % 23);
                    if (used1[pos])
                        continue;
                    used1[pos] = 1;
                    flip_codeword_bit(air, 24 + pos);
                    i++;
                }
                xombe_fec_init(&fec, XOMBE_DMR);
                xombe_fec_decode(&fec, air, back, &st);
                check(!memcmp(back, b, 7), "corrected data", trial);
                check(st.eps0 == e0 && st.eps1 == e1, "error counts", trial);
                check(st.bf == (e0 >= 2 && e0 + e1 >= 6), "bad-frame flag", trial);
                n++;
            }
        }
        {
            uint8_t b[7], air[9], back[7];
            int i;
            xombe_fec fec;
            xombe_fec_stats st;
            for (i = 0; i < 7; i++)
                b[i] = (uint8_t)rnd();
            b[6] &= 0x80;
            xombe_fec_encode(XOMBE_DMR, b, air);
            flip_codeword_bit(air, 23);
            xombe_fec_init(&fec, XOMBE_DMR);
            xombe_fec_decode(&fec, air, back, &st);
            check(!memcmp(back, b, 7) && st.eps0 == 1 && st.eps1 == 0, "parity-only error", trial);
            n++;
        }
    }
    printf("error injection: %d frames\n", n);
}

static void test_c0x4(const char *dir)
{
    char line[512], s49[64], a72[64], o49[64];
    FILE *f = open_vec(dir, "fec_decode_c0x4.txt");
    int n = 0;
    while (fgets(line, sizeof(line), f)) {
        uint8_t air[9], want[7], got[7];
        int e0, e1, bf;
        xombe_fec fec;
        xombe_fec_stats st;
        if (line[0] == '#' || sscanf(line, "%63s %63s %63s %d %d %d", s49, a72, o49, &e0, &e1, &bf) != 6)
            continue;
        hex_bytes(a72, air, 9);
        hex_bytes(o49, want, 7);
        xombe_fec_init(&fec, XOMBE_DMR);
        xombe_fec_decode(&fec, air, got, &st);
        check(!memcmp(got, want, 7), "c0x4 output bits", n);
        check(st.eps0 == e0 && st.eps1 == e1 && st.bf == bf, "c0x4 counts and BF", n);
        n++;
    }
    fclose(f);
    check(n == 200, "200 c0x4 frames", n);
    printf("fec_decode_c0x4: %d frames\n", n);
}

static void test_hb_onset(const char *dir)
{
    char line[512], a72[64], o49[64];
    FILE *f = open_vec(dir, "fec_hb_onset.txt");
    int n = 0, cur_seq = -1;
    xombe_fec fec;
    xombe_fec_init(&fec, XOMBE_DMR);
    while (fgets(line, sizeof(line), f)) {
        uint8_t air[9], want[7], got[7];
        int seq, fr, e0, e1, bf, hb;
        xombe_fec_stats st;
        if (line[0] == '#' ||
            sscanf(line, "%d %d %63s %63s %d %d %d %d", &seq, &fr, a72, o49, &e0, &e1, &bf, &hb) != 8)
            continue;
        if (seq != cur_seq) {
            xombe_fec_init(&fec, XOMBE_DMR);
            cur_seq = seq;
        }
        hex_bytes(a72, air, 9);
        hex_bytes(o49, want, 7);
        xombe_fec_decode(&fec, air, got, &st);
        check(!memcmp(got, want, 7), "hb_onset output bits", n);
        check(st.eps0 == e0 && st.eps1 == e1 && st.bf == bf && st.hb == hb, "hb_onset flags", n);
        n++;
    }
    fclose(f);
    printf("fec_hb_onset: %d frames\n", n);
}

static void test_fields(void)
{
    uint8_t d[49];
    int b[9], c, i, f;
    for (c = 0; c < 64; c++) {
        memset(d, 0, sizeof(d));
        xombe__put_bits(d, 0, 6, c);
        i = xombe__frame_class(d);
        check(i == (c < 60 ? XOMBE_CLASS_VOICE : c < 62 ? XOMBE_CLASS_ERASURE
                    : c == 62 ? XOMBE_CLASS_SILENCE : XOMBE_CLASS_TONE), "frame class", c);
    }
    for (f = 0; f < 9; f++) {
        static const int lsb[9] = {39, 35, 36, 43, 42, 44, 45, 46, 48};
        static const int msb[9] = {0, 4, 8, 12, 20, 24, 28, 31, 34};
        static const int width[9] = {7, 5, 5, 9, 7, 5, 4, 4, 3};
        memset(b, 0, sizeof(b));
        b[f] = 1;
        memset(d, 0, sizeof(d));
        xombe__put_fields(b, d);
        check(d[lsb[f]] == 1, "field LSB position", f);
        b[f] = 1 << (width[f] - 1);
        memset(d, 0, sizeof(d));
        xombe__put_fields(b, d);
        check(d[msb[f]] == 1, "field MSB position", f);
    }
    for (i = 0; i < 2000; i++) {
        int back[9];
        uint8_t packed[7], again[49];
        for (c = 0; c < 49; c++)
            d[c] = (uint8_t)(rnd() & 1);
        xombe__get_fields(d, b);
        memset(again, 0, sizeof(again));
        xombe__put_fields(b, again);
        check(!memcmp(again, d, 49), "field round trip", i);
        xombe__pack49(d, packed);
        xombe__unpack49(packed, again);
        check(!memcmp(again, d, 49), "pack round trip", i);
        xombe__get_fields(again, back);
        check(!memcmp(back, b, sizeof(b)), "field values", i);
    }
    printf("class and field map: ok\n");
}

/* Q24 running error rate against 0.95*eps_r + 0.05*eps_t/47 in double, over runs of error counts hovering
   around the HB threshold (eps_t = 4 or 5 in steady state gives eps_r 0.085 or 0.106) */
static void test_eps_r(void)
{
    xombe_fec f;
    double r = 0.0, worst = 0.0;
    int i, hb_diff = 0;
    xombe_fec_init(&f, XOMBE_DMR);
    for (i = 0; i < 200000; i++) {
        xombe_fec_stats st;
        int burst = (i / 1000) % 3;
        memset(&st, 0, sizeof(st));
        st.eps0 = (int)(rnd() % (burst == 0 ? 3u : burst == 1 ? 5u : 8u));
        st.eps1 = (int)(rnd() % 4u);
        if (st.eps0 > 4)
            st.eps0 = 4;
        xombe__fec_stats_update(&f, &st);
        r = 0.95 * r + 0.05 * st.eps_t / 47.0;
        if (fabs(st.eps_r / 16777216.0 - r) > worst)
            worst = fabs(st.eps_r / 16777216.0 - r);
        if (st.hb != (r > 0.096) && fabs(r - 0.096) > 1e-6)
            hb_diff++;
    }
    check(worst < 1e-6, "eps_r within 1e-6 of the double recursion", 0);
    check(hb_diff == 0, "HB as in the double recursion away from the threshold", hb_diff);
    printf("eps_r: worst %.2g against double\n", worst);
}

int main(int argc, char **argv)
{
    const char *dir = argc > 1 ? argv[1] : "vectors/fec";
    test_encode_vectors(dir);
    test_air_map(dir);
    test_injection();
    test_c0x4(dir);
    test_hb_onset(dir);
    test_fields();
    test_eps_r();
    printf("T1: %d/%d checks passed\n", checks - failures, checks);
    return failures ? 1 : 0;
}
