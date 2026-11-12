#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../src/internal.h"
#include "tone_fit.h"

#define FRAMES 30
#define FIT_START (12 * 160)

static void tone_bits(int id, int ad, uint8_t bytes[7])
{
    uint8_t d[49];
    xombe__tone_frame_bits(id, ad, d);
    xombe__pack49(d, bytes);
}

static void decode_tone(const uint8_t bytes[7], double *x)
{
    xombe_dec dec;
    int16_t pcm[160];
    int f, i;
    xombe_dec_set_noise(xombe_dec_init(&dec, XOMBE_DMR), XOMBE_NOISE_CONFORMANCE, 0);
    for (f = 0; f < FRAMES; f++) {
        xombe_dec_frame_bits(&dec, bytes, 0, pcm);
        for (i = 0; i < 160; i++)
            x[f * 160 + i] = pcm[i];
    }
}

static void compare_reference_table(const char *dir)
{
    static double x[FRAMES * 160];
    char path[4096], line[512], bits[64], kind[32];
    double max_df = 0.0, max_dl = 0.0;
    int rows = 0;
    FILE *f;
    snprintf(path, sizeof(path), "%s/tone_table.tsv", dir);
    f = fopen(path, "r");
    if (!f) {
        perror(path);
        return;
    }
    while (fgets(line, sizeof(line), f)) {
        double f1, f2 = 0.0, level, freq[2], amp[2];
        uint8_t d[49], bytes[7];
        int idx, i, nf, c;
        if (line[0] == '#')
            continue;
        nf = sscanf(line, "%d %31s %63s %lf %lf %lf", &idx, kind, bits, &f1, &f2, &level);
        if (nf == 5) {
            level = f2;
            nf = 1;
        } else if (nf == 6) {
            nf = 2;
        } else {
            continue;
        }
        for (i = 0; i < 49; i++)
            d[i] = (uint8_t)(bits[i] == '1');
        xombe__pack49(d, bytes);
        decode_tone(bytes, x);
        freq[0] = f1;
        freq[1] = f2;
        search_freqs(x + FIT_START, freq, nf, 0.01);
        sinusoid_fit(x + FIT_START, FIT_LEN, freq, nf, amp);
        for (c = 0; c < nf; c++) {
            double df = fabs(freq[c] - (c ? f2 : f1)), dl = fabs(20.0 * log10(amp[c] / 32768.0) - level);
            if (df > max_df)
                max_df = df;
            if (dl > max_dl)
                max_dl = dl;
        }
        rows++;
    }
    fclose(f);
    printf("reference tone table (informative): %d frames, max |frequency difference| %.2f Hz, "
           "max |peak level difference| %.2f dB\n", rows, max_df, max_dl);
}

int main(int argc, char **argv)
{
    static double x[FRAMES * 160];
    int id, fail = 0, cases = 0;
    double worst_f_single = 0.0, worst_f_dtmf = 0.0, worst_f_knox = 0.0, worst_f_progress = 0.0;
    double worst_l_single = 0.0, worst_l_dual = 0.0;
    for (id = 5; id <= 163; id++) {
        float flo, fhi;
        int kind = xombe__tone_id_accepted(id) ? xombe__tone_freqs(id, &flo, &fhi) : 0, ad;
        double freq[2], nominal[2], amp[2];
        int nf;
        uint8_t bytes[7];
        if (!kind)
            continue;
        nf = kind;
        nominal[0] = flo;
        nominal[1] = fhi;
        tone_bits(id, 90, bytes);
        decode_tone(bytes, x);
        freq[0] = nominal[0];
        freq[1] = nominal[1];
        search_freqs(x + FIT_START, freq, nf, kind == 1 ? 0.005 : 0.05);
        if (kind == 1) {
            double err = fabs(freq[0] - nominal[0]);
            if (err > worst_f_single)
                worst_f_single = err;
            if (err > 0.1) {
                fail++;
                printf("FAIL id %d: %.3f Hz, nominal %.3f\n", id, freq[0], nominal[0]);
            }
        } else {
            double tol = id >= 160 ? 0.04 : 0.01;
            int c;
            for (c = 0; c < 2; c++) {
                double rel = fabs(freq[c] - nominal[c]) / nominal[c];
                double *worst = id >= 160 ? &worst_f_progress : id >= 144 ? &worst_f_knox : &worst_f_dtmf;
                if (rel > *worst)
                    *worst = rel;
                if (rel > tol) {
                    fail++;
                    printf("FAIL id %d component %d: %.2f Hz, nominal %.0f\n", id, c, freq[c], nominal[c]);
                }
            }
        }
        for (ad = 30; ad <= (kind == 1 ? 120 : 115); ad++) {
            int c;
            tone_bits(id, ad, bytes);
            decode_tone(bytes, x);
            sinusoid_fit(x + FIT_START, FIT_LEN, freq, nf, amp);
            for (c = 0; c < nf; c++) {
                double err = 20.0 * log10(amp[c] / 32768.0 + 1e-30) - expected_db(ad, freq[c]);
                double *worst = kind == 1 ? &worst_l_single : &worst_l_dual;
                if (fabs(err) > fabs(*worst))
                    *worst = err;
                if (fabs(err) > 0.5) {
                    fail++;
                    printf("FAIL id %d AD %d component %d: level error %+.2f dB\n", id, ad, c, err);
                }
            }
            cases++;
        }
    }
    printf("T5: %d tone frames; worst single-tone frequency error %.4f Hz (limit 0.1)\n", cases, worst_f_single);
    printf("T5: worst frequency error DTMF %.2f%% (limit 1%%), Knox %.2f%% (limit 1%%), call progress %.2f%% "
           "(limit 4%%)\n",
           100.0 * worst_f_dtmf, 100.0 * worst_f_knox, 100.0 * worst_f_progress);
    printf("T5: worst level error single %+.3f dB, dual %+.3f dB (limit 0.5)\n", worst_l_single, worst_l_dual);
    if (argc > 1)
        compare_reference_table(argv[1]);
    printf("T5: %s\n", fail ? "FAIL" : "PASS");
    return fail ? 1 : 0;
}
