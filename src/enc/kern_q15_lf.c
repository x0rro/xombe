/* q15 line fit (kernels.h, XOMBE_KERNELS=q15): int16 block-scaled spectrum with exact prefix sums of |S|^2 and
   exact 64-bit MACs, so the residual comes from ee - |sw|^2/ww in one pass without float cancellation */
#include <assert.h>
#include <string.h>
#include "kern_q15.h"

/* line fit table runs past the float table's end on its clamped last value, interpolating equal tail values
   is exact */
XOMBE_STATIC_ASSERT(XOMBE_LF_WT_OFF >= XOMBE_WTAB_LEN, lf_table_padded);

/* weights (1 - f, f) in Q14 for the fraction f of Q16 position x */
static uint32_t lf_coef(int32_t x)
{
    uint32_t f = (uint32_t)(x & 0xffff) >> 2;
    return (f << 16) | (16384u - f);
}

/* n bins of re (im at XOMBE_LF_IM) against the window from entries w, w + XOMBE_WTAB_OVERSAMPLE, ..., weight
   pair c, two bins per MAC */
static void lf_run(const int16_t *re, int n, const uint32_t *w, uint32_t c, int64_t *psr, int64_t *psi,
                   int64_t *pww)
{
    int64_t sr = 0, si = 0, ww = 0;
    const int16_t *end = re + (n & ~1);
    for (; re != end; re += 2, w += 2 * XOMBE_WTAB_OVERSAMPLE) {
        uint32_t wp = q15_pack_lsl2((uint32_t)(q15_smuad(w[0], c) >> 14),
                                    (uint32_t)q15_smuad(w[XOMBE_WTAB_OVERSAMPLE], c));
        sr = q15_smlald(q15_ld2(re), wp, sr);
        si = q15_smlald(q15_ld2(re + XOMBE_LF_IM), wp, si);
        ww = q15_smlald(wp, wp, ww);
    }
    if (n & 1) {
        int32_t v = q15_smuad(w[0], c) >> 14;
        sr += (int32_t)re[0] * v;
        si += (int32_t)re[XOMBE_LF_IM] * v;
        ww += v * v;
    }
    *psr = sr;
    *psi = si;
    *pww = ww;
}

/* low 32 bits of x >> sh, 0 < sh < 32 */
static uint32_t lf_shr(uint64_t x, int sh)
{
    return (uint32_t)x >> sh | (uint32_t)(x >> 32) << (32 - sh);
}

float xombe__kern_line_fit(const xombe_enc *e, const xombe_lf_spec *sp, float period, int kmin, int kmax,
                           float band_resid[8], float band_energy[8])
{
    const int16_t *re = sp->S->v;
    const uint64_t *ep = e->kq.ep;
    float per_harm = XOMBE_FFT / period, total = 0.0f, inv_g2, g = sp->g;
    uint64_t be[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    /* period exact in Q19 for 16 <= period <= XOMBE_FFT (24-bit mantissa), so harmonic and band edges compare
       exactly in integers: bin j is past harmonic l when j*period >= (l + 1/2)*XOMBE_FFT, as xombe__harmonic_of
       (midway goes up), band b = floor(l*XOMBE_BANDS_PER_UNIT/period). H = per_harm in Q20 (< 2^27) from 2^49/pq,
       taken 64 low so 56..73 under the true value (float rounding <= 8), and (2l + 1)*73 < 2^21 for l <= 513, so
       the edge estimates are exact or one low */
    int32_t pq = (int32_t)(period * 524288.0f), H = (int32_t)(562949953421312.0f / (float)pq) - 64;
    int k = kmin < 1 ? 1 : kmin, l = xombe__harmonic_of(k, period);
    /* bins lie within per_harm/2 + 1 of their harmonic, inside the table's padding. callers' periods are >= 18
       samples: pitch refinement >= 20, voicing metric at table f0 <= 0.0524 with 2 % jitter. per_harm >= 1 gives
       every harmonic a bin */
    assert(per_harm <= 2.0f * (XOMBE_LF_PAD_BINS - 2) && per_harm >= 1.0f && kmax <= XOMBE_FFT / 2);
    memset(band_resid, 0, 8 * sizeof(float));
    memset(band_energy, 0, 8 * sizeof(float));
    if (g == 0.0f)
        return 0.0f; /* all-zero spectrum fits exactly */
    inv_g2 = 1.0f / (g * g);
    for (; k <= kmax; l++) {
        int64_t sr = 0, si = 0, ww = 0, t = (int64_t)(2 * l + 1) << 28; /* (l + 1/2)*XOMBE_FFT in Q19 */
        uint64_t ee;
        float resid, pos;
        int64_t lb = (int64_t)l << 23; /* l*XOMBE_BANDS_PER_UNIT in Q19 */
        int32_t xs;
        int kend = (int)(((2 * l + 1) * H + 2097151) >> 21), band;
        if ((int64_t)kend * pq < t)
            kend++;
        if (kend > kmax + 1)
            kend = kmax + 1;
        if (l == 0) {
            k = kend;
            continue;
        }
        ee = ep[kend] - ep[k];
        /* table position of the first bin, in bins, the rest lie within per_harm/2 + 1 bins. in Q16 entries it
           steps exactly XOMBE_WTAB_OVERSAMPLE per bin */
        pos = (float)k - (float)l * per_harm;
        xs = (int32_t)(pos * (XOMBE_WTAB_OVERSAMPLE * 65536.0f));
        lf_run(re + k, kend - k, &e->kq.wt[XOMBE_LF_WT_OFF + (xs >> 16)], lf_coef(xs), &sr, &si, &ww);
        k = kend;
        if (ww > 0) {
            /* sums to 31 bits so products fit in 64: |sr|, |si| <= sqrt(ee*ww) need no say in the shift.
               numerator >= 0 up to rounding */
            int32_t e3, w3, r3, i3;
            int64_t num;
            int sh = 33 - q15_clz64(ee | (uint64_t)ww);
            if (sh <= 0) {
                sh = 0;
                e3 = (int32_t)ee;
                w3 = (int32_t)ww;
                r3 = (int32_t)sr;
                i3 = (int32_t)si;
            } else { /* sh <= 32 - 22, sums below 2^41 */
                e3 = (int32_t)lf_shr(ee, sh);
                w3 = (int32_t)lf_shr((uint64_t)ww, sh);
                r3 = (int32_t)lf_shr((uint64_t)sr, sh);
                i3 = (int32_t)lf_shr((uint64_t)si, sh);
            }
            num = (int64_t)e3 * w3 - (int64_t)r3 * r3 - (int64_t)i3 * i3;
            resid = num > 0 ? i64_to_float(num) / (float)w3 * (float)(1u << sh) : 0.0f;
        } else {
            resid = i64_to_float((int64_t)ee);
        }
        resid *= inv_g2;
        band = (int)(((int64_t)l * H) >> 26);
        if ((int64_t)(band + 1) * pq <= lb)
            band++;
        if (band > 7)
            band = 7;
        band_resid[band] += resid;
        be[band] += ee;
        total += resid;
    }
    for (l = 0; l < 8; l++)
        band_energy[l] = i64_to_float((int64_t)be[l]) * inv_g2;
    return total;
}
