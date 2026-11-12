#ifndef XOMBE_KERN_STATE_H
#define XOMBE_KERN_STATE_H

/* kernel backend state, last member of struct xombe_enc, empty for float, and the analysis spectrum's layout,
   which only the backend reads (kernels.h) */
#if XOMBE_KERNELS_Q15
#define XOMBE_LF_BINS (XOMBE_FFT / 2 + 1)
#define XOMBE_LF_IM (XOMBE_LF_BINS + 1)
/* block floating point, bin k = (v[k], v[XOMBE_LF_IM + k]) * sc, components within +/-32767 and the largest
   within 2^-10 of it, sc 0 and v all zero for an all-zero spectrum. XOMBE_LF_IM even, so re[k] and im[k] share
   alignment for the line fit's pair loads */
typedef struct xombe_spec {
    int16_t v[XOMBE_LF_IM + XOMBE_LF_BINS + 1];
    float sc;
} xombe_spec;
/* q15 line fit window table and prepared spectrum, off the stack.
   wt: indexed by signed table position, out to XOMBE_LF_PAD_BINS bins either side of a harmonic
   ep: prefix sums of |S|^2 in the prepared spectrum's units
   spec: a prepared range shifted up, when its peak is below 16384 */
#define XOMBE_LF_PAD_BINS 64
#define XOMBE_LF_WT_OFF (XOMBE_LF_PAD_BINS * XOMBE_WTAB_OVERSAMPLE)
typedef struct {
    uint32_t wt[2 * XOMBE_LF_WT_OFF + 1];
    uint64_t ep[XOMBE_LF_BINS + 1];
    xombe_spec spec;
} xombe_kern_q15;
#define XOMBE_KERN_STATE xombe_kern_q15 kq;
#else
typedef struct xombe_spec {
    xombe_cpx b[XOMBE_SPEC_BINS];
} xombe_spec;
#define XOMBE_KERN_STATE
#endif

#endif
