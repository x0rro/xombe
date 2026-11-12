#ifndef XOMBE_PROF_H
#define XOMBE_PROF_H

/* encoder stage cycles for a Cortex-M bench (-DXOMBE_PROF_DWT): XOMBE_PF(k, stmt) runs stmt and adds its DWT
   cycles to xombe__pf[k], without XOMBE_PROF_DWT it runs stmt alone
   stages: 0 push, 1 analyze, 2 AMBE+2 voicing, 3 D-STAR voicing, 4 background, 5 AMBE+2 emit, 6 D-STAR emit,
   7 D-STAR tone push
   nested: 8 front end voicing and 9 voicing metric in 3, 10 AMBE+2 tone detector in 1, 11 detection and
   12 tone-in-noise pass in 7, 13 D-STAR quantizer, 14 D-STAR targets in 6
   15 sums L over the D-STAR quantizer calls */

#ifdef XOMBE_PROF_DWT
extern unsigned long xombe__pf[16];
#define XOMBE_PF_NOW() (*(volatile unsigned long *)0xE0001004u) /* DWT_CYCCNT, counter enabled by the caller */
#define XOMBE_PF(k, ...)                                                                                       \
    do {                                                                                                       \
        unsigned long xombe_pf_t0 = XOMBE_PF_NOW();                                                            \
        __VA_ARGS__;                                                                                           \
        xombe__pf[k] += XOMBE_PF_NOW() - xombe_pf_t0;                                                          \
    } while (0)
#define XOMBE_PF_ADD(k, n) (xombe__pf[k] += (unsigned long)(n))
#else
#define XOMBE_PF(k, ...)                                                                                       \
    do {                                                                                                       \
        __VA_ARGS__;                                                                                           \
    } while (0)
#define XOMBE_PF_ADD(k, n) ((void)0)
#endif

#endif
