#ifndef XOMBE_DSD_STATE_H
#define XOMBE_DSD_STATE_H

/* D-STAR frame logic on the shared decoder core (decoder.c), state xombe_dsd_logic in state.h */

#include "dstar_internal.h"

/* what a frame logged, a bad voice frame logs voice (or mute with hb set) but plays a repeat */
enum {
    XOMBE_DSD_LOG_VOICE,
    XOMBE_DSD_LOG_TONE,
    XOMBE_DSD_LOG_NULL,
    XOMBE_DSD_LOG_REPEAT,
    XOMBE_DSD_LOG_MUTE
};

/* sets up s and the D-STAR rules of a reset d */
void xombe__dsd_logic_init(xombe_dsd_logic *s, xombe_dec_core *d);
/* frame logic of one frame up to its sets. f from xombe__ds_frame_from_bits, NULL for a lost frame. flags: the
   error rule's XOMBE_FLAG_BF and XOMBE_FLAG_HB, eps_t and eps1 its error counts (xombe_fec_stats). sets the
   decision's e0 and e1 to -1 */
void xombe__dsd_start_frame(xombe_dsd_logic *s, xombe_dec_core *d, const xombe_ds_frame *f, unsigned flags,
                            int eps_t, int eps1);

/* xombe_dec's frame logic for D-STAR */
extern const struct xombe_dec_ops xombe__dsd_ops;

#endif
