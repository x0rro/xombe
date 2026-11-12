#include <string.h>
#include "../internal.h"

void xombe__tone_set(xombe_set *s, int id, int ad)
{
    float f0_hz, amp = xombe__tone_log2_amp(ad);
    int l, l1, l2, top;
    memset(s, 0, sizeof(*s));
    for (l = 0; l <= XOMBE_LMAX; l++)
        s->e[l] = XOMBE_AMP_OFF;
    xombe__tone_harmonics(id, &f0_hz, &l1, &l2);
    s->cls = XOMBE_CLASS_TONE;
    s->prepared = 1;
    s->f0 = XOMBE_HZ(f0_hz);
    top = l2 > l1 ? l2 : l1;
    s->L = top > 9 ? top : 9;
    if (l1 > 0)
        s->e[l1] = amp;
    if (l2 > 0)
        s->e[l2] = amp;
}
