#include <limits.h>
#include "../internal.h"

int xombe__conceal_decide(xombe_dec_core *d, int bad, int hb)
{
    int last_three_concealed = (d->hist & 7u) == 7u;
    int c = XOMBE_CONCEAL_NONE;
    if (hb)
        c = XOMBE_CONCEAL_MUTE;
    else if (bad)
        c = last_three_concealed ? XOMBE_CONCEAL_MUTE : XOMBE_CONCEAL_REPEAT;
    d->hist = ((d->hist << 1) | (c != XOMBE_CONCEAL_NONE ? 1u : 0u)) & 7u;
    d->run = c == XOMBE_CONCEAL_NONE ? 0 : d->run < INT_MAX ? d->run + 1 : INT_MAX;
    return c;
}
