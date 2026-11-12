#ifndef XOMBE_BIQUAD_H
#define XOMBE_BIQUAD_H

/* high-pass: first-order section (pole HPF_POLE1) and second-order (1 - HPF_A1 z^-1 + HPF_A2 z^-2),
   three zeros at DC */
#define HPF_POLE1 0.9428
#define HPF_A1 1.9394
#define HPF_A2 0.9428

#endif
