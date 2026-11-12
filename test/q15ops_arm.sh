#!/bin/sh
# q15ops_arm.sh CC QEMU SRC OUT: test_q15ops built for a Cortex-M33 (DSP extension) and run under qemu-arm
# with semihosting, checking q15ops.h's intrinsics as the host run checks the fallback
set -e
"$1" -mcpu=cortex-m33 -mthumb -O2 -std=c99 -Wall -Wextra -Werror --specs=rdimon.specs \
    -Wl,-Ttext-segment=0x100000 "$3/test/test_q15ops.c" -o "$4/test_q15ops_m33.elf" -lm
"$2" -cpu cortex-m33 "$4/test_q15ops_m33.elf"
