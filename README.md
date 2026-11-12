# xombe

C99 library implementing AMBE 2400 and AMBE+2 used by D-STAR and DMR. decoders,
a single input dual output encoder, parameter-space transcoders, and FEC.

## minimum usage

play a DMR call saved by dsd (`.amb`):

    xombe-dec call.amb - | aplay -q -f S16_LE -r 8000

encode a wav for DMR:

    sox talk.wav -t raw -r 8000 -c 1 -b 16 -e signed -L - | xombe-enc - talk.amb

to decode DMR air frames in C:

    #include <xombe.h>

    static xombe_dec_mem mem;
    xombe_dec *dec = xombe_dec_init(&mem, XOMBE_DMR);    /* or XOMBE_DSTAR */

    xombe_dec_frame(dec, frame, 0, pcm);    /* 9 bytes in, 160 samples out */

when a frame is lost set flags to `XOMBE_FLAG_LOST` (frame can be NULL). frames
with FEC applied (49 bits, `.amb`) are decoded by `xombe_dec_frame_bits`.

to encode:

    static xombe_enc_mem mem;
    xombe_enc *enc = xombe_enc_init(&mem, NULL);    /* NULL: default settings */

    if (xombe_enc_frame(enc, pcm, NULL, dmr, NULL) == 1)    /* 160 samples in */
        send(dmr);                                  /* lags input by 3 frames */

    for (i = 0; i < xombe_enc_delay(enc); i++) {    /* at end of input */
        xombe_enc_frame(enc, NULL, NULL, dmr, NULL);
        send(dmr);
    }

for D-STAR output, set `cfg.outputs` and pass a buffer as the last argument.
every call is documented in `include/xombe.h`.

## features

 - AMBE+2 and AMBE 2400 decoders with frame-loss concealment
 - encoders with noise suppression, DTX, and tone detection
 - send DTMF, Knox and call-progress tones from the encoder
 - DMR and D-STAR channel FEC, both directions
 - transcode DMR to D-STAR and back in parameter space (no PCM)
 - no malloc or globals: all state is caller-managed
 - fixed-point build for Cortex-M with the DSP extension
 - hardened against malformed and malicious inputs

benchmarks on a raspberry pi pico 2 (RP2350):

 - 150 MHz (stock), q15 encoder uses 41% of one core, decoder 6%
 - 200 MHz (no overvolt), three full-duplex channels fit across two cores

## building

needs cmake 3.13 and a C99 compiler.

    git clone https://github.com/x0rro/xombe && cd xombe
    make
    make check # optional
    sudo cmake --install build

pass build options to `make` as `CMAKE_FLAGS="-DXOMBE_KERNELS=q15 ..."`, or run
cmake yourself.

once installed, there are two main ways of using the library.

 - cmake: `find_package(xombe)` and link `xombe::xombe`
 - `pkg-config xombe` for the compile and link flags.

the library builds on all platforms. the tools and tests require POSIX.

### build options

    -DXOMBE_KERNELS=float|q15   q15 is fixed point for M33 DSP (default float)
    -DBUILD_SHARED_LIBS=ON      shared library (default static)

    -DXOMBE_ENCODER=OFF         decoder only
    -DXOMBE_DSTAR=OFF           DMR only

float and q15 builds encode slightly differently, and their decoded audio
differs by a few LSB. both pass the same tests.

### testing

    make check                         # release, q15 and sanitizer builds
    make conformance PYTHON=python3    # scores encoded and decoded speech
    make conformance-dstar PYTHON=python3

`make check` runs a few more tests if python 3 is installed. the conformance
targets need numpy, scipy, pesq, pystoi and visqol.

## file formats

    .raw     8 kHz signed 16-bit little-endian mono
    .amb     dsd's format: ".amb", then 8 bytes per frame
    .bits    49 bytes per frame, one bit per byte
    .air     72-bit DMR air frames (xombe does FEC and interleaving)
    .r33     same as .air
    .dstar   72-bit D-STAR voice frames (xombe does FEC and interleaving)

## xombe-dec

    xombe-dec [options] IN OUT.raw

    --bits              IN is .bits
    --air               IN is .air, decoded through the FEC
    --dstar             IN is .dstar
                        (default: .amb)

    --lost LIST         treat these frames as lost, from 0, e.g. 3,10-12
    --decisions         print each frame's class and concealment on stderr
    --help              list all options

## xombe-enc

    xombe-enc [options] IN.raw OUT

    --bits              write .bits
    --air               write .air
    --dstar             write .dstar only
    --dual DSTAR_OUT    write .dstar to DSTAR_OUT as well as AMBE+2 to OUT
                        (default: .amb)

    --ns on|off         noise suppression (default on)
    --ns-strength S     how hard to suppress (0-1, default 0.5)
    --dtx on|off        send silence frames for background noise (default on)
    --dstar-dtx on|off  same for D-STAR (default off)
    --tones on|off      detect tones and send as tone frames (default on)
    --vscale S          voicing decision scale, 0.6 voices less (default 0.35)
    --quant closed|open open is faster but sounds worse (default closed)

    --help              list all options

## xombe-xc

    xombe-xc [options] IN.r33 OUT.dstar             # DMR to D-STAR
    xombe-xc --reverse [options] IN.dstar OUT.r33   # D-STAR to DMR

    --amb               DMR side is .amb instead of .r33 (for D-STAR in dsd)
    --lost LIST         treat these frames as lost, from 0, e.g. 3,10-12
    --decisions         print each frame's class and concealment on stderr

## verifying a release

each release is a tag signed by x0rro (`x0rroRF@proton.me`), key fingerprint
`638C 5D53 640A CEFB 79B6 87AA 3E93 43F0 B819 4BFF`. the tag message has the
hashes of the release tarball and its buildinfo:

    git tag -v v0.1.0

### reproducing a release

release tarballs are reproducible. rebuild one on x86_64 in the debian container
named in `scripts/release_lock.sh`:

    git clone -b v0.1.0 https://github.com/x0rro/xombe && cd xombe && mkdir out
    podman run --rm -v "$PWD:/src" -v "$PWD/out:/out" \
        docker.io/library/debian@sha256:2dc419c6d7fe3dad0a2f9fd469c34cc90b0f9ba642d8fde837fe3d7d34a0e7f1 \
        sh -c 'cd /src && scripts/release_deps.sh && git config --global safe.directory /src &&
               scripts/release_tarball.sh HEAD /out && cd /out && sha256sum *'

the hashes it prints match the tag's.

## patent note

patents covering methods used in this code may still be in force in some
jurisdictions. if that applies to you, get a license from the patent holder
before using it.

## copyright

copyright (c) 2026 x0rro

released under MIT license, see LICENSE
