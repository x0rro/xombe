# changelog

## 0.1.0

first release

- decoder for AMBE 2400 (D-STAR) and AMBE+2 half-rate (DMR)
- encoder for D-STAR and DMR
- single encoder produces both outputs
- one API for both codecs, picked by `enum xombe_codec`
- transcoding between codecs in parameter space (without converting to PCM)
- q15 kernel for Cortex-M33 with DSP extension (RP2350)
- conformance and test suites
