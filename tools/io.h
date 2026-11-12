#ifndef XOMBE_TOOLS_IO_H
#define XOMBE_TOOLS_IO_H

/* file formats, checked I/O, and messages shared by the tools and test_stream */

#include <stdio.h>
#include <stdint.h>
#include "xombe.h"

enum io_format {
    IO_AMB,  /* ".amb" header, then 8 bytes per frame: a zero byte, d0..d47 MSB-first, d48 in the LSB */
    IO_BITS, /* 49 bytes per frame, one bit per byte */
    IO_AIR   /* 9 bytes per frame, the 72-bit frame */
};

#define IO_AMB_HEADER ".amb"
#define IO_AMB_RECORD 8
#define IO_BITS_RECORD 49
#define IO_AIR_RECORD XOMBE_FRAME_BYTES

/* bytes: d0..d48 MSB-first, d48 in the MSB of bytes[6] (the library's voice bits) */
void io_pack_bits(const uint8_t d[49], uint8_t bytes[XOMBE_BITS_BYTES]);
void io_unpack_bits(const uint8_t bytes[XOMBE_BITS_BYTES], uint8_t d[49]);
void io_amb_from_bytes(const uint8_t bytes[XOMBE_BITS_BYTES], uint8_t rec[IO_AMB_RECORD]);
void io_amb_to_bytes(const uint8_t rec[IO_AMB_RECORD], uint8_t bytes[XOMBE_BITS_BYTES]);

/* tool name for messages, each tool sets it first */
extern const char *io_prog;
/* usage error: prints "io_prog: <arg><msg>" (arg may be NULL) and exits 2 */
void io_fail(const char *msg, const char *arg);
/* names of enum xombe_class and enum xombe_conceal values, "?" out of range */
const char *io_class_name(int cls);
const char *io_conceal_name(int conceal);
/* --decisions line on stderr: "n class conceal run mute e0 e1 bf hb[ out]", mute the mute set's source frame,
   init, or - without a mute */
void io_print_decision(long n, const xombe_decision *d, const char *out);

/* "-" selects stdin or stdout, both report failures on stderr with the path and return NULL or 0 */
FILE *io_open(const char *path, const char *mode);
int io_close(FILE *f, const char *path);
int io_write(FILE *f, const void *p, size_t n, const char *path);
/* reads an n-byte record, returns 1, or 0 at the end of the input or on a read error (check ferror), or -1 for a
   partial final record, reported as an error */
int io_read_record(FILE *f, void *rec, size_t n, const char *path);
/* raw audio, signed 16-bit little-endian samples whatever the host byte order
   io_read_pcm returns the number of whole samples read */
size_t io_read_pcm(FILE *f, int16_t *pcm, size_t n);
int io_write_pcm(FILE *f, const int16_t *pcm, size_t n, const char *path);
/* parses a whole string as a number in [lo, hi], returns 0 if it is not one */
int io_parse_float(const char *s, float lo, float hi, float *out);

/* frame index set from "a,b-c,...", for --lost, io_frames_parse returns 0 for a malformed list */
typedef struct {
    long (*r)[2];
    int n;
} io_frames;
int io_frames_parse(const char *s, io_frames *out);
int io_frames_has(const io_frames *f, long i);

/* JSONL decoder streams, one {"init":...}, {"set":KEY,"value":N}, or half-frame object per line */
enum stream_event { STREAM_END, STREAM_HALF, STREAM_SET, STREAM_INIT, STREAM_ERROR };

typedef struct {
    FILE *f;
    long pid;
    const char *path;
    long lineno;
    char *line;
    /* latest "set" values, an "init" line resets them */
    int eps_t, eps_1;
    unsigned status_flags;
} stream_reader;

typedef struct {
    uint8_t bits[XOMBE_BITS_BYTES];
    uint8_t d[49];
    unsigned flags; /* the line's flags OR the latest status_flags */
    int second_half;
    int nsamples;
    const char *line; /* whole line, for further keys, valid until the next stream_next */
} stream_half;

/* a path ending in .gz is decompressed with gzip (run directly, without a shell) */
int stream_open(stream_reader *r, const char *path);
int stream_next(stream_reader *r, stream_half *h);
/* returns 0 when the input was read without error */
int stream_close(stream_reader *r);
/* text after "key": in s, or NULL */
const char *stream_key(const char *s, const char *key);
/* reads the integer after "key": in s, reports a missing or malformed value and returns 0 */
int stream_int(const stream_reader *r, const char *s, const char *key, long *out);

#endif
