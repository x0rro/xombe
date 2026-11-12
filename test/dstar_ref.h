#ifndef XOMBE_TEST_DSTAR_REF_H
#define XOMBE_TEST_DSTAR_REF_H

/* reader and writer for the D-STAR reference records, the writer for stage dumps test/dstar_cmp.py compares

   a record file is text, gzip-compressed if the path ends in .gz (read through gzip -dc):
     # xombe-dstar-ref 1          meta lines start with '#', "# key=value" lines are kept
     # stage=fec
     n:i in:x ea:i ... a:f*       column line, name:type, types i f s x b i* f*
     0 9e8d32... 0 ... 3 0.5 1 2  one record per line, whitespace separated
   i integer, f float (nan allowed), s token, x hex bytes, b string of 0/1, '-' is an empty s/x/b
   i* and f*: a count followed by that many values */

#include <stdint.h>
#include <stdio.h>

#define DSREF_MAX_COLS 64
#define DSREF_NAME_LEN 32
#define DSREF_MAX_META 16

typedef enum { DSREF_INT, DSREF_FLOAT, DSREF_STR, DSREF_HEX, DSREF_BITS, DSREF_IVEC, DSREF_FVEC } dsref_type;

typedef struct {
    char name[DSREF_NAME_LEN];
    dsref_type type;
    /* the current record's value */
    long i;
    double f;
    int len;          /* vector or token length, byte or bit count, 1 for i and f */
    const char *s;    /* s, x, and b: the token ("" for '-') */
    long *iv;
    double *fv;
    uint8_t *bytes;   /* x: decoded bytes, b: one 0/1 per bit */
    size_t cap;
} dsref_column;

typedef struct {
    const char *path;
    FILE *f;
    long pid;
    char *line;
    size_t line_cap;
    long lineno;
    long nrec;
    int failed;
    int ncol;
    dsref_column col[DSREF_MAX_COLS];
    int nmeta;
    char meta_key[DSREF_MAX_META][DSREF_NAME_LEN];
    char meta_val[DSREF_MAX_META][128];
} dsref;

/* opens path and reads the meta and column lines, returns 0 (reported) on failure */
int dsref_open(dsref *r, const char *path);
/* reads the next record: 1, 0 at the end, -1 on a malformed record (reported) */
int dsref_next(dsref *r);
/* returns 1 if the file was opened and read without error */
int dsref_close(dsref *r);
/* column index of name, or -1 */
int dsref_col(const dsref *r, const char *name);
/* column index of name with the given type, a missing column or another type is reported and fails the file (-1) */
int dsref_need(dsref *r, const char *name, dsref_type type);
/* "# key=value" meta value, or NULL */
const char *dsref_meta(const dsref *r, const char *key);

static inline long dsref_i(const dsref *r, int c) { return r->col[c].i; }
static inline double dsref_f(const dsref *r, int c) { return r->col[c].f; }
static inline int dsref_len(const dsref *r, int c) { return r->col[c].len; }
static inline const char *dsref_s(const dsref *r, int c) { return r->col[c].s; }
static inline const long *dsref_iv(const dsref *r, int c) { return r->col[c].iv; }
static inline const double *dsref_fv(const dsref *r, int c) { return r->col[c].fv; }
static inline const uint8_t *dsref_bytes(const dsref *r, int c) { return r->col[c].bytes; }

const char *dsref_type_name(dsref_type t);

/* whole binary files, PCM (int16 little endian) and 9-byte frame streams
   return malloc'd data (free() it) and the count, or NULL (reported) */
int16_t *dsref_read_pcm(const char *path, long *nsamples);
uint8_t *dsref_read_frames(const char *path, int frame_bytes, long *nframes);

typedef struct {
    FILE *f;
    int ncol, cur, prec;
    dsref_type type[DSREF_MAX_COLS];
    int failed;
} dsref_writer;

/* writes the meta and column lines, prec = significant digits of floats (9 if 0) */
int dsref_wopen(dsref_writer *w, const char *path, const char *stage, const char *cols, int prec);
void dsref_wi(dsref_writer *w, long v);
void dsref_wf(dsref_writer *w, double v);
void dsref_ws(dsref_writer *w, const char *s);
void dsref_wx(dsref_writer *w, const uint8_t *bytes, int n);
void dsref_wb(dsref_writer *w, const uint8_t *bits, int n);
void dsref_wiv(dsref_writer *w, const long *v, int n);
void dsref_wiv_int(dsref_writer *w, const int *v, int n);
void dsref_wfv(dsref_writer *w, const double *v, int n);
void dsref_wfv_float(dsref_writer *w, const float *v, int n);
/* ends the record, a wrong number of values or a type mismatch fails the writer */
void dsref_wend(dsref_writer *w);
/* returns 1 if every record was written correctly */
int dsref_wclose(dsref_writer *w);

#endif
