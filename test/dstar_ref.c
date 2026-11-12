#define _POSIX_C_SOURCE 200809L
#include "dstar_ref.h"
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

static void fail(dsref *r, const char *what)
{
    fprintf(stderr, "%s:%ld: %s\n", r->path, r->lineno, what);
    r->failed = 1;
}

static FILE *open_gzip(const char *path, long *pid_out)
{
    int fd[2];
    pid_t pid;
    FILE *f;
    if (pipe(fd))
        return NULL;
    pid = fork();
    if (pid < 0) {
        close(fd[0]);
        close(fd[1]);
        return NULL;
    }
    if (pid == 0) {
        dup2(fd[1], STDOUT_FILENO);
        close(fd[0]);
        close(fd[1]);
        execlp("gzip", "gzip", "-dc", "--", path, (char *)NULL);
        _exit(127);
    }
    close(fd[1]);
    f = fdopen(fd[0], "r");
    if (!f) {
        close(fd[0]);
        waitpid(pid, NULL, 0);
        return NULL;
    }
    *pid_out = (long)pid;
    return f;
}

static int read_line(dsref *r)
{
    ssize_t n = getline(&r->line, &r->line_cap, r->f);
    if (n < 0)
        return 0;
    r->lineno++;
    while (n > 0 && (r->line[n - 1] == '\n' || r->line[n - 1] == '\r'))
        r->line[--n] = 0;
    return 1;
}

static int parse_type(const char *s, dsref_type *t)
{
    static const char *names[] = {"i", "f", "s", "x", "b", "i*", "f*"};
    int k;
    for (k = 0; k < 7; k++)
        if (!strcmp(s, names[k])) {
            *t = (dsref_type)k;
            return 1;
        }
    return 0;
}

const char *dsref_type_name(dsref_type t)
{
    static const char *names[] = {"i", "f", "s", "x", "b", "i*", "f*"};
    return (unsigned)t < 7 ? names[t] : "?";
}

int dsref_open(dsref *r, const char *path)
{
    size_t n = strlen(path);
    memset(r, 0, sizeof(*r));
    r->path = path;
    r->pid = -1;
    if (n > 3 && !strcmp(path + n - 3, ".gz"))
        r->f = open_gzip(path, &r->pid);
    else
        r->f = fopen(path, "r");
    if (!r->f) {
        fprintf(stderr, "%s: cannot open: %s\n", path, strerror(errno));
        r->failed = 1;
        return 0;
    }
    for (;;) {
        char *tok, *save = NULL;
        if (!read_line(r)) {
            fail(r, "no column line");
            return 0;
        }
        if (r->line[0] == '#') {
            char *eq = strchr(r->line, '=');
            if (eq && r->nmeta < DSREF_MAX_META) {
                char *k = r->line + 1;
                while (*k == ' ')
                    k++;
                *eq = 0;
                snprintf(r->meta_key[r->nmeta], DSREF_NAME_LEN, "%s", k);
                snprintf(r->meta_val[r->nmeta], sizeof(r->meta_val[0]), "%s", eq + 1);
                r->nmeta++;
            }
            continue;
        }
        for (tok = strtok_r(r->line, " \t", &save); tok; tok = strtok_r(NULL, " \t", &save)) {
            char *colon = strrchr(tok, ':');
            dsref_column *c;
            if (r->ncol == DSREF_MAX_COLS || !colon || colon == tok || colon - tok >= DSREF_NAME_LEN) {
                fail(r, "bad column spec");
                return 0;
            }
            c = &r->col[r->ncol];
            memcpy(c->name, tok, (size_t)(colon - tok));
            c->name[colon - tok] = 0;
            if (!parse_type(colon + 1, &c->type)) {
                fail(r, "bad column type");
                return 0;
            }
            if (dsref_col(r, c->name) >= 0) {
                fail(r, "duplicate column");
                return 0;
            }
            r->ncol++;
        }
        if (!r->ncol) {
            fail(r, "empty column line");
            return 0;
        }
        return 1;
    }
}

int dsref_col(const dsref *r, const char *name)
{
    int k;
    for (k = 0; k < r->ncol; k++)
        if (!strcmp(r->col[k].name, name))
            return k;
    return -1;
}

int dsref_need(dsref *r, const char *name, dsref_type type)
{
    int c = dsref_col(r, name);
    if (c < 0 || r->col[c].type != type) {
        char msg[96];
        snprintf(msg, sizeof(msg), "column %s:%s %s", name, dsref_type_name(type), c < 0 ? "missing" : "has another type");
        fail(r, msg);
        return -1;
    }
    return c;
}

const char *dsref_meta(const dsref *r, const char *key)
{
    int k;
    for (k = 0; k < r->nmeta; k++)
        if (!strcmp(r->meta_key[k], key))
            return r->meta_val[k];
    return NULL;
}

static int grow(dsref_column *c, size_t n, size_t elem)
{
    void *p;
    if (n <= c->cap)
        return 1;
    if (c->type == DSREF_IVEC)
        p = realloc(c->iv, n * elem);
    else if (c->type == DSREF_FVEC)
        p = realloc(c->fv, n * elem);
    else
        p = realloc(c->bytes, n * elem);
    if (!p)
        return 0;
    if (c->type == DSREF_IVEC)
        c->iv = p;
    else if (c->type == DSREF_FVEC)
        c->fv = p;
    else
        c->bytes = p;
    c->cap = n;
    return 1;
}

static int to_long(const char *s, long *v)
{
    char *end;
    errno = 0;
    *v = strtol(s, &end, 10);
    return end != s && !*end && !errno;
}

static int to_double(const char *s, double *v)
{
    char *end;
    *v = strtod(s, &end);
    return end != s && !*end;
}

static int hexval(int ch)
{
    if (ch >= '0' && ch <= '9')
        return ch - '0';
    ch = tolower(ch);
    return ch >= 'a' && ch <= 'f' ? ch - 'a' + 10 : -1;
}

int dsref_next(dsref *r)
{
    char *save = NULL, *tok;
    int k;
    if (r->failed || !r->f)
        return -1;
    do {
        if (!read_line(r))
            return 0;
    } while (r->line[0] == '#' || r->line[0] == 0);
    tok = strtok_r(r->line, " \t", &save);
    for (k = 0; k < r->ncol; k++) {
        dsref_column *c = &r->col[k];
        long n, j;
        if (!tok) {
            fail(r, "too few values");
            return -1;
        }
        switch (c->type) {
        case DSREF_INT:
            if (!to_long(tok, &c->i)) {
                fail(r, "bad integer");
                return -1;
            }
            c->len = 1;
            break;
        case DSREF_FLOAT:
            if (!to_double(tok, &c->f)) {
                fail(r, "bad float");
                return -1;
            }
            c->len = 1;
            break;
        case DSREF_STR:
        case DSREF_HEX:
        case DSREF_BITS:
            c->s = strcmp(tok, "-") ? tok : "";
            c->len = (int)strlen(c->s);
            if (c->type == DSREF_HEX) {
                if (c->len % 2 || !grow(c, (size_t)c->len / 2 + 1, 1)) {
                    fail(r, "bad hex");
                    return -1;
                }
                for (j = 0; j < c->len / 2; j++) {
                    int hi = hexval(c->s[2 * j]), lo = hexval(c->s[2 * j + 1]);
                    if (hi < 0 || lo < 0) {
                        fail(r, "bad hex digit");
                        return -1;
                    }
                    c->bytes[j] = (uint8_t)(hi << 4 | lo);
                }
                c->len /= 2;
            } else if (c->type == DSREF_BITS) {
                if (!grow(c, (size_t)c->len + 1, 1)) {
                    fail(r, "out of memory");
                    return -1;
                }
                for (j = 0; j < c->len; j++) {
                    if (c->s[j] != '0' && c->s[j] != '1') {
                        fail(r, "bad bit string");
                        return -1;
                    }
                    c->bytes[j] = (uint8_t)(c->s[j] - '0');
                }
            }
            break;
        case DSREF_IVEC:
        case DSREF_FVEC:
            if (!to_long(tok, &n) || n < 0 || n > 100000 ||
                !grow(c, (size_t)n + 1, c->type == DSREF_IVEC ? sizeof(long) : sizeof(double))) {
                fail(r, "bad vector count");
                return -1;
            }
            for (j = 0; j < n; j++) {
                tok = strtok_r(NULL, " \t", &save);
                if (!tok || (c->type == DSREF_IVEC ? !to_long(tok, &c->iv[j]) : !to_double(tok, &c->fv[j]))) {
                    fail(r, tok ? "bad vector value" : "short vector");
                    return -1;
                }
            }
            c->len = (int)n;
            break;
        }
        tok = strtok_r(NULL, " \t", &save);
    }
    if (tok) {
        fail(r, "too many values");
        return -1;
    }
    r->nrec++;
    return 1;
}

int dsref_close(dsref *r)
{
    int k, ok = !r->failed;
    if (r->f) {
        if (ferror(r->f))
            ok = 0;
        fclose(r->f);
    }
    if (r->pid >= 0) {
        int status = 0;
        if (waitpid((pid_t)r->pid, &status, 0) < 0 || !WIFEXITED(status) || WEXITSTATUS(status)) {
            fprintf(stderr, "%s: gzip failed\n", r->path);
            ok = 0;
        }
    }
    for (k = 0; k < r->ncol; k++) {
        free(r->col[k].iv);
        free(r->col[k].fv);
        free(r->col[k].bytes);
    }
    free(r->line);
    memset(r, 0, sizeof(*r));
    r->pid = -1;
    return ok;
}

static void *read_all(const char *path, long *size)
{
    FILE *f = fopen(path, "rb");
    char *buf = NULL;
    long n = 0, cap = 0;
    size_t k;
    if (!f) {
        fprintf(stderr, "%s: cannot open: %s\n", path, strerror(errno));
        return NULL;
    }
    for (;;) {
        if (n + 65536 > cap) {
            char *p = realloc(buf, (size_t)(cap = 2 * cap + 65536));
            if (!p) {
                free(buf);
                fclose(f);
                return NULL;
            }
            buf = p;
        }
        k = fread(buf + n, 1, (size_t)(cap - n), f);
        n += (long)k;
        if (k == 0)
            break;
    }
    if (ferror(f)) {
        fprintf(stderr, "%s: read error\n", path);
        free(buf);
        buf = NULL;
    }
    fclose(f);
    *size = n;
    return buf;
}

int16_t *dsref_read_pcm(const char *path, long *nsamples)
{
    long n, i;
    uint8_t *b = read_all(path, &n);
    int16_t *x;
    if (!b)
        return NULL;
    if (n % 2) {
        fprintf(stderr, "%s: odd size %ld for int16 PCM\n", path, n);
        free(b);
        return NULL;
    }
    x = malloc((size_t)(n / 2 + 1) * sizeof(int16_t));
    if (x)
        for (i = 0; i < n / 2; i++)
            x[i] = (int16_t)(uint16_t)(b[2 * i] | b[2 * i + 1] << 8);
    free(b);
    *nsamples = n / 2;
    return x;
}

uint8_t *dsref_read_frames(const char *path, int frame_bytes, long *nframes)
{
    long n;
    uint8_t *b = read_all(path, &n);
    if (!b)
        return NULL;
    if (n % frame_bytes) {
        fprintf(stderr, "%s: size %ld is not a multiple of %d\n", path, n, frame_bytes);
        free(b);
        return NULL;
    }
    *nframes = n / frame_bytes;
    return b;
}

int dsref_wopen(dsref_writer *w, const char *path, const char *stage, const char *cols, int prec)
{
    char buf[1024], *tok, *save = NULL;
    memset(w, 0, sizeof(*w));
    w->prec = prec > 0 ? prec : 9;
    snprintf(buf, sizeof(buf), "%s", cols);
    for (tok = strtok_r(buf, " \t", &save); tok; tok = strtok_r(NULL, " \t", &save)) {
        char *colon = strrchr(tok, ':');
        if (w->ncol == DSREF_MAX_COLS || !colon || !parse_type(colon + 1, &w->type[w->ncol])) {
            fprintf(stderr, "%s: bad column spec '%s'\n", path, tok);
            return 0;
        }
        w->ncol++;
    }
    w->f = fopen(path, "w");
    if (!w->f) {
        fprintf(stderr, "%s: cannot create: %s\n", path, strerror(errno));
        return 0;
    }
    fprintf(w->f, "# xombe-dstar-ref 1\n# stage=%s\n# prec=%d\n# src=xombe\n%s\n", stage, w->prec, cols);
    return 1;
}

static int wcol(dsref_writer *w, dsref_type t)
{
    if (w->failed || !w->f)
        return 0;
    if (w->cur >= w->ncol || w->type[w->cur] != t) {
        fprintf(stderr, "dsref writer: value %d is not of type %s\n", w->cur, dsref_type_name(t));
        w->failed = 1;
        return 0;
    }
    if (w->cur++)
        fputc(' ', w->f);
    return 1;
}

static void put_f(dsref_writer *w, double v)
{
    if (v == 0.0)
        fputs("0", w->f);
    else if (isnan(v))
        fputs("nan", w->f);
    else
        fprintf(w->f, "%.*g", w->prec, v);
}

void dsref_wi(dsref_writer *w, long v)
{
    if (wcol(w, DSREF_INT))
        fprintf(w->f, "%ld", v);
}

void dsref_wf(dsref_writer *w, double v)
{
    if (wcol(w, DSREF_FLOAT))
        put_f(w, v);
}

void dsref_ws(dsref_writer *w, const char *s)
{
    if (wcol(w, DSREF_STR))
        fputs(*s ? s : "-", w->f);
}

void dsref_wx(dsref_writer *w, const uint8_t *bytes, int n)
{
    int k;
    if (!wcol(w, DSREF_HEX))
        return;
    if (!n)
        fputc('-', w->f);
    for (k = 0; k < n; k++)
        fprintf(w->f, "%02x", bytes[k]);
}

void dsref_wb(dsref_writer *w, const uint8_t *bits, int n)
{
    int k;
    if (!wcol(w, DSREF_BITS))
        return;
    if (!n)
        fputc('-', w->f);
    for (k = 0; k < n; k++)
        fputc(bits[k] ? '1' : '0', w->f);
}

void dsref_wiv(dsref_writer *w, const long *v, int n)
{
    int k;
    if (!wcol(w, DSREF_IVEC))
        return;
    fprintf(w->f, "%d", n);
    for (k = 0; k < n; k++)
        fprintf(w->f, " %ld", v[k]);
}

void dsref_wiv_int(dsref_writer *w, const int *v, int n)
{
    int k;
    if (!wcol(w, DSREF_IVEC))
        return;
    fprintf(w->f, "%d", n);
    for (k = 0; k < n; k++)
        fprintf(w->f, " %d", v[k]);
}

void dsref_wfv(dsref_writer *w, const double *v, int n)
{
    int k;
    if (!wcol(w, DSREF_FVEC))
        return;
    fprintf(w->f, "%d", n);
    for (k = 0; k < n; k++) {
        fputc(' ', w->f);
        put_f(w, v[k]);
    }
}

void dsref_wfv_float(dsref_writer *w, const float *v, int n)
{
    int k;
    if (!wcol(w, DSREF_FVEC))
        return;
    fprintf(w->f, "%d", n);
    for (k = 0; k < n; k++) {
        fputc(' ', w->f);
        put_f(w, v[k]);
    }
}

void dsref_wend(dsref_writer *w)
{
    if (!w->f || w->failed)
        return;
    if (w->cur != w->ncol) {
        fprintf(stderr, "dsref writer: record has %d values, %d columns\n", w->cur, w->ncol);
        w->failed = 1;
        return;
    }
    fputc('\n', w->f);
    w->cur = 0;
}

int dsref_wclose(dsref_writer *w)
{
    int ok = w->f && !w->failed && w->cur == 0;
    if (w->f) {
        if (ferror(w->f) || fclose(w->f))
            ok = 0;
    }
    w->f = NULL;
    return ok;
}
