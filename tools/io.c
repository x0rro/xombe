#define _POSIX_C_SOURCE 200809L
#include "io.h"
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#ifndef _WIN32
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#endif

#define LINE_MAX_LEN (1 << 14)

const char *io_prog = "xombe";

void io_fail(const char *msg, const char *arg)
{
    fprintf(stderr, "%s: %s%s\n", io_prog, arg ? arg : "", msg);
    exit(2);
}

const char *io_class_name(int cls)
{
    static const char *const names[] = {"voice", "erasure", "silence", "tone"};
    return cls >= 0 && cls < 4 ? names[cls] : "?";
}

const char *io_conceal_name(int conceal)
{
    static const char *const names[] = {"none", "repeat", "mute"};
    return conceal >= 0 && conceal < 3 ? names[conceal] : "?";
}

void io_print_decision(long n, const xombe_decision *d, const char *out)
{
    fprintf(stderr, "%ld %s %s %d ", n, io_class_name(d->cls), io_conceal_name(d->conceal), d->run);
    if (d->conceal != XOMBE_CONCEAL_MUTE)
        fputc('-', stderr);
    else if (d->mute_src < 0)
        fputs("init", stderr);
    else
        fprintf(stderr, "%lld", (long long)d->mute_src);
    fprintf(stderr, " %d %d %d %d", d->e0, d->e1, d->bf, d->hb);
    if (out)
        fprintf(stderr, " %s", out);
    fputc('\n', stderr);
}

void io_pack_bits(const uint8_t d[49], uint8_t bytes[XOMBE_BITS_BYTES])
{
    int i;
    memset(bytes, 0, XOMBE_BITS_BYTES);
    for (i = 0; i < 49; i++)
        if (d[i])
            bytes[i >> 3] |= (uint8_t)(0x80 >> (i & 7));
}

void io_unpack_bits(const uint8_t bytes[XOMBE_BITS_BYTES], uint8_t d[49])
{
    int i;
    for (i = 0; i < 49; i++)
        d[i] = (uint8_t)((bytes[i >> 3] >> (7 - (i & 7))) & 1);
}

void io_amb_from_bytes(const uint8_t bytes[XOMBE_BITS_BYTES], uint8_t rec[IO_AMB_RECORD])
{
    rec[0] = 0;
    memcpy(rec + 1, bytes, 6);
    rec[7] = (uint8_t)(bytes[6] >> 7);
}

void io_amb_to_bytes(const uint8_t rec[IO_AMB_RECORD], uint8_t bytes[XOMBE_BITS_BYTES])
{
    memcpy(bytes, rec + 1, 6);
    bytes[6] = (uint8_t)((rec[7] & 1) << 7);
}

FILE *io_open(const char *path, const char *mode)
{
    FILE *f;
    if (!strcmp(path, "-"))
        return mode[0] == 'r' ? stdin : stdout;
    f = fopen(path, mode);
    if (!f)
        fprintf(stderr, "%s: %s\n", path, strerror(errno));
    return f;
}

int io_close(FILE *f, const char *path)
{
    int bad = ferror(f);
    if (f == stdin || f == stdout)
        bad |= fflush(f);
    else
        bad |= fclose(f);
    if (bad)
        fprintf(stderr, "%s: I/O error\n", path);
    return !bad;
}

int io_write(FILE *f, const void *p, size_t n, const char *path)
{
    if (fwrite(p, 1, n, f) == n)
        return 1;
    fprintf(stderr, "%s: write failed: %s\n", path, strerror(errno));
    return 0;
}

int io_read_record(FILE *f, void *rec, size_t n, const char *path)
{
    size_t got = fread(rec, 1, n, f);
    if (got == n)
        return 1;
    if (got && !ferror(f)) {
        fprintf(stderr, "%s: error: partial final frame\n", path);
        return -1;
    }
    return 0;
}

size_t io_read_pcm(FILE *f, int16_t *pcm, size_t n)
{
    uint8_t b[2 * 256];
    size_t done = 0;
    while (done < n) {
        size_t want = n - done < 256 ? n - done : 256, got = fread(b, 2, want, f), i;
        for (i = 0; i < got; i++)
            pcm[done + i] = (int16_t)(uint16_t)(b[2 * i] | (b[2 * i + 1] << 8));
        done += got;
        if (got < want)
            break;
    }
    return done;
}

int io_write_pcm(FILE *f, const int16_t *pcm, size_t n, const char *path)
{
    uint8_t b[2 * 256];
    while (n) {
        size_t k = n < 256 ? n : 256, i;
        for (i = 0; i < k; i++) {
            uint16_t u = (uint16_t)pcm[i];
            b[2 * i] = (uint8_t)(u & 0xFF);
            b[2 * i + 1] = (uint8_t)(u >> 8);
        }
        if (!io_write(f, b, 2 * k, path))
            return 0;
        pcm += k;
        n -= k;
    }
    return 1;
}

int io_parse_float(const char *s, float lo, float hi, float *out)
{
    char *end;
    double v;
    errno = 0;
    v = strtod(s, &end);
    if (end == s || *end || errno || !(v >= lo && v <= hi))
        return 0;
    *out = (float)v;
    return 1;
}

int io_frames_parse(const char *s, io_frames *out)
{
    const char *p = s;
    int n = 1;
    for (; *p; p++)
        n += *p == ',';
    out->n = 0;
    out->r = malloc((size_t)n * sizeof(*out->r));
    if (!out->r)
        return 0;
    for (p = s;;) {
        char *end;
        long a, b;
        errno = 0;
        a = strtol(p, &end, 10);
        if (end == p || a < 0 || errno)
            break;
        b = a;
        if (*end == '-') {
            p = end + 1;
            b = strtol(p, &end, 10);
            if (end == p || b < a || errno)
                break;
        }
        out->r[out->n][0] = a;
        out->r[out->n][1] = b;
        out->n++;
        if (!*end)
            return 1;
        if (*end != ',')
            break;
        p = end + 1;
    }
    free(out->r);
    out->r = NULL;
    out->n = 0;
    return 0;
}

int io_frames_has(const io_frames *f, long i)
{
    int k;
    for (k = 0; k < f->n; k++)
        if (i >= f->r[k][0] && i <= f->r[k][1])
            return 1;
    return 0;
}

#ifndef _WIN32
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
#endif

int stream_open(stream_reader *r, const char *path)
{
    size_t n = strlen(path);
    memset(r, 0, sizeof(*r));
    r->path = path;
    r->pid = -1;
    r->line = malloc(LINE_MAX_LEN);
    if (!r->line)
        return 0;
    if (n > 3 && !strcmp(path + n - 3, ".gz")) {
#ifdef _WIN32
        fprintf(stderr, "%s: .gz streams are unsupported on Windows\n", path);
#else
        r->f = open_gzip(path, &r->pid);
        if (!r->f)
            fprintf(stderr, "%s: cannot run gzip: %s\n", path, strerror(errno));
#endif
    } else {
        r->f = io_open(path, "r");
    }
    if (!r->f) {
        free(r->line);
        r->line = NULL;
        return 0;
    }
    return 1;
}

int stream_close(stream_reader *r)
{
    int ok = 1;
    if (r->f) {
        ok = !ferror(r->f);
        if (r->f != stdin)
            fclose(r->f);
    }
#ifndef _WIN32
    if (r->pid >= 0) {
        int status = 0;
        if (waitpid((pid_t)r->pid, &status, 0) < 0 || !WIFEXITED(status) || WEXITSTATUS(status)) {
            fprintf(stderr, "%s: gzip failed\n", r->path);
            ok = 0;
        }
    }
#endif
    free(r->line);
    r->line = NULL;
    r->f = NULL;
    return ok;
}

const char *stream_key(const char *s, const char *key)
{
    char pat[64];
    const char *p;
    snprintf(pat, sizeof(pat), "\"%s\":", key);
    p = strstr(s, pat);
    return p ? p + strlen(pat) : NULL;
}

int stream_int(const stream_reader *r, const char *s, const char *key, long *out)
{
    const char *p = stream_key(s, key);
    char *end;
    if (p) {
        errno = 0;
        *out = strtol(p, &end, 10);
        if (end != p && !errno && (*end == ',' || *end == '}' || *end == ' '))
            return 1;
    }
    fprintf(stderr, "%s:%ld: %s \"%s\"\n", r->path, r->lineno, p ? "malformed" : "missing", key);
    return 0;
}

static int parse_bits(const stream_reader *r, const char *p, stream_half *h)
{
    int i;
    p += strspn(p, " ");
    if (*p++ != '"')
        goto bad;
    for (i = 0; i < 49; i++) {
        if (p[i] != '0' && p[i] != '1')
            goto bad;
        h->d[i] = (uint8_t)(p[i] == '1');
    }
    if (p[49] != '"')
        goto bad;
    io_pack_bits(h->d, h->bits);
    return 1;
bad:
    fprintf(stderr, "%s:%ld: \"bits\" is not a string of 49 binary digits\n", r->path, r->lineno);
    return 0;
}

int stream_next(stream_reader *r, stream_half *h)
{
    const char *p;
    long v;
    for (;;) {
        size_t n;
        if (!fgets(r->line, LINE_MAX_LEN, r->f))
            return ferror(r->f) ? STREAM_ERROR : STREAM_END;
        r->lineno++;
        n = strlen(r->line);
        if (n == LINE_MAX_LEN - 1 && r->line[n - 1] != '\n') {
            fprintf(stderr, "%s:%ld: line too long\n", r->path, r->lineno);
            return STREAM_ERROR;
        }
        if (r->line[strspn(r->line, " \t\r\n")])
            break;
    }
    if ((p = stream_key(r->line, "set"))) {
        p += strspn(p, " ");
        if (!stream_int(r, r->line, "value", &v))
            return STREAM_ERROR;
        if (!strncmp(p, "\"eps_T\"", 7))
            r->eps_t = (int)v;
        else if (!strncmp(p, "\"eps_1\"", 7))
            r->eps_1 = (int)v;
        else if (!strncmp(p, "\"status_flags\"", 14))
            r->status_flags = (unsigned)v;
        else {
            fprintf(stderr, "%s:%ld: unknown \"set\" key\n", r->path, r->lineno);
            return STREAM_ERROR;
        }
        return STREAM_SET;
    }
    if (stream_key(r->line, "init")) {
        r->eps_t = r->eps_1 = 0;
        r->status_flags = 0;
        return STREAM_INIT;
    }
    if ((p = stream_key(r->line, "bits"))) {
        long second, ns, flags;
        if (!parse_bits(r, p, h) || !stream_int(r, r->line, "second_half", &second) ||
            !stream_int(r, r->line, "nsamples", &ns) || !stream_int(r, r->line, "flags", &flags))
            return STREAM_ERROR;
        h->second_half = second != 0;
        h->nsamples = (int)ns;
        h->flags = (unsigned)flags | r->status_flags;
        h->line = r->line;
        return STREAM_HALF;
    }
    fprintf(stderr, "%s:%ld: unrecognised line\n", r->path, r->lineno);
    return STREAM_ERROR;
}
