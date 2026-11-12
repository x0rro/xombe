/* checks a D-STAR reference set is complete and well formed: every file in MANIFEST
   exists with its size, every record file parses, and the per-clip files agree on the frame counts
   tests the harness and the references, not the codec

     test_dstar_refs DIR */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dstar_ref.h"

#define MAX_FILES 512

static char names[MAX_FILES][128];
static int nfiles;
static int failures;

static void bad(const char *what, const char *name)
{
    fprintf(stderr, "FAIL %s: %s\n", name, what);
    failures++;
}

static int has(const char *name)
{
    int k;
    for (k = 0; k < nfiles; k++)
        if (!strcmp(names[k], name))
            return 1;
    return 0;
}

static int ends(const char *s, const char *suffix)
{
    size_t n = strlen(s), m = strlen(suffix);
    return n >= m && !strcmp(s + n - m, suffix);
}

static long file_size(const char *dir, const char *name)
{
    char path[1024];
    FILE *f;
    long n;
    if (snprintf(path, sizeof(path), "%s/%s", dir, name) >= (int)sizeof(path))
        return -1;
    f = fopen(path, "rb");
    if (!f)
        return -1;
    fseek(f, 0, SEEK_END);
    n = ftell(f);
    fclose(f);
    return n;
}

/* parses every record, returns the record count, or -1 */
static long records(const char *dir, const char *name)
{
    char path[1024];
    dsref r;
    long n = 0;
    int st;
    if (snprintf(path, sizeof(path), "%s/%s", dir, name) >= (int)sizeof(path) || !dsref_open(&r, path))
        return -1;
    if (!dsref_meta(&r, "stage"))
        fprintf(stderr, "%s: no stage meta line\n", path), r.failed = 1;
    while ((st = dsref_next(&r)) > 0)
        n++;
    if (!dsref_close(&r) || st < 0)
        return -1;
    return n;
}

static long count_of(const char *dir, const char *clip, const char *suffix, long unit)
{
    char name[160];
    long n;
    if (snprintf(name, sizeof(name), "%s%s", clip, suffix) >= (int)sizeof(name) || !has(name))
        return -2;
    if (ends(name, ".txt.gz"))
        return records(dir, name);
    n = file_size(dir, name);
    if (n < 0 || n % unit) {
        bad("size is not a multiple of the unit", name);
        return -1;
    }
    return n / unit;
}

static void expect(const char *dir, const char *clip, const char *suffix, long unit, long want)
{
    long n = count_of(dir, clip, suffix, unit);
    if (n == -2)
        return;
    if (n != want) {
        char msg[96], name[160];
        if (snprintf(name, sizeof(name), "%s%s", clip, suffix) >= (int)sizeof(name))
            name[0] = 0;
        snprintf(msg, sizeof(msg), "%ld units, expected %ld", n, want);
        bad(msg, name);
    }
}

/* writes a record of every type with the C writer and reads it back */
static void selftest(void)
{
    const char *path = "dstar_ref_selftest.txt";
    static const uint8_t bytes[3] = {0x9e, 0x8d, 0xff}, bits[5] = {1, 0, 1, 1, 0};
    static const double fv[3] = {-15.0, 0.125, 1e-7};
    static const long iv[2] = {-3, 123456};
    dsref_writer w;
    dsref r;
    int ok = dsref_wopen(&w, path, "selftest", "n:i x:f s:s h:x b:b iv:i* fv:f* e:f*", 0);
    dsref_wi(&w, 7);
    dsref_wf(&w, 0.1);
    dsref_ws(&w, "voice");
    dsref_wx(&w, bytes, 3);
    dsref_wb(&w, bits, 5);
    dsref_wiv(&w, iv, 2);
    dsref_wfv(&w, fv, 3);
    dsref_wfv(&w, fv, 0);
    dsref_wend(&w);
    ok = dsref_wclose(&w) && ok;
    if (ok && dsref_open(&r, path)) {
        int cn = dsref_need(&r, "n", DSREF_INT), cx = dsref_need(&r, "x", DSREF_FLOAT);
        int cs = dsref_need(&r, "s", DSREF_STR), ch = dsref_need(&r, "h", DSREF_HEX);
        int cb = dsref_need(&r, "b", DSREF_BITS), ci = dsref_need(&r, "iv", DSREF_IVEC);
        int cf = dsref_need(&r, "fv", DSREF_FVEC), ce = dsref_need(&r, "e", DSREF_FVEC);
        ok = cn >= 0 && cx >= 0 && cs >= 0 && ch >= 0 && cb >= 0 && ci >= 0 && cf >= 0 && ce >= 0 &&
             dsref_next(&r) == 1 && dsref_i(&r, cn) == 7 && dsref_f(&r, cx) == 0.1 &&
             !strcmp(dsref_s(&r, cs), "voice") && dsref_len(&r, ch) == 3 && !memcmp(dsref_bytes(&r, ch), bytes, 3) &&
             dsref_len(&r, cb) == 5 && !memcmp(dsref_bytes(&r, cb), bits, 5) && dsref_len(&r, ci) == 2 &&
             dsref_iv(&r, ci)[1] == 123456 && dsref_len(&r, cf) == 3 && dsref_fv(&r, cf)[0] == -15.0 &&
             dsref_fv(&r, cf)[2] == 1e-7 && dsref_len(&r, ce) == 0 && dsref_next(&r) == 0 &&
             !strcmp(dsref_meta(&r, "stage"), "selftest");
        ok = dsref_close(&r) && ok;
    } else {
        ok = 0;
    }
    if (!ok)
        bad("writer/reader round trip", path);
    remove(path);
}

int main(int argc, char **argv)
{
    const char *dir;
    char path[1024], line[512];
    FILE *mf;
    int k, nrecfiles = 0, nclips = 0;
    if (argc != 2) {
        fprintf(stderr, "usage: test_dstar_refs DIR\n");
        return 2;
    }
    dir = argv[1];
    selftest();
    if (snprintf(path, sizeof(path), "%s/MANIFEST", dir) >= (int)sizeof(path))
        return 2;
    mf = fopen(path, "r");
    if (!mf) {
        fprintf(stderr, "%s: cannot open\n", path);
        return 1;
    }
    while (fgets(line, sizeof(line), mf)) {
        char sha[80], name[128];
        long size;
        if (line[0] == '#')
            continue;
        if (sscanf(line, "%79s %ld %127s", sha, &size, name) != 3 || nfiles == MAX_FILES) {
            bad("malformed MANIFEST line", line);
            continue;
        }
        snprintf(names[nfiles++], sizeof(names[0]), "%s", name);
        if (file_size(dir, name) != size)
            bad("missing or size differs from MANIFEST", name);
        if (ends(name, ".txt.gz")) {
            nrecfiles++;
            if (records(dir, name) <= 0)
                bad("unreadable or empty record file", name);
        }
    }
    fclose(mf);

    /* per-clip consistency, keyed by the stream files */
    for (k = 0; k < nfiles; k++) {
        char clip[128];
        size_t n = strlen(names[k]);
        if (ends(names[k], ".dstar") && !ends(names[k], ".xc.dstar")) {
            long nf;
            memcpy(clip, names[k], n - 6);
            clip[n - 6] = 0;
            nf = count_of(dir, clip, ".dstar", 9);
            if (nf <= 0) {
                bad("empty D-Star stream", names[k]);
                continue;
            }
            nclips++;
            expect(dir, clip, ".fec.txt.gz", 1, nf);
            expect(dir, clip, ".deq.txt.gz", 1, nf);
            expect(dir, clip, ".dec.txt.gz", 1, nf);
            expect(dir, clip, ".half.txt.gz", 1, 2 * nf);
            expect(dir, clip, ".dec.raw", 2 * 160, nf);
            expect(dir, clip, ".xrev.txt.gz", 1, nf);
            expect(dir, clip, ".xrev.r33", 9, nf);
            expect(dir, clip, ".enc.txt.gz", 1, nf);
            if (has(strcat(strcpy(path, clip), ".in.raw"))) {
                long ns = file_size(dir, path) / 2;
                if ((ns + 159) / 160 != nf)
                    bad("encoder frame count is not ceil(samples / 160)", path);
            }
        } else if (ends(names[k], ".r33") && !ends(names[k], ".xrev.r33")) {
            long nf;
            memcpy(clip, names[k], n - 4);
            clip[n - 4] = 0;
            nf = count_of(dir, clip, ".r33", 9);
            nclips++;
            expect(dir, clip, ".xc.txt.gz", 1, nf);
            expect(dir, clip, ".xc.dstar", 9, nf);
        }
    }
    printf("%d files, %d record files, %d streams checked, %d failures\n", nfiles, nrecfiles, nclips, failures);
    if (!nrecfiles || !nclips)
        return 1;
    return failures ? 1 : 0;
}
