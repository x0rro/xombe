#!/usr/bin/env python3
"""compares a C stage dump against a D-STAR reference file, standard library only

records:
  dstar_cmp.py REF TEST [--key n,h] [--cols a,b,...] [--skip c,...] [--tol COL=ABS[:REL]]... [--min COL=RATE]...
               [--rows N] [--run CMD ARGS...]
rows are matched on the key columns (default n, plus h when both files have it, else row order)
compared columns are those in both files (or --cols), minus --skip. integers, tokens, hex and bit strings,
and integer vectors must be equal, floats and float vectors agree when |test - ref| <= ABS + REL*|ref| for
every element (default 1e-6:1e-6, --tol '*=ABS:REL' changes the default), nan only matches nan
a column passes when its agreement rate is >= its --min (default 1)

PCM:
  dstar_cmp.py --pcm REF.raw TEST.raw [--lsb N] [--max-rms R] [--min-within RATE] [--samples A:B]
               [--run CMD ARGS...]
int16 little endian (--samples: only samples A..B-1 of both, B empty for the end), reports the maximum and
rms difference in LSB and the fraction of samples within --lsb (default 1), fails when the rms difference
exceeds --max-rms (default 0) or the fraction is below --min-within (default 1)

--run (must come last) runs the producer first, e.g. the C dump program that writes TEST
exit status: 0 agree, 1 disagree, 2 usage or I/O error
"""
import gzip
import math
import os
import struct
import subprocess
import sys


def read(path):
    """-> (meta, [(name, type)], [row dict])"""
    with open(path, 'rb') as fh:
        gz = fh.read(2) == b'\x1f\x8b'
    op = gzip.open if gz else open
    meta, cols, rows = {}, None, []
    with op(path, 'rt') as fh:
        for ln, line in enumerate(fh, 1):
            line = line.rstrip('\n')
            if line.startswith('#'):
                if '=' in line:
                    k, v = line[1:].strip().split('=', 1)
                    meta[k] = v
                continue
            tok = line.split()
            if not tok:
                continue
            if cols is None:
                cols = [tuple(t.rsplit(':', 1)) for t in tok]
                continue
            r, i = {}, 0
            try:
                for name, t in cols:
                    if t == 'i':
                        r[name] = int(tok[i]); i += 1
                    elif t == 'f':
                        r[name] = float(tok[i]); i += 1
                    elif t in ('s', 'x', 'b'):
                        r[name] = '' if tok[i] == '-' else tok[i]; i += 1
                    elif t in ('i*', 'f*'):
                        n = int(tok[i]); i += 1
                        conv = int if t == 'i*' else float
                        r[name] = [conv(x) for x in tok[i:i + n]]
                        if len(r[name]) != n:
                            raise ValueError('short vector ' + name)
                        i += n
                    else:
                        raise ValueError('bad type ' + t)
            except (IndexError, ValueError) as e:
                raise ValueError('%s:%d: %s' % (path, ln, e))
            if i != len(tok):
                raise ValueError('%s:%d: %d values, %d used' % (path, ln, len(tok), i))
            rows.append(r)
    if cols is None:
        raise ValueError('%s: no column line' % path)
    return meta, cols, rows


def fclose(a, b, ab, rel):
    if math.isnan(a) or math.isnan(b):
        return math.isnan(a) and math.isnan(b), 0.0
    d = abs(a - b)
    return d <= ab + rel * abs(a), d


def cmp_records(ref, test, opt):
    _, rc, rr = read(ref)
    _, tc, tr = read(test)
    rtypes, ttypes = dict(rc), dict(tc)
    key = opt['key'] or (['n', 'h'] if 'h' in rtypes and 'h' in ttypes else ['n'])
    key = [k for k in key if k in rtypes and k in ttypes]
    cols = opt['cols'] or [c for c, _ in rc if c in ttypes]
    cols = [c for c in cols if c not in opt['skip'] and c not in key]
    for c in cols:
        if c not in rtypes or c not in ttypes:
            print('column %s missing from %s' % (c, 'REF' if c not in rtypes else 'TEST'))
            return 2
        if rtypes[c] != ttypes[c]:
            print('column %s: type %s in REF, %s in TEST' % (c, rtypes[c], ttypes[c]))
            return 2
    if key:
        tindex = {tuple(r[k] for k in key): r for r in tr}
        pairs = [(r, tindex.get(tuple(r[k] for k in key))) for r in rr]
    else:
        pairs = list(zip(rr, tr)) + [(r, None) for r in rr[len(tr):]]
    if opt['rows'] is not None:
        pairs = pairs[:opt['rows']]
    missing = sum(1 for _, t in pairs if t is None)
    print('%s vs %s: %d reference rows, %d test rows, %d unmatched; key %s' % (
        os.path.basename(ref), os.path.basename(test), len(rr), len(tr), missing, ','.join(key) or 'row order'))
    status = 0 if pairs else 1
    for c in cols:
        t = rtypes[c]
        ab, rel = opt['tol'].get(c, opt['tol'].get('*', (1e-6, 1e-6)))
        agree, total, worst, wrow = 0, 0, 0.0, None
        for r, s in pairs:
            total += 1
            if s is None:
                continue
            a, b = r[c], s[c]
            if t == 'f':
                ok, d = fclose(a, b, ab, rel)
            elif t == 'f*':
                ok, d = len(a) == len(b), (math.inf if len(a) != len(b) else 0.0)
                if ok:
                    for x, y in zip(a, b):
                        o, e = fclose(x, y, ab, rel)
                        ok, d = ok and o, max(d, e)
            else:
                ok, d = a == b, 0.0 if a == b else math.inf
            agree += ok
            if d > worst or (not ok and wrow is None):
                worst, wrow = max(d, worst), tuple(r[k] for k in key) if key else None
        rate = agree / total if total else 0.0
        need = opt['min'].get(c, opt['min'].get('*', 1.0))
        flag = 'ok' if rate >= need else 'FAIL'
        if flag == 'FAIL':
            status = 1
        tol = (' tol %g:%g' % (ab, rel)) if t in ('f', 'f*') else ''
        wtxt = ''
        if worst:
            wtxt = ' worst %s at %s' % ('mismatch' if worst == math.inf else '%.3g' % worst, wrow)
        print('  %-10s %-2s %6d/%-6d %7.3f%%%s%s  %s' % (c, t, agree, total, 100.0 * rate, tol, wtxt, flag))
    return status


def read_pcm(path):
    b = open(path, 'rb').read()
    if len(b) % 2:
        raise ValueError('%s: odd size' % path)
    return struct.unpack('<%dh' % (len(b) // 2), b)


def cmp_pcm(ref, test, opt):
    x, y = read_pcm(ref), read_pcm(test)
    if opt['samples']:
        a, b = opt['samples']
        x, y = x[a:b], y[a:b]
    n = min(len(x), len(y))
    if len(x) != len(y):
        print('length: REF %d samples, TEST %d' % (len(x), len(y)))
    d = [abs(a - b) for a, b in zip(x[:n], y[:n])]
    rms = math.sqrt(sum(v * v for v in d) / n) if n else 0.0
    sig = math.sqrt(sum(v * v for v in x[:n]) / n) if n else 0.0
    within = sum(1 for v in d if v <= opt['lsb']) / n if n else 0.0
    snr = 20 * math.log10(sig / rms) if rms > 0 and sig > 0 else math.inf
    print('%s vs %s: %d samples, max diff %d LSB, rms diff %.3f LSB, %.3f%% within %d LSB, SNR %.1f dB' % (
        os.path.basename(ref), os.path.basename(test), n, max(d) if d else 0, rms, 100 * within, opt['lsb'], snr))
    ok = len(x) == len(y) and rms <= opt['max_rms'] and within >= opt['min_within']
    print('  %s' % ('ok' if ok else 'FAIL'))
    return 0 if ok else 1


def parse(argv):
    opt = dict(key=None, cols=None, skip=set(), tol={}, min={}, rows=None, pcm=False, lsb=1, max_rms=0.0,
               min_within=1.0, run=None, samples=None)
    pos, i = [], 0
    while i < len(argv):
        a = argv[i]
        if a == '--run':
            opt['run'] = argv[i + 1:]
            break
        if a == '--pcm':
            opt['pcm'] = True
        elif a in ('--key', '--cols', '--skip', '--tol', '--min', '--rows', '--lsb', '--max-rms', '--min-within',
                   '--samples'):
            i += 1
            v = argv[i]
            if a == '--key':
                opt['key'] = [s for s in v.split(',') if s]
            elif a == '--cols':
                opt['cols'] = [s for s in v.split(',') if s]
            elif a == '--skip':
                opt['skip'] |= set(v.split(','))
            elif a == '--tol':
                c, t = v.split('=')
                ab, _, rel = t.partition(':')
                opt['tol'][c] = (float(ab), float(rel or 0.0))
            elif a == '--min':
                c, t = v.split('=')
                opt['min'][c] = float(t)
            elif a == '--rows':
                opt['rows'] = int(v)
            elif a == '--samples':
                a0, b0 = v.split(':')
                opt['samples'] = (int(a0), int(b0) if b0 else None)
            elif a == '--lsb':
                opt['lsb'] = int(v)
            elif a == '--max-rms':
                opt['max_rms'] = float(v)
            else:
                opt['min_within'] = float(v)
        elif a.startswith('--'):
            raise SystemExit('unknown option %s' % a)
        else:
            pos.append(a)
        i += 1
    if len(pos) != 2:
        raise SystemExit(__doc__)
    return pos, opt


def main(argv):
    try:
        (ref, test), opt = parse(argv)
    except SystemExit as e:
        print(e, file=sys.stderr)
        return 2
    if opt['run']:
        r = subprocess.run(opt['run'])
        if r.returncode:
            print('producer failed (%d): %s' % (r.returncode, ' '.join(opt['run'])))
            return 2
    try:
        return cmp_pcm(ref, test, opt) if opt['pcm'] else cmp_records(ref, test, opt)
    except (OSError, ValueError) as e:
        print(e)
        return 2


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
