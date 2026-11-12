#!/usr/bin/env python3
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (DEFAULT_DEC, DUAL, E7V, FS, MAX_JOBS, decode_many, frame_class, read_amb,
                    read_raw, val)
import argparse
import glob
import json
import math
import numpy as np

DTMF_NAMES = '0123456789ABCD*#'
LIM_LEVEL, LIM_ONSET, SEARCH = 1.0, 2, 0.05
PKG_CASE = dict(name='tones', file=os.path.join(E7V, 'tones.raw'), segments=[
    dict(start=2400, length=4800, id=32, freqs=[1000.0]),
    dict(start=9600, length=4800, id=133, freqs=[1336.0, 770.0]),
    dict(start=16800, length=4800, id=161, freqs=[480.0, 440.0])])


def gen(outdir):
    os.makedirs(outdir, exist_ok=True)
    cases = []
    specs = [('dtmf_%s' % ('star' if ch == '*' else 'hash' if ch == '#' else ch), 128 + i, DUAL[128 + i], -16.0)
             for i, ch in enumerate(DTMF_NAMES)]
    specs += [('knox_%d' % t, t, DUAL[t], -20.0) for t in range(144, 160)]
    specs += [('cp_%d' % t, t, DUAL[t], -20.0) for t in (160, 161, 162, 163)]
    specs += [('tone_%d_%d' % (f, -lev), int(round(f / 31.25)), (float(f),), float(lev))
              for f in (300, 1000, 2000, 3000) for lev in (-6, -16, -26, -36, -46)]
    rng = np.random.default_rng(20261019)
    for name, tid, freqs, lev in specs:
        n = np.arange(int(0.4 * FS))
        A = 32768 * 10 ** (lev / 20)
        tone = sum(A * np.cos(2 * np.pi * f / FS * n + rng.uniform(0, 2 * np.pi)) for f in freqs)
        x = np.concatenate([np.zeros(FS // 2), tone, np.zeros(int(0.3 * FS))])
        x = np.clip(np.round(x), -32768, 32767).astype('<i2')
        x.tofile(os.path.join(outdir, name + '.raw'))
        cases.append(dict(name=name, file=name + '.raw', segments=[
            dict(start=FS // 2, length=len(n), id=tid, freqs=[float(f) for f in freqs], level_dbfs=lev)]))
    json.dump(cases, open(os.path.join(outdir, 'cases.json'), 'w'), indent=1)
    print('%d E6 inputs written to %s' % (len(cases), outdir))


def fit_levels(seg, freqs):
    seg = np.asarray(seg, float)
    t = np.arange(len(seg))

    def fit(fs_):
        A = np.stack([g(2 * np.pi * f / FS * t) for f in fs_ for g in (np.cos, np.sin)], 1)
        c, *_ = np.linalg.lstsq(A, seg, rcond=None)
        return [math.hypot(c[2 * i], c[2 * i + 1]) for i in range(len(fs_))], float(np.sum((seg - A @ c) ** 2))
    best = list(freqs)
    for i in range(len(freqs)):
        cands = freqs[i] * (1 + np.linspace(-SEARCH, SEARCH, 401))
        errs = []
        for f in cands:
            fs_ = list(best)
            fs_[i] = f
            errs.append(fit(fs_)[1])
        best[i] = float(cands[int(np.argmin(errs))])
    amps, _ = fit(best)
    return [20 * math.log10(a / 32768 + 1e-12) for a in amps], best


def steady(x, lo, hi):
    env = np.convolve(np.asarray(x[lo:hi], float) ** 2, np.ones(160) / 160, 'same')
    if not len(env) or env.max() <= 0:
        return None
    on = np.where(env > 0.25 * env.max())[0]
    st = on[0] + 800
    if on[-1] - st <= 960:
        return None
    return lo + st, min(1600, on[-1] - st - 160)


def tone_ids(amb):
    return [val(d, range(12, 20)) if frame_class(d) == 'tone' else None for d in read_amb(amb)]


def check(cases, cand_dir, pattern, dec, cache, jobs):
    ambs = [os.path.join(cand_dir, pattern.format(name=c['name'])) for c in cases]
    for p in ambs:
        if not os.path.exists(p):
            sys.exit('missing: ' + p)
    pcm = decode_many(ambs, dec, cache, jobs)
    n_ok = n_all = 0
    rows = []
    for c, amb in zip(cases, ambs):
        x = read_raw(c['file'])
        ids = tone_ids(amb)
        y = pcm[amb]
        for s in c['segments']:
            n_all += 1
            on_f = s['start'] // 160
            end_f = min(len(ids), (s['start'] + s['length']) // 160 + 10)
            first = next((i for i in range(on_f, end_f) if ids[i] == s['id']), None)
            other = sorted({i for i in ids[on_f:end_f] if i is not None and i != s['id']})
            ok_id = first is not None and first - on_f <= LIM_ONSET and not other
            si = steady(x, s['start'] - 160, s['start'] + s['length'] + 160)
            so = steady(y, s['start'] - 160, min(len(y), s['start'] + s['length'] + 2400))
            if si and so:
                lin, _ = fit_levels(x[si[0]:si[0] + si[1]], s['freqs'])
                lout, fo = fit_levels(y[so[0]:so[0] + so[1]], s['freqs'])
                err = [o - i for o, i in zip(lout, lin)]
            else:
                lin, lout, err, fo = [None], [None], [None], []
            ok_lev = None not in err and all(abs(e) <= LIM_LEVEL for e in err)
            ok = ok_id and ok_lev
            n_ok += ok
            rows.append(dict(case=c['name'], id=s['id'], first=None if first is None else first - on_f,
                             other_ids=other, level_in=lin, level_out=lout, err=err, f_out=fo, ok=ok))
            print('%-14s ID %3d  first tone frame %s  other IDs %-8s  level in %s  err %s  %s' % (
                c['name'], s['id'], 'none' if first is None else 'onset+%d' % (first - on_f), other or '-',
                ' '.join('%6.2f' % v if v is not None else '  n/a' for v in lin),
                ' '.join('%+5.2f' % v if v is not None else ' n/a' for v in err), 'ok' if ok else 'FAIL'))
    print('E6 tone cases: %d/%d pass (ID within %d frames of onset, no other IDs, level within ±%.1f dB)'
          % (n_ok, n_all, LIM_ONSET, LIM_LEVEL))
    return n_ok == n_all, rows


def speech(paths, exclude=('tones',)):
    files = []
    for p in paths:
        files += sorted(glob.glob(os.path.join(p, '*.amb'))) if os.path.isdir(p) else [p]
    files = [f for f in files if not any(os.path.basename(f).startswith(x) for x in exclude if x)]
    bad = {}
    for f in files:
        n = sum(1 for i in tone_ids(f) if i is not None)
        if n:
            bad[f] = n
    for f, n in bad.items():
        print('  %s: %d tone frames' % (f, n))
    print('E6 speech check: %d files, %s' % (len(files), 'no tone frames: PASS' if not bad else
                                              '%d with tone frames: FAIL' % len(bad)))
    return not bad, bad


def main():
    ap = argparse.ArgumentParser(description='E6 tone frames')
    sub = ap.add_subparsers(dest='cmd', required=True)
    g = sub.add_parser('gen')
    g.add_argument('outdir')
    c = sub.add_parser('check')
    c.add_argument('cases', help="cases.json from 'gen', or 'pkg' for vectors/e7/tones.raw")
    c.add_argument('cand_dir', help="candidate directory ('ref' with 'pkg': the reference bitstream)")
    c.add_argument('--pattern', default='{name}.amb')
    c.add_argument('--decoder', default=DEFAULT_DEC)
    c.add_argument('--cache')
    c.add_argument('--jobs', type=int, default=MAX_JOBS)
    c.add_argument('--json')
    s = sub.add_parser('speech')
    s.add_argument('paths', nargs='+')
    s.add_argument('--exclude', default='tones')
    a = ap.parse_args()
    if a.cmd == 'gen':
        gen(a.outdir)
        return 0
    if a.cmd == 'speech':
        return 0 if speech(a.paths, a.exclude.split(','))[0] else 1
    if a.cases == 'pkg':
        cases = [PKG_CASE]
        if a.cand_dir == 'ref':
            a.cand_dir, a.pattern = E7V, '{name}.ref.amb'
    else:
        cases = json.load(open(a.cases))
        base = os.path.dirname(os.path.abspath(a.cases))
        for cs in cases:
            cs['file'] = os.path.join(base, cs['file'])
    ok, rows = check(cases, a.cand_dir, a.pattern, a.decoder, a.cache, a.jobs)
    if a.json:
        json.dump(dict(cases=rows, pass_=ok), open(a.json, 'w'), indent=1)
    print('E6 (tone cases): ' + ('PASS' if ok else 'FAIL'))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
