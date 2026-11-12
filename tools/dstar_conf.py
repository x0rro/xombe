#!/usr/bin/env python3
"""D-STAR conformance checks (make conformance-dstar) on in-repo data

usage: dstar_conf.py refs BUILD [--json OUT]       stage comparisons against vectors/dstar (ctest dstar_*)
       dstar_conf.py speech BUILD OUTDIR [--json]  encode/decode and both transcoders on vectors/, scored
       dstar_conf.py tones BUILD OUTDIR [--json]   D-STAR tone frames on the E6 tone inputs, and transcoded
speech and tones read the dithered E3 inputs from --inputs DIR (scripts/gen_inputs.py, default BUILD/e3in)
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import CLIPS, E5V, FS, MAX_JOBS, read_raw
import argparse
import json
import re
import subprocess
from concurrent.futures import ProcessPoolExecutor

# speech: per path, the mean over the clips and the worst clip. limits are the measured values (the lower of
# float and q15) less 0.05 on mean PESQ and ViSQOL, 0.10 on the worst clip, 0.01 STOI, and 0.10 dB LSD. the scores
# against the input are low for every path (AMBE+2 itself: PESQ 2.45), */src rows score a transcoded decode
# against the source codec's own decode
SPEECH_LIMITS = {
    # path: (mean PESQ, worst PESQ, mean ViSQOL, worst ViSQOL, mean STOI, mean LSD max)
    'ambe2': None,  # reference row
    'dstar': (2.39, 1.72, 3.14, 1.18, 0.81, 7.71),
    'dual': (2.38, 1.70, 3.13, 1.14, 0.81, 7.71),
    'a2ds': (2.25, 1.69, 2.80, 1.24, 0.74, 8.29),
    'dsa2': (2.21, 1.64, 2.88, 1.07, 0.77, 8.10),
    'a2ds/src': (4.14, 3.82, 3.87, 3.23, 0.95, 2.74),
    'dsa2/src': (3.99, 3.68, 4.04, 3.86, 0.94, 3.13),
}
# tones: frames from onset to the first tone frame and the level limit as E6. tones without a D-STAR code play as
# voice frames, their level limits are the measured worst case (float and q15) plus about 1 dB
TONE_ONSET, TONE_LEVEL_DB = 2, 1.0  # measured onset+0, 0.6 dB
VOICE_LEVEL_DB = 7.0     # measured 6.2 (knox_158)
XC_KNOX_LEVEL_DB = 10.5  # measured 9.7 (knox_146)
KNOX = range(144, 160)


def run(cmd, **kw):
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if r.returncode:
        sys.exit('failed (%d): %s\n%s' % (r.returncode, ' '.join(cmd), r.stderr))
    return r


def refs(build):
    """runs ctest's dstar_ tests verbosely and collects each comparison's agreement"""
    r = subprocess.run(['ctest', '-V', '-R', '^dstar_', '-j', str(os.cpu_count() or 4)], cwd=build,
                       capture_output=True, text=True)
    out = {}
    col = re.compile(r'^(\d+):   (\S+)\s+\S+\s+(\d+)/(\d+)\s+([\d.]+)%')
    pcm = re.compile(r'^(\d+): .* (\d+) samples, max diff (\d+) LSB, rms diff ([\d.]+) LSB, ([\d.]+)% within (\d+) '
                     r'LSB, SNR ([-\d.inf]+) dB')
    names, res = {}, {}
    for line in r.stdout.splitlines():
        m = re.match(r'^\s*\d+/\d+ Test\s+#(\d+): (\S+) \.+\s*(\*{3})?(\w+)', line)
        if m:
            names[m.group(1)] = m.group(2)
            res[m.group(1)] = m.group(4) == 'Passed'
            continue
        m = col.match(line)
        if m:
            t = out.setdefault(m.group(1), {'cols': {}})
            t['cols'][m.group(2)] = (int(m.group(3)), int(m.group(4)))
            continue
        m = pcm.match(line)
        if m:
            out.setdefault(m.group(1), {})['pcm'] = dict(samples=int(m.group(2)), max=int(m.group(3)),
                                                         rms=float(m.group(4)), within=float(m.group(5)),
                                                         lsb=int(m.group(6)), snr=m.group(7))
    rows = []
    for k in sorted(names, key=lambda k: names[k]):
        n = names[k]
        t = out.get(k, {})
        row = dict(test=n, passed=res[k])
        if t.get('cols'):
            cols = t['cols']
            worst = min(cols, key=lambda c: cols[c][0] / max(1, cols[c][1]))
            row['records'] = max(v[1] for v in cols.values())
            row['worst'] = (worst, 100.0 * cols[worst][0] / max(1, cols[worst][1]))
            bc = [c for c in ('out', 'enc', 'd49', 'b') if c in cols]
            if bc:
                row['bytes'] = (bc[0], 100.0 * cols[bc[0]][0] / max(1, cols[bc[0]][1]))
                row['bytes_n'] = cols[bc[0]]
        if t.get('pcm'):
            row['pcm'] = t['pcm']
        rows.append(row)
    print('%-34s %7s  %-24s %-18s %s' % ('test', 'records', 'lowest column', 'frames/bytes', ''))
    for w in rows:
        if 'pcm' in w:
            p = w['pcm']
            desc = '%d samples, rms %.2f LSB, %.2f%% within %d LSB, SNR %s dB' % (p['samples'], p['rms'],
                                                                                  p['within'], p['lsb'], p['snr'])
            print('%-34s %7s  %-43s %s' % (w['test'], '', desc, 'ok' if w['passed'] else 'FAIL'))
        elif 'records' in w:
            print('%-34s %7d  %-10s %8.3f%%     %-18s %s' % (
                w['test'], w['records'], w['worst'][0], w['worst'][1],
                '%s %.3f%%' % w['bytes'] if 'bytes' in w else '-', 'ok' if w['passed'] else 'FAIL'))
        else:
            print('%-34s %7s  %-43s %s' % (w['test'], '', 'exact file comparison', 'ok' if w['passed'] else 'FAIL'))
    npass = sum(w['passed'] for w in rows)
    agree = sum(w['bytes_n'][0] for w in rows if 'bytes_n' in w and 'selfcheck' not in w['test'])
    total = sum(w['bytes_n'][1] for w in rows if 'bytes_n' in w and 'selfcheck' not in w['test'])
    print('frames or records with equal output bytes/fields over the stage tests: %d/%d (%.3f%%)'
          % (agree, total, 100.0 * agree / max(1, total)))
    print('D-STAR refs: %d/%d tests pass' % (npass, len(rows)))
    ok = r.returncode == 0 and npass == len(rows) and rows
    return bool(ok), rows


def tool(build, name):
    return os.path.join(build, name)


def make_paths(build, raw, od):
    """encodes, transcodes, and decodes one clip on every path, returns {path: (decoded, source-codec decoded)}"""
    b = os.path.splitext(os.path.basename(raw))[0]
    p = lambda ext: os.path.join(od, b + ext)
    E, D, X = tool(build, 'xombe-enc'), tool(build, 'xombe-dec'), tool(build, 'xombe-xc')
    run([E, '--air', raw, p('.r33')])
    run([D, '--conformance-noise', '--air', p('.r33'), p('.a2.raw')])
    run([E, '--dstar', raw, p('.dstar')])
    run([D, '--conformance-noise', '--dstar', p('.dstar'), p('.ds.raw')])
    run([E, '--dual', p('.dual.dstar'), raw, p('.dual.amb')])
    run([D, '--conformance-noise', '--dstar', p('.dual.dstar'), p('.dual.raw')])
    run([X, p('.r33'), p('.a2ds.dstar')])
    run([D, '--conformance-noise', '--dstar', p('.a2ds.dstar'), p('.a2ds.raw')])
    run([X, '--reverse', p('.dstar'), p('.dsa2.r33')])
    run([D, '--conformance-noise', '--air', p('.dsa2.r33'), p('.dsa2.raw')])
    return {'ambe2': p('.a2.raw'), 'dstar': p('.ds.raw'), 'dual': p('.dual.raw'), 'a2ds': p('.a2ds.raw'),
            'dsa2': p('.dsa2.raw')}, {'a2ds': p('.a2.raw'), 'dsa2': p('.ds.raw')}


def score_one(args):
    ref, deg = args
    import score as S
    r, d, lag = S.align(read_raw(ref), read_raw(deg))
    o = S.score(r, d)
    o['lag'] = lag
    return o


def speech(build, od, clips):
    os.makedirs(od, exist_ok=True)
    jobs, keys = [], []
    for raw in clips:
        dec, src = make_paths(build, raw, od)
        for k, f in dec.items():
            jobs.append((raw, f))
            keys.append((os.path.basename(raw), k))
        for k, f in src.items():
            jobs.append((f, dec[k]))
            keys.append((os.path.basename(raw), k + '/src'))
    with ProcessPoolExecutor(MAX_JOBS) as ex:
        res = list(ex.map(score_one, jobs))
    table = {}
    for (clip, k), o in zip(keys, res):
        table.setdefault(k, {})[clip] = o
    ok = True
    summary = {}
    print('%-10s %6s %6s %7s %7s %6s %6s   %s' % ('path', 'PESQ', 'min', 'ViSQOL', 'min', 'STOI', 'LSD', ''))
    for k in ['ambe2', 'dstar', 'dual', 'a2ds', 'dsa2', 'a2ds/src', 'dsa2/src']:
        rows = list(table[k].values())
        if any(r.get('visqol_speech') is None for r in rows):
            print('%-10s ViSQOL failed' % k)
            ok = False
            continue
        mp = sum(r['pesq_nb'] for r in rows) / len(rows)
        np_ = min(r['pesq_nb'] for r in rows)
        mv = sum(r['visqol_speech'] for r in rows) / len(rows)
        nv = min(r['visqol_speech'] for r in rows)
        ms = sum(r['stoi'] for r in rows) / len(rows)
        ml = sum(r['lsd_db'] for r in rows) / len(rows)
        lim = SPEECH_LIMITS.get(k)
        good = lim is None or (mp >= lim[0] and np_ >= lim[1] and mv >= lim[2] and nv >= lim[3] and ms >= lim[4]
                               and ml <= lim[5])
        ok &= good
        summary[k] = dict(pesq=mp, pesq_min=np_, visqol=mv, visqol_min=nv, stoi=ms, lsd=ml)
        print('%-10s %6.3f %6.3f %7.3f %7.3f %6.3f %6.3f   %s' % (
            k, mp, np_, mv, nv, ms, ml, 'informative' if lim is None else 'ok' if good else 'FAIL'))
    print('paths: ambe2 AMBE+2 encode/decode, dstar D-STAR only, dual the D-STAR output of the dual encoder,')
    print('       a2ds AMBE+2 -> xombe-xc -> D-STAR decode, dsa2 D-STAR -> xombe-xc --reverse -> AMBE+2 decode,')
    print('       */src the transcoded decode against the source codec\'s decode')
    print('%d clips, scored against the encoder input (*/src: the source decode); limits (mean PESQ, worst PESQ, '
          'mean ViSQOL, worst ViSQOL, mean STOI, mean LSD):' % len(clips))
    for k, lim in SPEECH_LIMITS.items():
        if lim:
            print('  %-8s >= %.2f, >= %.2f, >= %.2f, >= %.2f, >= %.2f, <= %.2f' % ((k,) + lim))
    return ok, dict(summary=summary, clips=table)


def tone_ids(build, dstar):
    out = run([tool(build, 'dstar_dump_fec'), '--ids', dstar]).stdout.split('\n')
    return [int(l.split()[2]) if l.split()[1] == 'tone' else None for l in out if l.strip()]


def ds163_case(od):
    """D-STAR detects tone 163 as 340 + 550 Hz and plays AMBE+2's 350 + 490 Hz"""
    import numpy as np
    n = np.arange(int(0.4 * FS))
    A = 32768 * 10 ** (-20 / 20)
    rng = np.random.default_rng(20261024)
    tone = sum(A * np.cos(2 * np.pi * f / FS * n + rng.uniform(0, 2 * np.pi)) for f in (550.0, 340.0))
    x = np.concatenate([np.zeros(FS // 2), tone, np.zeros(int(0.3 * FS))])
    np.clip(np.round(x), -32768, 32767).astype('<i2').tofile(os.path.join(od, 'ds_163.raw'))
    return dict(name='ds_163', file='ds_163.raw', segments=[dict(start=FS // 2, length=len(n), id=163,
                                                                 freqs=[550.0, 340.0], play=[490.0, 350.0])])


def tones(build, od, clips):
    import e6
    os.makedirs(od, exist_ok=True)
    e6.gen(od)
    cases = json.load(open(os.path.join(od, 'cases.json'))) + [ds163_case(od)]
    E, D, X = tool(build, 'xombe-enc'), tool(build, 'xombe-dec'), tool(build, 'xombe-xc')
    ok = True
    n_ok = 0
    rows = []
    fmt = lambda e: ' '.join('%+5.2f' % v if v is not None else '  n/a' for v in e)
    print('%-14s %4s %-6s %-22s %-12s  %-12s  %s' % ('case', 'ID', 'kind', 'D-STAR frames', 'level dB', 'a2ds dB',
                                                     ''))
    for c in cases:
        name = c['name']
        raw = os.path.join(od, c['file'])
        s = c['segments'][0]
        kind = 'ds163' if name == 'ds_163' else 'voice' if s['id'] in KNOX or s['id'] == 163 else 'tone'
        ds, r33, xds = (os.path.join(od, name + e) for e in ('.dstar', '.r33', '.xc.dstar'))
        run([E, '--dstar', raw, ds])
        run([D, '--conformance-noise', '--dstar', ds, ds + '.raw'])
        x, y = read_raw(raw), read_raw(ds + '.raw')
        ids = tone_ids(build, ds)
        on_f, end_f = s['start'] // 160, min(len(ids), (s['start'] + s['length']) // 160 + 10)
        seen = sorted({i for i in ids[on_f:end_f] if i is not None})
        first = next((i for i in range(on_f, end_f) if ids[i] == s['id']), None)
        if kind == 'voice':
            ok_id = not seen
            desc = 'voice frames' if ok_id else 'tone IDs %s' % seen
        else:
            ok_id = first is not None and first - on_f <= TONE_ONSET and seen == [s['id']]
            desc = 'none' if first is None else 'onset+%d' % (first - on_f)
            if seen != [s['id']]:
                desc += ' IDs %s' % seen

        def lev(sig):
            """level errors of the decoded components"""
            si = e6.steady(x, s['start'] - 160, s['start'] + s['length'] + 160)
            so = e6.steady(sig, s['start'] - 160, min(len(sig), s['start'] + s['length'] + 2400))
            if not si or not so:
                return [None]
            lin, _ = e6.fit_levels(x[si[0]:si[0] + si[1]], s['freqs'])
            play = s.get('play', s['freqs'])
            lout, _ = e6.fit_levels(sig[so[0]:so[0] + so[1]], play)
            return [o - i for o, i in zip(lout, lin)]
        lim = VOICE_LEVEL_DB if kind == 'voice' else TONE_LEVEL_DB
        err = lev(y)
        good = ok_id and None not in err and all(abs(e) <= lim for e in err)
        xdesc, xerr = '-', []
        if kind != 'ds163':
            # from AMBE+2: tone frames keep the ID, Knox tones become voice frames
            run([E, '--air', raw, r33])
            run([X, r33, xds])
            run([D, '--conformance-noise', '--dstar', xds, xds + '.raw'])
            xids = tone_ids(build, xds)
            xseen = sorted({i for i in xids[on_f:end_f] if i is not None})
            xknox = s['id'] in KNOX
            ok_xid = (not xseen) if xknox else xseen == [s['id']]
            xerr = lev(read_raw(xds + '.raw'))
            xlim = XC_KNOX_LEVEL_DB if xknox else TONE_LEVEL_DB
            good &= ok_xid and None not in xerr and all(abs(e) <= xlim for e in xerr)
            xdesc = fmt(xerr) + ('' if ok_xid else ' IDs %s' % xseen)
        n_ok += good
        ok &= good
        rows.append(dict(case=name, id=s['id'], kind=kind, first=None if first is None else first - on_f, ids=seen,
                         err=err, xc_err=xerr, ok=good))
        print('%-14s %4d %-6s %-22s %-12s  %-12s  %s' % (name, s['id'], kind, desc, fmt(err), xdesc,
                                                         'ok' if good else 'FAIL'))
    print('D-STAR tone cases: %d/%d pass' % (n_ok, len(cases)))
    print('  tone: D-STAR tone frame of the ID within %d frames of onset, no other ID, level +/-%.1f dB, '
          'and the same from AMBE+2 (a2ds)' % (TONE_ONSET, TONE_LEVEL_DB))
    print('  voice: Knox 144..159 and AMBE+2\'s 163 (350 + 490 Hz) have no D-STAR tone code: voice frames, levels '
          '+/-%.1f dB, Knox from AMBE+2 +/-%.1f dB' % (VOICE_LEVEL_DB, XC_KNOX_LEVEL_DB))
    print('  ds163: 340 + 550 Hz is D-STAR\'s tone 163, played as 350 + 490 Hz')
    bad = []
    for raw in clips:
        ds = os.path.join(od, 'speech_' + os.path.basename(raw) + '.dstar')
        run([E, '--dstar', raw, ds])
        n = sum(i is not None for i in tone_ids(build, ds))
        if n:
            bad.append((os.path.basename(raw), n))
    print('D-STAR speech: %d clips, %s' % (len(clips), 'no tone frames' if not bad else 'tone frames in %s' % bad))
    ok &= not bad
    return ok, dict(cases=rows, speech_tone_frames=bad)


def main():
    ap = argparse.ArgumentParser(description='D-STAR conformance')
    ap.add_argument('cmd', choices=['refs', 'speech', 'tones'])
    ap.add_argument('build')
    ap.add_argument('outdir', nargs='?')
    ap.add_argument('--json')
    ap.add_argument('--inputs')
    a = ap.parse_args()
    build = os.path.abspath(a.build)
    inputs = a.inputs or os.path.join(build, 'e3in')
    clips = [os.path.join(inputs, c + '.d0.raw') for c in CLIPS] + \
        [os.path.join(E5V, c + '.raw') for c in ('low5', 'hp5')]
    if a.cmd == 'refs':
        ok, data = refs(build)
    elif a.cmd == 'speech':
        ok, data = speech(build, a.outdir, clips)
    else:
        ok, data = tones(build, a.outdir, clips)
    if a.json:
        json.dump(dict(ok=ok, data=data), open(a.json, 'w'), indent=1)
    print('D-STAR %s: %s' % (a.cmd, 'PASS' if ok else 'FAIL'))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
