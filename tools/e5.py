#!/usr/bin/env python3
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (DEFAULT_DEC, E3IN, E3V, E5V, FS, MAX_JOBS, decode_many, frame_class, p56, read_amb,
                    read_raw)
import argparse
import json
import math
import numpy as np

CLEAN_CLIPS = ['es1', 'es2', 'es3', 'ttsf', 'ttsc', 'ttsl']
LIM_CLIP = -0.30
LIM_CN_LEVEL, LIM_CN_SHAPE, CN_RANGE = 3.0, 4.0, 30.0
CLIP_MARGIN, CLIP_SKIP = 15.0, 3


def classes(amb):
    return [frame_class(d) for d in read_amb(amb)]


def cn_stats(pcm, cls, min_run=5, trim=2):
    from scipy.signal import welch
    first_v = cls.index('voice') if 'voice' in cls else len(cls)
    segs, i = [], first_v
    while i < len(cls):
        if cls[i] != 'silence':
            i += 1
            continue
        j = i
        while j < len(cls) and cls[j] == 'silence':
            j += 1
        if j - i >= min_run:
            segs.append(pcm[160 * (i + trim):160 * (j - trim)])
        i = j
    segs = [s for s in segs if len(s)]
    if not segs:
        return None
    x = np.concatenate(segs)
    f, P = welch(x, FS, nperseg=256)
    bands = np.array([10 * math.log10(np.sum(P[(f >= 500 * b) & (f < 500 * (b + 1))]) + 1e-12) for b in range(8)])
    return dict(frames=len(x) // 160, level=10 * math.log10(np.mean(x * x) + 1e-12), bands=bands)


def pesq_vs(orig, deg):
    import score as S
    from pesq import pesq
    r, d, _ = S.align(orig, deg)
    return float(pesq(FS, r / 32768, d / 32768, 'nb'))


def clipped_frames(x, cls, tol=1):
    al = p56(x)
    n = min(len(cls), len(x) // 160)
    fr = x[:n * 160].reshape(-1, 160)
    loud = 10 * np.log10(np.mean(fr ** 2, 1) + 1e-9) - 20 * math.log10(32768) > al - CLIP_MARGIN
    return al, [i for i in range(CLIP_SKIP, n) if cls[i] == 'silence' and loud[max(0, i - tol):i + tol + 1].all()]


def main():
    ap = argparse.ArgumentParser(description='E5 noise clips, comfort noise and clipped speech')
    ap.add_argument('e3_cand')
    ap.add_argument('--e5-cand', help='encodes of vectors/e5/low5.raw and hp5.raw')
    ap.add_argument('--cand-pattern', default='{clip}.d0.amb')
    ap.add_argument('--e5-pattern', default='{name}.amb')
    ap.add_argument('--inputs', default=E3IN, help='dithered E3 inputs from scripts/gen_inputs.py')
    ap.add_argument('--ref-mode', action='store_true')
    ap.add_argument('--decoder', default=DEFAULT_DEC)
    ap.add_argument('--cache')
    ap.add_argument('--jobs', type=int, default=MAX_JOBS)
    ap.add_argument('--json')
    a = ap.parse_args()
    if a.ref_mode:
        a.e3_cand, a.cand_pattern, a.e5_cand, a.e5_pattern = E3V, '{clip}.ref.d0.amb', E5V, '{name}.ref.amb'
    noise = [('sn10', os.path.join(a.e3_cand, a.cand_pattern.format(clip='sn10', k=0)),
              os.path.join(E3V, 'sn10.ref.d0.amb'), os.path.join(E3V, 'sn10_orig.raw'))]
    if a.e5_cand:
        for nm in ('low5', 'hp5'):
            noise.append((nm, os.path.join(a.e5_cand, a.e5_pattern.format(name=nm)),
                          os.path.join(E5V, nm + '.ref.amb'), os.path.join(E5V, nm + '.raw')))
    for _, c, r, o in noise:
        for p in (c, r, o):
            if not os.path.exists(p):
                sys.exit('missing: ' + p)
    pcm = decode_many([p for _, c, r, _ in noise for p in (c, r)], a.decoder, a.cache, a.jobs)
    clean = read_raw(os.path.join(E5V, 'es1_lead0.5s.raw'))
    out, fails = {'noise': {}, 'clipped': {}}, []
    print('E5: decoder = %s, g = 0 input' % a.decoder)
    print('%-5s %22s %22s %9s %9s %9s  %s' % ('clip', 'PESQ noisy c/r/Δ', 'PESQ clean c/r/Δ', 'CN frames',
                                             'CN level', 'shape', 'band diffs counted (dB)'))
    for nm, c, r, o in noise:
        noisy = read_raw(o)
        pn = (pesq_vs(noisy, pcm[c]), pesq_vs(noisy, pcm[r]))
        pc = (pesq_vs(clean, pcm[c]), pesq_vs(clean, pcm[r]))
        res = dict(pesq_noisy=[pn[0], pn[1], pn[0] - pn[1]], pesq_clean=[pc[0], pc[1], pc[0] - pc[1]])
        if pn[0] - pn[1] < LIM_CLIP:
            fails.append('%s ΔPESQ vs noisy %.3f' % (nm, pn[0] - pn[1]))
        cm, cr = cn_stats(pcm[c], classes(c)), cn_stats(pcm[r], classes(r))
        cn_txt = 'no comfort noise in %s' % ('candidate' if cr else 'reference') if not (cm and cr) else ''
        if cm and cr:
            keep = cr['bands'] >= cr['bands'].max() - CN_RANGE
            dl = cm['level'] - cr['level']
            shape = (cm['bands'] - cr['bands']) - dl
            sm = float(np.max(np.abs(shape[keep])))
            res.update(cn_frames=[cm['frames'], cr['frames']], cn_level_diff=dl, band_shape=shape.tolist(),
                       bands_counted=keep.tolist(), shape_dev_max=sm)
            if abs(dl) > LIM_CN_LEVEL:
                fails.append('%s comfort-noise level %+.2f dB' % (nm, dl))
            if sm > LIM_CN_SHAPE:
                fails.append('%s comfort-noise shape %.2f dB' % (nm, sm))
            cn_txt = '%4d/%-4d %+9.2f %9.2f  %s' % (cm['frames'], cr['frames'], dl, sm,
                                                   ' '.join(('%+.1f' % v) if k else '  .' for v, k in zip(shape, keep)))
        elif cr and not cm:
            fails.append('%s: candidate sends no comfort noise' % nm)
        out['noise'][nm] = res
        print('%-5s %6.3f/%5.3f/%+6.3f %6.3f/%5.3f/%+6.3f  %s' % (nm, pn[0], pn[1], pn[0] - pn[1], pc[0], pc[1],
                                                               pc[0] - pc[1], cn_txt))
    print('  limits: ΔPESQ vs noisy >= %.2f; CN level within ±%.0f dB; shape within ±%.0f dB over bands within'
          ' %.0f dB of the reference\'s strongest' % (LIM_CLIP, LIM_CN_LEVEL, LIM_CN_SHAPE, CN_RANGE))
    print('Clipped speech: silence frames n >= %d with input frames n-1..n+1 within %.0f dB of the active level'
          % (CLIP_SKIP, CLIP_MARGIN))
    tot_c = tot_r = 0
    for c in CLEAN_CLIPS:
        amb = os.path.join(a.e3_cand, a.cand_pattern.format(clip=c, k=0))
        if not os.path.exists(amb):
            sys.exit('missing: ' + amb)
        cls, cls_r = classes(amb), classes(os.path.join(E3V, c + '.ref.d0.amb'))
        x = read_raw(os.path.join(a.inputs, c + '.d0.raw'))
        al, bad = clipped_frames(x, cls)
        _, bad_r = clipped_frames(x, cls_r)
        tot_c += len(bad)
        tot_r += len(bad_r)
        out['clipped'][c] = dict(cand=bad, ref=bad_r)
        print('  %-5s active %6.1f dBov  silence frames %3d  clipped: cand %d %s  ref %d' %
              (c, al, cls.count('silence'), len(bad), bad[:10] if bad else '', len(bad_r)))
    print('  total clipped: cand %d, ref %d (limit: cand <= ref)' % (tot_c, tot_r))
    if tot_c > tot_r:
        fails.append('clipped speech frames %d > reference %d' % (tot_c, tot_r))
    print('E5: ' + ('PASS' if not fails else 'FAIL (' + '; '.join(fails) + ')') +
          ('' if a.e5_cand else '  [low5/hp5 not checked]'))
    if a.json:
        json.dump(dict(out, pass_=not fails, fails=fails), open(a.json, 'w'), indent=1)
    return 0 if not fails else 1


if __name__ == '__main__':
    sys.exit(main())
