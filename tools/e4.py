#!/usr/bin/env python3
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import CLIPS, DEFAULT_DEC, E3V, MAX_JOBS, decode_many, p56, read_raw
import argparse
import json

LIM_E4 = 1.5


def main():
    ap = argparse.ArgumentParser(description='E4 active speech level (ITU-T P.56 method B), g = 0 inputs')
    ap.add_argument('cand_dir')
    ap.add_argument('--clips', default=','.join(CLIPS))
    ap.add_argument('--cand-pattern', default='{clip}.d0.amb')
    ap.add_argument('--decoder', default=DEFAULT_DEC)
    ap.add_argument('--vectors', default=E3V)
    ap.add_argument('--cache')
    ap.add_argument('--jobs', type=int, default=MAX_JOBS)
    ap.add_argument('--json')
    a = ap.parse_args()
    clips = [c for c in a.clips.split(',') if c]
    pairs = []
    for c in clips:
        cand = os.path.join(a.cand_dir, a.cand_pattern.format(clip=c, k=0))
        ref = os.path.join(a.vectors, c + '.ref.d0.amb')
        for p in (cand, ref, os.path.join(a.vectors, c + '_orig.raw')):
            if not os.path.exists(p):
                sys.exit('missing: ' + p)
        pairs.append((c, cand, ref))
    pcm = decode_many([p for _, c_, r_ in pairs for p in (c_, r_)], a.decoder, a.cache, a.jobs)
    print('E4: active speech level, decoder = %s' % a.decoder)
    print('%-6s %9s %9s %9s %7s' % ('clip', 'cand', 'ref', 'orig', 'diff'))
    out, fails = {}, []
    for c, cand, ref in pairs:
        lc, lr = p56(pcm[cand]), p56(pcm[ref])
        lo = p56(read_raw(os.path.join(a.vectors, c + '_orig.raw')))
        ok = abs(lc - lr) <= LIM_E4
        if not ok:
            fails.append('%s %+.2f dB' % (c, lc - lr))
        out[c] = dict(cand=lc, ref=lr, orig=lo, diff=lc - lr, ok=ok)
        print('%-6s %9.2f %9.2f %9.2f %+7.2f %s' % (c, lc, lr, lo, lc - lr, '' if ok else 'FAIL'))
    print('levels in dBov; limit |cand - ref| <= %.1f dB' % LIM_E4)
    print('E4: ' + ('PASS' if not fails else 'FAIL (' + '; '.join(fails) + ')'))
    if a.json:
        json.dump(dict(clips=out, pass_=not fails, fails=fails), open(a.json, 'w'), indent=1)
    return 0 if not fails else 1


if __name__ == '__main__':
    sys.exit(main())
