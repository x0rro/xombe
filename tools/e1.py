#!/usr/bin/env python3
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import ROOT, TONE_IDS, read_amb, val

FEC_REF = os.path.join(ROOT, 'vectors', 'fec', 'fec_encode_vectors.txt')
B0 = [0, 1, 2, 3, 37, 38, 39]
B1 = [4, 5, 6, 7, 35]


def check_frame(d):
    c = val(d, range(6))
    bad = []
    if c in (60, 61):
        return 'erasure', None, None, ['erasure class (d0..d5 = %s) must not be emitted' % format(c, '06b')]
    if c == 63:
        ids = [val(d, range(s, s + 8)) for s in (12, 20, 28, 36)]
        if len(set(ids)) != 1:
            bad.append('tone ID copies differ %s' % ids)
        if ids[0] not in TONE_IDS:
            bad.append('tone ID %d not in {5..122, 128..143, 160..162}' % ids[0])
        if any(d[45:49]):
            bad.append('tone d45..d48 != 0')
        return 'tone', None, ids[0], bad
    b0, b1 = val(d, B0), val(d, B1)
    if c == 62:
        if b0 != 124:
            bad.append('silence b0 = %d (expected 124)' % b0)
        if b1 != 16:
            bad.append('silence b1 = %d (expected 16)' % b1)
        return 'silence', b0, b1, bad
    if b0 > 119:
        bad.append('voice b0 = %d > 119' % b0)
    if 18 <= b1 <= 31:
        bad.append('voice b1 = %d (pulse code)' % b1)
    return 'voice', b0, b1, bad


def check_fec(path):
    def load(p):
        m = {}
        for line in open(p):
            line = line.split('#')[0].split()
            if len(line) >= 2:
                m[line[0].upper()] = line[1].upper()
        return m
    ref, cand = load(FEC_REF), load(path)
    miss = [k for k in ref if k not in cand]
    wrong = [k for k in ref if k in cand and cand[k] != ref[k]]
    return len(ref), miss, wrong


def main():
    ap = argparse.ArgumentParser(description='E1 bitstream validity')
    ap.add_argument('files', nargs='*')
    ap.add_argument('--input', help='raw input: also check one frame per 160 samples')
    ap.add_argument('--fec', help='candidate FEC pairs in the fec_encode_vectors.txt format')
    ap.add_argument('--frames', action='store_true')
    ap.add_argument('--quiet', action='store_true')
    a = ap.parse_args()
    ok = True
    for f in a.files:
        frames = read_amb(f)
        counts, nbad, tones = {}, 0, {}
        for n, d in enumerate(frames):
            cls, _, x, bad = check_frame(d)
            counts[cls] = counts.get(cls, 0) + 1
            if cls == 'tone':
                tones[x] = tones.get(x, 0) + 1
            if bad:
                nbad += 1
                if not a.quiet and (nbad <= 20 or a.frames):
                    print('  %s frame %d: %s' % (os.path.basename(f), n, '; '.join(bad)))
        line = '%s: %d frames (%s)' % (f, len(frames), ', '.join('%s %d' % kv for kv in sorted(counts.items())))
        if tones:
            line += '; tone IDs ' + ', '.join('%d×%d' % kv for kv in sorted(tones.items()))
        if a.input:
            want = os.path.getsize(a.input) // 320
            if len(frames) != want:
                line += '; FRAME COUNT %d != %d' % (len(frames), want)
                nbad += 1
        print(line + ('; E1 PASS' if not nbad else '; E1 FAIL (%d frames)' % nbad))
        ok &= not nbad
    if a.fec:
        n, miss, wrong = check_fec(a.fec)
        good = not miss and not wrong
        print('FEC: %d reference pairs, %d missing, %d wrong; %s' % (n, len(miss), len(wrong), 'PASS' if good else 'FAIL'))
        ok &= good
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
