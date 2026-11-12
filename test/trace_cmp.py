#!/usr/bin/env python3
"""compares two stage dumps from test/trace_dump (or another implementation writing the same records) and
reports, per stage, how far the second is from the first

absolute differences: log-amplitude stages (DEQUANT, SMOOTH, SET, TARGET), VAD decision and metric (VAD),
and phases (PHASE, wrapped)
log2 ratios: linear amplitudes (ENHANCED, SHAPED), ignoring harmonics below --floor in the reference
SNR over the whole dump: SYNTH
records match by frame, half, and stage, a record whose length or f0 differs is a structural mismatch
exits 1 if any limit given on the command line is exceeded or any record is unmatched"""
import argparse
import math
import struct
import sys

STAGES = ['DEQUANT', 'SMOOTH', 'SET', 'ENHANCED', 'PHASE', 'SHAPED', 'SYNTH', 'TARGET', 'VAD']
LOG_STAGES = {'DEQUANT', 'SMOOTH', 'SET', 'TARGET', 'VAD'}
LIN_STAGES = {'ENHANCED', 'SHAPED'}


def read(path):
    recs = {}
    with open(path, 'rb') as f:
        data = f.read()
    pos = 0
    while pos < len(data):
        if pos + 12 > len(data):
            sys.exit('%s: truncated record header at byte %d' % (path, pos))
        frame, half, stage, n, f0 = struct.unpack_from('<IBBHf', data, pos)
        pos += 12
        if pos + 4 * n > len(data):
            sys.exit('%s: truncated record at byte %d' % (path, pos))
        v = struct.unpack_from('<%df' % n, data, pos)
        pos += 4 * n
        key = (frame, half, stage)
        # a stage may fire more than once per half, keep them in order
        recs.setdefault(key, []).append((f0, v))
    return recs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('ref')
    ap.add_argument('test')
    ap.add_argument('--floor', type=float, default=1e-3, help='linear amplitude below which ENHANCED/SHAPED '
                    'harmonics are not compared (default 1e-3, far below one output LSB)')
    ap.add_argument('--max-log', type=float, help='limit on the log2 stages and linear-stage log2 ratios')
    ap.add_argument('--max-phase', type=float, help='limit on PHASE differences, radians')
    ap.add_argument('--min-snr', type=float, help='limit on SYNTH SNR, dB')
    args = ap.parse_args()

    a, b = read(args.ref), read(args.test)
    stats = {s: {'n': 0, 'max': 0.0, 'where': None, 'sig': 0.0, 'err': 0.0, 'struct': 0} for s in STAGES}
    unmatched = len(set(a) ^ set(b))
    for key in sorted(set(a) & set(b)):
        name = STAGES[key[2]] if key[2] < len(STAGES) else 'stage%d' % key[2]
        st = stats.setdefault(name, {'n': 0, 'max': 0.0, 'where': None, 'sig': 0.0, 'err': 0.0, 'struct': 0})
        for (f0a, va), (f0b, vb) in zip(a[key], b[key]):
            st['n'] += 1
            if len(va) != len(vb) or abs(f0a - f0b) > 1e-6 * max(abs(f0a), 1e-3):
                st['struct'] += 1
                continue
            for i, (x, y) in enumerate(zip(va, vb)):
                if name in LOG_STAGES:
                    d = abs(x - y)
                elif name == 'PHASE':
                    d = abs(math.remainder(x - y, 2.0 * math.pi))
                elif name in LIN_STAGES:
                    if x < args.floor:
                        continue
                    d = abs(math.log2(max(y, 1e-30) / x))
                else:
                    st['sig'] += x * x
                    st['err'] += (x - y) * (x - y)
                    d = abs(x - y)
                if d > st['max']:
                    st['max'], st['where'] = d, (key[0], key[1], i)

    bad = unmatched > 0
    print('%-9s %8s %6s %12s  %s' % ('stage', 'records', 'struct', 'max diff', 'worst (frame, half, index)'))
    for name, st in stats.items():
        if not st['n']:
            continue
        line = '%-9s %8d %6d %12.3g  %s' % (name, st['n'], st['struct'], st['max'], st['where'])
        if name == 'SYNTH':
            snr = 10 * math.log10(st['sig'] / st['err']) if st['err'] > 0 else float('inf')
            line += '  SNR %.1f dB' % snr
            if args.min_snr is not None and snr < args.min_snr:
                bad = True
        elif name == 'PHASE':
            if args.max_phase is not None and st['max'] > args.max_phase:
                bad = True
        elif args.max_log is not None and st['max'] > args.max_log:
            bad = True
        bad |= st['struct'] > 0
        print(line)
    if unmatched:
        print('%d records present in only one dump' % unmatched)
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
