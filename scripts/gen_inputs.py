#!/usr/bin/env python3
import argparse
import hashlib
import os
import struct
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'tools'))
from common import CLIPS, DITHER_DB as GAINS_DB, E3IN, E3V


def dither(samples, gain_db):
    g = 10.0 ** (gain_db / 20.0)
    out = []
    for x in samples:
        y = round(g * x)
        out.append(32767 if y > 32767 else -32768 if y < -32768 else y)
    return out


def main():
    ap = argparse.ArgumentParser(description='Regenerate the dithered E3 inputs and check them against SHA256SUMS')
    ap.add_argument('--dir', default=E3V)
    ap.add_argument('--out', default=E3IN, help='output directory (default $XOMBE_BUILD/e3in)')
    ap.add_argument('--write-manifest', action='store_true')
    a = ap.parse_args()
    out_dir = a.out
    os.makedirs(out_dir, exist_ok=True)
    manifest_path = os.path.join(a.dir, 'SHA256SUMS')
    expected = {}
    if not a.write_manifest and not os.path.exists(manifest_path):
        print('%s is missing; regenerate it with --write-manifest' % manifest_path, file=sys.stderr)
        return 1
    if not a.write_manifest:
        for line in open(manifest_path):
            digest, name = line.split()
            expected[name] = digest
    lines, bad = [], 0
    for clip in CLIPS:
        raw = open(os.path.join(a.dir, clip + '_orig.raw'), 'rb').read()
        x = struct.unpack('<%dh' % (len(raw) // 2), raw)
        for k, g in enumerate(GAINS_DB):
            data = struct.pack('<%dh' % len(x), *dither(x, g))
            name = '%s.d%d.raw' % (clip, k)
            open(os.path.join(out_dir, name), 'wb').write(data)
            digest = hashlib.sha256(data).hexdigest()
            lines.append('%s  %s' % (digest, name))
            if not a.write_manifest and expected.get(name) != digest:
                print('MISMATCH %s' % name)
                bad += 1
    if a.write_manifest:
        open(manifest_path, 'w').write('\n'.join(lines) + '\n')
        print('wrote %s' % manifest_path)
    else:
        print('%d inputs regenerated, %d mismatches against SHA256SUMS' % (len(lines), bad))
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
