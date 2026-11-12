#!/usr/bin/env python3
import os
import sys
import argparse
import hashlib
import json
import shlex
import subprocess
import warnings
from concurrent.futures import ProcessPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from common import BUILD, CLIPS, DEFAULT_DEC, DITHER_DB as DB, E3V, MAX_JOBS, decoder_id

LIM_CLIP, LIM_MEAN, LIM_VQ = -0.30, -0.08, -0.20


def sha(*parts):
    h = hashlib.sha256()
    for p in parts:
        h.update(p if isinstance(p, bytes) else str(p).encode())
        h.update(b'\0')
    return h.hexdigest()[:24]


def score_pair(orig, deg, visqol):
    import score as S
    from pesq import pesq
    r, d, lag = S.align(S.load(orig), S.load(deg))
    o = {'pesq_nb': float(pesq(S.FS, r / 32768, d / 32768, 'nb')), 'lag': lag}
    if visqol:
        from scipy.signal import resample_poly
        from visqol import VisqolApi
        api = VisqolApi()
        api.create(mode='speech')
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            o['visqol_speech'] = float(api.measure_from_arrays(resample_poly(r, 2, 1), resample_poly(d, 2, 1),
                                                               16000).moslqo)
    return o


def job(args):
    amb, orig, visqol, dec_cmd, dec_id, cache = args
    key = sha(dec_id, open(amb, 'rb').read(), open(orig, 'rb').read(), visqol)
    cf = os.path.join(cache, key + '.json')
    if os.path.exists(cf):
        return json.load(open(cf))
    raw = os.path.join(cache, key + '.raw')
    if not os.path.exists(raw):
        tmp = raw + '.%d.tmp' % os.getpid()
        subprocess.run(dec_cmd.format(inp=shlex.quote(amb), out=shlex.quote(tmp)), shell=True, check=True,
                       capture_output=True)
        os.replace(tmp, raw)
    o = score_pair(orig, raw, visqol)
    tmp = cf + '.%d.tmp' % os.getpid()
    json.dump(o, open(tmp, 'w'))
    os.replace(tmp, cf)
    return o


def main():
    ap = argparse.ArgumentParser(description='E3 five-dither perceptual comparison against the reference encoder')
    ap.add_argument('cand_dir')
    ap.add_argument('--clips', default=','.join(CLIPS))
    ap.add_argument('--jobs', type=int, default=MAX_JOBS)
    ap.add_argument('--json')
    ap.add_argument('--cand-pattern', default='{clip}.d{k}.amb')
    ap.add_argument('--decoder', default=DEFAULT_DEC)
    ap.add_argument('--vectors', default=E3V)
    ap.add_argument('--cache', default=os.path.join(BUILD, 'e3cache'))
    ap.add_argument('--quiet', action='store_true')
    a = ap.parse_args()
    a.jobs = max(1, min(a.jobs, MAX_JOBS))
    clips = [c for c in a.clips.split(',') if c]
    os.makedirs(a.cache, exist_ok=True)
    dec_id = decoder_id(a.decoder)
    tasks = []
    for c in clips:
        orig = os.path.join(a.vectors, c + '_orig.raw')
        for k in range(len(DB)):
            ref = os.path.join(a.vectors, '%s.ref.d%d.amb' % (c, k))
            cand = os.path.join(a.cand_dir, a.cand_pattern.format(clip=c, k=k))
            for p in (orig, ref, cand):
                if not os.path.exists(p):
                    sys.exit('missing: ' + p)
            for who, amb in (('ref', ref), ('cand', cand)):
                tasks.append(((c, k, who), (os.path.abspath(amb), orig, k == 0, a.decoder, dec_id, a.cache)))
    res = {}
    with ProcessPoolExecutor(a.jobs) as ex:
        for (tag, _), o in zip(tasks, ex.map(job, [t[1] for t in tasks])):
            res[tag] = o
            if not a.quiet:
                print('  %-5s d%d %-4s PESQ %.3f%s' % (tag[0], tag[1], tag[2], o['pesq_nb'],
                      ('  ViSQOL %.3f' % o['visqol_speech']) if 'visqol_speech' in o else ''),
                      file=sys.stderr, flush=True)
    out = {'decoder': a.decoder, 'dithers_db': DB, 'clips': {}}
    print('E3: PESQ-NB vs undithered original, decoder = %s' % a.decoder)
    print('dither k:           ' + ' '.join('%+7.3f' % g for g in DB) + '  (dB)')
    print('%-5s %-4s ' % ('clip', 'enc') + ' '.join('    d%d ' % k for k in range(len(DB))) +
          '   mean   ΔPESQ  ViSQOL  ΔViSQOL')
    dps, dvs, fails = [], [], []
    for c in clips:
        row = {}
        for who in ('ref', 'cand'):
            p = [res[(c, k, who)]['pesq_nb'] for k in range(len(DB))]
            row[who] = {'pesq': p, 'pesq_mean': sum(p) / len(p), 'visqol': res[(c, 0, who)]['visqol_speech']}
        dp = row['cand']['pesq_mean'] - row['ref']['pesq_mean']
        dv = row['cand']['visqol'] - row['ref']['visqol']
        row['d_pesq'], row['d_visqol'] = dp, dv
        row['d_pesq_per_dither'] = [x - y for x, y in zip(row['cand']['pesq'], row['ref']['pesq'])]
        dps.append(dp)
        dvs.append(dv)
        ok = dp >= LIM_CLIP
        if not ok:
            fails.append('%s ΔPESQ %.3f < %.2f' % (c, dp, LIM_CLIP))
        for who in ('ref', 'cand'):
            r = row[who]
            print('%-5s %-4s ' % (c if who == 'ref' else '', who) + ' '.join('%7.3f' % x for x in r['pesq']) +
                  ' %7.3f' % r['pesq_mean'] + ('  %+.3f%s' % (dp, '' if ok else '!') if who == 'cand' else '        ') +
                  ' %7.3f' % r['visqol'] + ('  %+.3f' % dv if who == 'cand' else ''))
        print('%-5s %-4s ' % ('', 'Δ') + ' '.join('%+7.3f' % x for x in row['d_pesq_per_dither']))
        out['clips'][c] = row
    mp, mv = sum(dps) / len(dps), sum(dvs) / len(dvs)
    if mp < LIM_MEAN:
        fails.append('mean ΔPESQ %.3f < %.2f' % (mp, LIM_MEAN))
    if mv < LIM_VQ:
        fails.append('mean ΔViSQOL %.3f < %.2f' % (mv, LIM_VQ))
    rm = sum(out['clips'][c]['ref']['pesq_mean'] for c in clips) / len(clips)
    cm = sum(out['clips'][c]['cand']['pesq_mean'] for c in clips) / len(clips)
    print('mean over %d clips: ref %.3f  cand %.3f  ΔPESQ %+.3f  ΔViSQOL %+.3f' % (len(clips), rm, cm, mp, mv))
    print('thresholds: per clip ΔPESQ >= %.2f, mean ΔPESQ >= %.2f, mean ΔViSQOL >= %.2f' % (LIM_CLIP, LIM_MEAN, LIM_VQ))
    full = set(clips) == set(CLIPS)
    print('E3%s: ' % ('' if full else ' (informative, partial clip set)') +
          ('PASS' if not fails else 'FAIL (' + '; '.join(fails) + ')'))
    out.update(mean_d_pesq=mp, mean_d_visqol=mv, pass_=not fails and full, fails=fails)
    if a.json:
        json.dump(out, open(a.json, 'w'), indent=1)
    return 0 if not fails else 1


if __name__ == '__main__':
    sys.exit(main())
