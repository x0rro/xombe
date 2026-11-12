#!/usr/bin/env python3
import os, sys, json, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import FS, read_raw as load
import numpy as np
from scipy.signal import resample_poly, stft

def align(r, d, maxlag=800):
    n = min(len(r), len(d))
    c = np.correlate(d[:n], r[:n], 'full')[n - 1 - maxlag:n + maxlag]
    lag = int(np.argmax(np.abs(c))) - maxlag
    if lag > 0: d = d[lag:]
    elif lag < 0: d = np.concatenate([np.zeros(-lag), d])
    n = min(len(r), len(d))
    return r[:n], d[:n], lag

def lsd(r, d):
    _, _, R = stft(r, FS, nperseg=256, noverlap=192)
    _, _, D = stft(d, FS, nperseg=256, noverlap=192)
    pr, pd = np.abs(R) ** 2, np.abs(D) ** 2
    act = 10 * np.log10(pr.sum(0) + 1e-9) > 10 * np.log10(pr.sum(0).max() + 1e-9) - 40
    band = (np.arange(pr.shape[0]) * FS / 256 >= 100) & (np.arange(pr.shape[0]) * FS / 256 <= 3800)
    floor = pr.max() * 1e-6
    x = 10 * np.log10(pr[band][:, act] + floor) - 10 * np.log10(pd[band][:, act] + floor)
    return float(np.mean(np.sqrt(np.mean(x ** 2, 0))))

def segsnr(r, d, n=160):
    k = len(r) // n
    s = []
    for i in range(k):
        a, b = r[i * n:(i + 1) * n], d[i * n:(i + 1) * n]
        e = np.sum(a ** 2)
        if e < 1e-3 * n * 32768 ** 2 * 1e-6: continue
        s.append(np.clip(10 * np.log10(e / (np.sum((a - b) ** 2) + 1e-9)), -10, 35))
    return float(np.mean(s)) if s else float('nan')

def score(r, d):
    out = {}
    from pesq import pesq
    out['pesq_nb'] = float(pesq(FS, r / 32768, d / 32768, 'nb'))
    from pystoi import stoi
    out['stoi'] = float(stoi(r, d, FS, extended=False))
    try:
        from visqol import VisqolApi
        api = VisqolApi(); api.create(mode='speech')
        r16, d16 = resample_poly(r, 2, 1), resample_poly(d, 2, 1)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            out['visqol_speech'] = float(api.measure_from_arrays(r16, d16, 16000).moslqo)
    except Exception as e:
        out['visqol_speech'] = None
        out['visqol_error'] = str(e)
    out['lsd_db'] = lsd(r, d)
    out['segsnr_db'] = segsnr(r, d)
    return out

if __name__ == '__main__':
    a = [x for x in sys.argv[1:] if not x.startswith('--')]
    r, d = load(a[0]), load(a[1])
    lag = 0
    if '--no-align' not in sys.argv: r, d, lag = align(r, d)
    o = score(r, d); o['lag'] = lag; o['seconds'] = len(r) / FS
    print(json.dumps(o) if '--json' in sys.argv else '\n'.join(f'{k:14s} {v}' for k, v in o.items()))
    if 'visqol_error' in o:
        print('ViSQOL failed: ' + o['visqol_error'], file=sys.stderr)
        sys.exit(1)
