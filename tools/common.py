"""paths, clip lists, file formats, and helpers shared by the scoring and checking scripts

import before numpy, it pins BLAS/OpenMP threads (unless the environment sets them) for the scripts' own pools
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import hashlib
import math
import shlex
import shutil
import subprocess
import tempfile
import uuid
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
VEC = os.path.join(ROOT, 'vectors')
E3V = os.path.join(VEC, 'e3')
E5V = os.path.join(VEC, 'e5')
E7V = os.path.join(VEC, 'e7')
BUILD = os.environ.get('XOMBE_BUILD', os.path.join(ROOT, 'build'))
E3IN = os.path.join(BUILD, 'e3in')  # dithered E3 inputs, written by scripts/gen_inputs.py
FS = 8000
MAX_JOBS = 4
DEFAULT_DEC = shlex.quote(os.path.join(BUILD, 'xombe-dec')) + ' --conformance-noise {inp} {out}'

# E3 clips and the level dithers of their five inputs in dB, input k is round(10^(g/20)*x)
CLIPS = ['es1', 'es2', 'es3', 'sn10', 'ttsf', 'ttsc', 'ttsl']
DITHER_DB = [0.0, +0.025, -0.025, +0.05, -0.05]

# dual-tone IDs -> (high, low) Hz (DTMF, Knox 144..159, call progress 160..163), single tones are IDs 5..122 at
# 31.25 Hz x ID, must match src/common/tone_table.c
DUAL = {128: (1336, 941), 129: (1209, 697), 130: (1336, 697), 131: (1477, 697), 132: (1209, 770),
        133: (1336, 770), 134: (1477, 770), 135: (1209, 852), 136: (1336, 852), 137: (1477, 852),
        138: (1633, 697), 139: (1633, 770), 140: (1633, 852), 141: (1633, 941), 142: (1209, 941),
        143: (1477, 941), 144: (1162, 820), 145: (1052, 606), 146: (1162, 606), 147: (1297, 606),
        148: (1052, 672), 149: (1162, 672), 150: (1297, 672), 151: (1052, 743), 152: (1162, 743),
        153: (1297, 743), 154: (1430, 606), 155: (1430, 672), 156: (1430, 743), 157: (1430, 820),
        158: (1052, 820), 159: (1297, 820), 160: (440, 350), 161: (480, 440), 162: (620, 480), 163: (490, 350)}
TONE_IDS = set(range(5, 123)) | set(DUAL)


def read_amb(path):
    """frames of an .amb file as lists of 49 bits d0..d48"""
    data = open(path, 'rb').read()
    if data[:4] != b'.amb':
        raise ValueError('%s: missing .amb header' % path)
    if (len(data) - 4) % 8:
        raise ValueError('%s: length is not 4 + 8*n bytes' % path)
    frames = []
    for off in range(4, len(data), 8):
        fr = data[off:off + 8]
        frames.append([(fr[1 + i // 8] >> (7 - i % 8)) & 1 for i in range(48)] + [fr[7] & 1])
    return frames


def val(d, idx):
    """bits of d at positions idx, MSB first, as an integer"""
    v = 0
    for i in idx:
        v = (v << 1) | d[i]
    return v


def frame_class(d):
    c = val(d, range(6))
    return 'voice' if c < 60 else 'erasure' if c < 62 else 'silence' if c == 62 else 'tone'


def read_raw(path):
    import numpy as np
    return np.fromfile(path, '<i2').astype(float)


def decoder_id(cmd):
    """cache key for a decoder command: the command line and the executable's contents"""
    h = hashlib.sha256(cmd.encode())
    exe = shlex.split(cmd)[0]
    if os.path.exists(exe):
        h.update(open(exe, 'rb').read())
    return h.hexdigest()[:24]


def decode_many(ambs, dec_cmd=DEFAULT_DEC, cache=None, jobs=MAX_JOBS):
    """decodes each .amb with dec_cmd (a shell command with {inp} and {out}), returns {amb: pcm}"""
    tmpdir = None
    if cache is None:
        tmpdir = cache = tempfile.mkdtemp(prefix='xdec_')
    os.makedirs(cache, exist_ok=True)
    did = decoder_id(dec_cmd)

    def one(amb):
        key = hashlib.sha256((did + hashlib.sha256(open(amb, 'rb').read()).hexdigest()).encode()).hexdigest()[:24]
        out = os.path.join(cache, key + '.raw')
        if not os.path.exists(out):
            tmp = out + '.%s.tmp' % uuid.uuid4().hex
            subprocess.run(dec_cmd.format(inp=shlex.quote(os.path.abspath(amb)), out=shlex.quote(tmp)),
                           shell=True, check=True, capture_output=True)
            os.replace(tmp, out)
        return read_raw(out)
    try:
        with ThreadPoolExecutor(max(1, min(jobs, MAX_JOBS))) as ex:
            pcm = list(ex.map(one, ambs))
    finally:
        if tmpdir:
            shutil.rmtree(tmpdir, ignore_errors=True)
    return dict(zip(ambs, pcm))


def p56(x, fs=FS):
    """active speech level in dBov (ITU-T P.56 method B)"""
    import numpy as np
    from scipy.signal import lfilter
    x = np.asarray(x, float)
    sq = float(np.sum(x * x))
    if sq <= 0:
        return -200.0
    g = math.exp(-1 / (fs * 0.03))
    q = lfilter([1 - g], [1, -g], lfilter([1 - g], [1, -g], np.abs(x)))
    hang = int(0.2 * fs)
    A, C = [], []
    for j in range(16):
        c = 2.0 ** j
        cs = np.cumsum((q >= c).astype(float))
        win = cs - np.concatenate([np.zeros(hang + 1), cs[:-(hang + 1)]])
        n_act = int(np.sum(win > 0))
        A.append(10 * math.log10(sq / n_act) if n_act else -200.0)
        C.append(20 * math.log10(c))
    d = np.array(A) - np.array(C) - 15.9
    for j in range(1, 16):
        if d[j - 1] >= 0 > d[j]:
            t = d[j - 1] / (d[j - 1] - d[j])
            return float(A[j - 1] + t * (A[j] - A[j - 1]) - 20 * math.log10(32768))
    return float(A[0] - 20 * math.log10(32768))
