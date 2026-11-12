#!/bin/sh
# fetch_librispeech.sh [OUT]: four LibriSpeech dev-clean-2 utterances as 8 kHz raw at -26 dBov
# (OUT defaults to build/librispeech in the source tree), needs curl, sox, and a Python with numpy and scipy
set -eu
ROOT=$(cd "$(dirname "$0")/.." && pwd)
OUT=${1:-$ROOT/build/librispeech}
PY=${PYTHON:-python3}
URL=https://www.openslr.org/resources/31/dev-clean-2.tar.gz
SHA256=176ec501490eced2d6c1f89f4f0ddc7dfe799e649e5322f8ba49fe3ff50c8012
for t in curl sox; do
    command -v $t > /dev/null || { echo "$t not found" >&2; exit 1; }
done
if command -v sha256sum > /dev/null; then
    sha() { sha256sum "$1"; }
elif command -v shasum > /dev/null; then
    sha() { shasum -a 256 "$1"; }
else
    echo "sha256sum or shasum not found" >&2; exit 1
fi
"$PY" -c 'import numpy, scipy' 2> /dev/null || { echo "$PY has no numpy or scipy (set PYTHON)" >&2; exit 1; }
mkdir -p "$OUT"
cd "$OUT"
if [ ! -f dev-clean-2.tar.gz ]; then
    curl -L -o dev-clean-2.tar.gz.part "$URL"
    mv dev-clean-2.tar.gz.part dev-clean-2.tar.gz
fi
digest=$(sha dev-clean-2.tar.gz | cut -d' ' -f1)
[ "$digest" = "$SHA256" ] || { echo "unexpected SHA-256 $digest for dev-clean-2.tar.gz" >&2; exit 1; }
for spec in lsf1:2035-147960-0000 lsm1:6295-244435-0010 lsm2:1272-135031-0010 lsf2:84-121550-0001; do
    name=${spec%%:*}
    utt=${spec#*:}
    spk=$(echo "$utt" | cut -d- -f1)
    chap=$(echo "$utt" | cut -d- -f2)
    member=LibriSpeech/dev-clean-2/$spk/$chap/$utt.flac
    tar -xzf dev-clean-2.tar.gz "$member"
    sox "$member" -r 8000 -b 16 -e signed-integer -c 1 -t raw "$name.raw" rate -v
    "$PY" - "$ROOT/tools" "$name.raw" <<'PYEOF'
import sys
sys.path.insert(0, sys.argv[1])
import numpy as np
from common import p56
x = np.fromfile(sys.argv[2], '<i2').astype(float)
y = np.clip(np.round(x * 10 ** ((-26.0 - p56(x)) / 20)), -32768, 32767).astype('<i2')
y.tofile(sys.argv[2])
print('%s: %.1f dBov' % (sys.argv[2], p56(y)))
PYEOF
done
echo "wrote lsf1.raw lsm1.raw lsm2.raw lsf2.raw to $OUT"
