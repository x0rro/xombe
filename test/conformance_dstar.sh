#!/bin/sh
# conformance_dstar.sh [BUILD [PYTHON]]: D-STAR gates from the repository's own vectors (tools/dstar_conf.py)
set -u
BUILD=${1:-build}
PY=${2:-python3}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
case $BUILD in /*) ;; *) BUILD=$ROOT/$BUILD ;; esac
OUT=$BUILD/conformance-dstar
export XOMBE_BUILD=$BUILD
mkdir -p "$OUT"
status=0
summary=""

record() {
    if [ "$2" -eq 0 ]; then summary="$summary
$1: PASS"; else summary="$summary
$1: FAIL"; status=1; fi
}

for t in xombe-enc xombe-dec xombe-xc dstar_dump_fec; do
    if [ ! -x "$BUILD/$t" ]; then
        echo "$BUILD/$t missing: conformance-dstar needs XOMBE_DSTAR, XOMBE_ENCODER and the tests" >&2
        exit 1
    fi
done
echo "== inputs"
"$PY" "$ROOT/scripts/gen_inputs.py" --out "$BUILD/e3in" > /dev/null
record inputs $?

echo "== refs"
"$PY" "$ROOT/tools/dstar_conf.py" refs "$BUILD" --json "$OUT/refs.json"
record refs $?

echo "== speech"
"$PY" "$ROOT/tools/dstar_conf.py" speech "$BUILD" "$OUT/speech" --inputs "$BUILD/e3in" --json "$OUT/speech.json"
record speech $?

echo "== tones"
"$PY" "$ROOT/tools/dstar_conf.py" tones "$BUILD" "$OUT/tones" --inputs "$BUILD/e3in" --json "$OUT/tones.json"
record tones $?

echo "== summary$summary"
exit $status
