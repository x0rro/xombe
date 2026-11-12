#!/bin/sh
set -u
BUILD=${1:-build}
PY=${2:-python3}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
case $BUILD in /*) ;; *) BUILD=$ROOT/$BUILD ;; esac
VEC=$ROOT/vectors
OUT=$BUILD/conformance
IN=$BUILD/e3in  # dithered E3 inputs, from gen_inputs.py
JOBS=4
export XOMBE_BUILD=$BUILD
mkdir -p "$OUT/t3" "$OUT/e3" "$OUT/e5" "$OUT/e6in" "$OUT/e6" "$OUT/e7"
status=0
summary=""

record() {
    if [ "$2" -eq 0 ]; then summary="$summary
$1: PASS"; else summary="$summary
$1: FAIL"; status=1; fi
}

# runs a command, reports it on stderr and returns 1 if it fails
run() {
    "$@" || { echo "FAILED ($?): $*" >&2; return 1; }
}

echo "== inputs"
run "$PY" "$ROOT/scripts/gen_inputs.py" --out "$IN"
record inputs $?

echo "== T3"
t3=0
for s in es1_nsdtx_off es1_nsdtx_on es2 es3 es4; do
    run "$BUILD/xombe-dec" --conformance-noise "$VEC/perceptual/$s.amb" "$OUT/t3/$s.raw" || t3=1
done
# each scorer's stderr goes to its own log, printed if that scorer fails
for s in es1_nsdtx_off es1_nsdtx_on es2 es3 es4; do
    echo "$s"
done | xargs -P $JOBS -I{} sh -c "\"$PY\" \"$ROOT/tools/score.py\" \"$VEC/perceptual/{}_ref.raw\" \"$OUT/t3/{}.raw\" --json > \"$OUT/t3/{}.json\" 2> \"$OUT/t3/{}.log\" || { echo 'score.py failed for {}:' >&2; cat \"$OUT/t3/{}.log\" >&2; exit 1; }" || t3=1
"$PY" - "$OUT/t3" <<'EOF' || t3=1
import json, sys, os
d = sys.argv[1]
sets = ['es1_nsdtx_off', 'es1_nsdtx_on', 'es2', 'es3', 'es4']
keys = ['pesq_nb', 'visqol_speech', 'stoi', 'lsd_db']
ok = True
rows = []
print('%-14s %7s %7s %7s %7s' % ('set', 'PESQ', 'ViSQOL', 'STOI', 'LSD'))
for s in sets:
    try:
        r = json.load(open(os.path.join(d, s + '.json')))
        bad = [k for k in keys if not isinstance(r.get(k), (int, float))]
    except (OSError, ValueError) as e:
        r, bad = {}, ['unreadable scores (%s)' % e]
    if bad:
        print('%-14s missing or failed: %s  FAIL' % (s, ', '.join(bad)))
        ok = False
        continue
    good = r['pesq_nb'] >= 4.05 and r['visqol_speech'] >= 4.15 and r['stoi'] >= 0.98 and r['lsd_db'] <= 1.30
    ok &= good
    rows.append(r)
    print('%-14s %7.3f %7.3f %7.3f %7.3f %s' % (s, r['pesq_nb'], r['visqol_speech'], r['stoi'], r['lsd_db'],
                                               'ok' if good else 'FAIL'))
if len(rows) == len(sets):
    mp = sum(r['pesq_nb'] for r in rows) / len(rows)
    ml = sum(r['lsd_db'] for r in rows) / len(rows)
    mv = sum(r['visqol_speech'] for r in rows) / len(rows)
    ms = sum(r['stoi'] for r in rows) / len(rows)
    print('%-14s %7.3f %7.3f %7.3f %7.3f' % ('mean (5 sets)', mp, mv, ms, ml))
    ok &= mp >= 4.25 and ml <= 1.00
print('limits: each PESQ >= 4.05, ViSQOL >= 4.15, STOI >= 0.98, LSD <= 1.30; mean PESQ >= 4.25, mean LSD <= 1.00')
print('T3: ' + ('PASS' if ok else 'FAIL'))
sys.exit(0 if ok else 1)
EOF
record T3 $t3

echo "== T4"
run "$BUILD/test_stream" --metrics "$VEC"/stream/mixed_stream_0.jsonl.gz "$VEC"/stream/mixed_stream_1.jsonl.gz \
    "$VEC"/stream/mixed_stream_2.jsonl.gz "$VEC"/stream/mixed_stream_3.jsonl.gz
record T4 $?

echo "== T5"
run "$BUILD/test_tones" "$VEC/tones"
record T5 $?

echo "== encodes"
enc=0
for c in es1 es2 es3 sn10 ttsf ttsc ttsl; do
    for k in 0 1 2 3 4; do echo "$c.d$k"; done
done | xargs -P $JOBS -I{} sh -c "\"$BUILD/xombe-enc\" \"$IN/{}.raw\" \"$OUT/e3/{}.amb\" || { echo 'xombe-enc failed for {}' >&2; exit 1; }" || enc=1
for n in low5 hp5; do run "$BUILD/xombe-enc" "$VEC/e5/$n.raw" "$OUT/e5/$n.amb" || enc=1; done
for n in glide mix saw150 tones; do run "$BUILD/xombe-enc" "$VEC/e7/$n.raw" "$OUT/e7/$n.amb" || enc=1; done
record encodes $enc

echo "== E1"
e1=0
for f in "$OUT"/e3/*.amb; do
    b=$(basename "$f" .amb)
    "$PY" "$ROOT/tools/e1.py" "$f" --input "$IN/$b.raw" --quiet > /dev/null || { echo "E1 fail: $f"; e1=1; }
done
for d in e5 e7; do
    for f in "$OUT"/$d/*.amb; do
        b=$(basename "$f" .amb)
        "$PY" "$ROOT/tools/e1.py" "$f" --input "$VEC/$d/$b.raw" --quiet > /dev/null || { echo "E1 fail: $f"; e1=1; }
    done
done
echo "E1 (e1.py): $(ls "$OUT"/e3/*.amb "$OUT"/e5/*.amb "$OUT"/e7/*.amb | wc -l) streams checked"
record "E1 (streams)" $e1

echo "== E1/E2/E7"
run "$BUILD/test_mirror" "$VEC" --e7 --inputs "$IN" --stress 600
record "E1/E2/E7" $?

echo "== E3"
run "$PY" "$ROOT/tools/e3.py" "$OUT/e3" --jobs $JOBS --quiet --json "$OUT/e3.json"
record E3 $?

echo "== E4"
run "$PY" "$ROOT/tools/e4.py" "$OUT/e3" --json "$OUT/e4.json"
record E4 $?

echo "== E5"
run "$PY" "$ROOT/tools/e5.py" "$OUT/e3" --e5-cand "$OUT/e5" --inputs "$IN" --json "$OUT/e5.json"
record E5 $?

echo "== E6"
e6=0
run "$PY" "$ROOT/tools/e6.py" gen "$OUT/e6in" > /dev/null || e6=1
for f in "$OUT"/e6in/*.raw; do run "$BUILD/xombe-enc" "$f" "$OUT/e6/$(basename "$f" .raw).amb" || e6=1; done
run "$PY" "$ROOT/tools/e6.py" check "$OUT/e6in/cases.json" "$OUT/e6" > "$OUT/e6.txt" || e6=1
tail -2 "$OUT/e6.txt"
run "$PY" "$ROOT/tools/e6.py" check pkg "$OUT/e7" > "$OUT/e6pkg.txt" || e6=1
tail -2 "$OUT/e6pkg.txt"
run "$PY" "$ROOT/tools/e6.py" speech "$OUT/e3" "$OUT/e5" "$OUT/e7" || e6=1
record E6 $e6

echo "== summary$summary"
exit $status
