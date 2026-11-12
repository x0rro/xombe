#!/bin/sh
# bit-exactness snapshot of a xombe source tree, compare two with test/snapcmp.py
#
# usage: snapshot.sh [-q] [-j N] [-D CMAKE_DEF]... SRC OUT
#   SRC   xombe source tree (a worktree or an exported `git archive`)
#   OUT   output dir, removed first if it exists, the build goes to OUT.build (reused, reconfigured each run)
#   -q    build with -DXOMBE_KERNELS=q15
#   -j N  build parallelism (default nproc)
#   -D X  extra cmake definition, e.g. -D XOMBE_DSTAR=OFF (repeatable)
#
# OUT holds one file per run, named <mode>/<vector path, / -> _>[.ext]:
#   dec/    xombe-dec PCM for every .amb (default, --lost 3,10-12, --conformance-noise), .bits, .r33 (--air), .dstar
#           (--dstar, and --lost 3,10-12), stream .jsonl.gz (--stream)
#   enc/    xombe-enc for every .raw: .amb, --air, --dstar, --dual (amb + dstar), and --dstar-dtx on, each with
#           default settings (def) and with --ns off --dtx off (quiet)
#   xc/     xombe-xc forward on every .r33, --reverse on every .dstar
#   status  exit status of every run, a mode the build lacks shows as a nonzero status
# stderr of each run is kept as .err next to its output (it holds --decisions where used)
set -u

kern=
jobs=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
defs=
while getopts qj:D: o; do
    case $o in
    q) kern=-DXOMBE_KERNELS=q15 ;;
    j) jobs=$OPTARG ;;
    D) defs="$defs -D$OPTARG" ;;
    *) sed -n '4,9p' "$0"; exit 2 ;;
    esac
done
shift $((OPTIND - 1))
[ $# -eq 2 ] || { sed -n '4,9p' "$0"; exit 2; }
src=$(cd "$1" && pwd) || exit 2
out=$2
case $out in /*) ;; *) out=$(pwd)/$out ;; esac
bld=$out.build

# shellcheck disable=SC2086
cmake -S "$src" -B "$bld" -DCMAKE_BUILD_TYPE=Release -DXOMBE_BUILD_TOOLS=ON -DXOMBE_BUILD_TESTS=OFF $kern $defs \
    >"$bld.log" 2>&1 &&
    cmake --build "$bld" -j "$jobs" >>"$bld.log" 2>&1 || { echo "snapshot.sh: build failed, see $bld.log" >&2; exit 1; }

tdir=$(dirname "$(find "$bld" -name xombe-dec -type f -perm -u+x | head -1)")
[ -x "$tdir/xombe-dec" ] || { echo "snapshot.sh: no xombe-dec in $bld" >&2; exit 1; }
DEC=$tdir/xombe-dec ENC=$tdir/xombe-enc XC=$tdir/xombe-xc

rm -rf "$out"
mkdir -p "$out/dec" "$out/enc" "$out/xc"
st=$out/status
: >"$st"

# run NAME CMD...: CMD with stderr to NAME.err, its status recorded
run() {
    n=$1; shift
    "$@" 2>"$n.err"
    echo "$? ${n#"$out"/}" >>"$st"
}

cd "$src" || exit 2
key() { echo "$1" | sed 's|^vectors/||; s|/|_|g'; }

for f in $(find vectors -name '*.amb' | sort); do
    k=$out/dec/$(key "$f")
    run "$k.raw" "$DEC" --decisions "$f" "$k.raw"
    run "$k.lost.raw" "$DEC" --lost 3,10-12 "$f" "$k.lost.raw"
    run "$k.conf.raw" "$DEC" --conformance-noise "$f" "$k.conf.raw"
done
for f in $(find vectors -name '*.bits' | sort); do
    k=$out/dec/$(key "$f")
    run "$k.raw" "$DEC" --bits --decisions "$f" "$k.raw"
done
for f in $(find vectors -name '*.r33' | sort); do
    k=$out/dec/$(key "$f")
    run "$k.raw" "$DEC" --air --decisions "$f" "$k.raw"
    k=$out/xc/$(key "$f")
    run "$k.dstar" "$XC" --decisions "$f" "$k.dstar"
done
for f in $(find vectors -name '*.dstar' | sort); do
    k=$out/dec/$(key "$f")
    run "$k.raw" "$DEC" --dstar --decisions "$f" "$k.raw"
    run "$k.lost.raw" "$DEC" --dstar --lost 3,10-12 "$f" "$k.lost.raw"
    k=$out/xc/$(key "$f")
    run "$k.r33" "$XC" --reverse --decisions "$f" "$k.r33"
done
for f in $(find vectors -name '*.jsonl.gz' | sort); do
    k=$out/dec/$(key "$f")
    run "$k.raw" "$DEC" --stream "$f" "$k.raw"
done
for f in $(find vectors -name '*.raw' | sort); do
    for set in def quiet; do
        if [ $set = quiet ]; then o="--ns off --dtx off"; else o=; fi
        k=$out/enc/$(key "$f").$set
        # shellcheck disable=SC2086
        {
            run "$k.amb" "$ENC" $o "$f" "$k.amb"
            run "$k.air" "$ENC" $o --air "$f" "$k.air"
            run "$k.dstar" "$ENC" $o --dstar "$f" "$k.dstar"
            run "$k.dual.amb" "$ENC" $o --dual "$k.dual.dstar" "$f" "$k.dual.amb"
            run "$k.dsdtx.dstar" "$ENC" $o --dstar --dstar-dtx on "$f" "$k.dsdtx.dstar"
        }
    done
done
n=$(wc -l <"$st")
bad=$(awk '$1 != 0' "$st" | wc -l)
echo "snapshot.sh: $n runs into $out ($bad nonzero exits, see status)"
