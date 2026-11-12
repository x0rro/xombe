#!/bin/sh
# tools.sh BUILD VECTORS: command-line contracts of xombe-enc, xombe-dec, and xombe-xc
B=$1; V=$2; T=$B/tools_out; mkdir -p "$T"; fails=0
expect() { # expect STATUS DESC CMD...
    want=$1; what=$2; shift 2
    "$@" > "$T/out.txt" 2> "$T/err.txt"; got=$?
    if [ "$got" -ne "$want" ]; then echo "FAIL $what: exit $got, want $want"; cat "$T/err.txt"; fails=$((fails + 1)); fi
}
have() { [ -x "$B/$1" ]; }
partial() { # partial WHO: the last command reported the partial final frame
    grep -q "error: partial final frame" "$T/err.txt" || { echo "FAIL $1 partial frame message"; fails=$((fails + 1)); }
}
for t in xombe-enc xombe-dec xombe-xc; do
    have $t || continue
    expect 0 "$t --help" "$B/$t" --help
    grep -q "usage: $t" "$T/out.txt" || { echo "FAIL $t --help text"; fails=$((fails + 1)); }
    expect 2 "$t unknown option" "$B/$t" --nope a b
    expect 2 "$t missing paths" "$B/$t"
done
head -c 1001 "$V/dstar/sp.in.raw" > "$T/odd.raw"
head -c 100 "$V/dstar/sp.dstar" > "$T/odd.dstar"
if have xombe-enc; then
    expect 0 "enc partial final frame" "$B/xombe-enc" "$T/odd.raw" "$T/odd.amb"
    grep -q "partial final frame" "$T/err.txt" || { echo "FAIL enc partial frame warning"; fails=$((fails + 1)); }
    expect 2 "enc bad on/off" "$B/xombe-enc" --ns maybe "$T/odd.raw" "$T/x"
    expect 2 "enc bad number" "$B/xombe-enc" --vscale 0 "$T/odd.raw" "$T/x"
    expect 2 "enc bad quantizer" "$B/xombe-enc" --quant maybe "$T/odd.raw" "$T/x"
    expect 0 "enc open quantizer" "$B/xombe-enc" --quant open "$T/odd.raw" "$T/x"
    if "$B/xombe-enc" --help | grep -q "^  --dstar .*not in this build"; then
        expect 2 "enc --dstar without D-STAR" "$B/xombe-enc" --dstar "$T/odd.raw" "$T/x"
    else
        expect 2 "enc --dstar with --air" "$B/xombe-enc" --dstar --air "$T/odd.raw" "$T/x"
        expect 0 "enc --dual" "$B/xombe-enc" --dual "$T/dual.dstar" "$T/odd.raw" "$T/dual.amb"
        [ "$(wc -c < "$T/dual.dstar")" -eq 27 ] || { echo "FAIL enc --dual frame count"; fails=$((fails + 1)); }
    fi
fi
if have xombe-dec; then
    expect 2 "dec two input formats" "$B/xombe-dec" --bits --air "$T/odd.dstar" "$T/x"
    expect 2 "dec bad --lost" "$B/xombe-dec" --lost 3-1 "$T/odd.dstar" "$T/x"
    head -c 15 "$V/e3/es1.ref.d0.amb" > "$T/odd.amb"
    expect 1 "dec truncated .amb frame" "$B/xombe-dec" "$T/odd.amb" "$T/x"
    partial "dec"
    if "$B/xombe-dec" --help | grep -q "^  --dstar .*not in this build"; then
        expect 2 "dec --dstar without D-STAR" "$B/xombe-dec" --dstar "$T/odd.dstar" "$T/x"
    else
        expect 1 "dec truncated D-STAR frame" "$B/xombe-dec" --dstar "$T/odd.dstar" "$T/x"
        partial "dec --dstar"
        expect 0 "dec --lost" "$B/xombe-dec" --dstar --decisions --lost 2,5-6 "$V/dstar/sp.dstar" "$T/x"
        [ "$(grep -c " erasure " "$T/err.txt")" -eq 3 ] || { echo "FAIL dec --lost decisions"; fails=$((fails + 1)); }
    fi
fi
if have xombe-xc; then
    expect 1 "xc truncated frame" "$B/xombe-xc" "$T/odd.dstar" "$T/x"
    partial "xc"
    expect 0 "xc --lost" "$B/xombe-xc" --decisions --lost 1 "$V/dstar/sp.r33" "$T/x"
    [ "$(grep -c " erasure " "$T/err.txt")" -eq 1 ] || { echo "FAIL xc --lost decisions"; fails=$((fails + 1)); }
    if "$B/xombe-xc" --help | grep -q "^  --reverse .*not in this build"; then
        expect 2 "xc --reverse without the encoder" "$B/xombe-xc" --reverse "$V/dstar/sp.dstar" "$T/x"
    else
        expect 0 "xc --reverse --lost" "$B/xombe-xc" --reverse --decisions --lost 0 "$V/dstar/sp.dstar" "$T/x"
        [ "$(grep -c " erasure " "$T/err.txt")" -eq 1 ] || { echo "FAIL xc --reverse --lost"; fails=$((fails + 1)); }
        expect 0 "xc --reverse --amb" "$B/xombe-xc" --reverse --amb "$V/dstar/sper.dstar" "$T/rt.amb"
        expect 0 "xc --amb" "$B/xombe-xc" --amb "$T/rt.amb" "$T/rt_amb.dstar"
        "$B/xombe-xc" --reverse "$V/dstar/sper.dstar" "$T/rt.r33" && "$B/xombe-xc" "$T/rt.r33" "$T/rt_r33.dstar"
        cmp -s "$T/rt_amb.dstar" "$T/rt_r33.dstar" || { echo "FAIL xc --amb round trip differs from r33"; fails=$((fails + 1)); }
    fi
    expect 1 "xc --amb on a non-.amb file" "$B/xombe-xc" --amb "$V/dstar/sp.r33" "$T/x"
fi
[ $fails -eq 0 ] && echo "tools: ok"
exit $fails
