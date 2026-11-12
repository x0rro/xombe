#!/bin/sh
# size_table.sh [WORKDIR]: state sizes (test/sizes.c) and worst-case stack per public call for the four builds
# (float and q15 kernels, XOMBE_DSTAR off and on), as Markdown tables. stack is from -fcallgraph-info=su
# (GCC 10 or later) walked over the library's call graph, on the host and, when arm-none-eabi-gcc is installed,
# on a Cortex-M33 build. WORKDIR (default: a temporary directory, removed at the end) holds the builds.
# needs cmake and python3, CC picks the host compiler, ARM_CC the cross compiler
set -eu
ROOT=$(cd "$(dirname "$0")/.." && pwd)
PY=${PYTHON:-python3}
ARM_CC=${ARM_CC:-arm-none-eabi-gcc}
JOBS=${JOBS:-4}
if [ $# -gt 0 ]; then
    W=$1
    mkdir -p "$W"
else
    W=$(mktemp -d)
    trap 'rm -rf "$W"' EXIT
fi
CI="-fstack-usage -fcallgraph-info=su"
ARM_FLAGS="-mcpu=cortex-m33 -mthumb -mfloat-abi=hard -mfpu=fpv5-sp-d16"
CONFIGS="float:OFF float:ON q15:OFF q15:ON"
command -v "$ARM_CC" > /dev/null && arm=1 || arm=0

build() { # build DIR KERNELS DSTAR [cmake options]
    d=$1 k=$2 ds=$3
    shift 3
    cmake -S "$ROOT" -B "$d" -DCMAKE_BUILD_TYPE=Release -DXOMBE_KERNELS=$k -DXOMBE_DSTAR=$ds -DXOMBE_BUILD_TOOLS=OFF \
        -DXOMBE_INSTALL=OFF "$@" > "$d.log" 2>&1 && cmake --build "$d" -j "$JOBS" >> "$d.log" 2>&1 ||
        { cat "$d.log" >&2; exit 1; }
}

for c in $CONFIGS; do
    k=${c%:*} ds=${c#*:}
    build "$W/host-$k-$ds" $k $ds -DXOMBE_BUILD_TESTS=ON "-DCMAKE_C_FLAGS=$CI"
    if [ $arm = 1 ]; then
        build "$W/m33-$k-$ds" $k $ds -DXOMBE_BUILD_TESTS=OFF -DCMAKE_SYSTEM_NAME=Generic -DCMAKE_C_COMPILER="$ARM_CC" \
            -DCMAKE_TRY_COMPILE_TARGET_TYPE=STATIC_LIBRARY "-DCMAKE_C_FLAGS_RELEASE=-O2 -DNDEBUG" "-DCMAKE_C_FLAGS=$ARM_FLAGS $CI"
    fi
done

echo "### State (bytes) and encoder delay (frames), $(uname -m)"
echo
"$W/host-float-OFF/sizes" --header
for c in $CONFIGS; do
    "$W/host-${c%:*}-${c#*:}/sizes"
done
echo
targets="x86-64"
[ $arm = 1 ] && targets="$targets M33"
echo "### Worst-case stack per public call (bytes)"
echo
"$PY" - "$ROOT" "$W" "$targets" $CONFIGS << 'EOF'
import glob, os, re, sys

root, work, targets, configs = sys.argv[1], sys.argv[2], sys.argv[3].split(), sys.argv[4:]
node_re = re.compile(r'node: \{ title: "([^"]*)" label: "([^"]*)"')
edge_re = re.compile(r'edge: \{ sourcename: "([^"]*)" targetname: "([^"]*)"(?: label: "([^"]*)")?')
comment_re = re.compile(r'/\*.*?\*/|//[^\n]*', re.S)


def short(title):
    return title.rsplit(':', 1)[-1]


def address_taken(stack):
    """functions named in the library sources other than in a call or declaration: the targets an indirect call
    may reach, for a static function only in its own file"""
    text = {}
    for f in glob.glob(os.path.join(root, 'src', '**', '*.[ch]'), recursive=True):
        text[os.path.realpath(f)] = comment_re.sub(' ', open(f).read())
    out = set()
    for t in stack:
        name = short(t)
        files = [os.path.realpath(t.rsplit(':', 1)[0])] if ':' in t else list(text)
        pat = re.compile(r'(?<![.>\w])%s\b(?!\s*\()' % re.escape(name))
        if any(pat.search(text.get(f, '')) for f in files):
            out.add(t)
    return out


def load(build):
    """call graph of the library's objects: {fn: bytes}, {fn: [callees]}, flags, unresolved indirect calls"""
    stack, calls, flag, src_lines, where_def = {}, {}, {}, {}, {}
    for ci in glob.glob(os.path.join(build, 'CMakeFiles', '*', 'src', '**', '*.ci'), recursive=True):
        for line in open(ci):
            m = node_re.match(line)
            if m:
                parts = m.group(2).split('\\n')
                if len(parts) > 2:
                    n, kind = re.match(r'(\d+) bytes \((.*)\)', parts[2]).groups()
                    t = m.group(1)
                    stack[t] = max(stack.get(t, 0), int(n))
                    where_def[t] = parts[1].rsplit(':', 2)[0]
                    if kind == 'dynamic':
                        flag[t] = 'unbounded'
                    calls.setdefault(t, set())
                continue
            m = edge_re.match(line)
            if m:
                src, dst, where = m.groups()
                if dst == '__indirect_call':
                    dst = 'indirect@' + (where or '?')
                calls.setdefault(src, set()).add(dst)
    by_short = {}
    for t in stack:
        by_short.setdefault(short(t), []).append(t)
    taken = address_taken(stack)
    unresolved = set()
    for src in list(calls):
        out = set()
        for dst in calls[src]:
            if dst.startswith('indirect@'):
                f, l = dst[len('indirect@'):].rsplit(':', 2)[:2]
                text = src_lines.setdefault(f, open(f).read().split('\n') if os.path.exists(f) else [])
                text = text[int(l) - 1] if int(l) <= len(text) else ''
                if 'XOMBE_TAP' in text:
                    continue  # test tap callback, runs on the caller's own stack
                if taken:
                    out.update(taken)
                else:
                    unresolved.add('%s:%s' % (os.path.relpath(f, root), l))
            elif dst in stack:
                out.add(dst)
            elif dst in by_short and len(by_short[dst]) == 1:
                out.add(by_short[dst][0])
            # otherwise a library outside xombe (libc, libm), not counted
        calls[src] = out
    return stack, calls, flag, unresolved, taken, where_def


def worst(fn, stack, calls, memo, path):
    if fn in path:
        return None  # recursion
    if fn in memo:
        return memo[fn]
    path.add(fn)
    best = 0
    for c in calls.get(fn, ()):
        w = worst(c, stack, calls, memo, path)
        if w is None:
            path.discard(fn)
            return None
        best = max(best, w)
    path.discard(fn)
    memo[fn] = stack.get(fn, 0) + best
    return memo[fn]


cols, notes = [], []
for t in targets:
    for c in configs:
        k, ds = c.split(':')
        b = os.path.join(work, '%s-%s-%s' % ('host' if t == 'x86-64' else 'm33', k, ds))
        stack, calls, flag, unresolved, taken, where_def = load(b)
        memo, res = {}, {}
        for fn in stack:
            if fn.startswith('xombe_') and not fn.startswith('xombe__') and ':' not in fn:
                w = worst(fn, stack, calls, memo, set())
                res[fn] = 'recursive' if w is None else str(w)
        if ds == 'ON':
            # the DMR-only encoder never enters src/dstar: the same walk without those functions
            dsdir = os.path.join(os.path.realpath(root), 'src', 'dstar') + os.sep
            amb = {f: {c for c in cs if not os.path.realpath(where_def.get(c, '')).startswith(dsdir)}
                   for f, cs in calls.items()}
            for fn in ('xombe_enc_frame', 'xombe_enc_frame_bits'):
                if fn in stack:
                    w = worst(fn, stack, amb, {}, set())
                    res[fn + ', DMR output only'] = 'recursive' if w is None else str(w)
        cols.append(('%s %s D-STAR %s' % (t, k, ds.lower()), res))
        if unresolved:
            notes.append('%s %s D-STAR %s: unresolved indirect calls at %s' % (t, k, ds.lower(),
                                                                                  ', '.join(sorted(unresolved))))
        notes.append('%s %s D-STAR %s: indirect calls bounded by %s' % (t, k, ds.lower(),
                                                                         ', '.join(sorted(map(short, taken)))))
        if flag:
            notes.append('%s %s D-STAR %s: unbounded dynamic stack in %s' % (t, k, ds.lower(),
                                                                                ', '.join(sorted(flag))))
names = sorted(set().union(*[set(r) for _, r in cols]))
print('| call | ' + ' | '.join(c for c, _ in cols) + ' |')
print('|---|' + '---|' * len(cols))
for n in names:
    fn, _, what = n.partition(', ')
    print('| `%s`%s | ' % (fn, ' (%s)' % what if what else '') + ' | '.join(r.get(n, '-') for _, r in cols) + ' |')
print()
print('stack of libc and libm callees not included, test tap callbacks not included')
for n in notes:
    print(n)
EOF
