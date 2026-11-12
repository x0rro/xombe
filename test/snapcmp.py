#!/usr/bin/env python3
"""usage: snapcmp.py [-v] A B, compares two test/snapshot.sh outputs, exit 0 when identical

per output: frames changed for each bitstream (AMBE+2 .amb, air, D-STAR, dual, xc, xr) and samples changed with the
largest difference for each decoder output, then totals per output kind. -v also lists unchanged files"""
import os
import sys
from array import array

AMB_HEADER = b".amb"
AMB_RECORD = 8
AIR_RECORD = 9  # air frames, D-STAR frames, and rate-33 frames are all 72 bits


def kind(path):
    """output kind and frame size (0: 16-bit PCM, None: compared as text) of a snapshot file"""
    top, name = path.split("/", 1)
    if name.endswith(".err"):
        return "stderr", None
    if top == "dec":
        return "decoder", 0
    if top == "xc":
        return ("xc", AIR_RECORD) if name.endswith(".dstar") else ("xr", AIR_RECORD)
    for suffix, k, n in ((".dual.amb", "dual AMBE+2", AMB_RECORD), (".dual.dstar", "dual D-STAR", AIR_RECORD),
                         (".dsdtx.dstar", "D-STAR dtx", AIR_RECORD), (".dstar", "D-STAR", AIR_RECORD),
                         (".air", "AMBE+2 air", AIR_RECORD), (".amb", "AMBE+2", AMB_RECORD)):
        if name.endswith(suffix):
            return k, n
    return "other", None


def frames(data, n):
    if n == AMB_RECORD and data.startswith(AMB_HEADER):
        data = data[len(AMB_HEADER):]
    return [data[i:i + n] for i in range(0, len(data), n)]


def pcm(data):
    a = array("h")
    a.frombytes(data[:len(data) & ~1])
    if sys.byteorder != "little":
        a.byteswap()
    return a


def compare(a, b, n):
    """(changed, total, detail) for one file pair"""
    if n is None:
        la, lb = a.decode(errors="replace").splitlines(), b.decode(errors="replace").splitlines()
        changed = sum(x != y for x, y in zip(la, lb)) + abs(len(la) - len(lb))
        return changed, max(len(la), len(lb)), "lines"
    if n == 0:
        sa, sb = pcm(a), pcm(b)
        diffs = [abs(x - y) for x, y in zip(sa, sb) if x != y]
        extra = abs(len(sa) - len(sb))
        detail = "samples, max diff %d" % max(diffs) if diffs else "samples"
        if extra:
            detail += ", length %d vs %d" % (len(sa), len(sb))
        return len(diffs) + extra, max(len(sa), len(sb)), detail
    fa, fb = frames(a, n), frames(b, n)
    changed = sum(x != y for x, y in zip(fa, fb)) + abs(len(fa) - len(fb))
    detail = "frames" if len(fa) == len(fb) else "frames, length %d vs %d" % (len(fa), len(fb))
    return changed, max(len(fa), len(fb)), detail


def files(root):
    out = set()
    for d in ("dec", "enc", "xc"):
        top = os.path.join(root, d)
        if os.path.isdir(top):
            out.update(d + "/" + f for f in os.listdir(top))
    return out


def read(path):
    with open(path, "rb") as f:
        return f.read()


def main(argv):
    verbose = "-v" in argv
    args = [a for a in argv if a != "-v"]
    if len(args) != 2:
        sys.exit(__doc__)
    ra, rb = args
    totals = {}
    differ = 0
    for path in sorted(files(ra) | files(rb)):
        k, n = kind(path)
        pa, pb = os.path.join(ra, path), os.path.join(rb, path)
        if not os.path.exists(pa) or not os.path.exists(pb):
            print("%s: only in %s" % (path, ra if os.path.exists(pa) else rb))
            differ = 1
            continue
        a, b = read(pa), read(pb)
        t = totals.setdefault(k, [0, 0, 0, 0, 0])  # files, files changed, units changed, units, max diff
        t[0] += 1
        if a == b:
            if verbose:
                print("%s: same" % path)
            continue
        changed, total, detail = compare(a, b, n)
        t[1] += 1
        t[2] += changed
        t[3] += total
        if n == 0 and "max diff" in detail:
            t[4] = max(t[4], int(detail.split("max diff ")[1].split(",")[0]))
        print("%s: %d of %d %s" % (path, changed, total, detail))
        differ = 1
    sa, sb = read(os.path.join(ra, "status")), read(os.path.join(rb, "status"))
    if sa != sb:
        print("status differs")
        differ = 1
    print("summary:")
    for k in sorted(totals):
        f, fc, uc, u, md = totals[k]
        line = "  %-12s %3d files, %3d changed" % (k, f, fc)
        if fc:
            unit = "lines" if kind_unit(k) is None else ("samples" if kind_unit(k) == 0 else "frames")
            line += ", %d of %d %s in those" % (uc, u, unit)
            if kind_unit(k) == 0:
                line += ", max diff %d" % md
        print(line)
    print("identical" if not differ else "DIFFERENT")
    return differ


def kind_unit(k):
    return {"decoder": 0, "stderr": None, "other": None}.get(k, AIR_RECORD)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
