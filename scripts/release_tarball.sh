#!/bin/sh
# release_tarball.sh [REF [OUT]]: writes xombe-X.Y.Z.tar.xz (git archive of REF, default HEAD) and
# xombe-X.Y.Z.buildinfo to OUT (default dist). byte-identical across hosts in the release container
# (release_deps.sh), the version is REF's include/xombe.h
set -eu
ref=${1:-HEAD}
root=$(cd "$(dirname "$0")/.." && pwd)
out=${2:-$root/dist}
part() { git -C "$root" show "$ref:include/xombe.h" | sed -n "s/^#define XOMBE_VERSION_$1 \([0-9][0-9]*\)\$/\1/p"; }
ver=$(part MAJOR).$(part MINOR).$(part PATCH)
case $ver in
*[!0-9.]* | .* | *..* | *.) echo "no XOMBE_VERSION_MAJOR/MINOR/PATCH in $ref:include/xombe.h" >&2; exit 1 ;;
esac
name=xombe-$ver
mkdir -p "$out"
git -C "$root" archive --format=tar --prefix="$name/" "$ref" -- . ':!.gitignore' | xz -9 -T1 > "$out/$name.tar.xz"
"$root/scripts/release_buildinfo.sh" > "$out/$name.buildinfo"
