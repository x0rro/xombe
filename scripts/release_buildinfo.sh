#!/bin/sh
# release_buildinfo.sh: the tools behind a release tarball and, on Debian, every installed package
set -eu
git --version
xz --version | head -n 1
if command -v dpkg-query > /dev/null; then
    dpkg-query -W -f='${db:Status-Abbrev} ${Package}:${Architecture} ${Version}\n' |
        awk '$1 == "ii" { print $2, $3 }' | LC_ALL=C sort
fi
