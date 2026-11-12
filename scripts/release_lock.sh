#!/bin/sh
# release_lock.sh: prints the packages (name=version) that release_deps.sh's snapshot adds to the base image, the
# contents of release_packages.lock. CONTAINER picks the engine (podman)
set -eu
base=docker.io/library/debian@sha256:2dc419c6d7fe3dad0a2f9fd469c34cc90b0f9ba642d8fde837fe3d7d34a0e7f1
here=$(cd "$(dirname "$0")" && pwd)
engine=${CONTAINER:-podman}
snap=$(grep -o 'http://snapshot\.debian\.org/archive/debian/[0-9TZ]*' "$here/release_deps.sh")
"$engine" run --rm "$base" sh -euc '
    export DEBIAN_FRONTEND=noninteractive
    q() { dpkg-query -W -f="\${Package}=\${Version}\n" | LC_ALL=C sort; }
    q > /tmp/base
    rm -f /etc/apt/sources.list.d/*
    echo "deb [check-valid-until=no] $1 trixie main" > /etc/apt/sources.list
    apt-get -qq -o Acquire::Retries=5 update
    apt-get -qq install -y --no-install-recommends ca-certificates git xz-utils > /dev/null
    q | LC_ALL=C comm -13 /tmp/base -
' _ "$snap"
