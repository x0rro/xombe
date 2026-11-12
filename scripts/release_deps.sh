#!/bin/sh
# release_deps.sh: installs release_packages.lock into the release container (release_lock.sh's base), from
# deb.debian.org or else the snapshot
set -eu
here=$(cd "$(dirname "$0")" && pwd)
rm -f /etc/apt/sources.list.d/*
printf '%s\n' 'deb http://deb.debian.org/debian trixie main' \
    'deb [check-valid-until=no] http://snapshot.debian.org/archive/debian/20261111T000000Z trixie main' \
    > /etc/apt/sources.list
apt-get -o Acquire::Retries=5 update || true
# shellcheck disable=SC2046
DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends $(cat "$here/release_packages.lock")
rm -rf /var/lib/apt/lists/*
