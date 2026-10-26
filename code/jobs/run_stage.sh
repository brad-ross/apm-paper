#!/usr/bin/env bash
set -euo pipefail
root=${1:?Repository root required}
script=${2:?Scientific script required}
cd "$root"
export R_LIBS_USER="${REPLICATION_LIBRARY:-$root/.replication/library}"
exec Rscript --vanilla "$script"
