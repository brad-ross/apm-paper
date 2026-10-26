#!/usr/bin/env bash
set -euo pipefail
root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
cd "$root"
source code/init_server.sh "${APM_MAX_THREADS:-1}"
exec Rscript --vanilla run.R --mode submit "$@"
