#!/usr/bin/env bash
set -euo pipefail
root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$root"
export DATA_PATH="${DATA_PATH:-$root/data}"
export OUTPUT_PATH="${OUTPUT_PATH:-$root/output}"
export REPLICATION_LIBRARY="${REPLICATION_LIBRARY:-$root/.replication/library}"
export R_LIBS_USER="$REPLICATION_LIBRARY"
export APM_MAX_THREADS="${APM_MAX_THREADS:-1}"
export APM_SEED="${APM_SEED:-9001}"
state="$OUTPUT_PATH/.replication/stamps"
mkdir -p "$state"

begin_stage() {
  stage=$1
  rm -f "$state/$stage.done"
}

finish_stage() {
  touch "$state/$stage.done.tmp.$$"
  mv "$state/$stage.done.tmp.$$" "$state/$stage.done"
  rm -f "$state/$stage.started"
}
