#!/usr/bin/env bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/common.sh"
mode=${1:-resume}
export REPLICATION_EVAL_STARTED="$state/evals.started"
export REPLICATION_RESUME_EVALS=0
if [[ "$mode" != fresh && -f "$state/evals.started" &&
      -f "$state/panels.done" && ! "$state/panels.done" -nt "$state/evals.started" ]]; then
  export REPLICATION_RESUME_EVALS=1
fi
if [[ "$REPLICATION_RESUME_EVALS" == 0 ]]; then rm -f "$state/evals.started"; fi
begin_stage evals
"${R:-Rscript}" --vanilla "code/analysis/leave_cohort_outcome_out_evals.R"
finish_stage
