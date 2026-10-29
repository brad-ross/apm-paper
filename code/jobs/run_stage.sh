#!/usr/bin/env bash
set -euo pipefail
root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$root"
stage=${1:?Stage required}
mode=${2:-resume}
export DATA_PATH="${DATA_PATH:-$root/data}"
export OUTPUT_PATH="${OUTPUT_PATH:-$root/output}"
export REPLICATION_LIBRARY="${REPLICATION_LIBRARY:-$root/.replication/library}"
export R_LIBS_USER="$REPLICATION_LIBRARY"
export APM_MAX_THREADS="${APM_MAX_THREADS:-1}"
export APM_SEED="${APM_SEED:-9001}"
case "$stage" in
  conversion) script=data/raw_vwh_dta_to_parquet.R ;;
  panels) script=data/construct_clustered_panels.R ;;
  panel-reports) script=reporting/construct_clustered_panels_results.R ;;
  evals) script=analysis/leave_cohort_outcome_out_evals.R ;;
  eval-reports) script=reporting/leave_cohort_outcome_out_evals_results.R ;;
  decomposition) script=analysis/match_outcome_decomposition.R ;;
  decomp-reports) script=reporting/match_outcome_decomposition_results.R ;;
  prepare-figure) script=data/construct_driver_panel_grid.R ;;
  figure-1a) script=reporting/driver_panel_results.R ;;
  *) echo "Unknown stage: $stage" >&2; exit 1 ;;
esac
state="$OUTPUT_PATH/.replication/stamps"
mkdir -p "$state"
# Resume interrupted evaluations only if their upstream panel has not rerun.
# Explicit force starts fresh. No code/configuration/package hashing is done.
if [[ "$stage" == evals ]]; then
  export REPLICATION_EVAL_STARTED="$state/evals.started"
  export REPLICATION_RESUME_EVALS=0
  if [[ "$mode" != fresh && -f "$state/evals.started" &&
        -f "$state/panels.done" && ! "$state/panels.done" -nt "$state/evals.started" ]]; then
    export REPLICATION_RESUME_EVALS=1
  fi
  if [[ "$REPLICATION_RESUME_EVALS" == 0 ]]; then rm -f "$state/evals.started"; fi
fi
rm -f "$state/$stage.done"
if [[ "$stage" != evals ]]; then touch "$state/$stage.started"; fi
"${R:-Rscript}" --vanilla "code/$script"
# Never mark an unsuccessful script complete.
touch "$state/$stage.done.tmp.$$"
mv "$state/$stage.done.tmp.$$" "$state/$stage.done"
rm -f "$state/$stage.started"
