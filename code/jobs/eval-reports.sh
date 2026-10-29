#!/usr/bin/env bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/common.sh"
begin_stage eval-reports
touch "$state/$stage.started"
"${R:-Rscript}" --vanilla "code/reporting/leave_cohort_outcome_out_evals_results.R"
finish_stage
