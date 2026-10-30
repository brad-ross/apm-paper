#!/usr/bin/env bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/common.sh"
begin_stage decomp_reports
touch "$state/$stage.started"
"${R:-Rscript}" --vanilla "code/reporting/match_outcome_decomposition_results.R"
finish_stage
