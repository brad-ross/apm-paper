#!/usr/bin/env bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/common.sh"
begin_stage panel-reports
touch "$state/$stage.started"
"${R:-Rscript}" --vanilla "code/reporting/construct_clustered_panels_results.R"
finish_stage
