#!/usr/bin/env bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/common.sh"
begin_stage figure-1a
touch "$state/$stage.started"
"${R:-Rscript}" --vanilla "code/reporting/driver_panel_results.R"
finish_stage
