#!/usr/bin/env bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/common.sh"
begin_stage prepare-figure
touch "$state/$stage.started"
"${R:-Rscript}" --vanilla "code/data/construct_driver_panel_grid.R"
finish_stage
