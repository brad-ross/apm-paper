#!/usr/bin/env bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/common.sh"
begin_stage figure_1a_data_prep
touch "$state/$stage.started"
"${R:-Rscript}" --vanilla "code/data/construct_driver_panel_grid.R"
finish_stage
