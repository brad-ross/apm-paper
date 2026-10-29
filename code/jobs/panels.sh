#!/usr/bin/env bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/common.sh"
begin_stage panels
touch "$state/$stage.started"
"${R:-Rscript}" --vanilla "code/data/construct_clustered_panels.R"
finish_stage
