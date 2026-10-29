#!/usr/bin/env bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/common.sh"
begin_stage decomposition
touch "$state/$stage.started"
"${R:-Rscript}" --vanilla "code/analysis/match_outcome_decomposition.R"
finish_stage
