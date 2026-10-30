#!/usr/bin/env bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/common.sh"
begin_stage conversion
touch "$state/$stage.started"
"${R:-Rscript}" --vanilla "code/data/raw_vwh_dta_to_parquet.R"
finish_stage
