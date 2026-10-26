#!/usr/bin/env bash
#SBATCH --job-name=apm-panel-reports
#SBATCH --cpus-per-task=16
#SBATCH --mem=50G
#SBATCH --time=12:00:00
set -euo pipefail
root=${1:?Repository root required}
shift
cd "$root"
source code/init_server.sh "${SLURM_CPUS_PER_TASK:?}"
exec Rscript --vanilla run.R --mode stage "$@"
