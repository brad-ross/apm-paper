#!/usr/bin/env bash
#SBATCH --job-name=apm-prepare-figure
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
set -euo pipefail
root=${1:?Repository root required}
shift
cd "$root"
source code/init_server.sh "${SLURM_CPUS_PER_TASK:?}"
exec Rscript --vanilla run.R --mode stage "$@"
