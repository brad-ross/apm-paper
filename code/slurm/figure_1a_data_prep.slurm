#!/usr/bin/env bash
#SBATCH --job-name=apm-figure_1a_data_prep
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
set -euo pipefail
root=${1:?Repository root required}
shift
cd "$root"
source code/jobs/init_server.sh "${SLURM_CPUS_PER_TASK:?}"
exec bash code/jobs/figure_1a_data_prep.sh "${1:-resume}"
