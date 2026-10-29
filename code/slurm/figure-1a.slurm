#!/usr/bin/env bash
#SBATCH --job-name=apm-figure-1a
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
set -euo pipefail
root=${1:?Repository root required}
shift
cd "$root"
source code/jobs/init_server.sh "${SLURM_CPUS_PER_TASK:?}"
exec bash code/jobs/figure-1a.sh "${1:-resume}"
