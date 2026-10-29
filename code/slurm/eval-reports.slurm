#!/usr/bin/env bash
#SBATCH --job-name=apm-eval-reports
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
set -euo pipefail
root=${1:?Repository root required}
shift
cd "$root"
source code/jobs/init_server.sh "${SLURM_CPUS_PER_TASK:?}"
exec bash code/jobs/eval-reports.sh "${1:-resume}"
