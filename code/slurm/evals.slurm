#!/usr/bin/env bash
#SBATCH --job-name=apm-evals
#SBATCH --cpus-per-task=250
#SBATCH --mem=50G
#SBATCH --time=2-00:00:00
set -euo pipefail
root=${1:?Repository root required}
shift
cd "$root"
source code/jobs/init_server.sh "${SLURM_CPUS_PER_TASK:?}"
exec bash code/jobs/evals.sh "${1:-resume}"
