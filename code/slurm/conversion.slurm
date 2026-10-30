#!/usr/bin/env bash
#SBATCH --job-name=apm-conversion
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
set -euo pipefail
root=${1:?Repository root required}
shift
cd "$root"
source code/jobs/init_server.sh "${SLURM_CPUS_PER_TASK:?}"
exec bash code/jobs/conversion.sh "${1:-resume}"
