#!/usr/bin/env bash
#SBATCH --job-name=apm-panels
#SBATCH --cpus-per-task=40
#SBATCH --mem=50G
#SBATCH --time=24:00:00
set -euo pipefail
root=${1:?Repository root required}
shift
cd "$root"
source code/jobs/init_server.sh "${SLURM_CPUS_PER_TASK:?}"
exec bash code/jobs/panels.sh "${1:-resume}"
