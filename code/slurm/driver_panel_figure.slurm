#!/usr/bin/env bash
#SBATCH --job-name=apm-driver-figure
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/dev/null
# Direct one-off submission: sbatch --chdir=/path/to/replication code/slurm/driver_panel_figure.slurm
set -euo pipefail
log_dir="${OUTPUT_PATH:-output}/logs"
mkdir -p "$log_dir"
exec >"$log_dir/driver-panel-figure-${SLURM_JOB_ID:?}.log" 2>&1
if [[ ! -f Makefile ]]; then
  echo "Submit from the replication root or pass sbatch --chdir=/path/to/replication" >&2
  exit 1
fi
source code/jobs/init_server.sh "${SLURM_CPUS_PER_TASK:?}"
exec make private-figure "$@"
