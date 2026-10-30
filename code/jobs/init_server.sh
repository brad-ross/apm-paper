#!/usr/bin/env bash
# Cluster-only environment. Pass threads explicitly from the Slurm allocation.
module load "${REPLICATION_R_MODULE:-R/4.3.0}"
export APM_MAX_THREADS="${1:-40}"
export APM_SEED="${APM_SEED:-9001}"
export R_LIBS_USER="${REPLICATION_LIBRARY:-$PWD/.replication/library}"
