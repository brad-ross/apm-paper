#!/bin/bash

ml R/4.3.0

export DATA_PATH=~/projects/apm/replication/data
export OUTPUT_PATH=~/projects/apm/replication/output

export APM_MAX_THREADS=${1:-40}
export APM_SEED="${APM_SEED:-9001}"

export FIRST_YEAR=1998
export LAST_YEAR=2001
export MAX_RANK=2
export MIN_COHORT_SIZE=75
export MIN_K=1
export MAX_K=25
export CHOSEN_K=3
