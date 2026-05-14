#!/bin/bash

# Run LHS replicated design comparison in parallel.
#
# rep_id=0: full-data models (all subjects measured)
# rep_id=1..N_REPS: subsampled replicates
# All run in parallel via GNU Parallel.
#
# Usage:
#   ./gnu_parallel.sh          # uses N_REPS and N_CORES from lhs_config.R
#   ./gnu_parallel.sh 50       # override: 50 reps
#   ./gnu_parallel.sh 50 14    # override: 50 reps, 14 cores

cd "$(dirname "$0")"

# Read defaults from lhs_config.R (single source of truth)
# --vanilla skips .Rprofile (renv activation) to avoid stdout contamination
CONFIG_N_REPS=$(Rscript --vanilla -e "source('lhs_config.R'); cat(N_REPS)")
CONFIG_N_CORES=$(Rscript --vanilla -e "source('lhs_config.R'); cat(N_CORES)")

N_REPS=${1:-$CONFIG_N_REPS}
N_CORES=${2:-$CONFIG_N_CORES}

mkdir -p results
mkdir -p logs

echo "========================================"
echo "LHS Replicated Design Comparison"
echo "  Replicates: $N_REPS"
echo "  Cores:      $N_CORES"
echo "========================================"

# Run rep 0 (full-data) + reps 1..N_REPS (replicates) in parallel
echo "Starting full-data fit + $N_REPS replicates on $N_CORES cores..."

parallel -j "$N_CORES" --progress --joblog logs/parallel.log \
    'Rscript --no-save run_rep.R {1} > "logs/rep_{1}.log" 2>&1' \
    ::: $(seq 0 "$N_REPS")

echo "========================================"
echo "Done! Results in ./results/"
echo "Logs in ./logs/"
echo "========================================"
