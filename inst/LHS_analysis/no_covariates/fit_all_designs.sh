#!/bin/bash

# Fit all designs x distributions (no covariates spec).
# full and srs_no_imp run once (normal only) since all x are observed.
# 16 GNU Parallel jobs x 1 sequential Stan chain = 16 cores.
#
# Usage:
#   ./fit_all_designs.sh        # default 16 jobs
#   ./fit_all_designs.sh 4      # override job count

cd "$(dirname "$0")"

N_CORES=${1:-16}

mkdir -p cr_results
mkdir -p logs

# Designs that need imputation: run for all 3 distributions
IMP_DESIGNS="srs ods_intercept ods_slope ods_ellipse bds_intercept bds_slope bds_ellipse"
DISTRIBUTIONS="normal beta_binomial negative_binomial"

# Designs with all x observed: run once (normal only)
NO_IMP_DESIGNS="full srs_no_imp"

# Build job list
JOBFILE=$(mktemp)
for d in $NO_IMP_DESIGNS; do echo "$d normal" >> "$JOBFILE"; done
for d in $IMP_DESIGNS; do for dist in $DISTRIBUTIONS; do echo "$d $dist" >> "$JOBFILE"; done; done

N_TOTAL=$(wc -l < "$JOBFILE" | tr -d ' ')

echo "========================================"
echo "Credible Regions: No Covariates"
echo "  Jobs:  $N_TOTAL"
echo "  Cores: $N_CORES"
echo "========================================"

parallel -j "$N_CORES" --progress --joblog logs/cr_parallel.log \
    --colsep ' ' \
    'Rscript --no-save fit_design.R {1} {2} > "logs/cr_{1}_{2}.log" 2>&1' \
    :::: "$JOBFILE"

rm "$JOBFILE"

echo "========================================"
echo "Done! Results in ./cr_results/"
echo "========================================"
