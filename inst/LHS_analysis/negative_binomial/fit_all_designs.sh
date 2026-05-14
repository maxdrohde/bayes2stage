#!/bin/bash

# Fit all 8 designs in parallel (negative-binomial imputation).
#
# Usage:
#   ./fit_all_designs.sh        # 8 cores (one per design)
#   ./fit_all_designs.sh 4      # override core count

cd "$(dirname "$0")"

N_CORES=${1:-8}

DESIGNS="full srs ods_intercept ods_slope ods_ellipse bds_intercept bds_slope bds_ellipse"

mkdir -p cr_results
mkdir -p logs

echo "========================================"
echo "Credible Regions: Negative-Binomial Imputation"
echo "  Designs: $(echo $DESIGNS | wc -w | tr -d ' ')"
echo "  Cores:   $N_CORES"
echo "========================================"

parallel -j "$N_CORES" --progress --joblog logs/cr_parallel.log \
    'Rscript --no-save fit_design.R {1} > "logs/cr_{1}.log" 2>&1' \
    ::: $DESIGNS

echo "========================================"
echo "Done! Results in ./cr_results/"
echo "========================================"
