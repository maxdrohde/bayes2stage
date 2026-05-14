#!/bin/bash

# Fit all credible region models: normal, beta-binomial, negative-binomial.
# Each runs 8 designs in parallel (2 parallel Stan chains each = 16 cores).

cd "$(dirname "$0")"

rm -f cr_results/*.parquet logs/cr_*.log
rm -f beta_binomial/cr_results/*.parquet beta_binomial/logs/cr_*.log
rm -f negative_binomial/cr_results/*.parquet negative_binomial/logs/cr_*.log

echo "=== Normal imputation ==="
./fit_all_designs.sh

echo ""
echo "=== Beta-binomial imputation ==="
cd beta_binomial && ./fit_all_designs.sh && cd ..

echo ""
echo "=== Negative-binomial imputation ==="
cd negative_binomial && ./fit_all_designs.sh
