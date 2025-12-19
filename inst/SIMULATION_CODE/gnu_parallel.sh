#!/bin/bash

I_MAX=1
J_MIN=1
J_MAX=2000

mkdir -p status_reports
mkdir -p results

run_task() {
    i=$1
    j=$2
    echo "Running task i=$i j=$j"
    
    Rscript --no-save run_sim.R $i $j > "status_reports/output_${i}_${j}.log" 2>&1
}

export -f run_task

parallel -j 14 --progress --joblog status_reports/parallel.log \
    --header : run_task {i} {j} \
    ::: j $(seq $J_MIN $J_MAX) \
    ::: i $(seq 1 $I_MAX)