#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=01:00:00
#SBATCH --output=install_package_%j.log

lscpu

module load r/4.4.0

echo "Job ID: $SLURM_JOB_ID"
echo "Running on node: $(hostname)"
echo "Start time: $(date)"

R -e 'pak::pak("maxdrohde/bayes2stage")'

echo "End time: $(date)"
