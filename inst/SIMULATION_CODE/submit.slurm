#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=1:00:00
#SBATCH --output=status_reports/output_1_%a.log
#SBATCH --array=1-2000

module load r/4.4.0

# Create output directories
mkdir -p status_reports
mkdir -p results

lscpu

# Simulation setting (adjust if running multiple settings)
SIM_SETTING=1

echo "Job Array ID: $SLURM_ARRAY_JOB_ID"
echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "Running on node: $(hostname)"
echo "Simulation setting: $SIM_SETTING, iteration: $SLURM_ARRAY_TASK_ID"

# Record start time
START_TIME=$(date +%s)
echo "Start time: $(date)"

/usr/bin/time -v Rscript --no-save run_sim.R $SIM_SETTING $SLURM_ARRAY_TASK_ID

# Record end time and calculate elapsed
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo "End time: $(date)"
echo "Elapsed time: $((ELAPSED / 3600))h $((ELAPSED % 3600 / 60))m $((ELAPSED % 60))s"