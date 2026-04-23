#!/bin/bash
#SBATCH --job-name=Temp_Files
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=1:00:00
#SBATCH --array=1-19

module load miniforge
mamba activate SL_MA_QC

SL_version="7k"

if [[ -z "$ID" ]]; then
    echo "No ID variable provided."
    exit 1
fi

if grep -q "^${ID}$" 5K.txt; then
    echo "ID $ID is to be run on 5k only"
    SL_version="5k"
fi

TASK=$(sed -n "${SLURM_ARRAY_TASK_ID}p" Cohorts.txt)

python MA/main.py $TASK $ID $SL_version $mode