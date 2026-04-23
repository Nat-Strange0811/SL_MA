#!/bin/bash
#SBATCH --job-name=Clean_Up
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=1:00:00

module load miniforge
mamba activate SL_MA_QC

if [[ -z "$ID" ]]; then
    echo "No ID variable provided."
    exit 1
fi

python MA/cleanUp.py $ID $mode