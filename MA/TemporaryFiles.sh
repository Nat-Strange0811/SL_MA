#!/bin/bash
#$ -cwd
#$ -N Temp_Files
#$ -pe smp 1
#$ -l h_vmem=32G
#$ -l h_rt=1:00:00
#$ -j y
#$ -t 1-19

module load miniforge
mamba activate SL_MA_QC

if [[ -z "$ID" ]]; then
    echo "No ID variable provided."
    exit 1
fi

TASK=$(sed -n "${SGE_TASK_ID}p" MA/Cohorts.txt)

python MA/main.py $TASK $ID