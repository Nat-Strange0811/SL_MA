#!/bin/bash
#$ -cwd
#$ -N Comparison
#$ -pe smp 1
#$ -l h_vmem=20G
#$ -l h_rt=1:00:00
#$ -o MA_Comparison/Logs/
#$ -j y
#$ -t 1-300

TASK=$(sed -n "${SGE_TASK_ID}p" MA/IDs.txt)

if [[ -z "$TASK" ]]; then
    echo "No TASK variable found for SGE_TASK_ID=${SGE_TASK_ID}."
    exit 1
fi

echo "Running comparison for ID: $TASK"

module load miniforge
mamba activate SL_MA_QC

python MA_Comparison/main.py $TASK