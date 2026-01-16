#!/bin/bash
#$ -cwd
#$ -N soma_qc
#$ -pe smp 1
#$ -l h_vmem=64G
#$ -l h_rt=12:00:00
#$ -j y
#$ -t 1-21

module load miniforge
mamba activate SL_MA_QC

TASK=$(sed -n "${SGE_TASK_ID}p" cohorts.txt)

if [[ -z "$TASK" ]]; then
    echo "No task for SGE_TASK_ID=${SGE_TASK_ID}"
    exit 1
fi

echo "Running task ${SGE_TASK_ID}: $TASK"
python main.py $TASK
echo "Task ${SGE_TASK_ID} completed"
