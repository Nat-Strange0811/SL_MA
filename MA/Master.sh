#!/bin/bash
#SBATCH --job-name=Master
#SBATCH --time=1:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-7578

# Load task ID
TASK=$(sed -n "${SLURM_ARRAY_TASK_ID}p" Proteins.txt)

if [[ -z "$TASK" ]]; then
    echo "No TASK variable found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}."
    exit 1
fi

echo "Running task ${SLURM_ARRAY_TASK_ID} with ID: $TASK"

# Submit temporary files job
echo "Starting Temporary Files generation for ID: $TASK"
job_id=$(sbatch --parsable --export=ID=$TASK,mode=$mode \
    --output="LogsMA/${TASK}/TemporaryFilesOutput/%j.out" \
    --error="LogsMA/${TASK}/TemporaryFilesOutput/%j.out" \
    MA/TemporaryFiles.sh)
echo "Submitted Temporary Files job with Job ID: $job_id"

# Submit metal job, held until temporary files job completes
echo "Starting Metal Analysis for ID: $TASK"
metal_id=$(sbatch --parsable --export=ID=$TASK,mode=$mode \
    --output="LogsMA/${TASK}/MetalOutput/%j.out" \
    --error="LogsMA/${TASK}/MetalOutput/%j.out" \
    --dependency=afterok:$job_id \
    MA/MetalAnalysis.sh)
echo "Submitted Metal Analysis job with Job ID: $metal_id"

# Submit cleanup job, held until metal job completes
echo "Scheduling Cleanup for ID: $TASK"
sbatch --parsable --export=ID=$TASK,mode=$mode \
    --output="LogsMA/${TASK}/CleanUpOutput/%j.out" \
    --error="LogsMA/${TASK}/CleanUpOutput/%j.out" \
    --dependency=afterok:$metal_id \
    MA/CleanUp.sh
echo "Submitted Cleanup job for ID: $TASK"