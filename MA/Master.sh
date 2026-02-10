#!/bin/bash
#$ -cwd
#$ -N Master
#$ -pe smp 1
#$ -l h_vmem=20G
#$ -l h_rt=1:00:00
#$ -o LogsMA/
#$ -j y
#$ -t 1-300

# Load the list of IDs from the file and select the one corresponding to the current task ID
# This process is paralised by the SGE_TASK_ID variable
TASK=$(sed -n "${SGE_TASK_ID}p" MA/IDs.txt)

#Ensure that the TASK variable is not empty
if [[ -z "$TASK" ]]; then
    echo "No TASK variable found for SGE_TASK_ID=${SGE_TASK_ID}."
    exit 1
fi

#Reporting call
echo "Running task ${SGE_TASK_ID} with ID: $TASK"

#This call is used to run the scheduler script in python, this enables a reduction in the number
#of 'jobs' called as the script handles multi-core processing.
#python MA/scheduler.py $TASK

#Activate the temporary files script, generating temporary files for the given ID, this is a 
#prerequisite for the metal analysis.
echo "Starting Temporary Files generation for ID: $TASK"
#Capture the job ID to ensure that METAL does not start until temporary files are generated.
job_id=$(qsub -terse -v ID="$TASK" -o "LogsMA/${TASK}/TemporaryFilesOutput/" MA/TemporaryFiles.sh)
#Trim the job ID so in proper format for the hold_jid parameter in the subsequent qsub calls.
job_id=${job_id%%.*}
echo "Submitted Temporary Files job with Job ID: $job_id"

echo

#Activate the metal analysis script, will not start until temporary files job has completed
echo "Starting Metal Analysis for ID: $TASK"
#Hold the jobID for the cleanup task, ensuring files remain until the process has finished.
metal_id=$(qsub -terse -v ID="$TASK" -o "LogsMA/${TASK}/MetalOutput/" -hold_jid "$job_id" MA/MetalAnalysis.sh)
#Cleanup the job ID for the same reason as above.
metal_id=${metal_id%%.*}
echo "Submitted Metal Analysis job with Job ID: $metal_id"

echo

#Call the clean up script.
echo "Scheduling Cleanup for ID: $TASK"
qsub -terse -v ID="$TASK" -o "LogsMA/${TASK}/CleanUpOutput/" -hold_jid "$metal_id" MA/CleanUp.sh
echo "Submitted Cleanup job for ID: $TASK"