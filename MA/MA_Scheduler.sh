#!/bin/bash
#SBATCH --job-name=MA_Scheduler
#SBATCH --output=LogsMA/Scheduler_%A.out
#SBATCH --error=LogsMA/Scheduler_%A.out
#SBATCH --time=10-00:00:00
#SBATCH --mem=4G

mode="LOO"
numberIDs=7578
batchSize=150

echo "Running MA for mode: $mode"
cd $SLURM_SUBMIT_DIR

wait_for_queue() {
    echo "Waiting for jobs to finish..."
    waited=0
    while [ $(squeue -u $USER | grep -v $SLURM_JOB_ID | grep -c " ") -gt 1 ]; do
        sleep 300
        waited=$((waited + 5))
    done
    echo "All jobs finished after waited for $waited minutes."
    echo "Queue empty, continuing..."
}

for start in $(seq 7351 $batchSize $numberIDs); do
    end=$((start + batchSize - 1))
    # cap at numberIDs for the final batch
    end=$((end > numberIDs ? numberIDs : end))
    
    echo "Submitting array job for IDs: $start-$end"
    sbatch --export=mode=$mode --array=${start}-${end} --output=LogsMA/Master_${start}-${end}/%A_%a.out --error=LogsMA/Master_${start}-${end}/%A_%a.out MA/Master.sh
    
    wait_for_queue
done

echo "All jobs complete"
