#!/bin/bash
#SBATCH --job-name=Sentinel-Variant-Analysis
#SBATCH --output=LogsSV/%A_%a.out
#SBATCH --error=LogsSV/%A_%a.out
#SBATCH --time=1:00:00
#SBATCH --mem=24G
#SBATCH --array=1-7578

cd $SLURM_SUBMIT_DIR

module load miniforge
mamba activate SL_MA_QC

MODE="STANDARD"
PROTEIN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" Proteins.txt)

python Sentinel-Variant-Analysis/main.py $PROTEIN $MODE

if [ "$SLURM_ARRAY_TASK_ID" -eq 1 ]; then
    sbatch --dependency=afterany:$SLURM_ARRAY_JOB_ID --wrap="awk -F',' 'NR==1{header=\$0; print; next} {print; for(i=1;i<=NF;i++) sum[i]+=\$i; count++} END{for(i=1;i<=NF;i++) printf (i==1?\"Average\":sum[i]/count) (i==NF?\"\n\":\",\")}' Sentinel-Variant-Analysis/Results/Comparison/comparison_${MODE}.csv > TemporaryFiles/tmp_averaged.csv && mv TemporaryFiles/tmp_averaged.csv Sentinel-Variant-Analysis/Results/Comparison/comparison_${MODE}.csv"
fi