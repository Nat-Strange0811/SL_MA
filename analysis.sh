#!/bin/bash
#$ -cwd
#$ -N QC_Analysis
#$ -pe smp 1
#$ -l h_vmem=8G
#$ -l h_rt=1:0:00
#$ -j y

module load miniforge
mamba activate SL_MA_QC

python analysis.py

echo "QC Analysis completed"