#!/bin/bash
#$ -cwd
#$ -N Combine_Viking
#$ -pe smp 1
#$ -l h_vmem=8G
#$ -l h_rt=1:0:00
#$ -j y

module load miniforge
mamba activate SL_MA_QC

python combine.py