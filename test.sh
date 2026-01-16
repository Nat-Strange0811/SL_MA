#!/bin/bash
#$ -cwd
#$ -N Test
#$ -pe smp 1
#$ -l h_vmem=64G
#$ -l h_rt=1:00:00
#$ -j y

module load miniforge
mamba activate SL_MA_QC

python test.py 