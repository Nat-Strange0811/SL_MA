#!/bin/bash
#$ -cwd
#$ -N Clean_Up
#$ -pe smp 1
#$ -l h_vmem=8G
#$ -l h_rt=1:00:00
#$ -j y

#Activate the conda environment
module load miniforge
mamba activate SL_MA_QC

#Ensure that ID variable is provided
if [[ -z "$ID" ]]; then
    echo "No ID variable provided."
    exit 1
fi

#Activate clean up script
python MA/delete.py $ID