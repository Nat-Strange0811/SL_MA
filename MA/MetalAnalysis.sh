#!/bin/bash
#$ -cwd
#$ -N Metal_Analysis
#$ -pe smp 1
#$ -l h_vmem=20G
#$ -l h_rt=1:00:00
#$ -l highmem
#$ -j y

#Set initial directories, where the files, metal executable and output directories are located.
directory="TemporaryFiles/$ID"
metal_exe="/data/PHURI-Langenberg/programs/METAL/random-metal-0.1.0/executables/metal"
outfile="MA/Results/${ID}/MetaAnalysis_${ID}_"

#Creates a temporary file to store the commands for the metal executable, is dynamically written to
cmd_file=$(mktemp)

#Make the output directory for the given ID
mkdir -p "MA/Results/${ID}"

#Write to the temporary file created above, refer to METAL documentation for detail on commands used
cat <<EOT > $cmd_file
SCHEME STDERR
MARKER SNPID
ALLELE Effect_Allele Other_Allele
EFFECT Beta
PVALUE Pval
WEIGHT N
STDERR SE
SEPARATOR TAB
EOT

#Loop through all temporary files and include them as 'PROCESS' files for the metal executable,
#We use >> instead of > to append to the file ensuring that the initial commands are not overwritten.
for file in $directory/*.txt.gz; do
    echo $file
    echo "PROCESS $file" >> "$cmd_file"
done

#Include the final commands calling for analysis and specifiying the output file.
cat <<EOT >> $cmd_file
OUTFILE $outfile .tbl
ANALYZE HETEROGENEITY
QUIT
EOT

#Activate metal passing the temporary command file as an argument
$metal_exe $cmd_file

#Delete the temporary command file as it is no longer needed
rm $cmd_file

#gzip the output file to save space
gzip -f "MA/Results/${ID}/MetaAnalysis_${ID}_1.tbl"