#!/bin/bash
#SBATCH --job-name=Metal_Analysis
#SBATCH --cpus-per-task=1
#SBATCH --mem=14G
#SBATCH --time=1:00:00

# Set initial directories
directory="TemporaryFiles/$ID"
metal_exe="/data/PHURI-Langenberg/programs/METAL/random-metal-0.1.0/executables/metal"
outfile="MA/Results/${ID}/MetaAnalysis_${ID}_${mode}_"

# Create temporary command file
cmd_file=$(mktemp)

# Make output directory
mkdir -p "MA/Results"
mkdir -p "MA/Results/${ID}"

# Write initial METAL commands
cat <<EOT > $cmd_file
SCHEME STDERR
AVERAGEFREQ ON
MINMAXFREQ ON
CUSTOMVARIABLE TotalSampleSize
LABEL TotalSampleSize as N

MARKER SNPID
ALLELE Effect_Allele Other_Allele
EFFECT Beta
PVALUE Pval
WEIGHT N
FREQ EAF
STDERR SE
SEPARATOR TAB
EOT

# Loop through temporary files and append PROCESS commands
for file in $directory/*.txt.gz; do
    echo $file
    echo "PROCESS $file" >> "$cmd_file"
done

# Append final METAL commands
cat <<EOT >> $cmd_file
OUTFILE $outfile .tbl
ANALYZE HETEROGENEITY
QUIT
EOT

# Run METAL
$metal_exe $cmd_file

# Remove temporary command file
rm $cmd_file

# Compress output
gzip -f "MA/Results/${ID}/MetaAnalysis_${ID}_${mode}_1.tbl"