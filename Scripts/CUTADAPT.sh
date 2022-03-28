#!/bin/bash

# Define list of samples
SAMPLES='ERR1713356 ERR2592255 ERR2683233 ERR4682862'

# Activate the conda environment
conda activate QC_env

# Create output folder
mkdir TRIMMED

# Run Cutadapt
for SAMPLE in $SAMPLES; do 
  cutadapt ~/Share/RAWDATA/$SAMPLE.R1.fastq.gz \
           ~/Share/RAWDATA/$SAMPLE.R2.fastq.gz \
           -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
           -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
           -o TRIMMED/$SAMPLE.R1.fastq.gz \
           -p TRIMMED/$SAMPLE.R2.fastq.gz \
           --minimum-length 50 \
           --quality-cutoff 20 \
           --cores 4 > TRIMMED/$SAMPLE.cutadapt.log.txt
done
