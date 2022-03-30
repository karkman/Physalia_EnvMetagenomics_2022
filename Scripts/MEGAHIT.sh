#!/bin/bash

# Define list of samples
SAMPLES='ERR1713356 ERR2592255 ERR2683233 ERR4682862'

# Activate the conda environment
conda activate megahit_env

# Create output folder
mkdir ASSEMBLY

# Run MEGAHIT
for SAMPLE in $SAMPLES; do
  megahit -1 TRIMMED/$SAMPLE.R1.fastq.gz \
          -2 TRIMMED/$SAMPLE.R2.fastq.gz \
          --out-dir ASSEMBLY/$SAMPLE \
          --k-min 27 \
          --k-max 127 \
          --k-step 10 \
          --num-cpu-threads 4
done
