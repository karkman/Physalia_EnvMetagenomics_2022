#!/bin/bash

# Define list of samples
SAMPLES='ERR1713356 ERR2592255 ERR2683233 ERR4682862'

# Activate the conda environment
source activate read_based_env

# Create output folder
mkdir MEGAN

# Resample with seqtk
for SAMPLE in $SAMPLES; do
  seqtk sample -s100 TRIMMED/$SAMPLE.R1.fastq.gz 500000 | gzip > MEGAN/$SAMPLE.R1.fastq.gz
  seqtk sample -s100 TRIMMED/$SAMPLE.R2.fastq.gz 500000 | gzip > MEGAN/$SAMPLE.R2.fastq.gz
done

# Blastx with DIAMOND
for SAMPLE in $SAMPLES; do 
  diamond blastx --query MEGAN/$SAMPLE.R1.fastq.gz \
                 --out MEGAN/$SAMPLE.R1.blastx.txt \
                 --db nr \
                 --outfmt 0 \
                 --threads 4

  diamond blastx --query MEGAN/$SAMPLE.R2.fastq.gz \
                 --out MEGAN/$SAMPLE.R2.blastx.txt \
                 --db nr \
                 --outfmt 0 \
                 --threads 4

  cat MEGAN/$SAMPLE.R1.blastx.txt MEGAN/$SAMPLE.R2.blastx.txt > MEGAN/$SAMPLE.blastx.txt
done

# Run MEGAN
for SAMPLE in $SAMPLES; do 
  megan/tools/blast2rma --in MEGAN/$SAMPLE.blastx.txt \
                        --out MEGAN/$SAMPLE.rma6 \
                        --mapDB megan-map-Feb2022.db \
                        --format BlastText \
                        --threads 4
done
