# Day 3

| Time                | Activity                                 | Slides                                        | Hands-on                                   |
|---------------------|------------------------------------------|-----------------------------------------------|--------------------------------------------|
| Morning & afternoon | Genome-resolved metagenomics with anvi'o | [Link here](genome-resolved-metagenomics.pdf) | [Link here](#genome-resolved-metagenomics) |

## Genome-resolved metagenomics with anvi'o

**TO DO:** Add short intro about anvi'o

### Preparing the data
First login to the Amazon Cloud, `cd` to your working directory and pull the changes from Github.

```bash
cd ~/Physalia_EnvMetagenomics_2022
git pull
```

And let's create a directory to store the binning files and activate the conda environment:

```bash
mkdir ANVIO
conda activate anvio7_env
```

Because we have four assemblies, we will run each step below essentially four times.  
One way to do it is manually: you run the same command four times, each time changing by hand the name of the input and output files to match one sample at a time.  
A better way to this - especilly if you have tens, hundreds, thousands of samples - is by automating this process.  
One way of automating a script is by creating a `for loop`, as we have done e.g. when we assembled our metagenomes.  
So take a look again at the `MEGAHIT` script from yesterday to see how the syntax of the `for loop` works, and let's create loops to run the commands below to prepare the files for anvi'o.  

To make it things more organized, let's create a directory for each sample inside the `ANVIO` directory:

```bash
mkdir $SAMPLE
```

The first thing we will do is reformat the `FASTA` files so that the headers are clean.  
In the command below we will also remove contigs shorter than 1,000 bp, to make things easier during binning.  

```bash
anvi-script-reformat-fasta ASSEMBLY/$SAMPLE/final.contigs.fa \
                           --output-file ANVIO/$SAMPLE/CONTIGS_2500nt.fa \
                           --report-file ANVIO/$SAMPLE/CONTIGS_reformat.txt \
                           --prefix $SAMPLE \
                           --min-len 1000 \
                           --simplify-names
```

Next thing after some cleaning is to build a contigs database from our assembly.  
In this step anvi'o will call genes from each contig using prodigal, calculate kmer frequencies for each contig with k=4 and "soft" split each contig longer than 20 kbp for visual reasons.

```bash
anvi-gen-contigs-database --contigs-fasta ANVIO/$SAMPLE/CONTIGS_2500nt.fa \
                          --output-db-path ANVIO/$SAMPLE/CONTIGS.db \
                          --project-name $SAMPLE \
                          --num-threads 4
```

In this step anvi'o will identify different sets of bacterial, archaeal and eukaryotic marker genes with HMM models and use them to predict the completeness and contamination (redundancy) of bins.

```bash
anvi-run-hmms --contigs-db ANVIO/$SAMPLE/CONTIGS.db \
              --num-threads 4
```

This step is related to the previous one, here anvi'o looks for a set of single-copy core genes that it will use to predict the taxonomy of individual bins.  

```bash
anvi-run-scg-taxonomy --contigs-db ANVIO/$SAMPLE/CONTIGS.db \
                      --num-threads 4
```

After we have populated to contigs database with various data about the contigs, it is time to use the information in the raw sequencing reads.  
We will map (basically align) each read to the assembly and then store the information we get from this process to separate file (`.bam`) for each sample.  
The most important information being where each read aligns and how well.  

Pay attention that below there will be an inception thing going on: a `for loop` inside another `for loop`!

```bash
mkdir ANVIO/$SAMPLE/MAPPING

bowtie2-build ANVIO/$SAMPLE/CONTIGS_2500nt.fa \
              ANVIO/$SAMPLE/MAPPING/contigs

for FILE in $SAMPLES; do
  bowtie2 -1 TRIMMED/$FILE.R1.fastq.gz \
          -2 TRIMMED/$FILE.R2.fastq.gz \
          -S ANVIO/$SAMPLE/MAPPING/$FILE.sam \
          -x ANVIO/$SAMPLE/MAPPING/contigs \
          --threads 4 \
          --no-unal

  samtools view -F 4 -bS ANVIO/$SAMPLE/MAPPING/$FILE.sam | samtools sort > ANVIO/$SAMPLE/MAPPING/$FILE.bam
  samtools index ANVIO/$SAMPLE/MAPPING/$FILE.bam

  rm -f ANVIO/$SAMPLE/MAPPING/$FILE.sam
done
```

From the information stored in the bam files, we (read: __anvi'o__) can calculate the coverage for each base in each contig and also the sequence variation (SNPs) in the community for each nucleotide in each contig.  

```bash
mkdir ANVIO/$SAMPLE/PROFILE

for FILE in $SAMPLES; do
  anvi-profile --input-file ANVIO/$SAMPLE/MAPPING/$FILE.bam \
               --output-dir ANVIO/$SAMPLE/PROFILE/$FILE \
               --contigs-db ANVIO/$SAMPLE/CONTIGS.db \
               --num-threads 4
done
```

**NOTE:** `for loop` ends here.  

In the next step we just combine the individual profiles (one per sample) and do some clustering based on the detection and frequency of the contigs in each sample.  

```bash
anvi-merge ANVIO/$SAMPLE/PROFILE/*/PROFILE.db \
           --output-dir ANVIO/$SAMPLE/MERGED \
           --contigs-db ANVIO/$SAMPLE/CONTIGS.db
```

### Tunneling the interactive interface

Although you can install anvi'o on your own computer (and you're free to do so, but we won't have time to help in that), we will run anvi'o in the cloud and tunnel the interactive interface to your local computer.  
To be able to to do this, everyone needs to use a different port for tunneling and your port will be __8080 + you user number__. So `user1` will use port 8081.

#### Linux/Mac

```bash
ssh -i KEY.pem -L PORT:localhost:PORT USERX@IP-ADDRESS
```

#### Windows (MobaXterm)

<!-- Add instructions here -->

### Binning MAGs
