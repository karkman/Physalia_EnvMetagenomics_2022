# Day 3

| Time                | Activity                                 | Slides                                        | Hands-on                                   |
|---------------------|------------------------------------------|-----------------------------------------------|--------------------------------------------|
| Morning & afternoon | Genome-resolved metagenomics with anvi'o | [Link here](genome-resolved-metagenomics.pdf) | [Link here](#genome-resolved-metagenomics) |

## Genome-resolved metagenomics

For binning contigs into metagenome-assembled genomes (MAGs), we will use `anvi'o`.  
You should definitely take a look at their [website](https://anvio.org) and maybe even join their [slack channnel](https://join.slack.com/t/anvio/shared_invite/zt-ov46uj90-9p2woLJFcVCfv7cdhANXSA).  

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
**TIP:** Open a text editor (either in your own computer or directly in the cloud) to organize the pipeline.

To make it things more organized, let's create a directory for each sample inside the `ANVIO` directory:

```bash
mkdir ANVIO/$SAMPLE
```

The first thing we will do is reformat the `FASTA` files so that the headers are clean.  
In the command below we will also remove contigs shorter than 1,000 bp, to make things easier during binning.  

```bash
anvi-script-reformat-fasta ~/Share/ASSEMBLY/$SAMPLE/final.contigs.fa \
                           --output-file ANVIO/$SAMPLE/CONTIGS_2500nt.fa \
                           --report-file ANVIO/$SAMPLE/CONTIGS_reformat.txt \
                           --prefix $SAMPLE \
                           --min-len 1000 \
                           --simplify-names
```

Next thing after some cleaning is to build a **contigs database** from our assembly.  
In this step anvi'o will call genes from each contig using prodigal, calculate kmer frequencies for each contig with k=4 and "soft" split each contig longer than 20 kbp for visual reasons.

```bash
anvi-gen-contigs-database --contigs-fasta ANVIO/$SAMPLE/CONTIGS_2500nt.fa \
                          --output-db-path ANVIO/$SAMPLE/CONTIGS.db \
                          --project-name $SAMPLE \
                          --num-threads 4
```

In this step anvi'o will identify different sets of bacterial, archaeal and eukaryotic single copy genes with HMM models and use them to predict the completeness and contamination (redundancy) of bins, among other things.

```bash
anvi-run-hmms --contigs-db ANVIO/$SAMPLE/CONTIGS.db \
              --num-threads 4
```

This step is related to the previous one, but here anvi'o looks at single-copy genes and will try to predict the taxonomy of the bins.  

```bash
anvi-run-scg-taxonomy --contigs-db ANVIO/$SAMPLE/CONTIGS.db \
                      --num-threads 4
```

After we have populated to contigs database with various data about the contigs, it is time to gather information on the abundance of the contigs using the raw sequencing reads.  
We will map (basically align) the reads to the contigs and then store the information we get from this process to a separate file (`.bam`) for each sample.  
The most important information being where each read aligns and how well.  

Pay attention that below there will be an inception thing going on: a `for loop` inside another `for loop`!  
This is because, for each assembly, we are mapping the raw sequences from all the four samples, not only the sample that was used to make the assembly.
Not sure why? Don't despair: this will become more clear when we are actually binning the MAGs.

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
This is done in the next step, which will store all this info into another object called a **profile database**:

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

In the next step we just combine the individual profiles (one per sample) and do some clustering based on the detection and frequency of the contigs in each sample:

```bash
anvi-merge ANVIO/$SAMPLE/PROFILE/*/PROFILE.db \
           --output-dir ANVIO/$SAMPLE/MERGED \
           --contigs-db ANVIO/$SAMPLE/CONTIGS.db
```

### Tunneling the interactive interface

Although you can install anvi'o on your own computer (and you're free to do so, but we won't have time to help in that, but take a look [here](https://anvio.org/install/)), we will run anvi'o in the cloud and tunnel the interactive interface to your local computer.  
To be able to to do this, everyone needs to use a different port for tunneling and your port will be __8080 + you user number__.  
So `user1` will use port 8081, `user2` will use port 8082 and so on, `user10` will use port 8090, `user20` will use port 8100, and so on.

#### Linux/Mac

```bash
ssh -i KEY.pem -L PORT:localhost:PORT USERX@IP-ADDRESS
```

#### Windows (MobaXterm)

Let's take a look together with Carlo to see how to set this up.

### Binning MAGs

Now we are ready for some serious MAG binning!  
Let's take a look together on how to use `anvi'o` to bin MAGs, and then take the rest of the day to bin the four samples yourself.  
Remember to change `--port-number XXXX` to match the port number that you used to log into the cloud.

```bash
SAMPLES='ERR1713356 ERR2592255 ERR2683233 ERR4682862'

for SAMPLE in $SAMPLES; do
  anvi-interactive --profile-db ANVIO/$SAMPLE/MERGED/PROFILE.db \
                   --contigs-db ANVIO/$SAMPLE/CONTIGS.db \
                   --port-number XXXX \
                   --server-only
done
```

REFINING  

CALLING MAGS  

SUMMARISING  