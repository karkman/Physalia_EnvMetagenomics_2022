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

In the next step we just combine the individual profiles (one per sample) and do some clustering based on the detection and frequency of the contigs in each sample:

```bash
anvi-merge ANVIO/$SAMPLE/PROFILE/*/PROFILE.db \
           --output-dir ANVIO/$SAMPLE/MERGED \
           --contigs-db ANVIO/$SAMPLE/CONTIGS.db
```

### Tunneling the interactive interface

Although you can install anvi'o on your own computer (and you're free to do so, but we won't have time to help in that, but take a look [here](https://anvio.org/install/)), we will run anvi'o in the cloud and tunnel the interactive interface to your local computer.  
To be able to to do this, everyone needs to log into the cloud using a different port; your port will be __8080 + you user number__.  
So `user1` will use port 8081, `user2` will use port 8082 and so on, `user10` will use port 8090, `user20` will use port 8100, and so on.  
If you're

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
anvi-interactive --profile-db ANVIO/$SAMPLE/MERGED/PROFILE.db \
                 --contigs-db ANVIO/$SAMPLE/CONTIGS.db \
                 --port-number XXXX \
                 --server-only
```

When you are done with the first round of binning, you can summarise your bins:  
**NOTE:** If you have saved your collection with a different name, you have to change `--collection-name` below accordingly.

```bash
anvi-summarize --contigs-db ANVIO/$SAMPLE/CONTIGS.db \
               --pan-or-profile-db ANVIO/$SAMPLE/MERGED/PROFILE.db \
               --output-dir ANVIO/"$SAMPLE"_SUMMARY \
               --collection-name default \
               --quick-summary

```

Download the four `$SAMPLE_SUMMARY` folder to your computer and open the `.html` file in your browser.  
Do you see some bin(s) that you will like to refine, or just check?  
You can do so like this:

```bash
anvi-refine --profile-db ANVIO/$SAMPLE/MERGED/PROFILE.db \
            --contigs-db ANVIO/$SAMPLE/CONTIGS.db \
            --collection-name default \
            --bin-id $BIN \
            --port-number XXXX \
            --server-only
```

When you are happy with binning and refining, it's time to rename things a bit to make things easier downstream.  
Here we will decide that all bins that are ≥50% complete and ≤10% redundant ("contaminated") are good enough.  
You should use, however, the thresholds that you feel are more aligned with your research question.  
Below we will create a new collection called **final**, and those bins that are above the threshold we selected will be renamed as MAGs.  
If you look at the command below, you will see that the MAG names will have a prefix indicating the user who binned them and which sample it came from.

```bash
anvi-rename-bins --profile-db ANVIO/$SAMPLE/MERGED/PROFILE.db \
                 --contigs-db ANVIO/$SAMPLE/CONTIGS.db \
                 --collection-to-read default \
                 --collection-to-write final \
                 --prefix "$USER"_"$SAMPLE" \
                 --report-file ANVIO/$SAMPLE/renamed_MAGs.txt \
                 --call-MAGs \
                 --min-completion-for-MAG 50 \
                 --max-redundancy-for-MAG 10
```

And now let's summarise everything again:

```bash
anvi-summarize --contigs-db ANVIO/$SAMPLE/CONTIGS.db \
               --pan-or-profile-db ANVIO/$SAMPLE/MERGED/PROFILE.db \
               --output-dir ANVIO/"$SAMPLE"_SUMMARY_FINAL \
               --collection-name final
```

Of course, you can go back many times as you wish and redo everything from scratch.  
You can try different ways of binning and refining, basically ad infinitum.  
But no matter how much you try to avoid, at one point you have to get over binning the MAGs and actually start writing that paper...  
For now, let's take a well deserved rest :)

![](https://i.pinimg.com/originals/61/0e/2b/610e2b80ac1632e0e7a0f47f5effc08d.jpg)