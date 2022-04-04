# Day 2

| Time      | Activity            | Slides                               | Hands-on                          |
|-----------|---------------------|--------------------------------------|-----------------------------------|
| Morning   | Read-based analyses | [Link here](read-based-analyses.pdf) | [Link here](#read-based-analyses) |
| Afternoon | Metagenome assembly | [Link here](assembly-and-qc.pdf) | [Link here](#metagenome-assembly) |
| Afternoon | Assembly QC         |                                      | [Link here](#assembly-QC)         |

## Read-based analyses

### Running the script
First login to the Amazon Cloud and `cd` to your working directory.  
We migh have made changes to the GitHub repository, so let's pull those changes now:

```bash
git pull
```

For the read-based analyses, we will use [seqtk](https://github.com/lh3/seqtk), [DIAMOND](https://github.com/bbuchfink/diamond), and [MEGAN](https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/megan6/).  
Like yesterday, the script is provided and can be found from the `Scripts` folder.  
Let's take a look at the script using `less`:

```bash
less Scripts/MEGAN.sh
```

**NOTE:** You can scroll up and down using the arrow keys on your keyboard, or move one "page" at a time using the spacebar.  
**NOTE:** To quit `less`, hit the letter **q**.  

As you will see, we are first running `seqtk` to subsample the data to an even number of sequences per sample (500,000).  
Then we are running `DIAMOND` to annotate the reads against the NCBI nr database.  
Then we will use `MEGAN` to parse the annotations and get taxonomic and functional assignments.  

Now let's run the script:

```bash
bash Scripts/MEGAN.sh
```

Both the `DIAMOND` and `MEGAN` steps will actually take a while and it's very likely that the jobs won't finish in a reasonable time.
For the purposes of this activity, it is enough if 1) we understand what the script is doing and 2) we are able to submit the script without any errors.
So now let's stop the script by hitting `ctrl+c`.  

Luckily, we have a copy of the `MEGAN` results in the shared folder.
Let's take a look at the `MEGAN` folder inside `~/Share`, which contains the actual output from the script we tried to run before:  

```bash
ls ~/Share/MEGAN
```

For each sample, you should find:
- `$SAMPLE.R1.fastq.gz` and `$SAMPLE.R2.fastq.gz`: the subsampled data obtained with `seqtk`
- `$SAMPLE.R1.blastx.txt` and `$SAMPLE.R2.blastx.txt`: the output from `DIAMOND` for each of the forward and reverse reads
- `$SAMPLE.blastx.txt`: a file combining both `DIAMOND` outputs
- `$SAMPLE.rma6`: `MEGAN` output

The `.rma6` files are compressed binary files created by `MEGAN` (command-line version).  
These describe the taxonomic and functional composition of the samples based on the `DIAMOND` annotation against the NCBI nr database.  

`MEGAN` also has a powerful GUI version, that you have installed in your own computer.  
First let's copy the four `.rma6` files to your own computers using FileZilla.  
When that's done let's launch `MEGAN` and take a look together at one of the samples.  

Now, by using the `Compare` tool, let's try to find differences between the samples.  
On the slides for the first day ("Course outline and practical info") we saw that we have two winter and two summer samples.  
Can we see major differences in community structure between the two seasons? For example:
- Are samples from the same season more similar to each other than to the other samples? **HINT:** Try changing the view to the Genus Rank and then going to "Window" > "Cluster Analysis" and chosing "UPGMA Tree".
- What is the most abundant phylum in the summer samples? And in the winter samples?
- Now looking at the functional profiles (e.g. SEED), can you spot differences between the two seasons? Specially regarding energy and metabolism?

## Metagenome assembly

Now we will go through the metagenomic assembly part, but not run the actual assembly script.  
The assembly takes some time and needs more resources than we have on our instance.  
So the assemblies will be provided.  

### Short-read assembly with MEGAHIT
The short reads will be assembled using [MEGAHIT](https://github.com/voutcn/megahit).  
Although, we won't be running the actual assembly, `MEGAHIT` is installed on our instance.  

So have a look at the different options you can change in the assembly.  
You can find more information about `MEGAHIT` in their [wiki](https://github.com/voutcn/megahit/wiki).  
You don't need to understand each and every option, but some of them can be important.

```bash
conda activate megahit_env
megahit -h
```

#### Questions about MEGAHIT
1. __What do you think would be important? What would you change or set?__  
2. __What version of MEGAHIT have we installed? Is it the latest?__

After that have a look at the assembly script `Scripts/MEGAHIT.sh`.  
Open it with a text editor or print it on the screen with `less`.  

__Would you have changed something else and why?__

When we're satisfied with the assembly options, we would start the assembly and wait from few hours to several days depending on your data and computational resources.  
But we won't do it, since we don't have to time or the resources.  
Instead, you can use the assemblies we provide in the shared foler:

```bash
ls ~/Share/ASSEMBLY
```

If you look inside the folders for each of the samples, you will see several files.  
But the most important is the `final.contigs.fa` which cointains the final contigs as one might expect.
You'll find the assembly logs inside the assembly folder for each sample.  
Start by looking at the assembly logs with `less`.

#### Questions about the assembly
1. __Which version of megahit did we actually use for the assemblies?__
2. __How long did the assemblies take to finish?__
3. __Which sample gave the longest contig?__

### Long-read assembly

Short read assemblers are not optimal for long-read data, but there are already several long-read assemblers available.  
Some of the most used metagenomic assemblers include [metaflye](https://github.com/fenderglass/Flye), [hicanu](https://github.com/marbl/canu) and [hifiasm-meta](https://github.com/xfengnefx/hifiasm-meta).  

For this course will use a data set of 434 105 PacBio HiFi-reads from one sample with mean length of ~7 500 bp (~3 Gbp).  
And we have used metaflye for the assembly with the following parameters:

```
flye \
    --meta \
    --threads 12 \
    --pacbio-hifi  \
    --min-overlap 4000 \
    --out-dir INF3_assembly \
    INF3.fastq.gz
```

The final assembly can be found from `PATH/TO/ASSEMBLY/assembly.fasta`.

#### Questions about long-read assembly
1. __What the min-overlap option controls and why we might have used it?__
2. __How would you change the options if you had nanopore data to assemble?__

## Assembly QC

Now we have all the assemblies ready and we can use [MetaQUAST](http://bioinf.spbau.ru/metaquast) to check how the assemblies look like.  
Let's activate the `MetaQuast` environment:

```bash
conda activate quast_env
```

Let's now have a look at the different options `MetaQuast` has with `metaquast -h`.  
You should at least check the options we are using.  
We will run `MetaQuast` inside a screen using the command `screen`.  
This way you can do other things or log out while `MetaQuast` is running and it won't be interrupted.

Mini manual for `screen`:
* `screen -S NAME` - open a screen and give it a session name `NAME`
* `screen` - open new screen without specifying any name
* `screen -ls` - list all open sessions
* `ctrl` + `a` + `d` - to detach from a session (from inside the screen)
* `screen -r` - re-attach to a detached session
* `screen -rD` - re-attach to a attached session
* `exit` - close the screen and kill all processes running inside the screen (from inside the screen)

```bash
screen -S metaquast

metaquast.py ~/Share/ASSEMBLY/*/final.contigs.fa ~/Share/ASSEMBLY_HIFI/assembly.fasta \
             --output-dir ASSEMBLY_QC \
             --threads 4 \
             --max-ref-number 0 \
             --min-contig 0
```
Detach from the screen with `ctrl` + `a` + `d`.  
This will take ~5 min.  
You can re-attach with `screen -r metaquast` to check whther it has finished.  
After it is done, we will go through the report together.  
Download the file `ASSEMBLY_QC/report.html` to your computer using FileZilla and open it on your favourite browser.

#### Questions about the assembly QC

1. __Which assembly has the longest contig when also long reads assembly is included?__
2. __Which assembly had the most contigs?__
3. __Was the long read assembly different from the short read assemblies?__
4. __If yes, in what way?__
