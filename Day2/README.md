# Day 2

| Time      | Activity            | Slides                               | Hands-on                          |
|-----------|---------------------|--------------------------------------|-----------------------------------|
| Morning   | Read-based analyses | [Link here](read-based-analyses.pdf) | [Link here](#read-based-analyses) |
| Afternoon | Metagenome assembly | [Link here](metagenome-assembly.pdf) | [Link here](#metagenome-assembly) |

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
