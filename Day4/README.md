# Day 4

| Time      | Activity                                 | Slides                                                | Hands-on                                                    |
|-----------|------------------------------------------|-------------------------------------------------------|-------------------------------------------------------------|
| Morning   | Genome-resolved metagenomics with anvi'o | [Link here](../Day3/genome-resolved-metagenomics.pdf) | [Link here](../Day3/README.md#genome-resolved-metagenomics) |
| Afternoon | Downstream analyses                      |                                                       | [Link here](#downstream-analyses)                           |

## Downstream analyses

When you're happy with the bins you got and have called MAGs in the last `anvi-rename-bins`and `anvi-summarize` commands, copy all the MAGs you got into one folder.  
This might work for you:

```bash
cd ~/Physalia_EnvMetagenomics_2022
mkdir MAG_folder
cp ANVIO/*/bin_by_bin/*MAG*/*-contigs.fa MAG_folder/
```

After it's done, list the content of the folder to see how many MAGs you have in total.

### CheckM

Although anvi'o gives us estimates about the completeness and redundancy, we'll double check those by using [checkM](https://ecogenomics.github.io/CheckM/).  
Make sure you change the input folder and the output folder in the command below.  
Running checkM took about 10 min for 16 MAGs.

```bash
conda activate MAG_annotation_env

checkm lineage_wf -x fa \
                  MAG_folder \
                  checkM_out \
                  -t 4 \
                  --tmpdir /tmp
```

The output table is printed on to the screen. If you missed it, you can quickly reproduce it with:

```
checkm qa checkM_out/lineage.ms checkM_out
```

__Were the checkM results in line with the estimates from anvi'o?__

### GTDB-Tk

Next step would be to run [GTDB-Tk](https://github.com/Ecogenomics/GtdbTk) to get the phylogenetic assignment of our MAGs.  
However, this step takes too much memory, so we won't run it.

```
## DO NOT RUN ###
gtdbtk classify_wf -x fa \
                   --genome_dir MAG_folder \
                   --out_dir GTDB \
                   --cpus 4 \
                   --tmpdir /tmp
```

### dRep

Often we retrieve MAGs that are very closely related to each other, specially if you are binning several samples at the same time.  
In this case, you might be interested in removing the redundancy (i.e. deplicating) and keep only one MAG from each cluster of closely related MAGs.  
For this, we can use a tool called [dRep](https://drep.readthedocs.io/en/latest/), but guess what, you can do that in `anvio` as well...

```bash
source activate drep_env

dRep compare MAGs_drep \
             --genomes MAG_folder/*.fa \
             --processors 4
```


### Antibiotic resistance gene annotation

One of our research questions for this course was related to antibiotic resistance genes (ARGs). We wanted to study the genetic contexts of different ARGs and see whether  they would be found from various context. We will use a program called [abricate](https://github.com/tseemann/abricate) and the antibiotic resistance gene database called [ResFinder](https://doi.org/10.1093/jac/dkaa345). And we will annotate the filtered contig files that were used to construct the anvi'o  contigs database.  
After activating the annotation environment (`conda activate annotation_env`), you can check all possible options for abricate with `abricate -h`. Abricate has also other pre-installed databases, try to list them on the command line.   

After this we can do the actual annotation for all four assemblies. The results will be printed on the screen.
```
conda activate annotation_env
cd ~/Physalia_EnvMetagenomics_2022/ANVIO
abricate --db resfinder */CONTIGS_2500nt.fa
```

## Long-read data

### Genome-resolved metagenomics with long-reads

After preparing the anvi'o files yourself and exploring the interactive interface by manually binning all four short-read data sets, you should have already a good understanding of the basic of anvi'o. So with the long-read data, we have prepared everything ready for you, so you can go straight to the interactive part.  

The long-read assembly is from a single sample (`INF3`), but there is also long-read data from two additional samples  (`INF1` & `INF2`). So as with short-read data, we have used all samples in the mapping step, although the assembly is only from a single sample.  

In addition, there are extra annotations for the gene calls in the contigs database. You can either use the search tab to search for specific functions or you can inspect a single split and see the annotations for each gene call in the split.   

Below are all the annotations steps. Some steps were done with functions in anvi'o (e.g. COG annotation with `anvi-run-ncbi-cogs`) and some with other programs (e.g. gene calling and functional annotation of the gene calls with `prokka`and taxonomic annotation of gene calls with `centrifuge`).  
Not all of the these programs are installed on the amazon cloud, but they are freely available, so you can install them where you need them.

```
#############################################################################################
##                               DO NOT TRY RUN THIS                                       ##
#############################################################################################

prokka --outdir INF3_PROKKA --prefix INF3 --metagenome --cpus 6 --coverage 0.5 assembly.fasta

anvi-script-process-genbank -i INF3_PROKKA/INF3.gbf -O INF3 \
                            --annotation-source prodigal \
                            --annotation-version v2.6.3

anvi-gen-contigs-database -f INF3-contigs.fa --external-gene-calls INF3-external-gene-calls.txt \
                          -o hifi_contigs.db -n Influent-HiFi -T 6 \
                          --ignore-internal-stop-codons

anvi-import-functions -c hifi_contigs.db  -i INF3-external-functions.txt
anvi-run-hmms -c hifi_contigs.db -T 6
anvi-run-ncbi-cogs -c hifi_contigs.db -T 6
anvi-run-scg-taxonomy -c hifi_contigs.db -T 6
anvi-get-sequences-for-gene-calls -o gene-calls.fa -c hifi_contigs.db
centrifuge -f -x ~/projappl/centrifuge/p_compressed gene-calls.fa -S centrifuge_hits.tsv
anvi-import-taxonomy-for-genes -c hifi_contigs.db -i centrifuge_report.tsv centrifuge_hits.tsv -p centrifuge

blastn -query gene-calls.fa -subject PATH/TO/resfinder.fasta -outfmt 6 -evalue 1e-20 -out resfinder.out -max_target_seqs 1
printf "gene_callers_id\tsource\taccession\tfunction\te_value\n" > resfinder_functions.txt
awk '{print $1"\tResfinder\t\t"$2"\t"$11}' resfinder.out >> resfinder_functions.txt  
anvi-import-functions -c hifi_contigs.db -i resfinder_functions.txt

blastn -query gene-calls.fa -subject PATH/TO/intI1.fasta -outfmt 6 -evalue 1e-20 -out integrons.out -max_target_seqs 1
printf "gene_callers_id\tsource\taccession\tfunction\te_value\n" > integron_functions.txt
awk '{print $1"\tintegron\t\t"$2"\t"$11}' integrons.out >> integron_functions.txt  
anvi-import-functions -c hifi_contigs.db -i integron_functions.txt  
```

The rest of the anvi'o pipeline was done in similar way as with the short-reads (`minimap2` instead of `Bowtie2`).  
The final files can be found from `~/Share/HIFI_ANVIO`. Copy the whole folder to your own course folder.

```
cd ~/Physalia_EnvMetagenomics_2022
cp -r ~/Share/HIFI_ANVIO/ ./
cd HIFI_ANVIO
```

Then launch the interacive interface for some more manual binning. All steps are exactly the same as with the long-reads.  

```
conda activate anvio7_env
anvi-interactive  -c hifi_contigs.db -p SAMPLES-MERGED/PROFILE.db --server-only -P $PORT
```


### ARG annotation and visualisation of the flanking regions

As you might have realised, there are already ARG annotations in the contigs database. So we can just use the search function to search for few example genes.
Few ARGs found from multiple contexts are: `blaOXA-491`, `tet(39)`, `ermB` and two that are found next to each other, `msr(E)` and `mph(E)`.  
Explore the genetic context of these genes in the interactive view of anvi'o by inspecting the splits containing the genes of interest.  

We can also extract these regions from the contigs database with `anvi-export-locus` base on the ARG annotation.

```
conda activate anvio7_env
cd ~/Physalia_EnvMetagenomics_2022/HIFI_ANVIO
mkdir oxa-491
anvi-export-locus -c hifi_contigs.db --num-genes 5,5 --search-term "oxa-491" -O oxa_491 -o oxa-491
```

This will give us each context in a separate fasta file. We will then use prokka to annotate these contigs before visualizing them.  
```
conda activate annotation_env
cd oxa-491

for file in *.fa
do
    prokka $file -o ${file%.fa}_PROKKA --prefix ${file%.fa} --cpus 4
done
```

And visualise the contigs with [clinker](https://github.com/gamcil/clinker).    
When clinker is ready, get the resulting `.html` file using e.g. FileZilla.
```
clinker *_PROKKA/*.gbk -o clinker -p clinker.html
```

Then you can re-run the above steps for some other ARG, starting from the `grep` part.
