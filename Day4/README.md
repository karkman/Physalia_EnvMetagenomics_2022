# Day 4

| Time      | Activity                     | Slides                                        | Hands-on                                   |
|-----------|------------------------------|-----------------------------------------------|--------------------------------------------|
| Morning   |  | [Link here](Day3/genome-resolved-metagenomics.pdf) | [Link here](#genome-resolved-metagenomics) |
| Afternoon |  |                                               | [Link here](#genome-resolved-metagenomics) |

### CheckM

### GTDB-Tk

### dREP

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

Then launch the interacive interface for some more manual binning.

```
conda activate anvio7_env
anvi-interactive  -c hifi_contigs.db -p SAMPLES-MERGED/PROFILE.db --server-only -P $PORT
```
