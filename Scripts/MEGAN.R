library(tidyverse)
library(phyloseq)
library(microViz)

# Read metadata
metadata <- read.table("sample_info.txt", sep = "\t", row.names = 1, header = F)
colnames(metadata) <- c("Date", "Season", "Year")

# Read MEGAN taxonomy at the genus level
megan_genus <- import_biom("Comparison-Taxonomy.biom")
sample_data(megan_genus) <- sample_data(metadata)

# Read MEGAN COG functions
megan_COG <- import_biom("Comparison-EGGNOG.biom")
sample_data(megan_COG) <- sample_data(metadata)

# Read MEGAN SEED functions
megan_SEED <- import_biom("Comparison-SEED.biom")
sample_data(megan_SEED) <- sample_data(metadata)
