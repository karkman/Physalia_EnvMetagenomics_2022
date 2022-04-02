library(tidyverse)
library(phyloseq)


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

# Lets explore  the taxonomic annotations
megan_genus

# See the first few OTUs
otu_table(megan_genus) %>% head
tax_table(megan_genus) %>% head

# Take a look at the samples
sample_data(megan_genus)
sample_sums(megan_genus)

# Plot the number of reads within each sample
sample_sums(megan_genus) %>% barplot(las=3)

# See the top 10 OTUs (most abundant throughout all samples)
megan_abund <- taxa_sums(megan_genus) %>%
  sort(decreasing = TRUE) %>%
  head(10) %>%
  names()

# See taxonomy for these OTUs
tax_table(megan_genus)[megan_abund,]

# And their abundance in our samples
otu_table(megan_genus)[megan_abund,]

# heatmap
megan_top10 <- prune_taxa(megan_abund, megan_genus)
otu_table(megan_top10) %>%
                    t() %>%
                    sqrt() %>%
                    as.matrix() %>%
                    heatmap(col=rev(heat.colors(20)))

# Calculate and plot Shannon diversity
otu_table(megan_genus) %>%
                    t() %>%
                    diversity(index="shannon") %>%
                    barplot(ylab="Shannon diversity", las=3)

# Calculate and plot richness
otu_table(megan_genus) %>%
                    t() %>%
                    specnumber() %>%
                    barplot(ylab="Observed taxa", las=3)


## NEEDS TO BE MADE TIDY

# Calculate distance matrix and do ordination
megan_dist <- vegdist(t(otu_table(megan_genus)))
megan_ord <- cmdscale(megan_dist)
megan_ord_df <- data.frame(megan_ord, Season = sample_data(megan_genus)$Season)

# Plot ordination
ggplot(megan_ord_df, aes(x = X1, y = X2, color = Season)) +
  geom_point(size = 3) +
  scale_color_manual(values=c("firebrick", "royalblue")) +
  theme_classic() +
  labs(x = "Axis-1", y = "Axis-2") +
  geom_text(label = row.names(megan_ord_df), nudge_y = 0.03) +
  theme(legend.position = "bottom")


####################################################################
#
#  Now go on and explore the COG and SEED functions.
#
####################################################################

## Optional
# If you can install microViz on your own computer, you can use it for data exploration

# Installing from github requires the devtools package
install.packages("devtools")

# To install the latest "released" version of this package
devtools::install_github("david-barnett/microViz@0.9.0") # check 0.9.0 is the latest release

# load the package
library(microViz)

# open an interactive shinyapp for some basic data exploration
ord_explore(megan_genus)
