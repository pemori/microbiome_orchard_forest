# =============================================================================
# Script: 01_prepare_16S_bacterial_objects.R
# Project: From forest to fields: The role of soil microbiome spillover in
#          agroecosystem sustainability
# Author: Pedro Mondaca
# Contact: pedromondaca@outlook.com
# Submitted to: Science Advances
#
# Description:
# This script imports 16S ASV abundance, taxonomy, metadata, and a phylogenetic
# tree; builds phyloseq objects; filters rare ASVs; agglomerates taxa at genus
# and family levels; standardizes unresolved genus labels; prepares crop and
# forest subsets; and exports derived tables and objects for downstream analyses.
#
# Required files in the working directory:
#   - ASV_abundance_Table_16s.csv
#   - TAX_16S.csv
#   - BGQ5.xlsx
#   - tree_16S.nwk
#
# Output files:
#   - written to the folder "results_16S"
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(phyloseq)
  library(ape)
  library(stringr)
})

# =============================================================================
# Paths and setup
# =============================================================================

setwd("c:\\Users\\pedro\\Research\\microbiome_orchard_forest\\Sci Adv - v2\\GitHub")

input_abund <- "ASV_abundance_Table_16s.csv"
input_tax   <- "TAX_16S.csv"
input_map   <- "BGQ5.xlsx"
input_tree  <- "tree_16S.nwk"
outdir      <- "results_16S"

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

required_files <- c(input_abund, input_tax, input_map, input_tree)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required files: ", paste(missing_files, collapse = ", "))
}

# =============================================================================
# Load ASV abundance table
# =============================================================================

abunda <- read.csv(
  input_abund,
  sep = ";",
  header = FALSE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

abunda_t <- as.data.frame(t(abunda), stringsAsFactors = FALSE)

seqs <- abunda_t[1, ]
samples <- abunda_t[, 1]

samples <- samples[samples != "#OTU ID"]
seqs <- seqs[seqs != "#OTU ID"]

abund <- abunda_t[-1, -1, drop = FALSE]
rownames(abund) <- samples
colnames(abund) <- seqs

abund <- abund %>%
  mutate(across(everything(), ~ as.numeric(as.character(.x))))

abund <- as.data.frame(abund, check.names = FALSE)

# =============================================================================
# Load taxonomy table
# =============================================================================

tax <- read.csv(
  input_tax,
  header = FALSE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

tax_sep <- tax %>%
  separate(
    V1,
    into = c("seqs", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Specie"),
    sep = ";",
    fill = "right",
    extra = "drop"
  ) %>%
  mutate(seqs = str_remove(seqs, pattern = "\t"))

rownames(tax_sep) <- tax_sep$seqs
tax_num <- tax_sep[, -1, drop = FALSE]

# =============================================================================
# Load metadata (mapping file)
# =============================================================================

BGQ <- read_excel(input_map, sheet = "All")
bgq <- as.data.frame(BGQ[, -1, drop = FALSE])

map_idx <- match(samples, BGQ$IMR)

if (any(is.na(map_idx))) {
  stop(
    "Some sample names from the abundance table were not found in BGQ$IMR: ",
    paste(samples[is.na(map_idx)], collapse = ", ")
  )
}

mapping <- bgq[map_idx, , drop = FALSE]
rownames(mapping) <- samples
mapping <- mapping %>%
  mutate(across(where(is.character), as.factor))

# =============================================================================
# Load phylogenetic tree
# =============================================================================

tree <- read.tree(input_tree)

shared_asvs <- Reduce(
  intersect,
  list(colnames(abund), rownames(tax_num), tree$tip.label)
)

if (length(shared_asvs) == 0) {
  stop("No shared ASV IDs among abundance table, taxonomy table and phylogenetic tree.")
}

abund   <- abund[, shared_asvs, drop = FALSE]
tax_num <- tax_num[shared_asvs, , drop = FALSE]
tree    <- ape::keep.tip(tree, shared_asvs)

OTU <- otu_table(as.matrix(abund), taxa_are_rows = FALSE)
SAM <- sample_data(mapping, errorIfNULL = TRUE)
TAX <- tax_table(as.matrix(tax_num))
PHY <- phy_tree(tree)

data_phylo0 <- phyloseq(OTU, TAX, SAM, PHY)

data_phylo <- subset_taxa(data_phylo0, Domain %in% "d__Bacteria")

# =============================================================================
# Filter rare ASVs
# Criterion: ASVs with counts > 2 in more than 1.587% of samples
# =============================================================================

data_phylo_filt <- filter_taxa(
  data_phylo,
  function(x) sum(x > 2) > (0.01587 * length(x)),
  prune = TRUE
)

data_otu_filt <- data.frame(otu_table(data_phylo_filt))

cat("Number of bacterial ASVs before filtering:", ntaxa(data_phylo), "\n")
cat("Number of bacterial ASVs after filtering:", ntaxa(data_phylo_filt), "\n")
cat("Total reads before filtering:", sum(otu_table(data_phylo)), "\n")
cat("Total reads after filtering:", sum(otu_table(data_phylo_filt)), "\n")

# =============================================================================
# Agglomeration at genus level
# =============================================================================

data_phylo_filt_gen <- tax_glom(data_phylo_filt, taxrank = "Genus")

tax_gen <- as.data.frame(tax_table(data_phylo_filt_gen))

genus_values    <- as.character(tax_gen$Genus)
family_values   <- as.character(tax_gen$Family)
order_values    <- as.character(tax_gen$Order)
class_values    <- as.character(tax_gen$Class)
phylum_values   <- as.character(tax_gen$Phylum)
domain_values   <- as.character(tax_gen$Domain)

genus_values[genus_values == "g__uncultured"] <- family_values[genus_values == "g__uncultured"]
genus_values[genus_values == "f__uncultured"] <- order_values[genus_values == "f__uncultured"]
genus_values[genus_values == "o__uncultured"] <- class_values[genus_values == "o__uncultured"]
genus_values[genus_values == "c__uncultured"] <- phylum_values[genus_values == "c__uncultured"]
genus_values[genus_values == "p__uncultured"] <- domain_values[genus_values == "p__uncultured"]

genus_values <- sub("^g__", "", genus_values)
genus_values <- sub("^f__", "", genus_values)
genus_values <- sub("^o__", "", genus_values)
genus_values <- sub("^c__", "", genus_values)
genus_values <- sub("^p__", "", genus_values)
genus_values <- sub("^d__", "", genus_values)

genus_values[is.na(genus_values)] <- "Unclassified"
genus_values[genus_values == ""] <- "Unclassified"

genus_values <- make.unique(genus_values, sep = "_")

tax_table(data_phylo_filt_gen)[, "Genus"] <- genus_values
tax_table_gen <- as.data.frame(tax_table(data_phylo_filt_gen))

if ("Specie" %in% colnames(tax_table_gen)) {
  tax_table_filt_gen <- tax_table_gen[, colnames(tax_table_gen) != "Specie", drop = FALSE]
} else {
  tax_table_filt_gen <- tax_table_gen
}

rownames(tax_table_filt_gen) <- genus_values

data_otu_filt_gen <- data.frame(otu_table(data_phylo_filt_gen))
colnames(data_otu_filt_gen) <- genus_values

OTU_g <- otu_table(as.matrix(data_otu_filt_gen), taxa_are_rows = FALSE)
SAM_g <- sample_data(mapping, errorIfNULL = TRUE)
TAX_g <- tax_table(as.matrix(tax_table_filt_gen))

data_phylo_g_bac <- phyloseq(OTU_g, TAX_g, SAM_g)
data_otu_g_bac <- data.frame(otu_table(data_phylo_g_bac))

# =============================================================================
# Agglomeration at family level for niche breadth analyses
# =============================================================================

data_phylo_filt_fam_b <- tax_glom(data_phylo_filt, taxrank = "Family")
tax_fam_b <- as.data.frame(tax_table(data_phylo_filt_fam_b))
data_otu_filt_fam_b <- data.frame(otu_table(data_phylo_filt_fam_b))

family_names <- as.character(tax_fam_b$Family)
family_names[is.na(family_names) | family_names == ""] <- "Unclassified_Family"
family_names <- make.unique(family_names, sep = "_")

colnames(data_otu_filt_fam_b) <- family_names

data_fam_b_RA <- data_otu_filt_fam_b / rowSums(data_otu_filt_fam_b) * 100

cat("Row sums of family relative abundance table:\n")
print(rowSums(data_fam_b_RA))

forest_fam_b_RA <- data_fam_b_RA[grepl("N", rownames(data_fam_b_RA)), , drop = FALSE]
crops_fam_b_RA  <- data_fam_b_RA[!grepl("N", rownames(data_fam_b_RA)), , drop = FALSE]

# =============================================================================
# Prepare crop and forest subsets
# =============================================================================

phylo_crops_bac  <- subset_samples(data_phylo_filt, cat != "F")
phylo_forest_bac <- subset_samples(data_phylo_filt, cat == "F")

otu_crops_bac  <- data.frame(otu_table(phylo_crops_bac))
otu_forest_bac <- data.frame(otu_table(phylo_forest_bac))

phylo_crops_g_bac  <- subset_samples(data_phylo_g_bac, cat != "F")
phylo_forest_g_bac <- subset_samples(data_phylo_g_bac, cat == "F")

otu_crops_g_bac  <- data.frame(otu_table(phylo_crops_g_bac))
otu_forest_g_bac <- data.frame(otu_table(phylo_forest_g_bac))

# =============================================================================
# Export outputs
# =============================================================================

saveRDS(data_phylo_filt,   file.path(outdir, "phyloseq_16S_asv_filtered.rds"))
saveRDS(data_phylo_g_bac,  file.path(outdir, "phyloseq_16S_genus.rds"))

write.csv(data_otu_filt,        file.path(outdir, "otu_table_16S_asv_filtered.csv"))
write.csv(data_otu_filt_gen,    file.path(outdir, "otu_table_16S_genus.csv"))
write.csv(data_fam_b_RA,        file.path(outdir, "otu_table_16S_family_relative_abundance.csv"))
write.csv(tax_table_filt_gen,   file.path(outdir, "taxonomy_16S_genus.csv"))
write.csv(mapping,              file.path(outdir, "mapping_16S_ordered.csv"))

write.csv(otu_crops_bac,        file.path(outdir, "otu_crops_16S_asv.csv"))
write.csv(otu_forest_bac,       file.path(outdir, "otu_forest_16S_asv.csv"))
write.csv(otu_crops_g_bac,      file.path(outdir, "otu_crops_16S_genus.csv"))
write.csv(otu_forest_g_bac,     file.path(outdir, "otu_forest_16S_genus.csv"))

write.csv(forest_fam_b_RA,      file.path(outdir, "family_relative_abundance_forest.csv"))
write.csv(crops_fam_b_RA,       file.path(outdir, "family_relative_abundance_crops.csv"))

writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo_16S.txt"))
