# =============================================================================
# Script: 05_hypothesis_2_16S_composition.R
# Project: From forest to fields: The role of soil microbiome spillover in
#          agroecosystem sustainability
# Author: Pedro Mondaca
# Contact: pedromondaca@outlook.com
# Submitted to: Science Advances
#
# Description:
# This script evaluates Hypothesis 2 for the bacterial community:
# compositional differences between orchard and forest soils diminish with
# increasing surrounding forest cover. It computes Bray-Curtis and weighted
# UniFrac distances, ordinations, PERMANOVA, PERMDISP, MNTD, taxonomic
# composition patterns, Mantel tests, and exploratory multivariate analyses.
#
# Required files in the working directory:
#   - results_16S/phyloseq_16S_asv_filtered.rds
#   - results_16S/mapping_16S_ordered.csv
#   - results_16S_hyp1/orchard_summary_16S.csv
#   - results_16S_hyp1/forest_summary_16S.csv
#
# Output files:
#   - written to the folder "results_16S_hyp2"
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(phyloseq)
  library(vegan)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(sf)
  library(terra)
  library(geodata)
  library(picante)
  library(microeco)
  library(reshape2)
  library(grid)
})

# =============================================================================
# Paths and setup
# =============================================================================

setwd("c:\\Users\\pedro\\Research\\microbiome_orchard_forest\\Sci Adv - v2\\GitHub")

input_phyloseq <- "results_16S/phyloseq_16S_asv_filtered.rds"
input_mapping  <- "results_16S/mapping_16S_ordered.csv"
input_orchard  <- "results_16S_hyp1/orchard_summary_16S.csv"
input_forest   <- "results_16S_hyp1/forest_summary_16S.csv"
outdir         <- "results_16S_hyp2"

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

required_files <- c(input_phyloseq, input_mapping, input_orchard, input_forest)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required files: ", paste(missing_files, collapse = ", "))
}

# =============================================================================
# Load objects
# =============================================================================

data_phylo_filt   <- readRDS(input_phyloseq)
mapping           <- read.csv(input_mapping, row.names = 1, check.names = FALSE)
orchard_summary_df <- read.csv(input_orchard, check.names = FALSE)
forest_summary_df  <- read.csv(input_forest, check.names = FALSE)

mapping <- mapping[sample_names(data_phylo_filt), , drop = FALSE]

if (!identical(rownames(mapping), sample_names(data_phylo_filt))) {
  stop("Mapping file row names do not match phyloseq sample names.")
}

sample_data(data_phylo_filt) <- sample_data(mapping)

mapping_df <- data.frame(mapping, check.names = FALSE)
mapping_df$SampleID <- rownames(mapping_df)

phylo_crops_bac  <- subset_samples(data_phylo_filt, cat != "F")
phylo_forest_bac <- subset_samples(data_phylo_filt, cat == "F")

mapping_orchard <- mapping_df %>%
  dplyr::filter(Group == "Orchard") %>%
  droplevels()

mapping_forest <- mapping_df %>%
  dplyr::filter(Group == "Forest") %>%
  droplevels()

# =============================================================================
# Hypothesis 2
# Compositional differences between orchard and forest diminish with increasing
# forest cover
# =============================================================================

# =============================================================================
# Beta diversity
# =============================================================================

# Resolve polytomies before UniFrac
tree_all <- phy_tree(data_phylo_filt)
tree_all_bin <- ape::multi2di(tree_all, random = TRUE)
phy_tree(data_phylo_filt) <- tree_all_bin

tree_orch <- phy_tree(phylo_crops_bac)
tree_orch_bin <- ape::multi2di(tree_orch, random = TRUE)
phy_tree(phylo_crops_bac) <- tree_orch_bin

tree_forest <- phy_tree(phylo_forest_bac)
tree_forest_bin <- ape::multi2di(tree_forest, random = TRUE)
phy_tree(phylo_forest_bac) <- tree_forest_bin

bray_dist    <- phyloseq::distance(data_phylo_filt, method = "bray")
unifrac_dist <- phyloseq::UniFrac(data_phylo_filt, weighted = TRUE, normalized = TRUE)

bray_dist_orchard    <- phyloseq::distance(phylo_crops_bac, method = "bray")
unifrac_dist_orchard <- phyloseq::UniFrac(phylo_crops_bac, weighted = TRUE, normalized = TRUE)

bray_dist_forest    <- phyloseq::distance(phylo_forest_bac, method = "bray")
unifrac_dist_forest <- phyloseq::UniFrac(phylo_forest_bac, weighted = TRUE, normalized = TRUE)

saveRDS(bray_dist,            file.path(outdir, "bray_dist_all.rds"))
saveRDS(unifrac_dist,         file.path(outdir, "unifrac_dist_all.rds"))
saveRDS(bray_dist_orchard,    file.path(outdir, "bray_dist_orchard.rds"))
saveRDS(unifrac_dist_orchard, file.path(outdir, "unifrac_dist_orchard.rds"))
saveRDS(bray_dist_forest,     file.path(outdir, "bray_dist_forest.rds"))
saveRDS(unifrac_dist_forest,  file.path(outdir, "unifrac_dist_forest.rds"))

# =============================================================================
# PCoA
# =============================================================================

pcoa_bc <- ordinate(data_phylo_filt, method = "PCoA", distance = bray_dist)
pcoa_uf <- ordinate(data_phylo_filt, method = "PCoA", distance = unifrac_dist)

# Bray-Curtis
var_explained_bc <- pcoa_bc$values$Relative_eig * 100
x_label_bc <- paste0("PCoA 1 (", round(var_explained_bc[1], 1), "%)")
y_label_bc <- paste0("PCoA 2 (", round(var_explained_bc[2], 1), "%)")

ord_plot_data_bc <- plot_ordination(data_phylo_filt, pcoa_bc, justDF = TRUE)

meta_ord <- data.frame(sample_data(data_phylo_filt), check.names = FALSE)

ord_plot_data_bc$Group  <- meta_ord$Group
ord_plot_data_bc$P_2500 <- meta_ord$P_2500
ord_plot_data_bc$P_100  <- meta_ord$P_100

p_pcoa_bc_2500 <- ggplot(ord_plot_data_bc, aes(x = Axis.1, y = Axis.2)) +
  geom_point(
    data = subset(ord_plot_data_bc, Group == "Forest"),
    shape = 17, color = "darkgreen", size = 3
  ) +
  geom_point(
    data = subset(ord_plot_data_bc, Group == "Orchard"),
    aes(color = P_2500),
    shape = 16, size = 3
  ) +
  scale_color_gradient(low = "gold", high = "darkgreen", name = expression(P[2500])) +
  theme_bw(base_size = 12) +
  labs(x = x_label_bc, y = y_label_bc) +
  theme(legend.position = "right", panel.grid = element_blank())

p_pcoa_bc_100 <- ggplot(ord_plot_data_bc, aes(x = Axis.1, y = Axis.2)) +
  geom_point(
    data = subset(ord_plot_data_bc, Group == "Forest"),
    shape = 17, color = "darkgreen", size = 3
  ) +
  geom_point(
    data = subset(ord_plot_data_bc, Group == "Orchard"),
    aes(color = P_100),
    shape = 16, size = 3
  ) +
  scale_color_gradient(low = "gold", high = "darkgreen", name = expression(P[100])) +
  theme_bw(base_size = 12) +
  labs(x = x_label_bc, y = y_label_bc) +
  theme(legend.position = "right", panel.grid = element_blank())

ggsave(file.path(outdir, "pcoa_bray_P2500.png"), p_pcoa_bc_2500, width = 5.5, height = 5, dpi = 300)
ggsave(file.path(outdir, "pcoa_bray_P100.png"),  p_pcoa_bc_100,  width = 5.5, height = 5, dpi = 300)

# Weighted UniFrac
var_explained_uf <- pcoa_uf$values$Relative_eig * 100
x_label_uf <- paste0("PCoA 1 (", round(var_explained_uf[1], 1), "%)")
y_label_uf <- paste0("PCoA 2 (", round(var_explained_uf[2], 1), "%)")

ord_plot_data_uf <- plot_ordination(data_phylo_filt, pcoa_uf, justDF = TRUE)

ord_plot_data_uf$Group  <- meta_ord$Group
ord_plot_data_uf$P_2500 <- meta_ord$P_2500
ord_plot_data_uf$P_100  <- meta_ord$P_100

p_pcoa_uf_2500 <- ggplot(ord_plot_data_uf, aes(x = Axis.1, y = Axis.2)) +
  geom_point(
    data = subset(ord_plot_data_uf, Group == "Forest"),
    shape = 17, color = "darkgreen", size = 3
  ) +
  geom_point(
    data = subset(ord_plot_data_uf, Group == "Orchard"),
    aes(color = P_2500),
    shape = 16, size = 3
  ) +
  scale_color_gradient(low = "gold", high = "darkgreen", name = expression(P[2500])) +
  theme_bw(base_size = 12) +
  labs(x = x_label_uf, y = y_label_uf) +
  theme(legend.position = "right", panel.grid = element_blank())

p_pcoa_uf_100 <- ggplot(ord_plot_data_uf, aes(x = Axis.1, y = Axis.2)) +
  geom_point(
    data = subset(ord_plot_data_uf, Group == "Forest"),
    shape = 17, color = "darkgreen", size = 3
  ) +
  geom_point(
    data = subset(ord_plot_data_uf, Group == "Orchard"),
    aes(color = P_100),
    shape = 16, size = 3
  ) +
  scale_color_gradient(low = "gold", high = "darkgreen", name = expression(P[100])) +
  theme_bw(base_size = 12) +
  labs(x = x_label_uf, y = y_label_uf) +
  theme(legend.position = "right", panel.grid = element_blank())

ggsave(file.path(outdir, "pcoa_unifrac_P2500.png"), p_pcoa_uf_2500, width = 5.5, height = 5, dpi = 300)
ggsave(file.path(outdir, "pcoa_unifrac_P100.png"),  p_pcoa_uf_100,  width = 5.5, height = 5, dpi = 300)

saveRDS(pcoa_bc, file.path(outdir, "pcoa_bray_object.rds"))
saveRDS(pcoa_uf, file.path(outdir, "pcoa_unifrac_object.rds"))

# =============================================================================
# PERMANOVA
# =============================================================================

adonis_group_bray <- adonis2(bray_dist ~ Group, data = mapping_df, permutations = 9999)
adonis_group_uf   <- adonis2(unifrac_dist ~ Group, data = mapping_df, permutations = 9999)

adonis_block_bray <- adonis2(bray_dist ~ Block, data = mapping_df, permutations = 9999)
adonis_block_uf   <- adonis2(unifrac_dist ~ Block, data = mapping_df, permutations = 9999)

adonis_group_block_bray <- adonis2(bray_dist ~ Group * Block, data = mapping_df, permutations = 9999)
adonis_group_block_uf   <- adonis2(unifrac_dist ~ Group * Block, data = mapping_df, permutations = 9999)

adonis_xy_bray <- adonis2(bray_dist ~ Coo_X * Coo_Y, data = mapping_df, permutations = 9999)
adonis_xy_uf   <- adonis2(unifrac_dist ~ Coo_X * Coo_Y, data = mapping_df, permutations = 9999)

adonis_orchard_block_bray <- adonis2(bray_dist_orchard ~ Block, data = mapping_orchard, permutations = 9999)
adonis_orchard_block_uf   <- adonis2(unifrac_dist_orchard ~ Block, data = mapping_orchard, permutations = 9999)

adonis_orchard_xy_bray <- adonis2(bray_dist_orchard ~ Coo_X * Coo_Y, data = mapping_orchard, permutations = 9999)
adonis_orchard_xy_uf   <- adonis2(unifrac_dist_orchard ~ Coo_X * Coo_Y, data = mapping_orchard, permutations = 9999)

adonis_orchard_P100_bray  <- adonis2(bray_dist_orchard ~ P_100,  data = mapping_orchard, permutations = 9999)
adonis_orchard_P250_bray  <- adonis2(bray_dist_orchard ~ P_250,  data = mapping_orchard, permutations = 9999)
adonis_orchard_P500_bray  <- adonis2(bray_dist_orchard ~ P_500,  data = mapping_orchard, permutations = 9999)
adonis_orchard_P1000_bray <- adonis2(bray_dist_orchard ~ P_1000, data = mapping_orchard, permutations = 9999)
adonis_orchard_P2500_bray <- adonis2(bray_dist_orchard ~ P_2500, data = mapping_orchard, permutations = 9999)

adonis_orchard_P100_uf  <- adonis2(unifrac_dist_orchard ~ P_100,  data = mapping_orchard, permutations = 9999)
adonis_orchard_P250_uf  <- adonis2(unifrac_dist_orchard ~ P_250,  data = mapping_orchard, permutations = 9999)
adonis_orchard_P500_uf  <- adonis2(unifrac_dist_orchard ~ P_500,  data = mapping_orchard, permutations = 9999)
adonis_orchard_P1000_uf <- adonis2(unifrac_dist_orchard ~ P_1000, data = mapping_orchard, permutations = 9999)
adonis_orchard_P2500_uf <- adonis2(unifrac_dist_orchard ~ P_2500, data = mapping_orchard, permutations = 9999)

adonis_forest_block_bray <- adonis2(bray_dist_forest ~ Block, data = mapping_forest, permutations = 9999)
adonis_forest_block_uf   <- adonis2(unifrac_dist_forest ~ Block, data = mapping_forest, permutations = 9999)

adonis_forest_xy_bray <- adonis2(bray_dist_forest ~ Coo_X * Coo_Y, data = mapping_forest, permutations = 9999)
adonis_forest_xy_uf   <- adonis2(unifrac_dist_forest ~ Coo_X * Coo_Y, data = mapping_forest, permutations = 9999)

capture.output(adonis_group_bray,         file = file.path(outdir, "adonis_group_bray.txt"))
capture.output(adonis_group_uf,           file = file.path(outdir, "adonis_group_unifrac.txt"))
capture.output(adonis_block_bray,         file = file.path(outdir, "adonis_block_bray.txt"))
capture.output(adonis_block_uf,           file = file.path(outdir, "adonis_block_unifrac.txt"))
capture.output(adonis_group_block_bray,   file = file.path(outdir, "adonis_group_block_bray.txt"))
capture.output(adonis_group_block_uf,     file = file.path(outdir, "adonis_group_block_unifrac.txt"))
capture.output(adonis_xy_bray,            file = file.path(outdir, "adonis_xy_bray.txt"))
capture.output(adonis_xy_uf,              file = file.path(outdir, "adonis_xy_unifrac.txt"))
capture.output(adonis_orchard_block_bray, file = file.path(outdir, "adonis_orchard_block_bray.txt"))
capture.output(adonis_orchard_block_uf,   file = file.path(outdir, "adonis_orchard_block_unifrac.txt"))
capture.output(adonis_orchard_xy_bray,    file = file.path(outdir, "adonis_orchard_xy_bray.txt"))
capture.output(adonis_orchard_xy_uf,      file = file.path(outdir, "adonis_orchard_xy_unifrac.txt"))
capture.output(adonis_orchard_P100_bray,  file = file.path(outdir, "adonis_orchard_P100_bray.txt"))
capture.output(adonis_orchard_P250_bray,  file = file.path(outdir, "adonis_orchard_P250_bray.txt"))
capture.output(adonis_orchard_P500_bray,  file = file.path(outdir, "adonis_orchard_P500_bray.txt"))
capture.output(adonis_orchard_P1000_bray, file = file.path(outdir, "adonis_orchard_P1000_bray.txt"))
capture.output(adonis_orchard_P2500_bray, file = file.path(outdir, "adonis_orchard_P2500_bray.txt"))
capture.output(adonis_orchard_P100_uf,    file = file.path(outdir, "adonis_orchard_P100_unifrac.txt"))
capture.output(adonis_orchard_P250_uf,    file = file.path(outdir, "adonis_orchard_P250_unifrac.txt"))
capture.output(adonis_orchard_P500_uf,    file = file.path(outdir, "adonis_orchard_P500_unifrac.txt"))
capture.output(adonis_orchard_P1000_uf,   file = file.path(outdir, "adonis_orchard_P1000_unifrac.txt"))
capture.output(adonis_orchard_P2500_uf,   file = file.path(outdir, "adonis_orchard_P2500_unifrac.txt"))
capture.output(adonis_forest_block_bray,  file = file.path(outdir, "adonis_forest_block_bray.txt"))
capture.output(adonis_forest_block_uf,    file = file.path(outdir, "adonis_forest_block_unifrac.txt"))
capture.output(adonis_forest_xy_bray,     file = file.path(outdir, "adonis_forest_xy_bray.txt"))
capture.output(adonis_forest_xy_uf,       file = file.path(outdir, "adonis_forest_xy_unifrac.txt"))

# =============================================================================
# PERMDISP
# =============================================================================

permdisp_bray <- betadisper(bray_dist, mapping_df$Group)
permdisp_bray_test <- permutest(permdisp_bray, pairwise = TRUE, permutations = 9999)

permdisp_uf <- betadisper(unifrac_dist, mapping_df$Group)
permdisp_uf_test <- permutest(permdisp_uf, pairwise = TRUE, permutations = 9999)

capture.output(permdisp_bray_test, file = file.path(outdir, "permdisp_bray.txt"))
capture.output(permdisp_uf_test,   file = file.path(outdir, "permdisp_unifrac.txt"))


# =============================================================================
# MNTD
# =============================================================================

shared_asvs_tree <- intersect(taxa_names(data_phylo_filt), phy_tree(data_phylo_filt)$tip.label)
data_phylo_pruned_all <- prune_taxa(shared_asvs_tree, data_phylo_filt)

otu_all <- as(otu_table(data_phylo_pruned_all), "matrix")
if (taxa_are_rows(data_phylo_pruned_all)) otu_all <- t(otu_all)

dist_mat_all <- cophenetic(phy_tree(data_phylo_pruned_all))
mntd_values_all <- mntd(otu_all, dist_mat_all, abundance.weighted = TRUE)
names(mntd_values_all) <- rownames(otu_all)

metadata_all <- data.frame(sample_data(data_phylo_filt), check.names = FALSE)
metadata_all$SampleID <- rownames(metadata_all)
metadata_all$mntd <- mntd_values_all[metadata_all$SampleID]

write.csv(metadata_all, file.path(outdir, "summary_mntd_all_samples.csv"), row.names = FALSE)

shared_asvs_tree_orch <- intersect(taxa_names(phylo_crops_bac), phy_tree(phylo_crops_bac)$tip.label)
data_phylo_pruned_orch <- prune_taxa(shared_asvs_tree_orch, phylo_crops_bac)

otu_orch <- as(otu_table(data_phylo_pruned_orch), "matrix")
if (taxa_are_rows(data_phylo_pruned_orch)) otu_orch <- t(otu_orch)

dist_mat_orch <- cophenetic(phy_tree(data_phylo_pruned_orch))
mntd_values_orch <- mntd(otu_orch, dist_mat_orch, abundance.weighted = TRUE)
names(mntd_values_orch) <- rownames(otu_orch)

metadata_orch <- data.frame(sample_data(phylo_crops_bac), check.names = FALSE)
metadata_orch$SampleID <- rownames(metadata_orch)
metadata_orch$mntd <- mntd_values_orch[metadata_orch$SampleID]

if ("SampleID" %in% colnames(orchard_summary_df)) {
  orchard_summary_df$mntd <- metadata_orch$mntd[match(orchard_summary_df$SampleID, metadata_orch$SampleID)]
} else {
  orchard_summary_df$mntd <- metadata_orch$mntd[match(rownames(orchard_summary_df), metadata_orch$SampleID)]
}

write.csv(orchard_summary_df, file.path(outdir, "orchard_summary_with_mntd.csv"), row.names = FALSE)

cat("NAs in metadata_all$mntd:", sum(is.na(metadata_all$mntd)), "\n")
cat("NAs in metadata_orch$mntd:", sum(is.na(metadata_orch$mntd)), "\n")
head(names(mntd_values_all))
head(metadata_all$SampleID)

p_mntd_group <- ggplot(metadata_all, aes(x = Group, y = mntd, fill = Group)) +
  geom_boxplot(alpha = 0.7) +
  theme_bw() +
  labs(x = "Group", y = "MNTD") +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = "none")

ggsave(file.path(outdir, "mntd_by_group.png"), p_mntd_group, width = 5, height = 4.5, dpi = 300)

# =============================================================================
# Composition of orchard-only, forest-only, and shared ASVs
# =============================================================================

orchard_phy <- subset_samples(data_phylo_filt, Group == "Orchard")
forest_phy  <- subset_samples(data_phylo_filt, Group == "Forest")

orchard_mat <- as(otu_table(orchard_phy), "matrix")
forest_mat  <- as(otu_table(forest_phy), "matrix")

if (!taxa_are_rows(orchard_phy)) orchard_mat <- t(orchard_mat)
if (!taxa_are_rows(forest_phy))  forest_mat  <- t(forest_mat)

orchard_asvs <- rownames(orchard_mat)[rowSums(orchard_mat) > 0]
forest_asvs  <- rownames(forest_mat)[rowSums(forest_mat) > 0]

shared_asvs         <- intersect(orchard_asvs, forest_asvs)
orchard_unique_asvs <- setdiff(orchard_asvs, forest_asvs)
forest_unique_asvs  <- setdiff(forest_asvs, orchard_asvs)

all_asvs <- taxa_names(data_phylo_filt)

asv_group <- case_when(
  all_asvs %in% orchard_unique_asvs ~ "Orchard-only",
  all_asvs %in% forest_unique_asvs  ~ "Forest-only",
  all_asvs %in% shared_asvs         ~ "Shared",
  TRUE ~ NA_character_
)
names(asv_group) <- all_asvs

otu_mat <- as(otu_table(data_phylo_filt), "matrix")
if (!taxa_are_rows(data_phylo_filt)) otu_mat <- t(otu_mat)

otu_rel <- sweep(otu_mat, 2, colSums(otu_mat), FUN = "/")
otu_rel <- as.data.frame(otu_rel)
otu_rel$ASV <- rownames(otu_rel)

tax_df <- as.data.frame(tax_table(data_phylo_filt))
tax_df$ASV <- rownames(tax_df)

otu_rel <- left_join(otu_rel, tax_df, by = "ASV")
otu_rel$Group_ASV <- asv_group[otu_rel$ASV]

otu_long <- otu_rel %>%
  pivot_longer(
    cols = colnames(otu_mat),
    names_to = "Sample",
    values_to = "Rel_Abund"
  )

otu_long_clean <- otu_long %>%
  mutate(Phylum = gsub("^p__", "", Phylum))

phylum_group <- otu_long_clean %>%
  filter(!is.na(Group_ASV)) %>%
  group_by(Group_ASV, Phylum) %>%
  summarise(Abund = sum(Rel_Abund, na.rm = TRUE), .groups = "drop")

top_phyla <- phylum_group %>%
  group_by(Phylum) %>%
  summarise(Total = sum(Abund, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  slice_head(n = 9) %>%
  pull(Phylum)

phylum_group <- phylum_group %>%
  mutate(Phylum = ifelse(Phylum %in% top_phyla, Phylum, "Others"))

phylum_group_norm <- phylum_group %>%
  group_by(Group_ASV, Phylum) %>%
  summarise(Abund = sum(Abund, na.rm = TRUE), .groups = "drop") %>%
  group_by(Group_ASV) %>%
  mutate(Prop = Abund / sum(Abund, na.rm = TRUE)) %>%
  ungroup()

phylum_group_norm$Group_ASV <- factor(
  phylum_group_norm$Group_ASV,
  levels = c("Orchard-only", "Shared", "Forest-only")
)

phylum_order <- phylum_group_norm %>%
  group_by(Phylum) %>%
  summarise(Total = sum(Prop, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  pull(Phylum)

phylum_order <- c(setdiff(phylum_order, "Others"), "Others")
phylum_order <- rev(phylum_order)

phylum_group_norm$Phylum <- factor(phylum_group_norm$Phylum, levels = phylum_order)

phylum_colors <- c(
  "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
  "#E6AB02", "#A6761D", "#1F78B4", "#B2DF8A", "#CCCCCC"
)

p_phylum_group <- ggplot(phylum_group_norm, aes(x = Group_ASV, y = Prop, fill = Phylum)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = rev(phylum_colors)) +
  labs(x = NULL, y = "Relative abundance (%)", fill = "Phylum") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    legend.position = "right"
  )

ggsave(file.path(outdir, "phylum_composition_shared_unique.png"), p_phylum_group, width = 5.82, height = 5.08, dpi = 300)

# =============================================================================
# Mantel tests
# =============================================================================

env_matrix <- mapping_orchard %>%
  dplyr::select(pH, SOM, N, P, claysilt) %>%
  as.data.frame()

cover_matrix <- mapping_orchard %>%
  dplyr::select(P_2500, P_250) %>%
  as.data.frame()

geo_matrix <- mapping_orchard %>%
  dplyr::select(Coo_X, Coo_Y) %>%
  as.data.frame()

rownames(env_matrix)   <- mapping_orchard$SampleID
rownames(cover_matrix) <- mapping_orchard$SampleID
rownames(geo_matrix)   <- mapping_orchard$SampleID

scale_env  <- scale(env_matrix, center = TRUE, scale = TRUE)
dist_env   <- dist(scale_env, method = "euclidean")
dist_cover <- dist(cover_matrix, method = "euclidean")
dist_geo   <- dist(geo_matrix, method = "euclidean")

utm_sf <- st_as_sf(geo_matrix, coords = c("Coo_X", "Coo_Y"), crs = 32719)
latlon <- st_transform(utm_sf, 4326)

coord_matrix <- as.data.frame(st_coordinates(latlon))
colnames(coord_matrix) <- c("X", "Y")
rownames(coord_matrix) <- mapping_orchard$SampleID

coords_vect <- vect(coord_matrix, geom = c("X", "Y"), crs = "EPSG:4326")

if (exists("prec")) {
  prec_monthly <- terra::extract(prec, coords_vect)
  prec_matrix  <- prec_monthly[, -1, drop = FALSE]
  dist_prec    <- dist(prec_matrix, method = "euclidean")
}

if (exists("tavg")) {
  tavg_monthly <- terra::extract(tavg, coords_vect)
  tavg_matrix  <- tavg_monthly[, -1, drop = FALSE]
  dist_tavg    <- dist(tavg_matrix, method = "euclidean")
}

bray_env   <- mantel(bray_dist_orchard, dist_env,   method = "spearman", permutations = 9999, na.rm = TRUE)
bray_geo   <- mantel(bray_dist_orchard, dist_geo,   method = "spearman", permutations = 9999, na.rm = TRUE)
bray_cover <- mantel(bray_dist_orchard, dist_cover, method = "spearman", permutations = 9999, na.rm = TRUE)

capture.output(bray_env,   file = file.path(outdir, "mantel_bray_env.txt"))
capture.output(bray_geo,   file = file.path(outdir, "mantel_bray_geo.txt"))
capture.output(bray_cover, file = file.path(outdir, "mantel_bray_cover.txt"))

if (exists("dist_prec")) {
  bray_prec <- mantel(bray_dist_orchard, dist_prec, method = "spearman", permutations = 9999, na.rm = TRUE)
  capture.output(bray_prec, file = file.path(outdir, "mantel_bray_prec.txt"))
}

if (exists("dist_tavg")) {
  bray_tavg <- mantel(bray_dist_orchard, dist_tavg, method = "spearman", permutations = 9999, na.rm = TRUE)
  capture.output(bray_tavg, file = file.path(outdir, "mantel_bray_tavg.txt"))
}

# =============================================================================
# composition analysis
# =============================================================================

tax_data <- data.frame(tax_table(data_phylo_filt))
otu_data <- data.frame(otu_table(data_phylo_filt))
sam_data <- data.frame(sample_data(data_phylo_filt), check.names = FALSE)

otu_data_t <- as.data.frame(t(otu_data))
rownames(sam_data) <- colnames(otu_data_t)

dataset <- microtable$new(
  sample_table = sam_data,
  otu_table    = otu_data_t,
  tax_table    = tax_data
)

dataset$tax_table <- subset(dataset$tax_table, Domain == "d__Bacteria")
dataset$sample_table <- dataset$sample_table %>%
  dplyr::filter(Group %in% c("Forest", "Orchard")) %>%
  droplevels()

dataset$tidy_dataset()
dataset$cal_abund()
dataset$cal_alphadiv(PD = FALSE)
dataset$cal_betadiv(unifrac = FALSE)

dataset$sample_table$Group <- factor(dataset$sample_table$Group, levels = c("Forest", "Orchard"))

# Orchard-only dataset kept for downstream ordination analyses if needed
otu_c <- data.frame(t(otu_table(phylo_crops_bac)))
tax_c <- data.frame(tax_table(phylo_crops_bac))
sam_c <- data.frame(sample_data(phylo_crops_bac), check.names = FALSE)
rownames(sam_c) <- colnames(otu_c)

dataset_c <- microtable$new(
  sample_table = sam_c,
  otu_table    = otu_c,
  tax_table    = tax_c
)
dataset_c$tidy_dataset()

# Heatmap by Group
t_heat <- trans_abund$new(dataset = dataset, taxrank = "Family", ntaxa = 30)
p_heat <- t_heat$plot_heatmap(facet = c("Group"), xtext_keep = FALSE, withmargin = FALSE)
ggsave(file.path(outdir, "heatmap_family_group.png"), p_heat, width = 8.36, height = 6.8, dpi = 300)

# Phylum barplot by Group
t_bar1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 8)
p_bar1 <- t_bar1$plot_bar(
  others_color = "grey70",
  facet = c("Group"),
  xtext_keep = FALSE,
  legend_text_italic = FALSE
)
ggsave(file.path(outdir, "barplot_phylum_group.png"), p_bar1, width = 8, height = 4, dpi = 300)

# Group mean phylum barplot: Forest vs Orchard only
t_bar2 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 10, groupmean = "Group")
p_bar2 <- t_bar2$plot_bar(others_color = "grey70", legend_text_italic = FALSE) +
  theme_classic() +
  theme(axis.title.y = element_text(size = 18))
ggsave(file.path(outdir, "barplot_phylum_groupmean.png"), p_bar2, width = 5.78, height = 4.26, dpi = 300)

# =============================================================================
# Differential abundance within microeco
# =============================================================================

dataset$cal_abund()

t_diff <- trans_diff$new(
  dataset = dataset,
  method = "ALDEx2_kw",
  group = "Group",
  taxa_level = "Family"
)

p_diff <- t_diff$plot_diff_abund(use_number = 1:11, group_order = c("Forest", "Orchard"))
ggsave(file.path(outdir, "aldex2_family_diff.png"), p_diff, width = 7, height = 5, dpi = 300)

abund_family <- dataset$taxa_abund$Family

abund_long <- reshape2::melt(as.matrix(abund_family))
colnames(abund_long) <- c("Family", "SampleID", "Abundance")
abund_long$Family <- gsub(".*f__", "", abund_long$Family)

df_diff <- data.frame(
  Taxa = t_diff$res_diff$Taxa,
  P.adj = t_diff$res_diff$P.adj,
  stringsAsFactors = FALSE
)

df_diff <- df_diff %>%
  dplyr::filter(!grepl("f__NA", Taxa)) %>%
  mutate(Family = gsub(".*f__", "", Taxa)) %>%
  dplyr::filter(!is.na(Family), Family != "", Family != "NA")

top_families <- df_diff %>%
  dplyr::filter(P.adj < 0.05) %>%
  arrange(P.adj) %>%
  distinct(Family, .keep_all = TRUE) %>%
  slice_head(n = 10) %>%
  pull(Family)

abund_long$Family <- factor(abund_long$Family, levels = rev(top_families))

abund_plot <- abund_long %>%
  dplyr::filter(Family %in% top_families)

meta <- dataset$sample_table
meta$SampleID <- rownames(meta)

abund_plot <- abund_plot %>%
  left_join(meta %>% dplyr::select(SampleID, Group), by = "SampleID")

abund_summary <- abund_plot %>%
  group_by(Family, Group) %>%
  summarise(
    mean_abund = mean(Abundance, na.rm = TRUE),
    se_abund = sd(Abundance, na.rm = TRUE) / sqrt(sum(!is.na(Abundance))),
    .groups = "drop"
  )

abund_summary$Family <- factor(abund_summary$Family, levels = levels(abund_plot$Family))
abund_summary$Group <- factor(abund_summary$Group, levels = c("Forest", "Orchard"))

custom_colors <- c("Orchard" = "gold", "Forest" = "darkgreen")

p_family_diff <- ggplot(abund_summary, aes(x = Family, y = mean_abund, fill = Group)) +
  geom_col(position = position_dodge(width = 0.9), color = "black") +
  geom_errorbar(
    aes(ymin = mean_abund - se_abund, ymax = mean_abund + se_abund),
    position = position_dodge(width = 0.9),
    width = 0.2,
    color = "black"
  ) +
  scale_fill_manual(values = custom_colors) +
  coord_flip() +
  labs(x = NULL, y = "Relative abundance", fill = "Land use") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(face = "plain", color = "black"),
    axis.text.x = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    axis.title.x = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(outdir, "family_diff_barplot.png"), p_family_diff, width = 6, height = 7.5, dpi = 300)

write.csv(df_diff, file.path(outdir, "aldex2_family_diff_table.csv"), row.names = FALSE)
write.csv(abund_summary, file.path(outdir, "family_diff_barplot_table.csv"), row.names = FALSE)

# =============================================================================
# dbRDA and RDA
# =============================================================================

tab <- as.data.frame(dataset_c$sample_table, check.names = FALSE)

vars_needed <- c("pH", "SOM", "N", "P", "claysilt", "SOMTXT", "P_250", "P_1000", "P_2500")
missing_vars <- setdiff(vars_needed, colnames(tab))
if (length(missing_vars) > 0) {
  stop("Missing variables in dataset_c$sample_table: ", paste(missing_vars, collapse = ", "))
}

tab_diver5 <- tab %>%
  dplyr::select(pH, SOM, N, P, claysilt, SOMTXT, P_250, P_1000, P_2500) %>%
  rename("Clay+Silt" = claysilt) %>%
  rename("SOM/(Clay+Silt)" = SOMTXT) %>%
  rename("FC250" = P_250) %>%
  rename("FC1000" = P_1000) %>%
  rename("FC2500" = P_2500)

t_env <- trans_env$new(dataset = dataset_c, add_data = tab_diver5)

t_env$cal_ordination(method = "dbRDA", use_measure = "bray")
t_env$trans_ordination()
p_dbrda <- t_env$plot_ordination()
ggsave(file.path(outdir, "dbrda_bray.png"), p_dbrda, width = 4.49, height = 4.21, dpi = 300)

t_env$cal_ordination_anova()
t_env$cal_ordination_envfit()

capture.output(t_env$res_ordination_terms,  file = file.path(outdir, "dbrda_terms.txt"))
capture.output(t_env$res_ordination_envfit, file = file.path(outdir, "dbrda_envfit.txt"))

t_env$cal_ordination(method = "RDA", taxa_level = "Phylum")
t_env$trans_ordination(
  show_taxa = 15,
  adjust_arrow_length = TRUE,
  max_perc_env = 1.5,
  max_perc_tax = 1.5,
  min_perc_env = 0.2,
  min_perc_tax = 0.2
)

res_trans <- t_env$res_ordination_trans
eigval <- res_trans$eigval

res_trans$df_arrows_spe$label <- gsub("^p__", "", res_trans$df_arrows_spe$Phylum)

p_rda <- ggplot() +
  geom_point(data = res_trans$df_sites, aes(x = x, y = y), color = "#FFD700", size = 4.5) +
  geom_segment(
    data = res_trans$df_arrows,
    aes(x = 0, y = 0, xend = x, yend = y),
    arrow = arrow(length = unit(0.25, "cm")),
    color = "black"
  ) +
  geom_text_repel(
    data = res_trans$df_arrows,
    aes(x = x, y = y, label = rownames(res_trans$df_arrows)),
    color = "black", size = 4
  ) +
  geom_segment(
    data = res_trans$df_arrows_spe,
    aes(x = 0, y = 0, xend = x, yend = y),
    arrow = arrow(length = unit(0.25, "cm")),
    color = "red"
  ) +
  geom_text_repel(
    data = res_trans$df_arrows_spe,
    aes(x = x, y = y, label = label),
    color = "red", size = 4
  ) +
  labs(x = eigval[1], y = eigval[2]) +
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    text = element_text(size = 14)
  )

ggsave(file.path(outdir, "rda_phylum.png"), p_rda, width = 8, height = 7, dpi = 300)

t_env$cal_mantel(use_measure = "bray", method = "spearman", p_adjust_method = "fdr")
capture.output(t_env$res_mantel, file = file.path(outdir, "microeco_mantel.txt"))

t_env$cal_cor(use_data = "Phylum", p_adjust_method = "fdr", p_adjust_type = "Env")
capture.output(t_env$res_cor$Significance, file = file.path(outdir, "microeco_cor_significance.txt"))

p_cor <- t_env$plot_cor(filter_feature = c("*", "**", "***"))
ggsave(file.path(outdir, "microeco_cor_plot.png"), p_cor, width = 7, height = 7, dpi = 300)

# =============================================================================
# Save session information
# =============================================================================

writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo_hypothesis2_16S.txt"))
