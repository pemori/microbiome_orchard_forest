# =============================================================================
# Script: 08_hypothesis_3_ITS_assembly.R
# Project: From forest to fields: The role of soil microbiome spillover in
#          agroecosystem sustainability
# Author: Pedro Mondaca
# Contact: pedromondaca@outlook.com
# Submitted to: Science Advances
#
# Description:
# This script evaluates Hypothesis 3 for the fungal community in orchard
# soils across a gradient of surrounding forest cover. It integrates
# taxonomic, ecological, and assembly-oriented analyses to examine whether
# fungal genera differ in niche breadth, abundance structure, and inferred
# assembly processes under contrasting forest-cover contexts.
#
# Specifically, the script:
#   - compares mean fungal relative abundance between orchard and forest
#     soils at the family and genus levels;
#   - summarizes orchard edaphic variation using PCA;
#   - estimates genus-level niche breadth from pH, SOM, N, and clay+silt
#     using ecospat;
#   - classifies genera into Specialist, Opportunist, and Generalist groups
#     using the lower 25%, middle 50%, and upper 25% of the geometric mean
#     niche-breadth distribution;
#   - quantifies the relative abundance of each niche-breadth category in
#     orchard samples and evaluates how these groups vary across forest-cover
#     classes;
#   - relates category-level abundance patterns to surrounding forest cover
#     using frequentist and Bayesian models;
#   - calculates environmental distance from orchard soils to block-matched
#     forest reference conditions;
#   - identifies the most abundant genera within each niche-breadth category;
#   - infers community assembly processes across orchard soils using iCAMP,
#     including phylogenetic signal assessment, bin-based process importance,
#     deterministic vs stochastic summaries, bootstrapping comparisons among
#     forest-cover classes, and complementary QPEN, neutral taxa, tNST, and
#     pNST analyses;
#   - fits a Sloan neutral model to orchard fungal communities and
#     integrates Sloan-above taxa with ASV partitioning and niche-breadth
#     categories to evaluate how neutral-model deviations relate to forest
#     spillover and ecological strategy.
#
# Required files in the working directory:
#   - results_ITS/phyloseq_ITS_asv_filtered.rds
#   - results_ITS/phyloseq_ITS_genus.rds
#   - results_ITS/mapping_ITS_ordered.csv
#   - results_ITS/otu_crops_ITS_genus.csv
#   - results_ITS/otu_forest_ITS_genus.csv
#   - results_ITS/otu_table_ITS_family_relative_abundance.csv
#   - results_ITS_hyp2/orchard_summary_with_mntd.csv
#
# Output files:
#   - written to the folder "results_ITS_hyp3"
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(phyloseq)
  library(ecospat)
  library(FactoMineR)
  library(factoextra)
  library(picante)
  library(ape)
  library(tibble)
  library(car)
  library(emmeans)
  library(performance)
  library(betareg)
})

# =============================================================================
# Paths and setup
# =============================================================================

setwd("c:\\Users\\pedro\\Research\\microbiome_orchard_forest\\Sci Adv - v2\\GitHub")

input_phyloseq_asv   <- "results_ITS/phyloseq_ITS_asv_filtered.rds"
input_phyloseq_genus <- "results_ITS/phyloseq_ITS_genus.rds"
input_mapping        <- "results_ITS/mapping_ITS_ordered.csv"
input_otu_crops      <- "results_ITS/otu_crops_ITS_genus.csv"
input_otu_forest     <- "results_ITS/otu_forest_ITS_genus.csv"
input_family_ra      <- "results_ITS/otu_table_ITS_family_relative_abundance.csv"
input_orchard_hyp2   <- "results_ITS_hyp2/orchard_summary_with_mntd.csv"
outdir               <- "results_ITS_hyp3"

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

required_files <- c(
  input_phyloseq_asv, input_phyloseq_genus, input_mapping,
  input_otu_crops, input_otu_forest, input_family_ra, input_orchard_hyp2
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required files: ", paste(missing_files, collapse = ", "))
}

# =============================================================================
# Load objects
# =============================================================================

data_phylo_filt  <- readRDS(input_phyloseq_asv)
data_phylo_genus <- readRDS(input_phyloseq_genus)

mapping <- read.csv(input_mapping, row.names = 1, check.names = FALSE)
otu_crops_g_fun  <- read.csv(input_otu_crops, row.names = 1, check.names = FALSE)
otu_forest_g_fun <- read.csv(input_otu_forest, row.names = 1, check.names = FALSE)
data_fam_f_RA    <- read.csv(input_family_ra, row.names = 1, check.names = FALSE)

orchard_summary_df <- read.csv(
  input_orchard_hyp2,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!"SampleID" %in% colnames(orchard_summary_df)) {
  stop("The file ", input_orchard_hyp2, " must contain a column named 'SampleID'.")
}

stopifnot(!anyDuplicated(orchard_summary_df$SampleID))
rownames(orchard_summary_df) <- orchard_summary_df$SampleID

# align mapping to phyloseq sample names
mapping <- mapping[sample_names(data_phylo_filt), , drop = FALSE]

if (!identical(rownames(mapping), sample_names(data_phylo_filt))) {
  stop("Mapping file row names do not match phyloseq ASV sample names.")
}

sample_data(data_phylo_filt)  <- sample_data(mapping)
sample_data(data_phylo_genus) <- sample_data(mapping)

mapping_df <- as(sample_data(data_phylo_filt), "data.frame")
mapping_df$SampleID <- rownames(mapping_df)

phylo_crops_fun    <- subset_samples(data_phylo_filt, cat != "F")
phylo_forest_fun   <- subset_samples(data_phylo_filt, cat == "F")
phylo_crops_g_fun  <- subset_samples(data_phylo_genus, cat != "F")
phylo_forest_g_fun <- subset_samples(data_phylo_genus, cat == "F")

mapping_orchard <- as(sample_data(phylo_crops_fun), "data.frame")
mapping_orchard$SampleID <- rownames(mapping_orchard)

mapping_forest <- as(sample_data(phylo_forest_fun), "data.frame")
mapping_forest$SampleID <- rownames(mapping_forest)

# =============================================================================
# Hypothesis 3
# Compared microbial assembly process, niche breadth distributions, and
# abundance structure of fungal genera in orchard soils with contrasting
# surrounding forest cover
# =============================================================================

# =============================================================================
# Comparison of mean relative abundance between forest and orchard
# =============================================================================

forest_fam_f_RA <- data_fam_f_RA[grepl("N", rownames(data_fam_f_RA)), , drop = FALSE]
crops_fam_f_RA  <- data_fam_f_RA[!grepl("N", rownames(data_fam_f_RA)), , drop = FALSE]

promedio_abundancia_f_forest <- colMeans(forest_fam_f_RA, na.rm = TRUE)
promedio_abundancia_f_crops  <- colMeans(crops_fam_f_RA, na.rm = TRUE)

shared_families <- intersect(names(promedio_abundancia_f_forest), names(promedio_abundancia_f_crops))

df_comparacion_family <- data.frame(
  Family = shared_families,
  Mean_Forest = promedio_abundancia_f_forest[shared_families],
  Mean_Orchard = promedio_abundancia_f_crops[shared_families]
)

p_family_forest_crop <- ggplot(df_comparacion_family, aes(x = Mean_Forest, y = Mean_Orchard)) +
  geom_point(color = "blue") +
  labs(
    x = "Mean relative abundance in forest",
    y = "Mean relative abundance in orchard"
  ) +
  theme_minimal()

ggsave(
  file.path(outdir, "family_mean_abundance_forest_vs_orchard.png"),
  p_family_forest_crop, width = 6, height = 5, dpi = 300
)

# Genus level
otu_crops_RA  <- otu_crops_g_fun / rowSums(otu_crops_g_fun)
otu_forest_RA <- otu_forest_g_fun / rowSums(otu_forest_g_fun)

otu_crops_RA[is.na(otu_crops_RA)]   <- 0
otu_forest_RA[is.na(otu_forest_RA)] <- 0

shared_genera <- intersect(colnames(otu_crops_RA), colnames(otu_forest_RA))

promedio_abundancia_crops  <- colMeans(otu_crops_RA[, shared_genera, drop = FALSE], na.rm = TRUE)
promedio_abundancia_forest <- colMeans(otu_forest_RA[, shared_genera, drop = FALSE], na.rm = TRUE)

df_comparacion_genus <- data.frame(
  Genus = shared_genera,
  Mean_Orchard = promedio_abundancia_crops,
  Mean_Forest = promedio_abundancia_forest
)

p_genus_forest_crop <- ggplot(df_comparacion_genus, aes(x = Mean_Forest, y = Mean_Orchard)) +
  geom_point(color = "blue", alpha = 0.7) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  labs(
    x = "Mean relative abundance in forest",
    y = "Mean relative abundance in orchard"
  ) +
  theme_minimal()

ggsave(
  file.path(outdir, "genus_mean_abundance_forest_vs_orchard.png"),
  p_genus_forest_crop, width = 6, height = 5, dpi = 300
)

write.csv(df_comparacion_family, file.path(outdir, "family_mean_abundance_forest_vs_orchard.csv"), row.names = FALSE)
write.csv(df_comparacion_genus,  file.path(outdir, "genus_mean_abundance_forest_vs_orchard.csv"),  row.names = FALSE)

# =============================================================================
# PCA of environmental variables in orchard soils
# =============================================================================

env_pca_df <- as(sample_data(phylo_crops_fun), "data.frame")[,
                                                             c("pH", "SOM", "N", "C", "CN", "P", "claysilt", "SOMTXT"),
                                                             drop = FALSE
]

env_pca_df <- as.data.frame(env_pca_df)
env_pca_df[] <- lapply(env_pca_df, as.numeric)
env_pca_df <- env_pca_df[, colSums(!is.na(env_pca_df)) > 0, drop = FALSE]
env_pca_df_complete <- env_pca_df[complete.cases(env_pca_df), , drop = FALSE]

pca_fun <- PCA(env_pca_df_complete, graph = FALSE)

capture.output(pca_fun$eig, file = file.path(outdir, "pca_env_eigenvalues.txt"))
write.csv(pca_fun$var$coord, file.path(outdir, "pca_env_variable_coordinates.csv"))

p_pca_vars <- fviz_pca_var(
  pca_fun,
  col.var = "black",
  repel = TRUE
) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black", face = "bold")
  )

ggsave(
  file.path(outdir, "pca_env_variables.png"),
  p_pca_vars, width = 6, height = 5, dpi = 300
)

# =============================================================================
# Niche breadth analysis
# =============================================================================

tab <- as(sample_data(phylo_crops_fun), "data.frame")
ID <- rownames(tab)

env_fun2 <- data.frame(
  ID = ID,
  pH = as.numeric(tab$pH),
  SOM = as.numeric(tab$SOM),
  N = as.numeric(tab$N),
  `Clay+Silt` = as.numeric(tab$claysilt),
  Block_crop = tab$Block,
  disp_fun_250 = as.numeric(tab$P_250),
  disp_fun_1000 = as.numeric(tab$P_1000),
  disp_fun_2500 = as.numeric(tab$P_2500),
  check.names = FALSE
)
rownames(env_fun2) <- env_fun2$ID

stopifnot(all(env_fun2$ID %in% rownames(otu_crops_RA)))
otu_crops_RA <- otu_crops_RA[env_fun2$ID, , drop = FALSE]

niche_fun_crops <- cbind(
  data.frame(ID = env_fun2$ID, check.names = FALSE),
  otu_crops_RA,
  env_fun2[, c("pH", "SOM", "N", "Clay+Silt", "disp_fun_2500", "Block_crop"), drop = FALSE]
)

genus_cols <- colnames(otu_crops_RA)

niche_fun_crops[, c(genus_cols, "pH", "SOM", "N", "Clay+Silt", "disp_fun_2500")] <-
  lapply(
    niche_fun_crops[, c(genus_cols, "pH", "SOM", "N", "Clay+Silt", "disp_fun_2500"), drop = FALSE],
    as.numeric
  )

colfreq_idx <- which(colnames(niche_fun_crops) %in% genus_cols)
colvar_idx  <- which(colnames(niche_fun_crops) %in% c("pH", "SOM", "N", "Clay+Silt"))

POSNB <- ecospat.nichePOSNB(
  niche_fun_crops,
  colvar = colvar_idx,
  colfreq = colfreq_idx
)

avg_NB <- ecospat.nicheNBmean(POSNB, w = c(1, 1, 1, 1))
write.csv(as.data.frame(avg_NB), file.path(outdir, "average_niche_breadth.csv"))

POSNB <- as.data.frame(POSNB)
POSNB$Genus <- rownames(POSNB)

nb_genus_df <- POSNB %>%
  transmute(
    Genus = Genus,
    pH_nb = pH_nb,
    SOM_nb = SOM_nb,
    N_nb = N_nb,
    claysilt_nb = `Clay+Silt_nb`
  ) %>%
  filter(
    !is.na(pH_nb) & pH_nb > 0,
    !is.na(SOM_nb) & SOM_nb > 0,
    !is.na(N_nb) & N_nb > 0,
    !is.na(claysilt_nb) & claysilt_nb > 0
  ) %>%
  mutate(
    geometric_mean_nb = exp((log(pH_nb) + log(SOM_nb) + log(N_nb) + log(claysilt_nb)) / 4)
  )

write.csv(nb_genus_df, file.path(outdir, "genus_niche_breadth_table.csv"), row.names = FALSE)

p_nb_hist <- ggplot(nb_genus_df, aes(x = geometric_mean_nb)) +
  geom_histogram(bins = 30, color = "black", fill = "grey70", alpha = 0.8) +
  labs(x = "Niche breadth", y = "Count") +
  theme_minimal(base_size = 14)

ggsave(
  file.path(outdir, "niche_breadth_histogram_raw.png"),
  p_nb_hist, width = 7.64, height = 5.07, dpi = 300
)

# =============================================================================
# Quantile-based classification: Specialist / Opportunist / Generalist
# =============================================================================

q25 <- quantile(nb_genus_df$geometric_mean_nb, 0.25, na.rm = TRUE)
q75 <- quantile(nb_genus_df$geometric_mean_nb, 0.75, na.rm = TRUE)

nb_genus_df <- nb_genus_df %>%
  mutate(
    NB_Category = case_when(
      geometric_mean_nb <= q25 ~ "Specialist",
      geometric_mean_nb >= q75 ~ "Generalist",
      TRUE ~ "Opportunist"
    ),
    NB_Category = factor(NB_Category, levels = c("Generalist", "Opportunist", "Specialist"))
  )

write.csv(
  data.frame(q25 = q25, q75 = q75),
  file.path(outdir, "niche_breadth_quantile_thresholds.csv"),
  row.names = FALSE
)

write.csv(
  nb_genus_df,
  file.path(outdir, "genus_niche_breadth_with_categories.csv"),
  row.names = FALSE
)

capture.output(table(nb_genus_df$NB_Category), file = file.path(outdir, "niche_breadth_category_counts.txt"))

p_nb_hist_cat <- ggplot(nb_genus_df, aes(x = geometric_mean_nb, fill = NB_Category)) +
  geom_histogram(bins = 30, color = "black", alpha = 0.7) +
  labs(x = "Niche breadth", y = "Count") +
  scale_fill_manual(values = c(
    "Specialist" = "red",
    "Opportunist" = "blue",
    "Generalist" = "green"
  )) +
  theme_minimal(base_size = 14)

ggsave(
  file.path(outdir, "niche_breadth_histogram_categories.png"),
  p_nb_hist_cat, width = 7.64, height = 5.07, dpi = 300
)

# =============================================================================
# Long table for sample-level relative abundance summaries
# =============================================================================

fungi_pal_graph <- data.frame(
  ID = rownames(otu_crops_RA),
  otu_crops_RA,
  check.names = FALSE
) %>%
  left_join(env_fun2, by = "ID")

fungi_pal_graph_long <- fungi_pal_graph %>%
  pivot_longer(
    cols = all_of(genus_cols),
    names_to = "Genus",
    values_to = "Relative_Abundance"
  ) %>%
  left_join(
    nb_genus_df %>% dplyr::select(Genus, pH_nb, SOM_nb, N_nb, claysilt_nb, geometric_mean_nb, NB_Category),
    by = "Genus"
  ) %>%
  mutate(
    Forest_coverage = case_when(
      disp_fun_2500 <= 0.3584376 ~ "LFC",
      disp_fun_2500 <= 0.5224324 ~ "MFC",
      disp_fun_2500 > 0.5224324  ~ "HFC"
    ),
    Forest_coverage = factor(Forest_coverage, levels = c("LFC", "MFC", "HFC")),
    NB_Category = factor(NB_Category, levels = c("Generalist", "Opportunist", "Specialist"))
  )

write.csv(
  fungi_pal_graph_long,
  file.path(outdir, "fungal_niche_breadth_long_table.csv"),
  row.names = FALSE
)

fungi_check_sum <- fungi_pal_graph_long %>%
  group_by(ID) %>%
  summarise(total_RA = sum(Relative_Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_RA))

write.csv(fungi_check_sum, file.path(outdir, "check_sample_relative_abundance_sums.csv"), row.names = FALSE)

sample_count_fc <- fungi_pal_graph_long %>%
  distinct(ID, Forest_coverage) %>%
  count(Forest_coverage)

write.csv(sample_count_fc, file.path(outdir, "sample_count_by_forest_coverage.csv"), row.names = FALSE)

genera_count_nb <- nb_genus_df %>%
  count(NB_Category)

write.csv(genera_count_nb, file.path(outdir, "genera_count_by_nb_category.csv"), row.names = FALSE)

# =============================================================================
# Per-sample summed relative abundance by niche breadth category
# =============================================================================

fungi_sum_per_sample <- fungi_pal_graph_long %>%
  filter(!is.na(NB_Category)) %>%
  group_by(ID, Forest_coverage, NB_Category) %>%
  summarise(
    Sum_Relative_Abundance = sum(Relative_Abundance, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  fungi_sum_per_sample,
  file.path(outdir, "sample_relative_abundance_by_nb_category_long.csv"),
  row.names = FALSE
)

fungi_sum_wide <- fungi_sum_per_sample %>%
  dplyr::select(ID, NB_Category, Sum_Relative_Abundance) %>%
  pivot_wider(
    names_from = NB_Category,
    values_from = Sum_Relative_Abundance,
    values_fill = 0
  )

for (nm in c("Generalist", "Opportunist", "Specialist")) {
  if (!nm %in% colnames(fungi_sum_wide)) fungi_sum_wide[[nm]] <- 0
}

fungi_sum_wide <- fungi_sum_wide %>%
  mutate(total = Generalist + Opportunist + Specialist)

write.csv(fungi_sum_wide, file.path(outdir, "sample_relative_abundance_by_nb_category_wide.csv"), row.names = FALSE)

# =============================================================================
# Summary plots from sample-level sums
# =============================================================================

fungi_plot_summary <- fungi_sum_per_sample %>%
  group_by(NB_Category, Forest_coverage) %>%
  summarise(
    Relative_Abundance = mean(Sum_Relative_Abundance, na.rm = TRUE),
    SE = sd(Sum_Relative_Abundance, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    Forest_coverage = factor(Forest_coverage, levels = c("LFC", "MFC", "HFC")),
    NB_Category = factor(NB_Category, levels = c("Generalist", "Opportunist", "Specialist"))
  )

write.csv(fungi_plot_summary, file.path(outdir, "summary_mean_se_by_nb_and_forestcover.csv"), row.names = FALSE)

colors_dispersion <- c(
  "LFC" = "#1f78b4",
  "MFC" = "gold",
  "HFC" = "#33a02c"
)

p_nb_bar_facet <- ggplot(
  fungi_plot_summary,
  aes(x = Forest_coverage, y = Relative_Abundance, fill = Forest_coverage)
) +
  geom_bar(
    stat = "identity",
    width = 0.88,
    color = "black",
    linewidth = 0.4,
    alpha = 0.8
  ) +
  geom_errorbar(
    aes(ymin = Relative_Abundance - SE, ymax = Relative_Abundance + SE),
    width = 0.18,
    linewidth = 0.4,
    color = "black"
  ) +
  facet_wrap(~ NB_Category, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = colors_dispersion) +
  labs(x = "", y = "Relative abundance (%)") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )

ggsave(
  file.path(outdir, "niche_breadth_barplot_facet.png"),
  p_nb_bar_facet, width = 8.5, height = 3.8, dpi = 300
)

Generalist_summary  <- fungi_plot_summary %>% filter(NB_Category == "Generalist")
Opportunist_summary <- fungi_plot_summary %>% filter(NB_Category == "Opportunist")
Specialist_summary  <- fungi_plot_summary %>% filter(NB_Category == "Specialist")

p_generalist <- ggplot(Generalist_summary, aes(x = Forest_coverage, y = Relative_Abundance, fill = Forest_coverage)) +
  geom_bar(stat = "identity", width = 0.9, color = "black", alpha = 0.8) +
  geom_errorbar(aes(ymin = Relative_Abundance - SE, ymax = Relative_Abundance + SE), width = 0.2, color = "black") +
  scale_fill_manual(values = colors_dispersion) +
  labs(x = "", y = "Relative abundance (%)") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )

p_opportunist <- ggplot(Opportunist_summary, aes(x = Forest_coverage, y = Relative_Abundance, fill = Forest_coverage)) +
  geom_bar(stat = "identity", width = 0.9, color = "black", alpha = 0.8) +
  geom_errorbar(aes(ymin = Relative_Abundance - SE, ymax = Relative_Abundance + SE), width = 0.2, color = "black") +
  scale_fill_manual(values = colors_dispersion) +
  labs(x = "", y = "Mean relative abundance (%)") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )

p_specialist <- ggplot(Specialist_summary, aes(x = Forest_coverage, y = Relative_Abundance, fill = Forest_coverage)) +
  geom_bar(stat = "identity", width = 0.9, color = "black", alpha = 0.8) +
  geom_errorbar(aes(ymin = Relative_Abundance - SE, ymax = Relative_Abundance + SE), width = 0.2, color = "black") +
  scale_fill_manual(values = colors_dispersion) +
  labs(x = "", y = "Relative abundance (%)") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )

ggsave(file.path(outdir, "generalist_fun.png"),  p_generalist,  width = 2.5, height = 4.5, dpi = 300)
ggsave(file.path(outdir, "opportunist_fun.png"), p_opportunist, width = 2.5, height = 4.5, dpi = 300)
ggsave(file.path(outdir, "specialist_fun.png"),  p_specialist,  width = 2.5, height = 4.5, dpi = 300)

# =============================================================================
# Join niche-breadth category abundances to orchard_summary_df
# =============================================================================

RA_by_cat <- fungi_sum_per_sample %>%
  dplyr::select(ID, NB_Category, Sum_Relative_Abundance) %>%
  pivot_wider(
    names_from = NB_Category,
    values_from = Sum_Relative_Abundance,
    values_fill = 0
  )

for (nm in c("Generalist", "Opportunist", "Specialist")) {
  if (!nm %in% colnames(RA_by_cat)) RA_by_cat[[nm]] <- 0
}

RA_by_cat <- RA_by_cat %>%
  mutate(
    total = Generalist + Opportunist + Specialist,
    Generalist = ifelse(total > 0, Generalist / total, NA_real_),
    Opportunist = ifelse(total > 0, Opportunist / total, NA_real_),
    Specialist = ifelse(total > 0, Specialist / total, NA_real_)
  )

n_beta <- nrow(RA_by_cat)

RA_by_cat <- RA_by_cat %>%
  mutate(
    Generalist_beta  = (Generalist  * (n_beta - 1) + 0.5) / n_beta,
    Opportunist_beta = (Opportunist * (n_beta - 1) + 0.5) / n_beta,
    Specialist_beta  = (Specialist  * (n_beta - 1) + 0.5) / n_beta
  )

write.csv(RA_by_cat, file.path(outdir, "relative_abundance_by_nb_category_with_beta.csv"), row.names = FALSE)

stopifnot(all(RA_by_cat$ID %in% orchard_summary_df$SampleID))

orchard_summary_df <- orchard_summary_df %>%
  dplyr::select(-any_of(c(
    "Generalist", "Opportunist", "Specialist",
    "Generalist_beta", "Opportunist_beta", "Specialist_beta",
    "FC"
  ))) %>%
  left_join(
    RA_by_cat %>%
      dplyr::select(ID, Generalist, Opportunist, Specialist,
                    Generalist_beta, Opportunist_beta, Specialist_beta) %>%
      rename(SampleID = ID),
    by = "SampleID"
  ) %>%
  mutate(
    FC = case_when(
      P_2500 <= 0.3584376 ~ "LFC",
      P_2500 <= 0.5224324 ~ "MFC",
      P_2500 > 0.5224324  ~ "HFC",
      TRUE ~ NA_character_
    ),
    FC = factor(FC, levels = c("LFC", "MFC", "HFC"))
  )

rownames(orchard_summary_df) <- orchard_summary_df$SampleID

write.csv(
  orchard_summary_df,
  file.path(outdir, "orchard_summary_with_nb_categories.csv"),
  row.names = FALSE
)

# =============================================================================
# Environmental distance to forest reference within block
# =============================================================================

forest_df <- as(sample_data(phylo_forest_g_fun), "data.frame")
forest_df$SampleID <- rownames(forest_df)

orchard_df <- as(sample_data(phylo_crops_g_fun), "data.frame")
orchard_df$SampleID <- rownames(orchard_df)

global_sds <- forest_df %>%
  summarise(
    g_sd_pH = sd(pH, na.rm = TRUE),
    g_sd_SOM = sd(SOM, na.rm = TRUE),
    g_sd_N = sd(N, na.rm = TRUE),
    g_sd_claysilt = sd(claysilt, na.rm = TRUE)
  )

forest_ref <- forest_df %>%
  group_by(Block) %>%
  summarise(
    ref_pH = mean(pH, na.rm = TRUE),
    ref_SOM = mean(SOM, na.rm = TRUE),
    ref_N = mean(N, na.rm = TRUE),
    ref_claysilt = mean(claysilt, na.rm = TRUE),
    sd_pH = sd(pH, na.rm = TRUE),
    sd_SOM = sd(SOM, na.rm = TRUE),
    sd_N = sd(N, na.rm = TRUE),
    sd_claysilt = sd(claysilt, na.rm = TRUE),
    n_forest = dplyr::n(),
    .groups = "drop"
  ) %>%
  mutate(
    sd_pH = if_else(is.na(sd_pH) | sd_pH == 0, global_sds$g_sd_pH, sd_pH),
    sd_SOM = if_else(is.na(sd_SOM) | sd_SOM == 0, global_sds$g_sd_SOM, sd_SOM),
    sd_N = if_else(is.na(sd_N) | sd_N == 0, global_sds$g_sd_N, sd_N),
    sd_claysilt = if_else(is.na(sd_claysilt) | sd_claysilt == 0, global_sds$g_sd_claysilt, sd_claysilt)
  )

orchard_env_dist <- orchard_df %>%
  left_join(forest_ref, by = "Block") %>%
  mutate(
    z_pH = (pH - ref_pH) / sd_pH,
    z_SOM = (SOM - ref_SOM) / sd_SOM,
    z_N = (N - ref_N) / sd_N,
    z_claysilt = (claysilt - ref_claysilt) / sd_claysilt,
    env_dist = sqrt(z_pH^2 + z_SOM^2 + z_N^2 + z_claysilt^2),
    env_dist_L1 = abs(z_pH) + abs(z_SOM) + abs(z_N) + abs(z_claysilt)
  )

write.csv(
  orchard_env_dist,
  file.path(outdir, "orchard_environmental_distance_to_forest.csv"),
  row.names = FALSE
)

nb_fun_df <- orchard_summary_df %>%
  left_join(
    orchard_env_dist %>% select(SampleID, env_dist, env_dist_L1),
    by = "SampleID"
  )

write.csv(nb_fun_df, file.path(outdir, "nb_fungal_model_dataframe.csv"), row.names = FALSE)

# =============================================================================
# Top genera by niche breadth category
# =============================================================================

top25_by_nb <- fungi_pal_graph_long %>%
  filter(!is.na(NB_Category)) %>%
  group_by(NB_Category, Genus) %>%
  summarise(total_RA = sum(Relative_Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(NB_Category, desc(total_RA)) %>%
  group_by(NB_Category) %>%
  slice_head(n = 25)

write.csv(top25_by_nb, file.path(outdir, "top25_genera_by_nb_category.csv"), row.names = FALSE)

# =============================================================================
# Frequentist statistical analysis by forest cover category
# lm - inverse.gaussian - lognorm - betaregression
# =============================================================================

orchard_summary_df$FC <- factor(orchard_summary_df$FC, levels = c("LFC", "MFC", "HFC"))

# Generalist
mod_lm <- lm(Generalist ~ FC, data = orchard_summary_df)
mod_invgauss <- glm(Generalist ~ FC, data = orchard_summary_df, family = inverse.gaussian(link = "log"))
mod_lognorm <- lm(log(Generalist) ~ FC, data = orchard_summary_df)
mod_beta <- betareg(
  Generalist_beta ~ FC,
  data = orchard_summary_df,
  link = "logit"
)

compare_performance(mod_lm, mod_invgauss, mod_lognorm, mod_beta,
                    rank = TRUE, verbose = TRUE)
check_model(mod_invgauss)
emmeans(mod_invgauss, pairwise ~ FC, type = "response")

# Opportunist
mod_lm <- lm(Opportunist ~ FC, data = orchard_summary_df)
mod_invgauss <- glm(Opportunist ~ FC, data = orchard_summary_df, family = inverse.gaussian(link = "log"))
mod_lognorm <- lm(log(Opportunist) ~ FC, data = orchard_summary_df)
mod_beta <- betareg(
  Opportunist_beta ~ FC,
  data = orchard_summary_df,
  link = "logit"
)

compare_performance(mod_lm, mod_invgauss, mod_lognorm, mod_beta,
                    rank = TRUE, verbose = TRUE)
check_model(mod_beta)
emmeans(mod_beta, pairwise ~ FC, type = "response")

# Specialist
mod_lm <- lm(Specialist ~ FC, data = orchard_summary_df)
mod_invgauss <- glm(Specialist ~ FC, data = orchard_summary_df, family = inverse.gaussian(link = "log"))
mod_lognorm <- lm(log(Specialist) ~ FC, data = orchard_summary_df)
mod_beta <- betareg(
  Specialist_beta ~ FC,
  data = orchard_summary_df,
  link = "logit"
)

compare_performance(mod_lm, mod_invgauss, mod_lognorm, mod_beta,
                    rank = TRUE, verbose = TRUE)
check_model(mod_beta)
emmeans(mod_beta, pairwise ~ FC, type = "response")

# -----------------------------------------------------------------------------
# Bayesian analysis
# -----------------------------------------------------------------------------

orchard_env_dist2 <- orchard_env_dist %>%
  dplyr::select(SampleID, env_dist, env_dist_L1)

RA_by_cat2 <- RA_by_cat %>%
  dplyr::select(
    ID,
    Generalist, Opportunist, Specialist,
    Generalist_beta, Opportunist_beta, Specialist_beta
  ) %>%
  dplyr::rename(SampleID = ID)

pattern <- "^(env_dist(_L1)?|Generalist|Opportunist|Specialist|Generalist_beta|Opportunist_beta|Specialist_beta)(\\..*)?$"

nb_df <- orchard_summary_df %>%
  dplyr::select(-dplyr::any_of(grep(pattern, names(.), value = TRUE))) %>%
  dplyr::left_join(orchard_env_dist2, by = "SampleID") %>%
  dplyr::left_join(RA_by_cat2, by = "SampleID")

stopifnot("SampleID" %in% colnames(nb_df))
stopifnot(!anyDuplicated(nb_df$SampleID))
rownames(nb_df) <- nb_df$SampleID

write.csv(nb_df, file.path(outdir, "nb_fungal_model_dataframe.csv"), row.names = FALSE)

str(nb_df)

# -----------------------------------------------------------------------------
# Specialist models
# -----------------------------------------------------------------------------

spe_df <- nb_df %>%
  dplyr::select(SampleID, Block, P_2500, env_dist, pH, SOM, Specialist_beta) %>%
  filter(complete.cases(.))

model_Spe0 <- brm(
  Specialist_beta ~ 1 + (1 | Block),
  data = spe_df,
  family = Beta(),
  chains = 4,
  cores = min(12, parallel::detectCores()),
  iter = 4000,
  control = list(adapt_delta = 0.99)
)

model_Spe1 <- brm(
  Specialist_beta ~ t2(P_2500) + (1 | Block),
  data = spe_df,
  family = Beta(),
  chains = 4,
  cores = min(12, parallel::detectCores()),
  iter = 4000,
  control = list(adapt_delta = 0.99)
)

model_Spe2 <- brm(
  Specialist_beta ~ t2(P_2500, env_dist) + (1 | Block),
  data = spe_df,
  family = Beta(),
  chains = 4,
  cores = min(12, parallel::detectCores()),
  iter = 4000,
  control = list(adapt_delta = 0.99)
)

model_Spe4 <- brm(
  Specialist_beta ~ t2(P_2500, pH) + (1 | Block),
  data = spe_df,
  family = Beta(),
  chains = 4,
  cores = min(12, parallel::detectCores()),
  iter = 4000,
  control = list(adapt_delta = 0.99)
)

model_Spe5 <- brm(
  Specialist_beta ~ t2(P_2500, SOM) + (1 | Block),
  data = spe_df,
  family = Beta(),
  chains = 4,
  cores = min(12, parallel::detectCores()),
  iter = 4000,
  control = list(adapt_delta = 0.99)
)

saveRDS(model_Spe0, file.path(outdir, "model_Spe0.rds"))
saveRDS(model_Spe1, file.path(outdir, "model_Spe1.rds"))
saveRDS(model_Spe2, file.path(outdir, "model_Spe2.rds"))
saveRDS(model_Spe4, file.path(outdir, "model_Spe4.rds"))
saveRDS(model_Spe5, file.path(outdir, "model_Spe5.rds"))

capture.output(summary(model_Spe0), file = file.path(outdir, "summary_model_Spe0.txt"))
capture.output(summary(model_Spe1), file = file.path(outdir, "summary_model_Spe1.txt"))
capture.output(summary(model_Spe2), file = file.path(outdir, "summary_model_Spe2.txt"))
capture.output(summary(model_Spe4), file = file.path(outdir, "summary_model_Spe4.txt"))
capture.output(summary(model_Spe5), file = file.path(outdir, "summary_model_Spe5.txt"))

capture.output(bayes_R2(model_Spe0), file = file.path(outdir, "bayesR2_model_Spe0.txt"))
capture.output(bayes_R2(model_Spe1), file = file.path(outdir, "bayesR2_model_Spe1.txt"))
capture.output(bayes_R2(model_Spe2), file = file.path(outdir, "bayesR2_model_Spe2.txt"))
capture.output(bayes_R2(model_Spe4), file = file.path(outdir, "bayesR2_model_Spe4.txt"))
capture.output(bayes_R2(model_Spe5), file = file.path(outdir, "bayesR2_model_Spe5.txt"))

png(file.path(outdir, "ppcheck_model_Spe1.png"), width = 1800, height = 1400, res = 220)
print(pp_check(model_Spe1))
dev.off()

png(file.path(outdir, "ppcheck_model_Spe2.png"), width = 1800, height = 1400, res = 220)
print(pp_check(model_Spe2))
dev.off()

png(file.path(outdir, "marginal_effects_model_Spe4_surface.png"), width = 1800, height = 1400, res = 220)
print(plot(marginal_effects(model_Spe4, surface = TRUE), points = TRUE))
dev.off()

# -----------------------------------------------------------------------------
# Opportunist models
# -----------------------------------------------------------------------------

opp_df <- nb_df %>%
  dplyr::select(SampleID, Block, P_2500, env_dist, pH, SOM, Opportunist_beta) %>%
  filter(complete.cases(.))

model_Opp0 <- brm(
  Opportunist_beta ~ 1 + (1 | Block),
  data = opp_df,
  family = Beta(),
  chains = 4,
  cores = min(12, parallel::detectCores()),
  iter = 4000,
  control = list(adapt_delta = 0.99)
)

model_Opp1 <- brm(
  Opportunist_beta ~ t2(P_2500) + (1 | Block),
  data = opp_df,
  family = Beta(),
  chains = 4,
  cores = min(12, parallel::detectCores()),
  iter = 4000,
  control = list(adapt_delta = 0.99)
)

model_Opp2 <- brm(
  Opportunist_beta ~ t2(P_2500, env_dist) + (1 | Block),
  data = opp_df,
  family = Beta(),
  chains = 4,
  cores = min(12, parallel::detectCores()),
  iter = 4000,
  control = list(adapt_delta = 0.99)
)

model_Opp4 <- brm(
  Opportunist_beta ~ t2(P_2500, pH) + (1 | Block),
  data = opp_df,
  family = Beta(),
  chains = 4,
  cores = min(12, parallel::detectCores()),
  iter = 4000,
  control = list(adapt_delta = 0.99)
)

model_Opp5 <- brm(
  Opportunist_beta ~ t2(P_2500, SOM) + (1 | Block),
  data = opp_df,
  family = Beta(),
  chains = 4,
  cores = min(12, parallel::detectCores()),
  iter = 4000,
  control = list(adapt_delta = 0.99)
)

saveRDS(model_Opp0, file.path(outdir, "model_Opp0.rds"))
saveRDS(model_Opp1, file.path(outdir, "model_Opp1.rds"))
saveRDS(model_Opp2, file.path(outdir, "model_Opp2.rds"))
saveRDS(model_Opp4, file.path(outdir, "model_Opp4.rds"))
saveRDS(model_Opp5, file.path(outdir, "model_Opp5.rds"))

capture.output(summary(model_Opp0), file = file.path(outdir, "summary_model_Opp0.txt"))
capture.output(summary(model_Opp1), file = file.path(outdir, "summary_model_Opp1.txt"))
capture.output(summary(model_Opp2), file = file.path(outdir, "summary_model_Opp2.txt"))
capture.output(summary(model_Opp4), file = file.path(outdir, "summary_model_Opp4.txt"))
capture.output(summary(model_Opp5), file = file.path(outdir, "summary_model_Opp5.txt"))

capture.output(bayes_R2(model_Opp0), file = file.path(outdir, "bayesR2_model_Opp0.txt"))
capture.output(bayes_R2(model_Opp1), file = file.path(outdir, "bayesR2_model_Opp1.txt"))
capture.output(bayes_R2(model_Opp2), file = file.path(outdir, "bayesR2_model_Opp2.txt"))
capture.output(bayes_R2(model_Opp4), file = file.path(outdir, "bayesR2_model_Opp4.txt"))
capture.output(bayes_R2(model_Opp5), file = file.path(outdir, "bayesR2_model_Opp5.txt"))

png(file.path(outdir, "ppcheck_model_Opp1.png"), width = 1800, height = 1400, res = 220)
print(pp_check(model_Opp1))
dev.off()

png(file.path(outdir, "ppcheck_model_Opp2.png"), width = 1800, height = 1400, res = 220)
print(pp_check(model_Opp2))
dev.off()

png(file.path(outdir, "marginal_effects_model_Opp4_surface.png"), width = 1800, height = 1400, res = 220)
print(plot(marginal_effects(model_Opp4, surface = TRUE), points = TRUE))
dev.off()

# -----------------------------------------------------------------------------
# Generalist models
# -----------------------------------------------------------------------------

gen_df <- nb_df %>%
  dplyr::select(SampleID, Block, P_2500, env_dist, pH, SOM, Generalist_beta) %>%
  filter(complete.cases(.))

model_Gen0 <- brm(
  Generalist_beta ~ 1 + (1 | Block),
  data = gen_df,
  family = Beta(),
  chains = 4,
  cores = min(12, parallel::detectCores()),
  iter = 4000,
  control = list(adapt_delta = 0.99)
)

model_Gen1 <- brm(
  Generalist_beta ~ t2(P_2500) + (1 | Block),
  data = gen_df,
  family = Beta(),
  chains = 4,
  cores = min(12, parallel::detectCores()),
  iter = 4000,
  control = list(adapt_delta = 0.99)
)

model_Gen2 <- brm(
  Generalist_beta ~ t2(P_2500, env_dist) + (1 | Block),
  data = gen_df,
  family = Beta(),
  chains = 4,
  cores = min(12, parallel::detectCores()),
  iter = 4000,
  control = list(adapt_delta = 0.99)
)

model_Gen4 <- brm(
  Generalist_beta ~ t2(P_2500, pH) + (1 | Block),
  data = gen_df,
  family = Beta(),
  chains = 4,
  cores = min(12, parallel::detectCores()),
  iter = 4000,
  control = list(adapt_delta = 0.99)
)

model_Gen5 <- brm(
  Generalist_beta ~ t2(P_2500, SOM) + (1 | Block),
  data = gen_df,
  family = Beta(),
  chains = 4,
  cores = min(12, parallel::detectCores()),
  iter = 4000,
  control = list(adapt_delta = 0.99)
)

saveRDS(model_Gen0, file.path(outdir, "model_Gen0.rds"))
saveRDS(model_Gen1, file.path(outdir, "model_Gen1.rds"))
saveRDS(model_Gen2, file.path(outdir, "model_Gen2.rds"))
saveRDS(model_Gen4, file.path(outdir, "model_Gen4.rds"))
saveRDS(model_Gen5, file.path(outdir, "model_Gen5.rds"))

capture.output(summary(model_Gen0), file = file.path(outdir, "summary_model_Gen0.txt"))
capture.output(summary(model_Gen1), file = file.path(outdir, "summary_model_Gen1.txt"))
capture.output(summary(model_Gen2), file = file.path(outdir, "summary_model_Gen2.txt"))
capture.output(summary(model_Gen4), file = file.path(outdir, "summary_model_Gen4.txt"))
capture.output(summary(model_Gen5), file = file.path(outdir, "summary_model_Gen5.txt"))

capture.output(bayes_R2(model_Gen0), file = file.path(outdir, "bayesR2_model_Gen0.txt"))
capture.output(bayes_R2(model_Gen1), file = file.path(outdir, "bayesR2_model_Gen1.txt"))
capture.output(bayes_R2(model_Gen2), file = file.path(outdir, "bayesR2_model_Gen2.txt"))
capture.output(bayes_R2(model_Gen4), file = file.path(outdir, "bayesR2_model_Gen4.txt"))
capture.output(bayes_R2(model_Gen5), file = file.path(outdir, "bayesR2_model_Gen5.txt"))

png(file.path(outdir, "ppcheck_model_Gen1.png"), width = 1800, height = 1400, res = 220)
print(pp_check(model_Gen1))
dev.off()

png(file.path(outdir, "ppcheck_model_Gen2.png"), width = 1800, height = 1400, res = 220)
print(pp_check(model_Gen2))
dev.off()

png(file.path(outdir, "marginal_effects_model_Gen4_surface.png"), width = 1800, height = 1400, res = 220)
print(plot(marginal_effects(model_Gen4, surface = TRUE), points = TRUE))
dev.off()

# -----------------------------------------------------------------------------
# LOO comparison
# -----------------------------------------------------------------------------

loo_Gen0 <- loo(model_Gen0)
loo_Gen1 <- loo(model_Gen1)
loo_Gen2 <- loo(model_Gen2)
loo_Gen4 <- loo(model_Gen4)
loo_Gen5 <- loo(model_Gen5)

loo_Opp0 <- loo(model_Opp0)
loo_Opp1 <- loo(model_Opp1)
loo_Opp2 <- loo(model_Opp2)
loo_Opp4 <- loo(model_Opp4)
loo_Opp5 <- loo(model_Opp5)

loo_Spe0 <- loo(model_Spe0)
loo_Spe1 <- loo(model_Spe1)
loo_Spe2 <- loo(model_Spe2)
loo_Spe4 <- loo(model_Spe4)
loo_Spe5 <- loo(model_Spe5)

capture.output(loo_Gen0, file = file.path(outdir, "loo_Gen0.txt"))
capture.output(loo_Gen1, file = file.path(outdir, "loo_Gen1.txt"))
capture.output(loo_Gen2, file = file.path(outdir, "loo_Gen2.txt"))
capture.output(loo_Gen4, file = file.path(outdir, "loo_Gen4.txt"))
capture.output(loo_Gen5, file = file.path(outdir, "loo_Gen5.txt"))

capture.output(loo_Opp0, file = file.path(outdir, "loo_Opp0.txt"))
capture.output(loo_Opp1, file = file.path(outdir, "loo_Opp1.txt"))
capture.output(loo_Opp2, file = file.path(outdir, "loo_Opp2.txt"))
capture.output(loo_Opp4, file = file.path(outdir, "loo_Opp4.txt"))
capture.output(loo_Opp5, file = file.path(outdir, "loo_Opp5.txt"))

capture.output(loo_Spe0, file = file.path(outdir, "loo_Spe0.txt"))
capture.output(loo_Spe1, file = file.path(outdir, "loo_Spe1.txt"))
capture.output(loo_Spe2, file = file.path(outdir, "loo_Spe2.txt"))
capture.output(loo_Spe4, file = file.path(outdir, "loo_Spe4.txt"))
capture.output(loo_Spe5, file = file.path(outdir, "loo_Spe5.txt"))

capture.output(
  loo_compare(loo_Gen0, loo_Gen1, loo_Gen2, loo_Gen4, loo_Gen5),
  file = file.path(outdir, "loo_compare_generalist.txt")
)

capture.output(
  loo_compare(loo_Opp0, loo_Opp1, loo_Opp2, loo_Opp4, loo_Opp5),
  file = file.path(outdir, "loo_compare_opportunist.txt")
)

capture.output(
  loo_compare(loo_Spe0, loo_Spe1, loo_Spe2, loo_Spe4, loo_Spe5),
  file = file.path(outdir, "loo_compare_specialist.txt")
)

print(loo_compare(loo_Gen0, loo_Gen1, loo_Gen2, loo_Gen4, loo_Gen5))
print(loo_compare(loo_Opp0, loo_Opp1, loo_Opp2, loo_Opp4, loo_Opp5))
print(loo_compare(loo_Spe0, loo_Spe1, loo_Spe2, loo_Spe4, loo_Spe5))









# =============================================================================
# iCAMP analyses for Hypothesis 3
# Assembly processes across orchard soils with contrasting forest cover
# =============================================================================

suppressPackageStartupMessages({
  library(iCAMP)
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ape)
  library(ggplot2)
  library(writexl)
  library(NST)
})

# =============================================================================
# iCAMP paths and safe output directory
# =============================================================================

save.wd <- normalizePath(file.path(getwd(), outdir, "iCAMP_ITS"), winslash = "/", mustWork = FALSE)

if (!dir.exists(save.wd)) {
  dir.create(save.wd, recursive = TRUE, showWarnings = FALSE)
}

if (!dir.exists(save.wd)) {
  stop("Could not create iCAMP output directory: ", save.wd)
}

prefix <- "orchard_ITS_hyp3"
rand.time <- 100
nworker <- min(12, parallel::detectCores())
memory.G <- 50

cat("iCAMP output directory:\n", save.wd, "\n")
setwd(save.wd)

# remove old partial files only if you want a clean rerun
old_files <- c("pd.bin", "pd.desc", "pd.taxon.name.csv", "path.rda")
old_files <- file.path(save.wd, old_files)
old_files_exist <- old_files[file.exists(old_files)]
if (length(old_files_exist) > 0) {
  file.remove(old_files_exist)
}

# -----------------------------------------------------------------------------
# Prepare orchard phyloseq object and metadata
# -----------------------------------------------------------------------------

orchard_phy <- phylo_crops_fun

tree <- phy_tree(orchard_phy)
if (is.null(tree)) stop("phylo_crops_fun has no phylogenetic tree.")
if (!ape::is.rooted(tree)) {
  tree <- ape::root(tree, outgroup = tree$tip.label[1], resolve.root = TRUE)
}

orchard_otu <- as(otu_table(orchard_phy), "matrix")
if (taxa_are_rows(orchard_phy)) orchard_otu <- t(orchard_otu)
comm <- as.data.frame(orchard_otu, check.names = FALSE)

orchard_tax <- as(tax_table(orchard_phy), "matrix")
if (is.null(rownames(orchard_tax))) rownames(orchard_tax) <- colnames(comm)
clas <- as.data.frame(orchard_tax, check.names = FALSE)

orchard_sam <- as(sample_data(orchard_phy), "data.frame")
orchard_sam$SampleID <- rownames(orchard_sam)

# -----------------------------------------------------------------------------
# Define treatment groups based on P_2500 classes
# Use the same cutpoints already used in the main script
# -----------------------------------------------------------------------------

orchard_sam <- orchard_sam %>%
  mutate(
    P_2500 = as.numeric(P_2500),
    FC_cat = case_when(
      P_2500 <= 0.3584376 ~ "L",
      P_2500 <= 0.5224324 ~ "M",
      P_2500 > 0.5224324  ~ "H",
      TRUE ~ NA_character_
    )
  )

treat <- orchard_sam %>%
  dplyr::select(FC_cat) %>%
  rename(cat = FC_cat)

rownames(treat) <- orchard_sam$SampleID

env <- orchard_sam %>%
  dplyr::select(pH, SOM, N, claysilt) %>%
  mutate(across(everything(), as.numeric))

rownames(env) <- orchard_sam$SampleID

# -----------------------------------------------------------------------------
# Match sample IDs across comm / treat / env
# -----------------------------------------------------------------------------

sampid.check <- match.name(rn.list = list(comm = comm, treat = treat, env = env))
comm <- sampid.check$comm
treat <- sampid.check$treat
env <- sampid.check$env

comm <- comm[, colSums(comm) > 0, drop = FALSE]

# -----------------------------------------------------------------------------
# Match taxa IDs across comm / clas / tree
# -----------------------------------------------------------------------------

spid.check <- match.name(
  cn.list = list(comm = comm),
  rn.list = list(class = clas),
  tree.list = list(tree = tree)
)

comm <- spid.check$comm
clas <- spid.check$class
tree <- spid.check$tree

write.csv(
  data.frame(SampleID = rownames(comm), treat, env, check.names = FALSE),
  file.path(save.wd, paste0(prefix, "_sample_table_for_icamp.csv")),
  row.names = FALSE
)

write.csv(
  data.frame(ASV = colnames(comm), clas[colnames(comm), , drop = FALSE], check.names = FALSE),
  file.path(save.wd, paste0(prefix, "_taxonomy_for_icamp.csv")),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Phylogenetic distance matrix
# -----------------------------------------------------------------------------

setwd(save.wd)

if (!file.exists(file.path(save.wd, "pd.desc"))) {
  pd.big <- iCAMP::pdist.big(
    tree = tree,
    wd = save.wd,
    nworker = nworker,
    memory.G = memory.G
  )
} else {
  pd.big <- list()
  pd.big$tip.label <- read.csv(
    file.path(save.wd, "pd.taxon.name.csv"),
    row.names = 1,
    stringsAsFactors = FALSE
  )[, 1]
  pd.big$pd.wd <- save.wd
  pd.big$pd.file <- "pd.desc"
  pd.big$pd.name.file <- "pd.taxon.name.csv"
}

# -----------------------------------------------------------------------------
# Niche difference matrix
# -----------------------------------------------------------------------------

niche.dif <- iCAMP::dniche(
  env = env,
  comm = comm,
  method = "niche.value",
  nworker = nworker,
  out.dist = FALSE,
  bigmemo = TRUE,
  nd.wd = save.wd
)

# -----------------------------------------------------------------------------
# Within-bin phylogenetic signal assessment
# -----------------------------------------------------------------------------

ds <- 0.2
bin.size.limit <- 24

if (!ape::is.rooted(tree)) {
  tree.rt <- iCAMP::midpoint.root.big(
    tree = tree,
    pd.desc = pd.big$pd.file,
    pd.spname = pd.big$tip.label,
    pd.wd = pd.big$pd.wd,
    nworker = nworker
  )
  tree <- tree.rt$tree
}

phylobin <- taxa.binphy.big(
  tree = tree,
  pd.desc = pd.big$pd.file,
  pd.spname = pd.big$tip.label,
  pd.wd = pd.big$pd.wd,
  ds = ds,
  bin.size.limit = bin.size.limit,
  nworker = nworker
)

sp.bin <- phylobin$sp.bin[, 3, drop = FALSE]

sp.ra <- colMeans(comm / rowSums(comm))
abcut <- 3
commc <- comm[, colSums(comm) >= abcut, drop = FALSE]
spname.use <- colnames(commc)

sp.bin <- sp.bin[rownames(sp.bin) %in% spname.use, , drop = FALSE]
sp.ra <- sp.ra[names(sp.ra) %in% spname.use]
sp.bin <- sp.bin[spname.use, , drop = FALSE]
sp.ra <- sp.ra[spname.use]

pd.spname.use <- pd.big$tip.label[pd.big$tip.label %in% spname.use]
nd.spname.use <- niche.dif$names[niche.dif$names %in% spname.use]

binps <- iCAMP::ps.bin(
  sp.bin = sp.bin,
  sp.ra = sp.ra,
  spname.use = spname.use,
  pd.desc = pd.big$pd.file,
  pd.spname = pd.spname.use,
  pd.wd = pd.big$pd.wd,
  nd.list = niche.dif$nd,
  nd.spname = nd.spname.use,
  ndbig.wd = niche.dif$nd.wd,
  cor.method = "spearman",
  r.cut = 0.1,
  p.cut = 0.05,
  min.spn = 5
)

write.csv(
  data.frame(ds = ds, n.min = bin.size.limit, binps$Index),
  file.path(save.wd, paste0(prefix, "_PhyloSignalSummary.csv")),
  row.names = FALSE
)

write.csv(
  data.frame(ds = ds, n.min = bin.size.limit, binID = rownames(binps$detail), binps$detail),
  file.path(save.wd, paste0(prefix, "_PhyloSignalDetail.csv")),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Main iCAMP analysis
# -----------------------------------------------------------------------------

sig.index <- "Confidence"

icres <- iCAMP::icamp.big(
  comm = comm,
  pd.desc = pd.big$pd.file,
  pd.spname = pd.big$tip.label,
  pd.wd = pd.big$pd.wd,
  rand = rand.time,
  tree = tree,
  prefix = prefix,
  ds = ds,
  pd.cut = NA,
  sp.check = TRUE,
  phylo.rand.scale = "within.bin",
  taxa.rand.scale = "across.all",
  phylo.metric = "bMPD",
  sig.index = sig.index,
  bin.size.limit = bin.size.limit,
  nworker = nworker,
  memory.G = memory.G,
  rtree.save = FALSE,
  detail.save = TRUE,
  qp.save = FALSE,
  detail.null = FALSE,
  ignore.zero = TRUE,
  output.wd = save.wd,
  correct.special = TRUE,
  unit.sum = rowSums(comm),
  special.method = "depend",
  ses.cut = 1.96,
  rc.cut = 0.95,
  conf.cut = 0.975,
  omit.option = "no",
  meta.ab = NULL
)

save(icres, file = file.path(save.wd, paste0(prefix, "_iCAMP_detail.rda")))

# -----------------------------------------------------------------------------
# Check of null-value normality using SES.RC
# -----------------------------------------------------------------------------

icres_ses <- iCAMP::icamp.big(
  comm = comm,
  pd.desc = pd.big$pd.file,
  pd.spname = pd.big$tip.label,
  pd.wd = pd.big$pd.wd,
  rand = rand.time,
  tree = tree,
  prefix = paste0(prefix, "_SESRC"),
  ds = ds,
  pd.cut = NA,
  sp.check = TRUE,
  phylo.rand.scale = "within.bin",
  taxa.rand.scale = "across.all",
  phylo.metric = "bMPD",
  sig.index = "SES.RC",
  bin.size.limit = bin.size.limit,
  nworker = nworker,
  memory.G = memory.G,
  rtree.save = FALSE,
  detail.save = TRUE,
  qp.save = FALSE,
  detail.null = TRUE,
  ignore.zero = TRUE,
  output.wd = save.wd,
  correct.special = TRUE,
  unit.sum = rowSums(comm),
  special.method = "depend",
  ses.cut = 1.96,
  rc.cut = 0.95,
  conf.cut = 0.975,
  omit.option = "no",
  meta.ab = NULL
)

nntest <- iCAMP::null.norm(
  icamp.output = icres_ses,
  p.norm.cut = 0.05,
  detail.out = FALSE
)

capture.output(
  nntest,
  file = file.path(save.wd, paste0(prefix, "_null_normality_test.txt"))
)

# -----------------------------------------------------------------------------
# Bin-level summaries
# -----------------------------------------------------------------------------

icbin <- iCAMP::icamp.bins(
  icamp.detail = icres$detail,
  treat = treat,
  clas = clas,
  silent = FALSE,
  boot = TRUE,
  rand.time = rand.time,
  between.group = TRUE
)

save(icbin, file = file.path(save.wd, paste0(prefix, "_iCAMP_Summary.rda")))

write.csv(icbin$Pt, file = file.path(save.wd, paste0(prefix, "_ProcessImportance_EachGroup.csv")), row.names = FALSE)
write.csv(icbin$Ptk, file = file.path(save.wd, paste0(prefix, "_ProcessImportance_EachBin_EachGroup.csv")), row.names = FALSE)
write.csv(icbin$Ptuv, file = file.path(save.wd, paste0(prefix, "_ProcessImportance_EachTurnover.csv")), row.names = FALSE)
write.csv(icbin$BPtk, file = file.path(save.wd, paste0(prefix, "_BinContributeToProcess_EachGroup.csv")), row.names = FALSE)
write.csv(
  data.frame(ID = rownames(icbin$Class.Bin), icbin$Class.Bin, stringsAsFactors = FALSE),
  file = file.path(save.wd, paste0(prefix, "_Taxon_Bin.csv")),
  row.names = FALSE
)
write.csv(icbin$Bin.TopClass, file = file.path(save.wd, paste0(prefix, "_Bin_TopTaxon.csv")), row.names = FALSE)

# -----------------------------------------------------------------------------
# Summary tables and figures
# -----------------------------------------------------------------------------

process_group <- icbin$Pt

process_group_clean <- process_group %>%
  mutate(
    Group = as.character(Group),
    HeS = as.numeric(HeS),
    HoS = as.numeric(HoS),
    DL  = as.numeric(DL),
    HD  = as.numeric(HD),
    DR  = as.numeric(DR)
  )

write.csv(
  process_group_clean,
  file.path(save.wd, paste0(prefix, "_ProcessImportance_EachGroup_clean.csv")),
  row.names = FALSE
)

process_group_long <- process_group_clean %>%
  pivot_longer(
    cols = c(HeS, HoS, DL, HD, DR),
    names_to = "Process",
    values_to = "Relative_Importance"
  ) %>%
  filter(!is.na(Group), !is.na(Relative_Importance))

process_group_main <- process_group_long %>%
  filter(Group %in% c("L", "M", "H"))

write.csv(
  process_group_main,
  file.path(save.wd, paste0(prefix, "_ProcessImportance_long_main_groups.csv")),
  row.names = FALSE
)

process_group_pie <- process_group_main %>%
  group_by(Group) %>%
  mutate(
    Relative_Percent = Relative_Importance / sum(Relative_Importance, na.rm = TRUE) * 100
  ) %>%
  ungroup()

write.csv(
  process_group_pie,
  file.path(save.wd, paste0(prefix, "_ProcessImportance_pie_table.csv")),
  row.names = FALSE
)

p_icamp_group_bar <- ggplot(
  process_group_pie,
  aes(x = Group, y = Relative_Percent, fill = Process)
) +
  geom_col(color = "black", width = 0.8) +
  theme_minimal(base_size = 13) +
  labs(x = "Forest cover class", y = "Relative importance (%)", fill = "Process") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(face = "bold")
  )

ggsave(
  file.path(save.wd, paste0(prefix, "_ProcessImportance_stackedbar.png")),
  p_icamp_group_bar,
  width = 5.5,
  height = 4.2,
  dpi = 300
)

p_icamp_group_pie <- ggplot(
  process_group_pie,
  aes(x = "", y = Relative_Percent, fill = Process)
) +
  geom_col(color = "black", width = 1) +
  coord_polar(theta = "y") +
  facet_wrap(~ Group, nrow = 1) +
  theme_void(base_size = 13) +
  labs(fill = "Process") +
  theme(strip.text = element_text(face = "bold"))

ggsave(
  file.path(save.wd, paste0(prefix, "_ProcessImportance_pieplot.png")),
  p_icamp_group_pie,
  width = 9,
  height = 3.8,
  dpi = 300
)

det_sto_summary <- process_group_main %>%
  mutate(
    Assembly_Class = case_when(
      Process %in% c("HeS", "HoS") ~ "Deterministic",
      Process %in% c("DL", "HD", "DR") ~ "Stochastic",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(Group, Assembly_Class) %>%
  summarise(
    Relative_Importance = sum(Relative_Importance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Group) %>%
  mutate(
    Relative_Percent = Relative_Importance / sum(Relative_Importance, na.rm = TRUE) * 100
  ) %>%
  ungroup()

write.csv(
  det_sto_summary,
  file.path(save.wd, paste0(prefix, "_Deterministic_vs_Stochastic_summary.csv")),
  row.names = FALSE
)

p_det_sto <- ggplot(
  det_sto_summary,
  aes(x = Group, y = Relative_Percent, fill = Assembly_Class)
) +
  geom_col(color = "black", width = 0.8) +
  theme_minimal(base_size = 13) +
  labs(x = "Forest cover class", y = "Relative importance (%)", fill = "") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(face = "bold")
  )

ggsave(
  file.path(save.wd, paste0(prefix, "_Deterministic_vs_Stochastic_barplot.png")),
  p_det_sto,
  width = 5.5,
  height = 4.2,
  dpi = 300
)

process_group_main %>%
  group_by(Group) %>%
  summarise(total = sum(Relative_Importance, na.rm = TRUE))

# -----------------------------------------------------------------------------
# Bootstrapping comparison among groups
# -----------------------------------------------------------------------------

treat.use <- treat[, 1, drop = FALSE]
icamp.result <- icres$CbMPDiCBraya

icboot <- iCAMP::icamp.boot(
  icamp.result = icamp.result,
  treat = treat.use,
  rand.time = rand.time,
  compare = TRUE,
  silent = FALSE,
  between.group = TRUE,
  ST.estimation = TRUE
)

save(icboot, file = file.path(save.wd, paste0(prefix, "_iCAMP_Boot_cat.rda")))
write.csv(icboot$summary, file = file.path(save.wd, paste0(prefix, "_iCAMP_BootSummary_cat.csv")), row.names = FALSE)
write.csv(icboot$compare, file = file.path(save.wd, paste0(prefix, "_iCAMP_Compare_cat.csv")), row.names = FALSE)

# -----------------------------------------------------------------------------
# Summary of deterministic vs stochastic assembly
# -----------------------------------------------------------------------------

if ("Process" %in% colnames(process_group_long)) {
  deterministic_terms <- c("HoS", "HeS")
  stochastic_terms <- c("DL", "HD", "DR")
  
  det_sto_summary <- process_group_long %>%
    filter(Group %in% c("L", "M", "H")) %>%
    mutate(
      Assembly_Class = case_when(
        Process %in% deterministic_terms ~ "Deterministic",
        Process %in% stochastic_terms ~ "Stochastic",
        TRUE ~ "Other"
      )
    ) %>%
    group_by(Group, Assembly_Class) %>%
    summarise(Relative_Importance = sum(Relative_Importance, na.rm = TRUE), .groups = "drop")
  
  write.csv(
    det_sto_summary,
    file.path(save.wd, paste0(prefix, "_Deterministic_vs_Stochastic_summary.csv")),
    row.names = FALSE
  )
  
  p_det_sto <- ggplot(det_sto_summary, aes(x = Group, y = Relative_Importance, fill = Assembly_Class)) +
    geom_col(color = "black", width = 0.8) +
    theme_minimal(base_size = 13) +
    labs(x = "Forest cover class", y = "Relative importance", fill = "") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(face = "bold")
    )
  
  ggsave(
    file.path(save.wd, paste0(prefix, "_Deterministic_vs_Stochastic_barplot.png")),
    p_det_sto,
    width = 5.5,
    height = 4.2,
    dpi = 300
  )
}

# -----------------------------------------------------------------------------
# Session info
# -----------------------------------------------------------------------------

writeLines(
  capture.output(sessionInfo()),
  file.path(save.wd, paste0(prefix, "_sessionInfo_iCAMP.txt"))
)











# =============================================================================
# Sloan neutral model and niche-breadth integration
# =============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(broom)
  library(stringr)
})

sloan_outdir <- file.path(outdir, "sloan_ITS")
if (!dir.exists(sloan_outdir)) dir.create(sloan_outdir, recursive = TRUE)

# -----------------------------------------------------------------------------
# Base objects
# -----------------------------------------------------------------------------

ps_base <- phylo_crops_fun
ps_all  <- data_phylo_filt

if (is.null(sample_data(ps_base, errorIfNULL = FALSE))) {
  stop("phylo_crops_fun must contain sample_data.")
}
if (is.null(tax_table(ps_base, errorIfNULL = FALSE))) {
  stop("phylo_crops_fun must contain a tax_table.")
}

# -----------------------------------------------------------------------------
# Orchard community matrix for Sloan model
# Samples in rows, taxa in columns
# -----------------------------------------------------------------------------

otu_orch <- as(otu_table(ps_base), "matrix")
if (taxa_are_rows(ps_base)) otu_orch <- t(otu_orch)

comm_sloan <- as.data.frame(otu_orch, check.names = FALSE)
comm_sloan <- comm_sloan[rowSums(comm_sloan) > 0, colSums(comm_sloan) > 0, drop = FALSE]

write.csv(
  data.frame(SampleID = rownames(comm_sloan), comm_sloan, check.names = FALSE),
  file.path(sloan_outdir, "orchard_asv_table_for_sloan.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Sloan model inputs: mean relative abundance and occurrence frequency
# -----------------------------------------------------------------------------

otu_rel <- sweep(otu_orch, 1, rowSums(otu_orch), "/")
otu_rel[is.na(otu_rel)] <- 0

mean_abund <- colMeans(otu_rel, na.rm = TRUE)
occurrence <- colSums(otu_orch > 0, na.rm = TRUE) / nrow(otu_orch)

df_neutral <- data.frame(
  Taxon = names(mean_abund),
  mean_abund = as.numeric(mean_abund),
  occurrence = as.numeric(occurrence),
  stringsAsFactors = FALSE
) %>%
  filter(mean_abund > 0, occurrence > 0)

N_sloan <- nrow(otu_orch)

write.csv(
  df_neutral,
  file.path(sloan_outdir, "sloan_input_table.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Fit Sloan neutral model
# -----------------------------------------------------------------------------

neutral_model <- nls(
  occurrence ~ 1 - exp(-m * N_sloan * mean_abund),
  data = df_neutral,
  start = list(m = 0.1),
  control = nls.control(maxiter = 500, warnOnly = TRUE)
)

capture.output(
  summary(neutral_model),
  file = file.path(sloan_outdir, "sloan_model_summary.txt")
)

m_est <- coef(neutral_model)[["m"]]

df_neutral <- df_neutral %>%
  mutate(
    predicted = 1 - exp(-m_est * N_sloan * mean_abund),
    se_pred = sqrt(predicted * (1 - predicted) / N_sloan),
    lower_ci = pmax(predicted - 1.96 * se_pred, 0),
    upper_ci = pmin(predicted + 1.96 * se_pred, 1),
    Sloan_Class = case_when(
      occurrence > upper_ci ~ "Above",
      occurrence < lower_ci ~ "Below",
      TRUE ~ "Neutral"
    )
  )

write.csv(
  df_neutral,
  file.path(sloan_outdir, "sloan_fitted_table.csv"),
  row.names = FALSE
)

capture.output(
  table(df_neutral$Sloan_Class),
  file = file.path(sloan_outdir, "sloan_class_counts.txt")
)

# -----------------------------------------------------------------------------
# Sloan model plots
# -----------------------------------------------------------------------------

p_sloan_basic <- ggplot(df_neutral, aes(x = mean_abund, y = occurrence)) +
  geom_point(alpha = 0.65, size = 2) +
  geom_line(aes(y = predicted), color = "black", linewidth = 0.9) +
  geom_line(aes(y = lower_ci), color = "grey40", linetype = "dashed", linewidth = 0.6) +
  geom_line(aes(y = upper_ci), color = "grey40", linetype = "dashed", linewidth = 0.6) +
  scale_x_log10() +
  theme_bw(base_size = 13) +
  labs(
    x = "Mean relative abundance (log scale)",
    y = "Occurrence frequency",
    title = "Sloan neutral model"
  ) +
  theme(panel.grid = element_blank())

ggsave(
  file.path(sloan_outdir, "sloan_neutral_model_basic.png"),
  p_sloan_basic,
  width = 6.2,
  height = 4.8,
  dpi = 300
)

p_sloan_class <- ggplot(df_neutral, aes(x = mean_abund, y = occurrence, color = Sloan_Class)) +
  geom_point(alpha = 0.75, size = 2) +
  geom_line(aes(x = mean_abund, y = predicted), color = "black", linewidth = 0.9, inherit.aes = FALSE) +
  geom_line(aes(x = mean_abund, y = lower_ci), color = "grey40", linetype = "dashed", linewidth = 0.6, inherit.aes = FALSE) +
  geom_line(aes(x = mean_abund, y = upper_ci), color = "grey40", linetype = "dashed", linewidth = 0.6, inherit.aes = FALSE) +
  scale_x_log10() +
  scale_color_manual(values = c("Above" = "#D55E00", "Neutral" = "grey50", "Below" = "#0072B2")) +
  theme_bw(base_size = 13) +
  labs(
    x = "Mean relative abundance (log scale)",
    y = "Occurrence frequency",
    color = "Neutral class",
    title = "Sloan neutral model with confidence envelope"
  ) +
  theme(panel.grid = element_blank())

ggsave(
  file.path(sloan_outdir, "sloan_neutral_model_by_class.png"),
  p_sloan_class,
  width = 6.4,
  height = 4.9,
  dpi = 300
)

# -----------------------------------------------------------------------------
# Shared, orchard-only, and forest-only ASVs
# -----------------------------------------------------------------------------

ps_orchard_all <- subset_samples(ps_all, Group == "Orchard")
ps_forest_all  <- subset_samples(ps_all, Group == "Forest")

otu_orch_all <- as(otu_table(ps_orchard_all), "matrix")
otu_forest_all <- as(otu_table(ps_forest_all), "matrix")

if (taxa_are_rows(ps_orchard_all)) otu_orch_all <- t(otu_orch_all)
if (taxa_are_rows(ps_forest_all)) otu_forest_all <- t(otu_forest_all)

orchard_asvs <- colnames(otu_orch_all)[colSums(otu_orch_all) > 0]
forest_asvs  <- colnames(otu_forest_all)[colSums(otu_forest_all) > 0]

shared_asvs <- intersect(orchard_asvs, forest_asvs)
orchard_unique_asvs <- setdiff(orchard_asvs, forest_asvs)
forest_unique_asvs  <- setdiff(forest_asvs, orchard_asvs)

asv_partition <- data.frame(
  Metric = c("Shared_ASVs", "Orchard_only_ASVs", "Forest_only_ASVs"),
  Count = c(length(shared_asvs), length(orchard_unique_asvs), length(forest_unique_asvs))
)

write.csv(
  asv_partition,
  file.path(sloan_outdir, "asv_partition_summary.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Cross-tabulation of Sloan classes against ASV partition
# -----------------------------------------------------------------------------

df_neutral <- df_neutral %>%
  mutate(
    ASV_Set = case_when(
      Taxon %in% shared_asvs ~ "Shared",
      Taxon %in% orchard_unique_asvs ~ "Orchard_only",
      Taxon %in% forest_unique_asvs ~ "Forest_only",
      TRUE ~ "Other"
    )
  )

write.csv(
  df_neutral,
  file.path(sloan_outdir, "sloan_fitted_table_with_asv_partition.csv"),
  row.names = FALSE
)

sloan_partition_summary <- df_neutral %>%
  count(ASV_Set, Sloan_Class) %>%
  tidyr::pivot_wider(names_from = Sloan_Class, values_from = n, values_fill = 0)

write.csv(
  sloan_partition_summary,
  file.path(sloan_outdir, "sloan_by_asv_partition_summary.csv"),
  row.names = FALSE
)

shared_types <- df_neutral$Sloan_Class[df_neutral$Taxon %in% shared_asvs]
orchard_types <- df_neutral$Sloan_Class[df_neutral$Taxon %in% orchard_unique_asvs]
forest_types <- df_neutral$Sloan_Class[df_neutral$Taxon %in% forest_unique_asvs]

capture.output(
  table(shared_types),
  file = file.path(sloan_outdir, "sloan_shared_class_counts.txt")
)
capture.output(
  table(orchard_types),
  file = file.path(sloan_outdir, "sloan_orchard_unique_class_counts.txt")
)
capture.output(
  table(forest_types),
  file = file.path(sloan_outdir, "sloan_forest_unique_class_counts.txt")
)

if (length(shared_types) > 0) {
  capture.output(
    chisq.test(table(shared_types)),
    file = file.path(sloan_outdir, "chisq_shared_sloan_classes.txt")
  )
}

# -----------------------------------------------------------------------------
# Taxonomy lookup and NB category mapping
# -----------------------------------------------------------------------------

tax_df <- as.data.frame(tax_table(ps_base), stringsAsFactors = FALSE) %>%
  rownames_to_column("Taxon")

if (!"Genus" %in% colnames(tax_df)) {
  genus_col <- grep("Genus", colnames(tax_df), value = TRUE)[1]
  if (is.na(genus_col)) stop("No Genus column found in tax_table(ps_base).")
  tax_df$Genus <- tax_df[[genus_col]]
}

tax_df <- tax_df %>%
  mutate(
    Genus = as.character(Genus),
    Genus_clean = Genus %>%
      str_remove("^g__") %>%
      str_remove("^g_") %>%
      str_replace_all("\\s+", "")
  )

if (!"NB_Category" %in% colnames(tax_df)) {
  if (!exists("nb_genus_df")) {
    stop("NB_Category is not present in tax_table(ps_base), and nb_genus_df is not available.")
  }
  
  genus2nb <- nb_genus_df %>%
    mutate(
      Genus_clean = as.character(Genus) %>%
        str_remove("^g__") %>%
        str_remove("^g_") %>%
        str_replace_all("\\s+", "")
    ) %>%
    dplyr::select(Genus_clean, NB_Category) %>%
    distinct()
  
  tax_df <- tax_df %>%
    left_join(genus2nb, by = "Genus_clean")
}

tax_df <- tax_df %>%
  mutate(
    NB_Category = ifelse(is.na(NB_Category) | NB_Category == "", "unassigned", as.character(NB_Category))
  ) %>%
  dplyr::select(Taxon, Genus, NB_Category)

# -----------------------------------------------------------------------------
# Sloan-above taxa and their partition into shared vs orchard-only
# -----------------------------------------------------------------------------

sloan_above_fun <- df_neutral %>% filter(Sloan_Class == "Above")
sloan_neutral_fun <- df_neutral %>% filter(Sloan_Class == "Neutral")
sloan_below_fun <- df_neutral %>% filter(Sloan_Class == "Below")

write.csv(sloan_above_fun,   file.path(sloan_outdir, "sloan_above_taxa.csv"),   row.names = FALSE)
write.csv(sloan_neutral_fun, file.path(sloan_outdir, "sloan_neutral_taxa.csv"), row.names = FALSE)
write.csv(sloan_below_fun,   file.path(sloan_outdir, "sloan_below_taxa.csv"),   row.names = FALSE)

ids_sloan_above_shared <- intersect(sloan_above_fun$Taxon, shared_asvs)
ids_sloan_above_unique <- intersect(sloan_above_fun$Taxon, orchard_unique_asvs)

df_sloan_shared <- data.frame(Taxon = ids_sloan_above_shared, stringsAsFactors = FALSE) %>%
  left_join(sloan_above_fun, by = "Taxon") %>%
  left_join(tax_df, by = "Taxon")

df_sloan_unique <- data.frame(Taxon = ids_sloan_above_unique, stringsAsFactors = FALSE) %>%
  left_join(sloan_above_fun, by = "Taxon") %>%
  left_join(tax_df, by = "Taxon")

write.csv(df_sloan_shared, file.path(sloan_outdir, "sloan_above_shared_taxa.csv"), row.names = FALSE)
write.csv(df_sloan_unique, file.path(sloan_outdir, "sloan_above_orchard_unique_taxa.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# Relative abundance of Sloan-above taxa across orchard samples
# Denominator is always the full orchard community (phylo_crops_fun)
# -----------------------------------------------------------------------------

otu_base <- as(otu_table(ps_base), "matrix")
if (taxa_are_rows(ps_base)) otu_base <- t(otu_base)

otu_base_rel <- sweep(otu_base, 2, colSums(otu_base), "/")
otu_base_rel[is.na(otu_base_rel)] <- 0

meta_df <- as(sample_data(ps_base), "data.frame") %>%
  rownames_to_column("Sample")

meta_df$P_2500 <- as.numeric(as.character(meta_df$P_2500))

# Shared Sloan-above taxa
shared_taxa_ok <- intersect(df_sloan_shared$Taxon, colnames(otu_base_rel))
shared_nb_lookup <- tax_df %>% filter(Taxon %in% shared_taxa_ok)

df_shared_nb <- data.frame()
if (length(shared_taxa_ok) > 0) {
  for (nb in unique(shared_nb_lookup$NB_Category)) {
    taxa_nb <- shared_nb_lookup$Taxon[shared_nb_lookup$NB_Category == nb]
    ra_nb <- rowSums(otu_base_rel[, taxa_nb, drop = FALSE], na.rm = TRUE)
    
    df_tmp <- data.frame(
      Sample = names(ra_nb),
      NB_Category = nb,
      RA_subset = as.numeric(ra_nb),
      Percent_subset = 100 * as.numeric(ra_nb),
      Set = "shared",
      stringsAsFactors = FALSE
    )
    
    df_shared_nb <- bind_rows(df_shared_nb, df_tmp)
  }
}

# Orchard-unique Sloan-above taxa
unique_taxa_ok <- intersect(df_sloan_unique$Taxon, colnames(otu_base_rel))
unique_nb_lookup <- tax_df %>% filter(Taxon %in% unique_taxa_ok)

df_unique_nb <- data.frame()
if (length(unique_taxa_ok) > 0) {
  for (nb in unique(unique_nb_lookup$NB_Category)) {
    taxa_nb <- unique_nb_lookup$Taxon[unique_nb_lookup$NB_Category == nb]
    ra_nb <- rowSums(otu_base_rel[, taxa_nb, drop = FALSE], na.rm = TRUE)
    
    df_tmp <- data.frame(
      Sample = names(ra_nb),
      NB_Category = nb,
      RA_subset = as.numeric(ra_nb),
      Percent_subset = 100 * as.numeric(ra_nb),
      Set = "unique",
      stringsAsFactors = FALSE
    )
    
    df_unique_nb <- bind_rows(df_unique_nb, df_tmp)
  }
}

df_nb <- bind_rows(df_shared_nb, df_unique_nb) %>%
  left_join(meta_df %>% dplyr::select(Sample, P_2500), by = "Sample") %>%
  filter(!is.na(P_2500))

write.csv(
  df_nb,
  file.path(sloan_outdir, "sloan_above_nb_abundance_by_sample.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Plot Sloan-above abundance vs P_2500 by NB category and set
# -----------------------------------------------------------------------------

deg <- 2

p_sloan_nb_shared <- ggplot(df_nb %>% filter(Set == "shared"),
                            aes(x = P_2500, y = Percent_subset, color = NB_Category)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", formula = y ~ poly(x, deg, raw = TRUE), se = TRUE) +
  facet_wrap(~ NB_Category, scales = "free_y") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank(), legend.position = "none") +
  labs(
    x = expression(P[2500]),
    y = "Sloan-above taxa (% of total orchard community)",
    title = "Shared Sloan-above taxa"
  )

ggsave(
  file.path(sloan_outdir, "sloan_above_shared_nb_vs_P2500.png"),
  p_sloan_nb_shared,
  width = 6.6,
  height = 4.8,
  dpi = 300
)

p_sloan_nb_unique <- ggplot(df_nb %>% filter(Set == "unique"),
                            aes(x = P_2500, y = Percent_subset, color = NB_Category)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", formula = y ~ poly(x, deg, raw = TRUE), se = TRUE) +
  facet_wrap(~ NB_Category, scales = "free_y") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank(), legend.position = "none") +
  labs(
    x = expression(P[2500]),
    y = "Sloan-above taxa (% of total orchard community)",
    title = "Orchard-unique Sloan-above taxa"
  )

ggsave(
  file.path(sloan_outdir, "sloan_above_orchard_unique_nb_vs_P2500.png"),
  p_sloan_nb_unique,
  width = 6.6,
  height = 4.8,
  dpi = 300
)

# -----------------------------------------------------------------------------
# Polynomial models for Sloan-above abundance vs P_2500
# -----------------------------------------------------------------------------

mods <- df_nb %>%
  group_by(Set, NB_Category) %>%
  nest() %>%
  mutate(
    n = purrr::map_int(data, nrow),
    model = purrr::map(
      data,
      ~ {
        dat <- .x
        if (nrow(dat) < (deg + 2)) return(NULL)
        if (sd(dat$Percent_subset, na.rm = TRUE) == 0) return(NULL)
        if (sd(dat$P_2500, na.rm = TRUE) == 0) return(NULL)
        lm(Percent_subset ~ poly(P_2500, deg, raw = TRUE), data = dat)
      }
    ),
    model_summary = purrr::map(
      model,
      ~ {
        if (is.null(.x)) {
          tibble(
            r.squared = NA_real_,
            adj.r.squared = NA_real_,
            statistic = NA_real_,
            df = NA_real_,
            df.residual = NA_real_,
            p.value = NA_real_,
            AIC = NA_real_
          )
        } else {
          broom::glance(.x) %>%
            dplyr::select(r.squared, adj.r.squared, statistic, df, df.residual, p.value, AIC)
        }
      }
    ),
    coefficients = purrr::map(
      model,
      ~ {
        if (is.null(.x)) {
          tibble(
            term = character(),
            estimate = numeric(),
            std.error = numeric(),
            statistic = numeric(),
            p.value = numeric()
          )
        } else {
          broom::tidy(.x)
        }
      }
    )
  ) %>%
  unnest(model_summary)

mods_eq <- mods %>%
  transmute(
    Set,
    NB_Category,
    n,
    r2 = r.squared,
    adj_r2 = adj.r.squared,
    F = statistic,
    df1 = df,
    df2 = df.residual,
    p_model = p.value,
    AIC = AIC
  )

coef_tables <- mods %>%
  dplyr::select(Set, NB_Category, coefficients) %>%
  unnest(coefficients)

write.csv(
  mods_eq,
  file.path(sloan_outdir, "sloan_above_nb_polynomial_model_summary.csv"),
  row.names = FALSE
)

write.csv(
  coef_tables,
  file.path(sloan_outdir, "sloan_above_nb_polynomial_model_coefficients.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Summary table for Sloan-above taxa by set and NB category
# -----------------------------------------------------------------------------

sloan_nb_summary <- df_nb %>%
  group_by(Set, NB_Category) %>%
  summarise(
    n_samples = n(),
    mean_percent = mean(Percent_subset, na.rm = TRUE),
    sd_percent = sd(Percent_subset, na.rm = TRUE),
    min_percent = min(Percent_subset, na.rm = TRUE),
    max_percent = max(Percent_subset, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  sloan_nb_summary,
  file.path(sloan_outdir, "sloan_above_nb_summary.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Session information
# -----------------------------------------------------------------------------

writeLines(
  capture.output(sessionInfo()),
  file.path(sloan_outdir, "sessionInfo_sloan_ITS.txt")
)
