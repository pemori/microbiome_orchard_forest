# =============================================================================
# Script: 04_hypothesis_3_bac.R
# Project: From forest to fields: The role of soil microbiome spillover in
#          agroecosystem sustainability
# Author: Pedro Mondaca
# Contact: pedromondaca@outlook.com
# Submitted to: Science Advances
#
# Description:
# This script evaluates hypothesis 3 for the bacterial community:
# microbial assembly-related patterns, niche breadth distributions, and
# co-occurrence-related abundance structure in orchard soils with contrasting
# surrounding forest cover. It compares forest and orchard relative abundance,
# estimates niche breadth categories, summarizes relative abundance of
# specialist/opportunist/generalist taxa, and fits frequentist, beta-regression,
# and Bayesian models.
#
# Required files in the working directory:
#   - results_16S/phyloseq_16S_asv_filtered.rds
#   - results_16S/phyloseq_16S_genus.rds
#   - results_16S/otu_crops_16S_genus.csv
#   - results_16S/otu_forest_16S_genus.csv
#   - results_16S/otu_table_16S_family_relative_abundance.csv
#   - results_16S/mapping_16S_ordered.csv
#   - results_hypothesis1_bac/orchard_summary_bact.csv
#   - results_hypothesis2_bac/orchard_summary_with_mntd.csv
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(phyloseq)
  library(ecospat)
  library(FactoMineR)
  library(factoextra)
  library(brms)
  library(ape)
  library(igraph)
  library(tibble)
  library(performance)
  library(emmeans)
  library(betareg)
})

if (!dir.exists("results_hypothesis3_bac")) dir.create("results_hypothesis3_bac")

# =============================================================================
# Load objects
# =============================================================================

data_phylo_filt <- readRDS("results_16S/phyloseq_16S_asv_filtered.rds")
data_phylo_g_bac <- readRDS("results_16S/phyloseq_16S_genus.rds")

mapping <- read.csv("results_16S/mapping_16S_ordered.csv", row.names = 1, check.names = FALSE)
otu_crops_g_bac <- read.csv("results_16S/otu_crops_16S_genus.csv", row.names = 1, check.names = FALSE)
otu_forest_g_bac <- read.csv("results_16S/otu_forest_16S_genus.csv", row.names = 1, check.names = FALSE)
data_fam_b_RA <- read.csv("results_16S/otu_table_16S_family_relative_abundance.csv", row.names = 1, check.names = FALSE)

orchard_summary_df <- read.csv("results_hypothesis2_bac/orchard_summary_with_mntd.csv", row.names = 1, check.names = FALSE)

sample_data(data_phylo_filt) <- sample_data(mapping)
sample_data(data_phylo_g_bac) <- sample_data(mapping)

phylo_crops_bac <- subset_samples(data_phylo_filt, cat != "F")
phylo_forest_bac <- subset_samples(data_phylo_filt, cat == "F")

phylo_crops_g_bac <- subset_samples(data_phylo_g_bac, cat != "F")
phylo_forest_g_bac <- subset_samples(data_phylo_g_bac, cat == "F")

forest_fam_b_RA <- data_fam_b_RA[grepl("N", rownames(data_fam_b_RA)), , drop = FALSE]
crops_fam_b_RA <- data_fam_b_RA[!grepl("N", rownames(data_fam_b_RA)), , drop = FALSE]

tab <- as.data.frame(sample_data(phylo_crops_bac))
env_bact <- tab[, c("pH", "SOM", "N", "claysilt"), drop = FALSE]
disp_bact <- tab[, c("P_250", "P_1000", "P_2500"), drop = FALSE]

# =============================================================================
# Hypothesis 3
# Compared microbial assembly process, niche breadth distributions, and
# co-occurrence patterns of the microbial community from orchard soils with
# <40% and >=40% forest cover
# =============================================================================

# =============================================================================
# Comparison of mean relative abundance between forest and orchard
# =============================================================================

promedio_abundancia_b_forest <- colMeans(forest_fam_b_RA, na.rm = TRUE)
promedio_abundancia_b_crops <- colMeans(crops_fam_b_RA, na.rm = TRUE)

df_comparacion <- data.frame(
  Microorganismo = names(promedio_abundancia_b_forest),
  Promedio_df1 = promedio_abundancia_b_forest,
  Promedio_df2 = promedio_abundancia_b_crops
)

p_family_forest_crop <- ggplot(df_comparacion, aes(x = Promedio_df1, y = Promedio_df2)) +
  geom_point(color = "blue") +
  labs(
    x = "Mean relative abundance in forest",
    y = "Mean relative abundance in orchard"
  ) +
  theme_minimal() +
  xlim(0, 5) + ylim(0, 5)

ggsave("results_hypothesis3_bac/family_mean_abundance_forest_vs_orchard.png",
       p_family_forest_crop, width = 6, height = 5, dpi = 300)

# =============================================================================
# Genus-level comparison
# =============================================================================

otu_crops_RA <- otu_crops_g_bac / rowSums(otu_crops_g_bac)
otu_forest_RA <- otu_forest_g_bac / rowSums(otu_forest_g_bac)

otu_crops_RA[is.na(otu_crops_RA)] <- 0
otu_forest_RA[is.na(otu_forest_RA)] <- 0

promedio_abundancia_crops <- colMeans(otu_crops_RA, na.rm = TRUE)
promedio_abundancia_forest <- colMeans(otu_forest_RA, na.rm = TRUE)

df_comparacion_genero <- data.frame(
  Genero = colnames(otu_crops_g_bac),
  Promedio_Cultivos = promedio_abundancia_crops,
  Promedio_Bosques = promedio_abundancia_forest
)

p_genus_forest_crop <- ggplot(df_comparacion_genero, aes(x = Promedio_Bosques, y = Promedio_Cultivos)) +
  geom_point(color = "blue", alpha = 0.7) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  labs(
    x = "Mean relative abundance in forest",
    y = "Mean relative abundance in orchard"
  ) +
  theme_minimal() +
  xlim(0, max(df_comparacion_genero$Promedio_Bosques, na.rm = TRUE)) +
  ylim(0, max(df_comparacion_genero$Promedio_Cultivos, na.rm = TRUE))

ggsave("results_hypothesis3_bac/genus_mean_abundance_forest_vs_orchard.png",
       p_genus_forest_crop, width = 6, height = 5, dpi = 300)

# =============================================================================
# PCA of environmental variables
# =============================================================================
tab
env_bact <- as.matrix(tab[,c(9,10,11,12,13,14,16,18)])
env_bact
pca_bact <- PCA(env_bact, graph = FALSE)
pca_bact
capture.output(pca_bact$eig, file = "results_hypothesis3_bac/pca_env_eigenvalues.txt")
write.csv(pca_bact$var$coord, "results_hypothesis3_bac/pca_env_variable_coordinates.csv")

p_pca_vars <- fviz_pca_var(
  pca_bact,
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

p_pca_vars

ggsave(
  "results_hypothesis3_bac/pca_env_variables.png",
  p_pca_vars,
  width = 6,
  height = 5,
  dpi = 300
)

# =============================================================================
# Niche breadth analysis
# =============================================================================

Block_crop <- tab$Block
ID <- rownames(tab)

env_bact2 <- data.frame(
  pH = as.numeric(env_bact[, "pH"]),
  SOM = as.numeric(env_bact[, "SOM"]),
  N = as.numeric(env_bact[, "N"]),
  `Clay+Silt` = as.numeric(env_bact[, "claysilt"]),
  disp_bact_250 = as.numeric(disp_bact$P_250),
  disp_bact_1000 = as.numeric(disp_bact$P_1000),
  disp_bact_2500 = as.numeric(disp_bact$P_2500),
  Block_crop = Block_crop,
  row.names = ID,
  check.names = FALSE
)

niche_bact_crops <- data.frame(
  otu_crops_RA,
  env_bact2[, c("pH", "SOM", "N", "Clay+Silt", "disp_bact_2500", "Block_crop")],
  check.names = FALSE
)

str(niche_bact_crops[, 500:513])

colnames(niche_bact_crops)[colnames(niche_bact_crops) == "claysilt"] <- "Clay+Silt"
dim(niche_bact_crops)
niche_bact_crops[1:5,500:513]
POSNB <- ecospat.nichePOSNB(
  niche_bact_crops,
  colvar = c(508:511),
  colfreq = c(1:507)
)
avg_NB <- ecospat.nicheNBmean(POSNB, w = c(1, 1, 1, 1))
write.csv(as.data.frame(avg_NB), "results_hypothesis3_bac/average_niche_breadth.csv")

POSNB <- as.data.frame(POSNB)
pH_nb <- POSNB$pH_nb
SOM_nb <- POSNB$SOM_nb
N_nb <- POSNB$N_nb
ClaySilt_nb <- POSNB$`Clay+Silt_nb`

bacti_pal_graph <- niche_bact_crops %>%
  tibble::rownames_to_column("ID") %>%
  mutate(disp_bact_2500 = as.double(disp_bact_2500))

bacti_pal_graph_long <- bacti_pal_graph %>%
  pivot_longer(
    cols = -c("ID", "pH", "SOM", "N", "Clay+Silt", "disp_bact_2500", "Block_crop"),
    names_to = "Genus",
    values_to = "Relative_Abundance"
  )

bacti_pal_graph_long$pH_nb <- rep(pH_nb, length.out = nrow(bacti_pal_graph_long))
bacti_pal_graph_long$SOM_nb <- rep(SOM_nb, length.out = nrow(bacti_pal_graph_long))
bacti_pal_graph_long$N_nb <- rep(N_nb, length.out = nrow(bacti_pal_graph_long))
bacti_pal_graph_long$claysilt_nb <- rep(ClaySilt_nb, length.out = nrow(bacti_pal_graph_long))

bacti_pal_graph_long <- bacti_pal_graph_long %>%
  filter(
    !is.na(pH_nb) & pH_nb > 0,
    !is.na(SOM_nb) & SOM_nb > 0,
    !is.na(N_nb) & N_nb > 0,
    !is.na(claysilt_nb) & claysilt_nb > 0
  ) %>%
  mutate(
    geometric_mean_nb = exp((log(pH_nb) + log(SOM_nb) + log(N_nb) + log(claysilt_nb)) / 4)
  )

# =============================================================================
# Forest cover categories and per-sample summed relative abundance
# =============================================================================

bacti_pal_graph_long <- bacti_pal_graph_long %>%
  mutate(
    Forest_coverage = case_when(
      disp_bact_2500 <= 0.3584376 ~ "LFC",
      disp_bact_2500 <= 0.5224324 ~ "MFC",
      disp_bact_2500 > 0.5224324 ~ "HFC"
    ),
    Forest_coverage = factor(Forest_coverage, levels = c("LFC", "MFC", "HFC")),
    NB_Category = factor(NB_Category, levels = c("Generalist", "Opportunist", "Specialist"))
  )

# Sum relative abundance within each soil sample for each niche-breadth category
bacti_sum_per_sample <- bacti_pal_graph_long %>%
  group_by(ID, Forest_coverage, NB_Category) %>%
  summarise(
    Sum_Relative_Abundance = sum(Relative_Abundance, na.rm = TRUE),
    .groups = "drop"
  )

# Summary across soil samples for plotting
bacti_pal_graph_summary <- bacti_sum_per_sample %>%
  group_by(NB_Category, Forest_coverage) %>%
  summarise(
    Mean_Sum_RA = mean(Sum_Relative_Abundance, na.rm = TRUE),
    SE = sd(Sum_Relative_Abundance, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    Forest_coverage = factor(Forest_coverage, levels = c("LFC", "MFC", "HFC")),
    NB_Category = factor(NB_Category, levels = c("Generalist", "Opportunist", "Specialist"))
  )

colors_dispersion <- c(
  "LFC" = "#4C8DBC",
  "MFC" = "#F0D62F",
  "HFC" = "#63B255"
)

p_nb_bar_facet <- ggplot(
  bacti_pal_graph_summary,
  aes(x = Forest_coverage, y = Mean_Sum_RA, fill = Forest_coverage)
) +
  geom_col(
    width = 0.88,
    color = "black",
    linewidth = 0.4
  ) +
  geom_errorbar(
    aes(ymin = Mean_Sum_RA - SE, ymax = Mean_Sum_RA + SE),
    width = 0.18,
    linewidth = 0.4
  ) +
  facet_wrap(~ NB_Category, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = colors_dispersion) +
  labs(
    x = "",
    y = "Relative abundance (%)"
  ) +
  theme_classic(base_size = 16) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 13),
    axis.text.x = element_text(face = "bold", color = "black"),
    axis.text.y = element_text(face = "bold", color = "black"),
    axis.title.y = element_text(face = "bold"),
    legend.position = "none",
    panel.spacing = unit(1.2, "lines")
  )

p_nb_bar_facet

ggsave(
  "results_hypothesis3_bac/niche_breadth_barplot_facet.png",
  p_nb_bar_facet,
  width = 8.5,
  height = 3.8,
  dpi = 300
)

# Optional barplot based on sample-level sums
p_nb_bar <- ggplot(
  bacti_pal_graph_summary,
  aes(x = NB_Category, y = Mean_Sum_RA, fill = Forest_coverage)
) +
  geom_col(
    position = position_dodge(width = 0.7),
    width = 0.6,
    color = "black",
    alpha = 0.85
  ) +
  geom_errorbar(
    aes(ymin = Mean_Sum_RA - SE, ymax = Mean_Sum_RA + SE),
    position = position_dodge(width = 0.7),
    width = 0.2,
    color = "black"
  ) +
  scale_fill_manual(values = colors_dispersion) +
  labs(
    x = "",
    y = "Summed relative abundance",
    fill = "Forest cover"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.position = "right"
  )

ggsave(
  "results_hypothesis3_bac/niche_breadth_barplot.png",
  p_nb_bar,
  width = 6.7,
  height = 4,
  dpi = 300
)

REVISAR NICHE BREADTH... USÉ ESAS 4 VARIABLES? HICE LA PARTICIÓN CON K-MEANS??

# =============================================================================
# Prepare per-sample niche-breadth abundances for downstream models
# =============================================================================

RA_by_cat_wide <- bacti_sum_per_sample %>%
  dplyr::select(ID, NB_Category, Sum_Relative_Abundance) %>%
  pivot_wider(
    names_from = NB_Category,
    values_from = Sum_Relative_Abundance,
    values_fill = 0
  )

orchard_summary_df <- orchard_summary_df %>%
  tibble::rownames_to_column("SampleID") %>%
  left_join(
    RA_by_cat_wide %>% rename(SampleID = ID),
    by = "SampleID"
  ) %>%
  tibble::column_to_rownames("SampleID")

# Keep the three forest-cover classes for downstream models
orchard_summary_df <- orchard_summary_df %>%
  mutate(
    FC = case_when(
      P_2500 <= 0.3584376 ~ "LFC",
      P_2500 <= 0.5224324 ~ "MFC",
      P_2500 > 0.5224324 ~ "HFC"
    )
  )

orchard_summary_df$FC <- factor(orchard_summary_df$FC, levels = c("LFC", "MFC", "HFC"))

# Beta-transformed versions for beta regression if needed
n <- nrow(orchard_summary_df)

orchard_summary_df <- orchard_summary_df %>%
  mutate(
    Generalist_beta = (Generalist * (n - 1) + 0.5) / n,
    Opportunist_beta = (Opportunist * (n - 1) + 0.5) / n,
    Specialist_beta = (Specialist * (n - 1) + 0.5) / n
  )

write.csv(
  orchard_summary_df,
  "results_hypothesis3_bac/orchard_summary_with_nb_categories.csv"
)

# =============================================================================
# Beta regression and summary plots
# =============================================================================

# Generalist
mod_lm_Gen <- lm(Generalist ~ FC, data = orchard_summary_df)
mod_invgauss_Gen <- glm(Generalist ~ FC, data = orchard_summary_df, family = inverse.gaussian(link = "log"))
mod_lognorm_Gen <- lm(log(Generalist) ~ FC, data = orchard_summary_df)
mod_beta_Gen <- betareg(Generalist_beta ~ FC, data = orchard_summary_df, link = "logit")

capture.output(compare_performance(mod_lm_Gen, mod_invgauss_Gen, mod_lognorm_Gen, mod_beta_Gen, rank = TRUE, verbose = TRUE),
               file = "results_hypothesis3_bac/compare_models_generalist.txt")
capture.output(emmeans(mod_beta_Gen, pairwise ~ FC, type = "response"),
               file = "results_hypothesis3_bac/posthoc_generalist.txt")

Generalist_summary <- orchard_summary_df %>%
  group_by(FC) %>%
  summarise(
    Relative_Abundance = mean(Generalist, na.rm = TRUE),
    SE = sd(Generalist, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

p_generalist <- ggplot(Generalist_summary, aes(x = FC, y = Relative_Abundance, fill = FC)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.1), width = 0.9, color = "black", alpha = 0.8) +
  geom_errorbar(aes(ymin = Relative_Abundance - SE, ymax = Relative_Abundance + SE),
                position = position_dodge(width = 0.1), width = 0.2, color = "black") +
  scale_fill_manual(values = c("LFC" = "#1f78b4", "MFC" = "gold", "HFC" = "#33a02c")) +
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

ggsave("results_hypothesis3_bac/generalist_fc_barplot.png", p_generalist, width = 2.5, height = 4.5, dpi = 300)

# Opportunist
mod_lm_Opp <- lm(Opportunist ~ FC, data = orchard_summary_df)
mod_invgauss_Opp <- glm(Opportunist ~ FC, data = orchard_summary_df, family = inverse.gaussian(link = "log"))
mod_lognorm_Opp <- lm(log(Opportunist) ~ FC, data = orchard_summary_df)
mod_beta_Opp <- betareg(Opportunist_beta ~ FC, data = orchard_summary_df, link = "logit")

capture.output(compare_performance(mod_lm_Opp, mod_invgauss_Opp, mod_lognorm_Opp, mod_beta_Opp, rank = TRUE, verbose = TRUE),
               file = "results_hypothesis3_bac/compare_models_opportunist.txt")
capture.output(emmeans(mod_beta_Opp, pairwise ~ FC, type = "response"),
               file = "results_hypothesis3_bac/posthoc_opportunist.txt")

Opportunist_summary <- orchard_summary_df %>%
  group_by(FC) %>%
  summarise(
    Relative_Abundance = mean(Opportunist, na.rm = TRUE),
    SE = sd(Opportunist, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

p_opportunist <- ggplot(Opportunist_summary, aes(x = FC, y = Relative_Abundance, fill = FC)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.1), width = 0.9, color = "black", alpha = 0.8) +
  geom_errorbar(aes(ymin = Relative_Abundance - SE, ymax = Relative_Abundance + SE),
                position = position_dodge(width = 0.1), width = 0.2, color = "black") +
  scale_fill_manual(values = c("LFC" = "#1f78b4", "MFC" = "gold", "HFC" = "#33a02c")) +
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

ggsave("results_hypothesis3_bac/opportunist_fc_barplot.png", p_opportunist, width = 2.5, height = 4.5, dpi = 300)

# Specialist
mod_lm_Spe <- lm(Specialist ~ FC, data = orchard_summary_df)
mod_invgauss_Spe <- glm(Specialist ~ FC, data = orchard_summary_df, family = inverse.gaussian(link = "log"))
mod_lognorm_Spe <- lm(log(Specialist) ~ FC, data = orchard_summary_df)
mod_beta_Spe <- betareg(Specialist_beta ~ FC, data = orchard_summary_df, link = "logit")

capture.output(compare_performance(mod_lm_Spe, mod_invgauss_Spe, mod_lognorm_Spe, mod_beta_Spe, rank = TRUE, verbose = TRUE),
               file = "results_hypothesis3_bac/compare_models_specialist.txt")
capture.output(emmeans(mod_beta_Spe, pairwise ~ FC, type = "response"),
               file = "results_hypothesis3_bac/posthoc_specialist.txt")

Specialist_summary <- orchard_summary_df %>%
  group_by(FC) %>%
  summarise(
    Relative_Abundance = mean(Specialist, na.rm = TRUE),
    SE = sd(Specialist, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

p_specialist <- ggplot(Specialist_summary, aes(x = FC, y = Relative_Abundance, fill = FC)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.1), width = 0.9, color = "black", alpha = 0.8) +
  geom_errorbar(aes(ymin = Relative_Abundance - SE, ymax = Relative_Abundance + SE),
                position = position_dodge(width = 0.1), width = 0.2, color = "black") +
  scale_fill_manual(values = c("LFC" = "#1f78b4", "MFC" = "gold", "HFC" = "#33a02c")) +
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

ggsave("results_hypothesis3_bac/specialist_fc_barplot.png", p_specialist, width = 2.5, height = 4.5, dpi = 300)

# =============================================================================
# Environmental distance to forest reference within block
# =============================================================================

forest_df <- as(sample_data(phylo_forest_g_bac), "data.frame") %>%
  rownames_to_column("SampleID") %>%
  filter(Group == "Forest")

orchard_df <- as(sample_data(phylo_crops_g_bac), "data.frame") %>%
  rownames_to_column("SampleID") %>%
  filter(Group == "Orchard")

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

# =============================================================================
# Top genera by niche breadth category
# =============================================================================

top25_by_nb <- bacti_pal_graph_long %>%
  group_by(NB_Category, Genus) %>%
  summarise(total_RA = sum(Relative_Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(NB_Category, desc(total_RA)) %>%
  group_by(NB_Category) %>%
  slice_head(n = 25)

write.csv(top25_by_nb, "results_hypothesis3_bac/top25_genera_by_nb_category.csv")

# =============================================================================
# Clean joined dataframe for Bayesian models
# =============================================================================

pattern <- "^(env_dist(_L1)?|Generalist|Opportunist|Specialist)(\\..*)?$"

nb_df <- orchard_summary_df %>%
  tibble::rownames_to_column("SampleID") %>%
  { dplyr::select(., -dplyr::any_of(grep(pattern, names(.), value = TRUE))) } %>%
  left_join(
    orchard_env_dist %>% dplyr::select(SampleID, env_dist, env_dist_L1),
    by = "SampleID"
  ) %>%
  left_join(
    RA_by_cat_wide %>% rename(SampleID = ID),
    by = "SampleID"
  ) %>%
  tibble::column_to_rownames("SampleID")

nb_bac_df <- nb_df
write.csv(nb_bac_df, "results_hypothesis3_bac/nb_bacterial_model_dataframe.csv")

# =============================================================================
# Bayesian models
# =============================================================================

model_Spe <- brm(
  Specialist ~ t2(P_2500) + (1 | Block),
  data = nb_df,
  family = Beta(),
  chains = 4,
  cores = 12,
  iter = 4000,
  control = list(adapt_delta = 0.99)
)

model_Spe2 <- brm(
  Specialist ~ t2(P_2500, env_dist) + (1 | Block),
  data = nb_df,
  family = Beta(),
  chains = 4,
  cores = 12,
  iter = 4000,
  control = list(adapt_delta = 0.99)
)

model_Spe4 <- brm(
  Specialist ~ t2(P_2500, pH) + (1 | Block),
  data = nb_df,
  family = Beta(),
  chains = 4,
  cores = 12,
  iter = 4000,
  control = list(adapt_delta = 0.99)
)

model_Spe5 <- brm(
  Specialist ~ t2(P_2500, SOM) + (1 | Block),
  data = nb_df,
  family = Beta(),
  chains = 4,
  cores = 12,
  iter = 4000,
  control = list(adapt_delta = 0.99)
)

model_Spe6 <- brm(
  Specialist ~ t2(Block),
  data = nb_df,
  family = Beta(),
  chains = 4,
  cores = 12,
  iter = 4000,
  control = list(adapt_delta = 0.99)
)

model_Opp <- brm(
  Opportunist ~ t2(P_2500) + (1 | Block),
  data = nb_df,
  family = Beta(),
  chains = 4,
  cores = 12,
  iter = 4000,
  control = list(adapt_delta = 0.99)
)

model_Opp2 <- brm(
  Opportunist ~ t2(P_2500, env_dist) + (1 | Block),
  data = nb_df,
  family = Beta(),
  chains = 4,
  cores = 12,
  iter = 4000,
  control = list(adapt_delta = 0.99)
)

model_Opp4 <- brm(
  Opportunist ~ t2(P_2500, pH) + (1 | Block),
  data = nb_df,
  family = Beta(),
  chains = 4,
  cores = 12,
  iter = 4000,
  control = list(adapt_delta = 0.99)
)

model_Opp5 <- brm(
  Opportunist ~ t2(P_2500, SOM) + (1 | Block),
  data = nb_df,
  family = Beta(),
  chains = 4,
  cores = 12,
  iter = 4000,
  control = list(adapt_delta = 0.99)
)

model_Opp6 <- brm(
  Opportunist ~ t2(Block),
  data = nb_df,
  family = Beta(),
  chains = 4,
  cores = 12,
  iter = 4000,
  control = list(adapt_delta = 0.99)
)

model_Gen <- brm(
  Generalist ~ t2(P_2500) + (1 | Block),
  data = nb_df,
  family = Beta(),
  chains = 4,
  cores = 12,
  iter = 4000,
  control = list(adapt_delta = 0.99)
)

model_Gen2 <- brm(
  Generalist ~ t2(P_2500, env_dist) + (1 | Block),
  data = nb_df,
  family = Beta(),
  chains = 4,
  cores = 12,
  iter = 4000,
  control = list(adapt_delta = 0.99)
)

model_Gen4 <- brm(
  Generalist ~ t2(P_2500, pH) + (1 | Block),
  data = nb_df,
  family = Beta(),
  chains = 4,
  cores = 12,
  iter = 4000,
  control = list(adapt_delta = 0.99)
)

model_Gen5 <- brm(
  Generalist ~ t2(P_2500, SOM) + (1 | Block),
  data = nb_df,
  family = Beta(),
  chains = 4,
  cores = 12,
  iter = 4000,
  control = list(adapt_delta = 0.99)
)

model_Gen6 <- brm(
  Generalist ~ t2(Block),
  data = nb_df,
  family = Beta(),
  chains = 4,
  cores = 12,
  iter = 4000,
  control = list(adapt_delta = 0.99)
)

saveRDS(model_Spe,  "results_hypothesis3_bac/model_Spe.rds")
saveRDS(model_Spe2, "results_hypothesis3_bac/model_Spe2.rds")
saveRDS(model_Spe4, "results_hypothesis3_bac/model_Spe4.rds")
saveRDS(model_Spe5, "results_hypothesis3_bac/model_Spe5.rds")
saveRDS(model_Spe6, "results_hypothesis3_bac/model_Spe6.rds")

saveRDS(model_Opp,  "results_hypothesis3_bac/model_Opp.rds")
saveRDS(model_Opp2, "results_hypothesis3_bac/model_Opp2.rds")
saveRDS(model_Opp4, "results_hypothesis3_bac/model_Opp4.rds")
saveRDS(model_Opp5, "results_hypothesis3_bac/model_Opp5.rds")
saveRDS(model_Opp6, "results_hypothesis3_bac/model_Opp6.rds")

saveRDS(model_Gen,  "results_hypothesis3_bac/model_Gen.rds")
saveRDS(model_Gen2, "results_hypothesis3_bac/model_Gen2.rds")
saveRDS(model_Gen4, "results_hypothesis3_bac/model_Gen4.rds")
saveRDS(model_Gen5, "results_hypothesis3_bac/model_Gen5.rds")
saveRDS(model_Gen6, "results_hypothesis3_bac/model_Gen6.rds")

# =============================================================================
# LOO comparison
# =============================================================================

loo_1 <- loo(model_Gen6)
loo_2 <- loo(model_Gen)
loo_3 <- loo(model_Gen2)
loo_4 <- loo(model_Gen4)
loo_5 <- loo(model_Gen5)

loo_6 <- loo(model_Opp6)
loo_7 <- loo(model_Opp)
loo_8 <- loo(model_Opp2)
loo_9 <- loo(model_Opp4)
loo_10 <- loo(model_Opp5)

loo_11 <- loo(model_Spe6)
loo_12 <- loo(model_Spe)
loo_13 <- loo(model_Spe2)
loo_14 <- loo(model_Spe4)
loo_15 <- loo(model_Spe5)

capture.output(loo_1, file = "results_hypothesis3_bac/loo_Gen6.txt")
capture.output(loo_2, file = "results_hypothesis3_bac/loo_Gen.txt")
capture.output(loo_3, file = "results_hypothesis3_bac/loo_Gen2.txt")
capture.output(loo_4, file = "results_hypothesis3_bac/loo_Gen4.txt")
capture.output(loo_5, file = "results_hypothesis3_bac/loo_Gen5.txt")

capture.output(loo_6, file = "results_hypothesis3_bac/loo_Opp6.txt")
capture.output(loo_7, file = "results_hypothesis3_bac/loo_Opp.txt")
capture.output(loo_8, file = "results_hypothesis3_bac/loo_Opp2.txt")
capture.output(loo_9, file = "results_hypothesis3_bac/loo_Opp4.txt")
capture.output(loo_10, file = "results_hypothesis3_bac/loo_Opp5.txt")

capture.output(loo_11, file = "results_hypothesis3_bac/loo_Spe6.txt")
capture.output(loo_12, file = "results_hypothesis3_bac/loo_Spe.txt")
capture.output(loo_13, file = "results_hypothesis3_bac/loo_Spe2.txt")
capture.output(loo_14, file = "results_hypothesis3_bac/loo_Spe4.txt")
capture.output(loo_15, file = "results_hypothesis3_bac/loo_Spe5.txt")

# =============================================================================
# Save session information
# =============================================================================

writeLines(
  capture.output(sessionInfo()),
  "results_hypothesis3_bac/sessionInfo_hypothesis3_bac.txt"
)

































































#

#ver cuanto representan estos generos según categoría en forest
bacti_pal_graph_long$Genus
bacti_pal_graph_long$NB_Category
sample_data(phylo_crops_g_bac)
tax_table(phylo_forest_g_bac)

library(dplyr)
library(ggplot2)

# Resumen de abundancia relativa total por categoría NB y bloque
niche_summary_forest <- bacti_pal_graph_long %>%
  group_by(Block_crop, NB_Category) %>%
  summarise(Total_RA = sum(RelAb_beta, na.rm = TRUE)) %>%
  ungroup()

# Ver tabla
print(niche_summary_forest)

#  Gráfico de barras por bloque
ggplot(niche_summary_forest, aes(x = factor(Block_crop), y = Total_RA, fill = NB_Category)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  labs(
    x = "Block (forest site)",
    y = "Total relative abundance (%)",
    fill = "Niche breadth category"
  ) +
  theme_minimal(base_size = 14)

#bacterias compartidos desde forest 583*509

# Convertir la tabla de conteos a abundancias relativas
phylo_forest_g_bac_RA <- transform_sample_counts(phylo_forest_g_bac, function(x) x / sum(x))

# Luego repetir los pasos
df_genus_RA <- psmelt(tax_glom(phylo_forest_g_bac_RA, taxrank = "Genus"))

top10_forest_RA <- df_genus_RA %>%
  group_by(Genus) %>%
  summarise(Total_RA = sum(Abundance, na.rm = TRUE)) %>%
  arrange(desc(Total_RA)) %>%
  slice_head(n = 13)

# Visualizar
print(top10_forest_RA)





#Ahora iCAMP


#Ahora... puedo sacar q es deterministico y q estocástico?

#wd="d:\\Users\\pedro\\Documents\\microbiome_orchard_forest"
# the folder to save the output. please change to a new folder even if you are just testing the example data.
#setwd("C:\\Users\\pedro\\OneDrive - Universidad Técnica Federico Santa María\\Documents\\microbiome_orchard_forest")
#save.wd="d:\\Users\\pedro\\Documents\\microbiome_orchard_forest\\assembly_ITS"
save.wd="C:\\Users\\pedro\\Research\\microbiome_orchard_forest\\Sci Adv - v2\\figures\\assembly_16S"
if(!dir.exists(save.wd)){dir.create(save.wd)}

#load previous progress (unsaved environment)
load("env.RData", envir = nuevo_env)


# 2 # key parameter setting
prefix="orchard_forest"  # prefix of the output file names. usually use a project ID.
rand.time=100  # randomization time, 1000 is usually enough. For example test, you may set as 100 or less to save time.
nworker=12 # nworker is thread number for parallel computing, which depends on the CPU core number of your computer.
memory.G=50 # to set the memory size as you need (but should be less than the available space in your hard disk), so that calculation of large tree will not be limited by physical memory. unit is Gb.

tree_16S <- read.tree("tree_16S.nwk")
#tree_ITS <- read.tree("tree_ITS.nwk")
#plot(tree_16S)
sample_data(orchard_phy)
#árbol filogenético
#tree_16S # ready pero no está rarefied.. es lo q hay #son 4970 ASVs
tree = phy_tree(orchard_phy)
tree
#tree<-tree_16S

#comm
dim(abund) #130 muestras x 4970 ASVs
data_phylo_filt
orchard_phy
orchard_asv <- as(otu_table(orchard_phy), "matrix")
orchard_asv <- as.data.frame(orchard_asv) %>%
  tibble::rownames_to_column("ASV")
rownames(orchard_asv) <- orchard_asv$ASV 
orchard_asv[1:5,1:5] #cumple con el formato
comm<-orchard_asv[,-1]
comm[1:5,1:5]
dim(comm)
#
orchard_tax <- as(tax_table(orchard_phy), "matrix")
orchard_tax <- as.data.frame(orchard_tax) 
orchard_tax[1:5,1:7]
dim(orchard_tax)
class<-orchard_tax
class
# class <- as.data.frame(lapply(orchard_tax, function(x) gsub("^.__", "", x)))
# rownames(class)
# rownames(class) <- colnames(comm)
# class[1:2,1:5]
# colnames(comm)


#treat
orchard_sam <- as(sample_data(orchard_phy), "matrix")
orchard_sam <- as.data.frame(orchard_sam) %>%
  tibble::rownames_to_column("ASV")
orchard_sam$P_2500

# Cortes por terciles
P2500 <- orchard_sam$P_2500

# convertir a numérico de forma segura
P2500_num <- as.numeric(P2500)

# (opcional) ver cuántos NA aparecen por conversión
sum(is.na(P2500_num))

#cuts <- quantile(P2500_num, probs = c(1/3, 2/3), na.rm = TRUE)
#cuts
quantiles
orchard_sam$Categoria <- cut(
  mapping_filtrado$P_2500,
  breaks = c(-Inf, quantiles[2], quantiles[3], Inf),
  labels = c("L", "M", "H"),
  include.lowest = TRUE
)
orchard_sam$Categoria
orchard_sam2<- orchard_sam[,c("pH","SOM", "N", "claysilt","Categoria")]
rownames(orchard_sam2)<-orchard_asv$ASV
orchard_sam2
sample_data(orchard_phy)
env0<-as.data.frame(orchard_sam2[,-5])
str(env0)
env <- env0[, c("pH","SOM","N","claysilt"), drop = FALSE]
env <- as.data.frame(lapply(env0, function(x) as.numeric(as.character(x))))
rownames(env) <- rownames (env0)
env
treat<-as.data.frame(orchard_sam2[,c(5)])
rownames(treat)<-rownames(env)
colnames(treat) <- "cat"
treat

#listo los archivos
comm [1:20,1:50] #bien
tree #bien
class[1:2,1:7] #bien
treat[1:10,] # bien
env[1:4,1:4] #bien

# 4 # match sample IDs in OTU table and treatment information table
sampid.check=match.name(rn.list=list(comm=comm,treat=treat,env=env)) #All match very well.

# sampid.check=match.name(rn.list=list(comm=comm,treat=treat)) # if you do not have env.file
# for the example data, the output should be "All match very well".
# for your data files, if you have not matched their IDs, the unmatched samples will be removed.
treat=sampid.check$treat
comm=sampid.check$comm
comm=comm[,colSums(comm)>0,drop=FALSE] # if some unmatched samples were removed, some OTUs may become ghosts, then you may use this line to remove them if necessary.
env=sampid.check$env # skip this if you do not have env.file

# 5 # match OTU IDs in OTU table and tree file
spid.check=match.name(cn.list=list(comm=comm),rn.list=list(class=class),tree.list=list(tree=tree))
# for the example data, the output should be "All match very well".
# for your data files, if you have not matched the IDs before, the unmatched OTUs will be removed.
comm=spid.check$comm
clas=spid.check$class
tree=spid.check$tree

# 6 # calculate pairwise phylogenetic distance matrix.
# since microbial community data usually has a large number of species (OTUs or ASVs), we use "big.matrix" in R package "bigmemory" to handle the large phylogenetic distance matrix. 
setwd(save.wd)
if(!file.exists("pd.desc")) 
{
  pd.big=iCAMP::pdist.big(tree = tree, wd=save.wd, nworker = nworker, memory.G = memory.G)
  # output files:
  # path.rda: a R object to list all the nodes and  edge lengthes from root to every tip. saved in R data format. an intermediate output when claculating phylogenetic distance matrix.
  # pd.bin: BIN file (backingfile) generated by function big.matrix in R package bigmemory. This is the big matrix storing pairwise phylogenetic distance values. By using this bigmemory format file, we will not need memory but hard disk when calling big matrix for calculation.
  # pd.desc: the DESC file (descriptorfile) to hold the backingfile (pd.bin) description.
  # pd.taxon.name.csv: comma delimited csv file storing the IDs of tree tips (OTUs), serving as the row/column names of the big phylogenetic distance matrix.
}else{
  # if you already calculated the phylogenetic distance matrix in a previous run
  pd.big=list()
  pd.big$tip.label=read.csv(paste0(save.wd,"/pd.taxon.name.csv"),row.names = 1,stringsAsFactors = FALSE)[,1]
  pd.big$pd.wd=save.wd
  pd.big$pd.file="pd.desc"
  pd.big$pd.name.file="pd.taxon.name.csv"
}

# 7 # assess niche preference difference between species
# env is required for this step.
# since microbial community data usually has a large number of species (OTUs or ASVs), we use "big.matrix" in R package "bigmemory" to handle the large niche difference matrix. 
setwd(save.wd)
niche.dif=iCAMP::dniche(env = env,comm = comm,method = "niche.value",
                        nworker = nworker,out.dist=FALSE,bigmemo=TRUE,
                        nd.wd=save.wd)

env
# 8 # within-bin phylogenetic signal assessment.
# For real data, you may try several different settings of binning, and choose the one leading to the best within-bin phylogenetic signal.
# env is required for this step.
# 8.1 # try phylogenetic binning using current setttings.
ds = 0.2 # setting can be changed to explore the best choice
bin.size.limit = 24 # setting can be changed to explore the best choice. # here set as 5 just for the small example dataset. For real data, usually try 12 to 48.

# The tree for taxa.binphy.big must be a rooted tree.
if(!ape::is.rooted(tree))
{
  tree.rt=iCAMP::midpoint.root.big(tree = tree, pd.desc = pd.big$pd.file,
                                   pd.spname = pd.big$tip.label,pd.wd = pd.big$pd.wd,
                                   nworker = nworker)
  tree=tree.rt$tree
}
phylobin=taxa.binphy.big(tree = tree, pd.desc = pd.big$pd.file,pd.spname = pd.big$tip.label,
                         pd.wd = pd.big$pd.wd, ds = ds, bin.size.limit = bin.size.limit,
                         nworker = nworker)
# 8.2 # test within-bin phylogenetic signal.
sp.bin=phylobin$sp.bin[,3,drop=FALSE]
dim(sp.bin)
comm[1:2,1:2]
head(rowSums(comm))
sp.ra=colMeans(comm/rowSums(comm))
length(sp.ra)
abcut=3 # you may remove some species, if they are too rare to perform reliable correlation test.
commc=comm[,colSums(comm)>=abcut,drop=FALSE]
dim(commc) 
spname.use=colnames(commc)

#el filtrado me dejó distintas dimensiones así que agrego estos pasos:
sp.bin <- sp.bin[rownames(sp.bin) %in% colnames(commc), , drop = FALSE]
sp.ra <- sp.ra[names(sp.ra) %in% colnames(commc)] #names xq es vector

dim(sp.bin)    # Debe ser [3957, 1]
length(sp.ra)  # Debe ser 3957

#length(niche.dif$names)

#ahora revisión de orden
all(spname.use == rownames(sp.bin))   # Debe devolver TRUE
all(spname.use == names(sp.ra))       # Debe devolver TRUE
all(spname.use == rownames(pd.big))   # Debe devolver TRUE
#sp.bin dio falso, así que
sp.bin <- sp.bin[spname.use, , drop = FALSE]

#revisión de NAs
anyNA(sp.bin)  # Revisa sp.bin
anyNA(sp.ra)   # Revisa sp.ra
anyNA(spname.use)  # Revisa spname.use

class(sp.bin)  # Debe ser "matrix" o "data.frame"
class(sp.ra)   # Debe ser "numeric" o "vector"
class(spname.use) # Debe ser "character" o "factor"

dim(sp.bin) #3957
length(sp.ra) #3957
length(spname.use) #3957
length(pd.big$tip.label) #8226, el árbol no está filtrado

# Filtrar el árbol para incluir solo especies presentes en spname.use
pd.big$tip.label <- pd.big$tip.label[pd.big$tip.label %in% spname.use]
length(pd.big$tip.label)#3957
niche.dif$names <- niche.dif$names[niche.dif$names %in% spname.use]
length(niche.dif$names) #8047
#ps.bin = Use Mantel test to evaluate phylogenetic signal within each bin, i.e. correlation between phylogenetic distance and niche difference.
binps=iCAMP::ps.bin(sp.bin = sp.bin,sp.ra = sp.ra,spname.use = spname.use,
                    pd.desc = pd.big$pd.file, pd.spname = pd.big$tip.label, pd.wd = pd.big$pd.wd,
                    nd.list = niche.dif$nd,nd.spname = niche.dif$names,ndbig.wd = niche.dif$nd.wd,
                    cor.method = "spearman",r.cut = 0.1, p.cut = 0.05, min.spn = 5)

if(file.exists(paste0(prefix,".PhyloSignalSummary.csv"))){appendy=TRUE;col.namesy=FALSE}else{appendy=FALSE;col.namesy=TRUE}
write.table(data.frame(ds=ds,n.min=bin.size.limit,binps$Index),file = paste0(prefix,".PhyloSignalSummary.csv"),
            append = appendy, quote=FALSE, sep=",", row.names = FALSE,col.names = col.namesy)
if(file.exists(paste0(prefix,".PhyloSignalDetail.csv"))){appendy2=TRUE;col.namesy2=FALSE}else{appendy2=FALSE;col.namesy2=TRUE}
write.table(data.frame(ds=ds,n.min=bin.size.limit,binID=rownames(binps$detail),binps$detail),file = paste0(prefix,".PhyloSignalDetail.csv"),
            append = appendy2, quote = FALSE, sep = ",", row.names = FALSE, col.names = col.namesy2)
# since this example small data is randomly generated, the correlation should be very weak.
# usually, you are looking for a binning setting lead to higher RAsig.abj (relative abundance of bins with significant phylogenetic signal) and relative high meanR (mean correlation coefficient across bins).
# see help document of the function "ps.bin" for the meaning of output.

#################### De aquí en adelante usa tree entonces cargo de nuevo 8.2
# 9 # iCAMP analysis
# 9.1 # without omitting small bins.
# commonly use # set sig.index as Confidence instead of SES.RC (betaNRI/NTI + RCbray)
bin.size.limit = 24 # For real data, usually use a proper number according to phylogenetic signal test or try some settings then choose the reasonable stochasticity level. our experience is 12, or 24, or 48. but for this example dataset which is too small, have to use 5.
sig.index="Confidence" # see other options in help document of icamp.big.
icres=iCAMP::icamp.big(comm=comm, pd.desc = pd.big$pd.file, pd.spname=pd.big$tip.label,
                       pd.wd = pd.big$pd.wd, rand = rand.time, tree=tree,
                       prefix = prefix, ds = 0.2, pd.cut = NA, sp.check = TRUE,
                       phylo.rand.scale = "within.bin", taxa.rand.scale = "across.all",
                       phylo.metric = "bMPD", sig.index=sig.index, bin.size.limit = bin.size.limit, 
                       nworker = nworker, memory.G = memory.G, rtree.save = FALSE, detail.save = TRUE, 
                       qp.save = FALSE, detail.null = FALSE, ignore.zero = TRUE, output.wd = save.wd, 
                       correct.special = TRUE, unit.sum = rowSums(comm), special.method = "depend",
                       ses.cut = 1.96, rc.cut = 0.95, conf.cut=0.975, omit.option = "no",meta.ab = NULL)
# there are quite a few parameters in this function, please check the help document of "icamp.big".
# output files:
# Test.iCAMP.detail.rda: the object "icres" saved in R data format. it is a list object. The first element bNRIiRCa is the result of relative importance of each assembly process in each pairwise comparison (each turnover). The second element "detail" including binning information (named taxabin), phylogenetic and taxonomic metrics results in each bin (named like bNRIi, RCa, etc.), relative abundance of each bin (bin.weight), relative importance of each process in each turnover between communities (processes), input settings (setting), and input community data matrix (comm). See help document of the function icamp.big for more details.
str(icres)
icres$CbMPDiCBraya
icres$detail$bin.weight #no se el output, seguiré

############################
# 9.2 to 9.4 are some optional special settings you may explore.
# 9.2 # explore different ways for null model significance test.
# 9.2.1 # set detail.null=TRUE, output all null values, to facilitate normality test and switch between different options
detail.null=TRUE
bin.size.limit = 24 
sig.index="SES.RC" # this is traditional way, with assumption that null values of phylogenetic metrics follow normal distribution. 
prefixb="TestB"

icres2=iCAMP::icamp.big(comm=comm, pd.desc = pd.big$pd.file, pd.spname=pd.big$tip.label,
                        pd.wd = pd.big$pd.wd, rand = rand.time, tree=tree,
                        prefix = prefixb, ds = 0.2, pd.cut = NA, sp.check = TRUE,
                        phylo.rand.scale = "within.bin", taxa.rand.scale = "across.all",
                        phylo.metric = "bMPD", sig.index=sig.index, bin.size.limit = bin.size.limit, 
                        nworker = nworker, memory.G = memory.G, rtree.save = FALSE, detail.save = TRUE, 
                        qp.save = FALSE, detail.null = detail.null, ignore.zero = TRUE, output.wd = save.wd, 
                        correct.special = TRUE, unit.sum = rowSums(comm), special.method = "depend",
                        ses.cut = 1.96, rc.cut = 0.95, conf.cut=0.975, omit.option = "no",meta.ab = NULL)

# 9.2.2 # normality test
nntest=iCAMP::null.norm(icamp.output=icres2, p.norm.cut=0.05, detail.out=FALSE)
# output shows non-normal distribution ratio in each bin, i.e. the proportion of turnovers which have null values significantly deviated from normal distribution.
# if some ratio values are very high, may need to change to use "Confidence" as sig.index.
icres2$bNRIiRCa

# 9.2.3 # change sig.index to "Confidence".
icres3=iCAMP::change.sigindex(icamp.output = icres2, sig.index = "Confidence", detail.save = TRUE, detail.null = FALSE, conf.cut = 0.975)
head(icres3$CbMPDiCBraya)

# 9.2.4 # change sig.index to "RC" for both phylogenetic and taxonomic metrics.
icres4=iCAMP::change.sigindex(icamp.output = icres2, sig.index = "RC", detail.save = TRUE, detail.null = FALSE, rc.cut = 0.95)
head(icres4$RCbMPDiRCbraya)

# 9.2.5 # the function can also change the significance threshold.
icres5=iCAMP::change.sigindex(icamp.output = icres2, sig.index = "SES.RC", detail.save = TRUE, detail.null = FALSE, ses.cut = 1.64, rc.cut = 0.9)
head(icres5$bNRIiRCbraya)


# 9.3 # Special settings for regional pool(s) (the metacommunity)
# 9.3.1 # you may specify the relative abundance of each species in the regional pool, if it is not the same as the average relative abundance from the "comm" you input.
meta.ab=rep(1,ncol(comm)) # here i assume all the species actuall have the same relative abundance in the regional pool.
prefix2=paste0(prefix,".MetaCrct")
sig.index="Confidence" # see other options in help document of icamp.big.
dim(comm)
length(pd.big$tip.label) #lo habia cortado pero volví a 6 para volver al original sin filtrar
tree
icres.meta=iCAMP::icamp.big(comm=comm, pd.desc = pd.big$pd.file, pd.spname=pd.big$tip.label,
                            pd.wd = pd.big$pd.wd, rand = rand.time, tree=tree,
                            prefix = prefix2, ds = 0.2, pd.cut = NA, sp.check = TRUE,
                            phylo.rand.scale = "within.bin", taxa.rand.scale = "across.all",
                            phylo.metric = "bMPD", sig.index=sig.index, bin.size.limit = bin.size.limit, 
                            nworker = nworker, memory.G = memory.G, rtree.save = FALSE, detail.save = TRUE, 
                            qp.save = FALSE, detail.null = FALSE, ignore.zero = TRUE, output.wd = save.wd, 
                            correct.special = TRUE, unit.sum = rowSums(comm), special.method = "depend",
                            ses.cut = 1.96, rc.cut = 0.95, conf.cut=0.975, omit.option = "no",meta.ab=meta.ab)


# 9.3.2 # if samples are from different regional pools (metacommunities)
# e.g., in this example, samples in 'south' are from one regional pool and those in 'north' are from another regional pool
sp_comm <- rownames(icres$detail$comm)
sp_comm
# Especies en meta.frequence
sp_meta <- rownames(icres$detail$meta.frequence)
sp_meta

names(icbin$detail)


treat
meta.group=treat[,1,drop=FALSE]
prefix3=paste0(prefix,".MetaGrp")
sig.index="Confidence"
icres.meta2=iCAMP::icamp.cm(comm=comm,tree=tree,meta.group=meta.group,
                            pd.desc=pd.big$pd.file, pd.spname=pd.big$tip.label, pd.wd=pd.big$pd.wd,
                            rand=rand.time,prefix=prefix3,ds=0.2,pd.cut=NA,
                            phylo.rand.scale="within.bin",taxa.rand.scale="across.all",
                            phylo.metric="bMPD",sig.index=sig.index,
                            bin.size.limit=bin.size.limit,nworker=nworker,memory.G=memory.G,
                            rtree.save=FALSE,detail.save=TRUE,qp.save=FALSE,detail.null=FALSE,
                            ignore.zero=TRUE,output.wd=save.wd,
                            correct.special=TRUE,unit.sum=rowSums(comm),
                            special.method="depend", ses.cut = 1.96,rc.cut = 0.95,conf.cut=0.975,
                            omit.option="no")
icres.meta2$CbMPDiCBraya

# 9.4 # consider to omit small bins
# 9.4.1 # if you would like to omit small bins rather than merging them to nearest relatives, set omit.option as "test" to check what will be omitted.
omit.option = "test"
icres.omit=iCAMP::icamp.big(comm=comm, pd.desc = pd.big$pd.file, pd.spname=pd.big$tip.label,
                            pd.wd = pd.big$pd.wd, rand = rand.time, tree=tree,
                            prefix = prefix, ds = 0.2, pd.cut = NA, sp.check = TRUE,
                            phylo.rand.scale = "within.bin", taxa.rand.scale = "across.all",
                            phylo.metric = "bMPD", sig.index=sig.index, bin.size.limit = bin.size.limit,
                            nworker = nworker, memory.G = memory.G, rtree.save = FALSE, detail.save = TRUE,
                            qp.save = FALSE, detail.null=FALSE, ignore.zero = TRUE, output.wd = save.wd,
                            correct.special = TRUE, unit.sum = rowSums(comm), special.method = "depend",
                            ses.cut = 1.96, rc.cut = 0.95, conf.cut=0.975, omit.option = omit.option)
# # "test" will return a detailed table of omitted species.

# 9.4.2 # then set it as "omit" to omit the small bins. # DE AQUI EN ADELANTE NO ME CORRIÓ
omit.option = "omit"
icres.omit2=iCAMP::icamp.big(comm=comm, pd.desc = pd.big$pd.file, pd.spname=pd.big$tip.label,
                             pd.wd = pd.big$pd.wd, rand = rand.time, tree=tree,
                             prefix = prefix, ds = 0.2, pd.cut = NA, sp.check = TRUE,
                             phylo.rand.scale = "within.bin", taxa.rand.scale = "across.all",
                             phylo.metric = "bMPD", sig.index=sig.index, bin.size.limit = bin.size.limit,
                             nworker = nworker, memory.G = memory.G, rtree.save = FALSE, detail.save = TRUE,
                             qp.save = FALSE, detail.null=FALSE, ignore.zero = TRUE, output.wd = save.wd,
                             correct.special = TRUE, unit.sum = rowSums(comm), special.method = "depend",
                             ses.cut = 1.96, rc.cut = 0.95, conf.cut=0.975, omit.option = omit.option)
# In this simple example, since all bins are small, "omit" should return an error. In real data, this will go ahead to do iCAMP analysis with the strict bins which are big enough (>bin.size.limit).


# 9.5 # input community matrix as relative abundances (values < 1) rather than counts
comra=comm/rowSums(comm)
anyNA(comra)  # Devuelve TRUE si hay NA en la matriz
comra[1:10,1:5]
prefixra=paste0(prefix,"RA")
bin.size.limit = 24 # For real data, usually use a proper number according to phylogenetic signal test or try some settings then choose the reasonable stochasticity level. our experience is 12, or 24, or 48. but for this example dataset which is too small, have to use 5.
dim(comra)
length(tree$tip.label)
icres6=iCAMP::icamp.big(comm=comra,tree=tree,pd.desc=pd.big$pd.file, pd.spname=pd.big$tip.label, pd.wd=pd.big$pd.wd,
                        rand=rand.time,prefix=prefixra,ds=0.2,pd.cut=NA,sp.check=TRUE,
                        phylo.rand.scale="within.bin",taxa.rand.scale="across.all",
                        phylo.metric="bMPD",sig.index="Confidence",
                        bin.size.limit=bin.size.limit,nworker=nworker,memory.G=memory.G,
                        rtree.save=FALSE,detail.save=TRUE,qp.save=FALSE,detail.null=FALSE,
                        ignore.zero=TRUE,output.wd=save.wd,correct.special=TRUE,unit.sum=rowSums(comra),
                        special.method="depend",ses.cut = 1.96,rc.cut = 0.95,conf.cut=0.975,
                        omit.option="no",meta.ab=NULL, taxo.metric="bray", transform.method=NULL,
                        logbase=2, dirichlet=TRUE)


# 9.6 # community data transformation and taxonomic dissimilarity index change
taxo.metric='euclidean'
transform.method='hellinger'
prefixtran=paste0(prefix,"Hel")
bin.size.limit = 5 # For real data, usually use a proper number according to phylogenetic signal test or try some settings then choose the reasonable stochasticity level. our experience is 12, or 24, or 48. but for this example dataset which is too small, have to use 5.
icres7=iCAMP::icamp.big(comm=comm,tree=tree,pd.desc=pd.big$pd.file, pd.spname=pd.big$tip.label, pd.wd=pd.big$pd.wd,
                        rand=rand.time,prefix=prefixtran,ds=0.2,pd.cut=NA,sp.check=TRUE,
                        phylo.rand.scale="within.bin",taxa.rand.scale="across.all",
                        phylo.metric="bMPD",sig.index="Confidence",
                        bin.size.limit=bin.size.limit,nworker=nworker,memory.G=memory.G,
                        rtree.save=FALSE,detail.save=TRUE,qp.save=FALSE,detail.null=FALSE,
                        ignore.zero=TRUE,output.wd=save.wd,correct.special=TRUE,unit.sum=rowSums(comra),
                        special.method="depend",ses.cut = 1.96,rc.cut = 0.95,conf.cut=0.975,
                        omit.option="no",meta.ab=NULL, taxo.metric=taxo.metric, transform.method=transform.method,
                        logbase=2, dirichlet=FALSE)

###############################
# 10 # iCAMP bin level statistics


# Cuáles están en comm pero no en meta.frequence
setdiff(sp_comm, sp_meta)


# Check if treatments align with samples
icbin=iCAMP::icamp.bins(icamp.detail = icres$detail,treat = treat,
                        clas=clas,silent=FALSE, boot = TRUE,
                        rand.time = rand.time,between.group = TRUE)
save(icbin,file = paste0(prefix,".iCAMP.Summary.rda")) # just to archive the result. rda file is automatically compressed, and easy to load into R.
write.csv(icbin$Pt,file = paste0(prefix,".ProcessImportance_EachGroup.csv"),row.names = FALSE)
write.csv(icbin$Ptk,file = paste0(prefix,".ProcessImportance_EachBin_EachGroup.csv"),row.names = FALSE)
write.csv(icbin$Ptuv,file = paste0(prefix,".ProcessImportance_EachTurnover.csv"),row.names = FALSE)
write.csv(icbin$BPtk,file = paste0(prefix,".BinContributeToProcess_EachGroup.csv"),row.names = FALSE)
write.csv(data.frame(ID=rownames(icbin$Class.Bin),icbin$Class.Bin,stringsAsFactors = FALSE),
          file = paste0(prefix,".Taxon_Bin.csv"),row.names = FALSE)
write.csv(icbin$Bin.TopClass,file = paste0(prefix,".Bin_TopTaxon.csv"),row.names = FALSE)

icbin$BPtk

# output files:
# Test.iCAMP.Summary.rda: the object "icbin" saved in R data format. see help document of the function icamp.bins for description of each element in the object.
# Test.ProcessImportance_EachGroup.csv: Relative importance of each process in governing the turnovers in a group of samples.
# Test.ProcessImportance_EachBin_EachGroup.csv: Relative importance of each process in governing the turnovers of each bin among a group of samples.
# Test.ProcessImportance_EachTurnover.csv: Relative importance of each process in governing the turnovers between each pair of communities (samples).
# Test.BinContributeToProcess_EachGroup.csv: Bin contribution to each process, measuring the contribution of each bin to the relative importance of each process in the assembly of a group of communities.
# Test.Taxon_Bin.csv: a matrix showing the bin ID and classification information for each taxon.
# Test.Bin_TopTaxon.csv: a matrix showing the bin relative abundance; the top taxon ID, percentage in bin, and classification; the most abundant name at each phylogeny level in the bin.


#pal graph
# Example data (replace this with icbin$BPtk)
data <- icbin$BPtk
data
library(dplyr)

# Calculate the average of all bins for each row
data_with_avg <- icbin$BPtk %>%
  mutate(Average_Bins = rowMeans(dplyr::select(., starts_with("bin")), na.rm = TRUE))%>%
  filter(Group != c("H_vs_L", "M_vs_H", "M_vs_L"))

# View the resulting data
head(data_with_avg)

# Plot the data
library(ggplot2)
ggplot(data_with_avg, aes(x = Process, y = Average_Bins, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(
    title = "Comparison of Processes Between Low (L) and High(H)",
    x = "Process",
    y = "Average Value",
    fill = "Group"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#MEJOR UN PIE PLOT

library(dplyr)
library(ggplot2)

# Step 1: Calculate relative abundance for each process in each group
data_pie <- data_with_avg %>%
  group_by(Group) %>%  # Group by Group (L and H)
  mutate(Relative_Abundance = Average_Bins / sum(Average_Bins) * 100) %>%  # Calculate relative abundance (%)
  ungroup()
library(writexl)
write_xlsx(data_pie, "data_pie_16S.xlsx")



# 11 # Bootstrapping test
# please specify which column in the treatment information table.
i=1
treat.use=treat[,i,drop=FALSE]
icamp.result=icres$CbMPDiCBraya
icboot=iCAMP::icamp.boot(icamp.result = icamp.result,treat = treat.use,rand.time = rand.time,
                         compare = TRUE,silent = FALSE,between.group = TRUE,ST.estimation = TRUE)
save(icboot,file=paste0(prefix,".iCAMP.Boot.",colnames(treat)[i],".rda"))
write.csv(icboot$summary,file = paste0(prefix,".iCAMP.BootSummary.",colnames(treat)[i],".csv"),row.names = FALSE)
write.csv(icboot$compare,file = paste0(prefix,".iCAMP.Compare.",colnames(treat)[i],".csv"),row.names = FALSE)
icboot$summary
icboot$compare
# output files:
# Test.iCAMP.Boot.Management.rda: the object "icboot" saved in R data format. see help document of the function icamp.boot for description of each element in the object.
# Test.BootSummary.Management.csv: a table to summarize bootstrapping results. see help document of the function icamp.boot for description of the output element "summary".
# Test.Compare.Management.csv: a table to summarize comparison index, effect size, and significance between each two groups. see help document of the function icamp.boot for description of the output element "compare".

# 12 # Other approach: QPEN (quantifying community assembly processes based on entire-community null model analysis)
# 12.1 # QPEN calculation
qpout=iCAMP::qpen(comm=comm,pd=pd.big$pd.file,pd.big.wd=pd.big$pd.wd,
                  pd.big.spname=pd.big$tip.label,ab.weight=TRUE,
                  rand.time=rand.time, nworker=nworker,project=prefix,
                  wd=save.wd, save.bNTIRC=TRUE)
# 12.2 # significance test
qptest=qpen.test(qpen.result = qpout,treat = treat,rand.time = rand.time,
                 between.group = TRUE,out.detail=TRUE,silent=FALSE)
write.csv(qptest$obs.summary,file = paste0(prefix,".QPEN.Index.Obs.Summary.csv"),row.names = FALSE)
write.csv(qptest$boot.summary,file = paste0(prefix,".QPEN.Bootstrapping.Summary.csv"),row.names = FALSE)
write.csv(qptest$compare,file = paste0(prefix,".QPEN.Comparison.Summary.csv"),row.names = FALSE)
save(qptest,file = paste0(prefix,".QPEN.bootstrap.rda"))

# 13 # Other approach: Neutral taxa percentage
snmout=iCAMP::snm.comm(comm = comm, treat = treat, 
                       rand=rand.time, alpha=0.05)
write.csv(snmout$stats,file = paste0(prefix,".NeutralModel.Stats.csv"))
write.csv(snmout$ratio.summary,file = paste0(prefix,".NeutralModel.TypeRatio.csv"))

# 14 # Other approach: tNST and pNST (taxonomic and phylogenetic normalized stochasticity ratio)
# need to install package NST if not yet
if(!("NST" %in% installed.packages()[,"Package"])){install.packages("NST")}
library(NST)
i=1
treat.use=treat[,i,drop=FALSE]

# 14.1a # tNST
tnstout=NST::tNST(comm=comm, group=treat.use, dist.method="bray", 
                  abundance.weighted=TRUE, rand=rand.time,  
                  nworker=nworker, null.model="PF", output.rand = TRUE,
                  SES = TRUE, RC = TRUE)
write.csv(tnstout$index.grp,file = paste0(prefix,".tNST.summary.",colnames(treat)[i],".csv"))
write.csv(tnstout$index.pair.grp,file = paste0(prefix,".tNST.pairwise.",colnames(treat)[i],".csv"))

# 14.1b # bootstrapping test for tNST
tnst.bt=NST::nst.boot(nst.result=tnstout, group=treat.use,
                      rand=rand.time, nworker=nworker)
write.csv(tnst.bt$NST.summary,file = paste0(prefix,".tNST.bootstr.",colnames(treat)[i],".csv"))
write.csv(tnst.bt$NST.compare,file = paste0(prefix,".tNST.compare.",colnames(treat)[i],".csv"))

# 14.2a # pNST
pnstout=NST::pNST(comm=comm, pd.desc=pd.big$pd.file, pd.wd=pd.big$pd.wd, 
                  pd.spname=pd.big$tip.label, group=treat.use, abundance.weighted=TRUE,
                  rand=rand.time, phylo.shuffle=TRUE, nworker=nworker,
                  output.rand = TRUE, SES=FALSE, RC=FALSE)
write.csv(pnstout$index.grp,file = paste0(prefix,".pNST.summary.",colnames(treat)[i],".csv"))
write.csv(pnstout$index.pair.grp,file = paste0(prefix,".pNST.pairwise.",colnames(treat)[i],".csv"))

pnst.bt=NST::nst.boot(nst.result=pnstout, group=treat.use,
                      rand=rand.time, nworker=nworker)
write.csv(pnst.bt$NST.summary,file = paste0(prefix,".pNST.bootstr.",colnames(treat)[i],".csv"))
write.csv(pnst.bt$NST.compare,file = paste0(prefix,".pNST.compare.",colnames(treat)[i],".csv"))

# 15 # summarize core, rare, and other taxa
# 15.1 # define the types of different taxa in category.txt
setwd(wd)
cate.file="category.txt"
cate=read.table(cate.file, header = TRUE, sep = "\t", row.names = 1,
                as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                check.names = FALSE)
cate=cate[which(rownames(cate) %in% colnames(comm)),,drop=FALSE] # remove unmatched taxa.
setwd(save.wd)

# 15.2
iccate=icamp.cate(icamp.bins.result = icbin,comm = comm,cate = cate,
                  treat = treat, silent = FALSE,between.group = TRUE)
write.csv(iccate$Ptuvx,file = paste0(prefix,".iCAMP.Process_EachTurnover_EachCategory.csv"))
write.csv(iccate$Ptx,file = paste0(prefix,".iCAMP.Process_EachGroup_EachCategory.csv"))


#MIC ENV MOD

# # Instalar (si no lo tienes)
if (!requireNamespace("MicEnvMod", quietly = TRUE)) {
  devtools::install_github("jdonhauser/MicEnvMod")
}

# Cargar el paquete
library(MicEnvMod)
#trait_table debe hacerse con los archivos fasta




#Sloan
library(devtools)
install_github("Russel88/MicEco")

library(MicEco)
library(dplyr)

# Asegúrate de tener ASVs como columnas y muestras como filas
comm <- as.data.frame(comm)
comm <- comm[rowSums(comm) > 0, colSums(comm) > 0] # Filtra ceros totales si hay

# Crear lista para guardar resultados por muestra
sloan_results <- list()

# Iterar por cada muestra
for (i in 1:nrow(comm)) {
  focal_sample <- rownames(comm)[i]
  focal_abund <- as.numeric(comm[focal_sample, ])
  
  # Pool comunitario: todas las demás muestras
  pool_abund <- colSums(comm[-i, ])
  
  # Filtrar para ASVs presentes en al menos una muestra
  use_asv <- which(pool_abund > 0)
  if (length(use_asv) < 10) next # saltar si hay pocos ASVs
  
  # Fit modelo de Sloan
  library(tidyverse)
  library(phyloseq)
  library(MASS)  # Para ajuste nls
  
  # Si no está instalado:
  # install.packages(c("phyloseq", "tidyverse", "MASS"))
  
  
  phylo_crops_bac
  # Matriz de abundancia
  otu_rel <- transform_sample_counts(phylo_crops_bac, function(x) x / sum(x))
  otu_table_rel <- as(t(otu_table(otu_rel)), "matrix")
  otu_table_rel
  # Abundancia media por taxón
  mean_abund <- rowMeans(otu_table_rel)
  mean_abund
  # Frecuencia de ocurrencia por taxón (número de muestras con presencia)
  occurrence <- apply(as(t(otu_table(phylo_crops_bac)),"matrix"), 1, function(x) sum(x > 0) / length(x))
  occurrence
  # Filtramos ASVs presentes al menos en una muestra
  use_asvs <- mean_abund > 0 & occurrence > 0
  
  # Vector con nombres de ASVs seleccionados
  asv_ids <- names(mean_abund[use_asvs])
  
  # Datos
  df_neutral <- data.frame(
    mean_abund = mean_abund[use_asvs],
    occurrence = occurrence[use_asvs]
  )
  
  N <- phyloseq::nsamples(phylo_crops_bac)
  
  # Ajuste del modelo
  neutral_model <- nls(
    occurrence ~ 1 - exp(-m * N * mean_abund),
    data = df_neutral,
    start = list(m = 0.1),
    control = nls.control(maxiter = 100, warnOnly = TRUE),
    algorithm = "default",
    trace = TRUE,
    na.action = na.exclude,
      )
  
  summary(neutral_model)
  
  neutral_model <- nls(
    occurrence ~ 1 - exp(-m * N * mean_abund),
    data = df_neutral,
    start = list(m = 0.1),
    control = nls.control(maxiter = 100, warnOnly = TRUE)
  )
  otu_table_rel2 <- as.data.frame(t(otu_table_rel))
  neutral.fit(otu_table_rel2)
  
  
  df_neutral$predicted <- 1 - exp(-coef(neutral_model)[["m"]] * N * df_neutral$mean_abund)
  m <- coef(neutral_model)[["m"]]
  
  df_neutral$predicted <- 1 - exp(-m * N * df_neutral$mean_abund)
  
  ggplot(df_neutral, aes(x = mean_abund, y = occurrence)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_line(aes(y = predicted), color = "blue", size = 1) +
    scale_x_log10() +
    labs(
      x = "Mean relative abundance (log scale)",
      y = "Occurrence frequency",
      title = "Neutral model of community assembly (Sloan)"
    ) +
    theme_minimal()

  df_neutral$type <- ifelse(
    df_neutral$occurrence > df_neutral$predicted + 1.96 * sqrt(df_neutral$predicted * (1 - df_neutral$predicted) / N),
    "Above",
    ifelse(df_neutral$occurrence < df_neutral$predicted - 1.96 * sqrt(df_neutral$predicted * (1 - df_neutral$predicted) / N),
           "Below", "Neutral")
  )
df_neutral$type
head(df_neutral)
sloan_above_bac <- df_neutral %>% dplyr::filter(type == "Above")
sloan_neutral_bac <- df_neutral %>% dplyr::filter(type == "Neutral")
sloan_below_bac <- df_neutral %>% dplyr::filter(type == "Below")
# sloan_above_bac
# sloan_neutral_bac
# sloan_below_bac
# sab <- rownames(sloan_above_bac)
# snb <- rownames(sloan_neutral_bac)
# sbb <- rownames(sloan_below_bac)
# 
# length(sab)
# sab_forest <- setdiff(sab,forest_unique_asvs) 
# sab_forest
# 
# sab_orchard <- setdiff(sab,orchard_unique_asvs) 
# sab_orchard
# sab_shared <- setdiff(sab, shared_asvs)
# sab_shared  
# length(sab_shared)
# length(sab_orchard)
# length(sab_forest)
# orchard_unique_asvs <- setdiff(orchard_asvs, forest_asvs)
# forest_unique_asvs <- setdiff(forest_asvs, orchard_asvs)
# str(orchard_unique_asvs)

ggplot(df_neutral, aes(x = mean_abund, y = occurrence)) +
  geom_point(aes(color = type), alpha = 0.7) +
  geom_line(aes(y = predicted), color = "black") +
  scale_x_log10() +
  theme_minimal()


shared_types <- df_neutral$type[rownames(df_neutral) %in% shared_asvs]
table(
  Sloan = shared_types,
  Shared = rep("Yes", length(shared_types))
)

orchard_types <- df_neutral$type[rownames(df_neutral) %in% orchard_unique_asvs]
table(
  Sloan = orchard_types,
  Shared = rep("Yes", length(orchard_types))
)

forest_types <- df_neutral$type[rownames(df_neutral) %in% forest_unique_asvs]
table(
  Sloan = forest_types,
  Shared = rep("Yes", length(forest_types))
)


shared_types
chisq.test(table(shared_types))


data_phylo_filt_bac<-data_phylo_filt

#a nivel de genero
# Subset of samples by group
orchard_samples <- rownames(sample_data(data_phylo_filt_bac)[sample_data(data_phylo_filt_bac)$Group == "Orchard", ])
forest_samples <- rownames(sample_data(data_phylo_filt_bac)[sample_data(data_phylo_filt_bac)$Group == "Forest", ])

orchard_otu <- otu_table(data_phylo_filt_bac)[orchard_samples,] 
forest_otu <- otu_table(data_phylo_filt_bac)[forest_samples,]

orchard_asvs <- taxa_names(prune_taxa(taxa_sums(orchard_otu) > 0, orchard_otu))
forest_asvs <- taxa_names(prune_taxa(taxa_sums(forest_otu) > 0, forest_otu))

shared_asvs <- intersect(orchard_asvs, forest_asvs)
shared_asvs
orchard_unique_asvs <- setdiff(orchard_asvs, forest_asvs)
forest_unique_asvs <- setdiff(forest_asvs, orchard_asvs)
str(orchard_unique_asvs)

# Print the results
cat("Shared ASVs:", length(shared_asvs), "\n")
cat("Orchard's unique ASVs:", length(orchard_unique_asvs), "\n")
cat("Forest's unique ASVs:", length(forest_unique_asvs), "\n")

# Extract the taxonomic table 
tax_data <- as.data.frame(tax_table(data_phylo_filt_bac))
tax_data$ASV <- rownames(tax_table(data_phylo_filt_bac))  # ASVs' names as column 
tax_data$ASV

# # Remove 'g__' prefix from Genus names
# tax_data$Genus <- gsub("^g__", "", tax_data$Genus)
# 
# # Filter genera for shared ASVs
# shared_genera <- tax_data %>%
#   filter(ASV %in% shared_asvs) %>%
#   select(Genus) %>%
#   distinct() %>%
#   pull(Genus)
# shared_genera

# View the shared ASVs
shared_asvs
orchard_phy <- subset_samples(data_phylo_filt_bac, Group == "Orchard")
forest_phy <- subset_samples(data_phylo_filt_bac, Group == "Forest")

#View(otu_table(orchard_phy))
str(shared_asvs)
# Subset phyloseq object for shared ASVs
shared_phy <- prune_taxa(taxa_names(data_phylo_filt_bac) %in% shared_asvs, data_phylo_filt_bac)
orchard_unique_phy <- prune_taxa(taxa_names(data_phylo_filt_bac) %in% orchard_unique_asvs, data_phylo_filt_bac)
forest_unique_phy <- prune_taxa(taxa_names(data_phylo_filt_bac) %in% forest_unique_asvs, data_phylo_filt_bac)

# Convert OTU table to data frames for shared and total ASVs
shared_otu_df <- as.data.frame(otu_table(shared_phy))  # Shared ASVs
total_otu_df <- as.data.frame(otu_table(data_phylo_filt_bac))         # Total ASVs

# OTUs should be in rows and samples in columns
if (!taxa_are_rows(shared_phy)) {
  shared_otu_df <- t(shared_otu_df)
}
if (!taxa_are_rows(data_phylo_filt_bac)) {
  total_otu_df <- t(total_otu_df)
}

# Analyze the abundance of Shared ASVs as a function of forest cover
# Calculate the number of shared and total ASVs per sample
num_shared_asvs <- colSums(shared_otu_df > 0)  # Non-zero counts for shared ASVs
num_total_asvs <- colSums(total_otu_df > 0)   # Non-zero counts for all ASVs

# Calculate the total reads for shared and total ASVs per sample
reads_shared_asvs <- colSums(shared_otu_df)
reads_total_asvs <- colSums(total_otu_df)

# Calculate percentages
percentage_shared_asvs <- (num_shared_asvs / num_total_asvs) * 100
percentage_shared_reads <- (reads_shared_asvs / reads_total_asvs) * 100

# Extraer metadata como data.frame
metadata <- as.data.frame(as(sample_data(phylo_crops_bac), "data.frame"))
metadata


shared_phy
tax_table(orchard_phy)
tax_table(orchard_unique_phy)
forest_unique_phy
dim(sloan_above_bac)
dim(sloan_neutral_bac)
dim(sloan_below_bac)

#ver cuantos hay en cada ps
sloan_ids <- rownames(sloan_above_bac)
sloan_above_bac
count_sloan_in_phy <- function(ps, ids = sloan_ids) {
  present <- intersect(taxa_names(ps), ids)
  list(
    n_present = length(present),
    n_total_sloan = length(ids),
    prop_present = length(present) / length(ids),
    present_ids = present
  )
}

res_orchard        <- count_sloan_in_phy(orchard_phy)
res_orchard_unique <- count_sloan_in_phy(orchard_unique_phy)
res_shared         <- count_sloan_in_phy(orchard_shared_phy)
res_forest_unique  <- count_sloan_in_phy(forest_unique_phy)

summary_df <- data.frame(
  phyloseq_object = c("orchard_phy", "orchard_unique_phy", "orchard_shared_phy","forest_unique_phy"),
  n_sloan_above_present = c(res_orchard$n_present, res_orchard_unique$n_present,res_shared$n_present, res_forest_unique$n_present),
  n_sloan_above_total = c(res_orchard$n_total_sloan, res_orchard_unique$n_total_sloan, res_shared$n_total_sloan,res_forest_unique$n_total_sloan),
  prop_present = c(res_orchard$prop_present, res_orchard_unique$prop_present,res_shared$prop_present, res_forest_unique$prop_present)
)

summary_df




# ver sloan_above as fx de FC2500
phylo_crops_bac
str(shared_asvs)
str(orchard_unique_asvs)
head(rownames(sloan_above_bac))
head(sloan_above_bac)


# ----------------------------
# 1) Intersecciones
# ----------------------------
ids_sloan_shared <- intersect(sloan_ids, shared_asvs)
ids_sloan_unique <- intersect(sloan_ids, orchard_unique_asvs)

# ----------------------------
# 2) Data frames solicitados
#    (solo los IDs que están en ambos)
# ----------------------------
df_sloan_shared <- data.frame(Taxon = ids_sloan_shared, stringsAsFactors = FALSE)
df_sloan_unique <- data.frame(Taxon = ids_sloan_unique, stringsAsFactors = FALSE)

# (opcional) agregar info de sloan_above_bac al df, por si quieres después
df_sloan_shared <- df_sloan_shared %>%
  left_join(
    sloan_above_bac %>%
      tibble::rownames_to_column("Taxon"),
    by = "Taxon"
  )

df_sloan_unique <- df_sloan_unique %>%
  left_join(
    sloan_above_bac %>%
      tibble::rownames_to_column("Taxon"),
    by = "Taxon"
  )

# chequeos rápidos
nrow(df_sloan_shared)
nrow(df_sloan_unique)
head(df_sloan_shared)
head(df_sloan_unique)

# ----------------------------
# 3) Podar phylo_crops_bac -> ps_shared y ps_unique
# ----------------------------
ps_shared <- prune_taxa(shared_asvs, phylo_crops_bac)
ps_unique <- prune_taxa(orchard_unique_asvs, phylo_crops_bac)

# (recomendado) eliminar taxa que queden en 0 y samples en 0
ps_shared <- prune_taxa(taxa_sums(ps_shared) > 0, ps_shared)
ps_shared <- prune_samples(sample_sums(ps_shared) > 0, ps_shared)

ps_unique <- prune_taxa(taxa_sums(ps_unique) > 0, ps_unique)
ps_unique <- prune_samples(sample_sums(ps_unique) > 0, ps_unique)

# chequeos
ps_shared
ps_unique
length(taxa_names(ps_shared))
length(taxa_names(ps_unique))
dim(df_sloan_shared)
dim(df_sloan_unique)

#veamos con shared first...















library(phyloseq)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(purrr)
library(broom)
library(stringr)

# ---------------------------
# Base phyloseq (DENOMINATOR)
# ---------------------------
ps_base <- phylo_crops_bac   # <- denominator is ALWAYS this (100%)

# ---------------------------
# Helper: sum RA of a taxon subset, where RA is computed vs ps_base totals
# ---------------------------
sum_subset_RA_over_base <- function(ps, subset_taxa, meta_var = "P_2500") {
  
  otu <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) otu <- t(otu)
  
  # Keep only taxa that truly exist in ps (prevents subscript out of bounds)
  subset_taxa_ok <- intersect(subset_taxa, rownames(otu))
  if (length(subset_taxa_ok) == 0) stop("None of the subset taxa are present in ps_base.")
  
  # Relative abundance matrix vs ps_base total counts per sample
  otu_rel <- sweep(otu, 2, colSums(otu), "/")
  
  # Sum RA of subset taxa per sample (this is your numerator; denominator is 1 == 100%)
  RA_subset <- colSums(otu_rel[subset_taxa_ok, , drop = FALSE], na.rm = TRUE)
  
  df <- tibble(
    Sample = names(RA_subset),
    RA_subset = as.numeric(RA_subset),
    Percent_subset = 100 * as.numeric(RA_subset)
  )
  
  # Add metadata if present
  md <- data.frame(sample_data(ps)) %>% rownames_to_column("Sample")
  if (meta_var %in% colnames(md)) {
    md[[meta_var]] <- as.numeric(as.character(md[[meta_var]]))
    df <- left_join(df, md %>% dplyr::select(Sample, all_of(meta_var)), by = "Sample")
  }
  
  # Also compute global % using counts (same denominator, but pooled across samples)
  sub_counts   <- colSums(otu[subset_taxa_ok, , drop = FALSE])
  total_counts <- colSums(otu)
  global_RA <- sum(sub_counts) / sum(total_counts)
  
  list(
    df = df,
    subset_taxa_used = subset_taxa_ok,
    n_used = length(subset_taxa_ok),
    n_requested = length(subset_taxa),
    global_RA = global_RA,
    global_percent = 100 * global_RA
  )
}

# ------------------------------------------------------------
# 1) SHARED case: taxa in df_sloan_shared over ps_base totals
# ------------------------------------------------------------
res_shared <- sum_subset_RA_over_base(
  ps = ps_base,
  subset_taxa = df_sloan_shared$Taxon,
  meta_var = "P_2500"
)

df_shared <- res_shared$df %>% mutate(Set = "shared")
res_shared$n_used
res_shared$global_percent

# ------------------------------------------------------------
# 2) UNIQUE case: taxa in df_sloan_unique over ps_base totals
# ------------------------------------------------------------
res_unique <- sum_subset_RA_over_base(
  ps = ps_base,
  subset_taxa = df_sloan_unique$Taxon,
  meta_var = "P_2500"
)

df_unique <- res_unique$df %>% mutate(Set = "unique")
res_unique$n_used
res_unique$global_percent

# ------------------------------------------------------------
# Combine and keep only samples with numeric P_2500
# ------------------------------------------------------------
df_both <- bind_rows(df_shared, df_unique) %>%
  filter(!is.na(P_2500))

# ------------------------------------------------------------
# Plot (percent of total phylo_crops_bac) vs P_2500
# ------------------------------------------------------------
ggplot(df_both, aes(x = P_2500, y = Percent_subset, color = Set)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~ Set, scales = "free_y") +
  theme_bw() +
  labs(
    x = "P_2500",
    y = "Sloan-above taxa (% of total phylo_crops_bac)",
    color = "Subset"
  )
#634*369
# ------------------------------------------------------------
# Polynomial fit per Set + equations + model stats
# ------------------------------------------------------------
deg <- 2  # cubic; change to 2 for quadratic, etc.

# Polynomial plot per Set
ggplot(df_both, aes(x = P_2500, y = Percent_subset, color = Set)) +
  geom_point(alpha = 0.7) +
  geom_smooth(
    method = "lm",
    formula = y ~ poly(x, deg, raw = TRUE),
    se = TRUE
  ) +
  facet_wrap(~ Set, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(
    x = "P_2500",
    y = "Sloan-above taxa (% of total phylo_crops_bac)"
  )

# Build equation string from model coefficients
coef_to_eq <- function(model, deg = 2, yname = "y", xname = "x") {
  b <- coef(model)
  b0 <- b["(Intercept)"]
  # remaining are in order: x, x^2, x^3 (because raw=TRUE)
  bs <- b[names(b) != "(Intercept)"]
  
  eq <- paste0(yname, " = ", format(b0, digits = 4, scientific = TRUE))
  for (i in seq_along(bs)) {
    bi <- bs[i]
    sign_str <- if (bi >= 0) " + " else " - "
    term <- if (i == 1) xname else paste0(xname, "^", i)
    eq <- paste0(eq, sign_str, format(abs(bi), digits = 4, scientific = TRUE), " * ", term)
  }
  eq
}

mods <- df_both %>%
  group_by(Set) %>%
  nest() %>%
  mutate(
    model  = map(data, ~ lm(Percent_subset ~ poly(P_2500, deg, raw = TRUE), data = .x)),
    glance = map(model, broom::glance),
    equation = map_chr(model, ~ coef_to_eq(.x, deg = deg, yname = "Percent", xname = "P_2500")),
    n = map_int(data, nrow)
  )

mods_eq <- mods %>%
  unnest(glance) %>%
  transmute(
    Set,
    n,
    equation,
    r2 = r.squared,
    adj_r2 = adj.r.squared,
    F = statistic,
    df1 = df,
    df2 = df.residual,
    p_model = p.value,
    AIC = AIC
  )

mods_eq



#now lets attach the NB info









library(phyloseq)
library(dplyr)
library(stringr)
library(tibble)

# ---------------------------
# INPUTS
# ---------------------------
ps_base <- phylo_crops_bac   # <- your base phyloseq (the 100% denominator)
# needs: bacti_pal_graph_long with columns ID, NB_Category, Genus, RelAb_beta (as you showed)

# ------------------------------------------------------------
# 0) SAFETY: make tax_table column names unique (prevents join errors)
#    keep the LAST occurrence of each duplicated name
# ------------------------------------------------------------
tax_mat0 <- as(tax_table(ps_base), "matrix")
tax_mat0 <- tax_mat0[, !duplicated(colnames(tax_mat0), fromLast = TRUE), drop = FALSE]
tax_table(ps_base) <- tax_table(tax_mat0)

# ------------------------------------------------------------
# 1) Build your genus->NB dictionary from bacti_pal_graph_long
# ------------------------------------------------------------
bacti_pal_graph_small <- bacti_pal_graph_long %>%
  group_by(ID, NB_Category, Genus) %>%
  summarise(relab_beta = sum(RelAb_beta, na.rm = TRUE), .groups = "drop")

clean_genus <- function(x) {
  x <- as.character(x)
  x <- str_remove(x, "^g__")
  x <- str_remove(x, "^g_")
  x <- str_replace_all(x, "\\s+", "")
  x
}

genus2nb <- bacti_pal_graph_small %>%
  mutate(Genus_clean = clean_genus(Genus)) %>%
  count(Genus_clean, NB_Category, name = "n") %>%
  group_by(Genus_clean) %>%
  slice_max(order_by = n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(Genus_clean, NB_Category)

# ------------------------------------------------------------
# 2) Extract tax_table from ps_base, preserve taxa_id, detect Genus column
# ------------------------------------------------------------
tax <- as.data.frame(tax_table(ps_base), stringsAsFactors = FALSE)
tax$taxa_id <- rownames(tax)

# find Genus column robustly
if (!"Genus" %in% colnames(tax)) {
  gcol <- grep("^Genus$|Genus", colnames(tax), value = TRUE)[1]
  if (is.na(gcol)) stop("No Genus column found in tax_table(ps_base).")
  tax$Genus <- tax[[gcol]]
}

tax$Genus_clean <- clean_genus(tax$Genus)

# ------------------------------------------------------------
# 3) Remove any existing NB_Category to avoid NB_Category.x/y
# ------------------------------------------------------------
if ("NB_Category" %in% colnames(tax)) {
  tax$NB_Category <- NULL
}

# ------------------------------------------------------------
# 4) Join NB from dictionary and fill unmatched
# ------------------------------------------------------------
tax2 <- tax %>%
  left_join(genus2nb, by = "Genus_clean") %>%
  mutate(
    NB_Category = ifelse(is.na(NB_Category) | NB_Category == "", "unassigned", as.character(NB_Category))
  ) %>%
  dplyr::select(-Genus_clean)

# ------------------------------------------------------------
# 5) CRITICAL: reorder to EXACT taxa_names(ps_base) before writing back
# ------------------------------------------------------------
tax2 <- tax2[match(taxa_names(ps_base), tax2$taxa_id), ]
stopifnot(all(tax2$taxa_id == taxa_names(ps_base)))

tax_mat <- as.matrix(tax2 %>% dplyr::select(-taxa_id))
rownames(tax_mat) <- taxa_names(ps_base)

tax_table(ps_base) <- tax_table(tax_mat)

# quick sanity check
table(tax_table(ps_base)[, "NB_Category"])


library(phyloseq)
library(dplyr)
library(tibble)
library(ggplot2)
library(purrr)
library(broom)

# ---------------------------
# SETTINGS
# ---------------------------
deg <- 2                    # polynomial degree (2 = cub)
meta_var <- "P_2500"         # x-axis variable

# ------------------------------------------------------------
# 0) SAFETY: make tax_table column names unique (prevents join errors)
# ------------------------------------------------------------
tax_mat0 <- as(tax_table(ps_base), "matrix")
tax_mat0 <- tax_mat0[, !duplicated(colnames(tax_mat0), fromLast = TRUE), drop = FALSE]
tax_table(ps_base) <- tax_table(tax_mat0)

# ------------------------------------------------------------
# 1) Taxonomy lookup: Taxon -> NB_Category (from ps_base)
# ------------------------------------------------------------
tax_df <- as.data.frame(tax_table(ps_base), stringsAsFactors = FALSE) %>%
  tibble::rownames_to_column("Taxon")

if (!"NB_Category" %in% colnames(tax_df)) {
  stop("NB_Category not found in tax_table(ps_base). Add/map it first.")
}

tax_df <- tax_df %>%
  mutate(NB_Category = ifelse(is.na(NB_Category) | NB_Category == "", "unassigned", as.character(NB_Category))) %>%
  dplyr::select(Taxon, NB_Category)

# ------------------------------------------------------------
# 2) Metadata: Sample -> P_2500 (numeric)
# ------------------------------------------------------------
meta_df <- data.frame(sample_data(ps_base)) %>%
  tibble::rownames_to_column("Sample")

if (!meta_var %in% colnames(meta_df)) stop(paste0(meta_var, " not found in sample_data(ps_base)."))

meta_df[[meta_var]] <- as.numeric(as.character(meta_df[[meta_var]]))

# ------------------------------------------------------------
# 3) OTU relative abundance matrix vs FULL ps_base totals
# ------------------------------------------------------------
otu <- as(otu_table(ps_base), "matrix")
if (!taxa_are_rows(ps_base)) otu <- t(otu)

otu_rel <- sweep(otu, 2, colSums(otu), "/")   # taxa x sample

# ------------------------------------------------------------
# 4) Function: compute NB-wise % abundance for one Sloan set
# ------------------------------------------------------------
compute_nb_percent <- function(set_name, df_sloan, otu_rel, tax_df, meta_df, meta_var) {
  subset_taxa <- df_sloan$Taxon
  subset_taxa_ok <- intersect(subset_taxa, rownames(otu_rel))
  
  if (length(subset_taxa_ok) == 0) stop(paste0("No taxa from ", set_name, " found in ps_base."))
  
  # NB categories for subset taxa
  sub_tax <- tax_df %>%
    filter(Taxon %in% subset_taxa_ok) %>%
    mutate(NB_Category = ifelse(is.na(NB_Category) | NB_Category == "", "unassigned", NB_Category))
  
  # taxa split by NB category
  taxa_by_nb <- split(sub_tax$Taxon, sub_tax$NB_Category)
  
  # sum RA per sample for each NB (matrix-efficient)
  df_list <- lapply(names(taxa_by_nb), function(nb) {
    tx <- taxa_by_nb[[nb]]
    ra <- colSums(otu_rel[tx, , drop = FALSE], na.rm = TRUE)
    tibble(
      Sample = names(ra),
      NB_Category = nb,
      RA_subset = as.numeric(ra),
      Percent_subset = 100 * as.numeric(ra),
      Set = set_name
    )
  })
  
  df <- bind_rows(df_list) %>%
    left_join(meta_df %>% dplyr::select(Sample, all_of(meta_var)), by = "Sample") %>%
    filter(!is.na(.data[[meta_var]]))
  
  # global pooled percent (counts-based), still vs full ps_base total
  sub_counts <- colSums(otu[subset_taxa_ok, , drop = FALSE])
  total_counts <- colSums(otu)
  global_percent <- 100 * (sum(sub_counts) / sum(total_counts))
  
  list(
    df = df,
    n_used = length(subset_taxa_ok),
    n_requested = length(subset_taxa),
    global_percent = global_percent
  )
}

# ------------------------------------------------------------
# 5) Build df for SHARED and UNIQUE
# ------------------------------------------------------------
res_shared <- compute_nb_percent("shared", df_sloan_shared, otu_rel, tax_df, meta_df, meta_var)
res_unique <- compute_nb_percent("unique", df_sloan_unique, otu_rel, tax_df, meta_df, meta_var)

df_nb <- bind_rows(res_shared$df, res_unique$df)

res_shared$global_percent
res_unique$global_percent

# ------------------------------------------------------------
# 6) PLOTS
# ------------------------------------------------------------

# (A) One big faceted plot: NB rows x Set columns
ggplot(df_nb, aes(x = .data[[meta_var]], y = Percent_subset)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", formula = y ~ poly(x, deg, raw = TRUE), se = TRUE) +
  facet_grid(NB_Category ~ Set, scales = "free_y") +
  theme_bw() +
  labs(
    x = meta_var,
    y = "Sloan-above taxa (% of total phylo_crops_bac)",
    title = paste0("NB-wise Sloan-above abundance vs ", meta_var, " (poly deg=", deg, ")")
  )

# (B) Separate plots per Set (shared vs unique)
ggplot(df_nb %>% filter(Set == "shared"),
       aes(x = .data[[meta_var]], y = Percent_subset, color = NB_Category)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", formula = y ~ poly(x, deg, raw = TRUE), se = TRUE) +
  facet_wrap(~ NB_Category, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = meta_var, y = "Percent of total", title = "Shared: NB-wise Sloan-above (%)")

ggplot(df_nb %>% filter(Set == "unique"),
       aes(x = .data[[meta_var]], y = Percent_subset, color = NB_Category)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", formula = y ~ poly(x, deg, raw = TRUE), se = TRUE) +
  facet_wrap(~ NB_Category, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = meta_var, y = "Percent of total", title = "Unique: NB-wise Sloan-above (%)")



library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(stringr)
library(tibble)

# set these if not already defined
# deg <- 3
# meta_var <- "P_2500"

coef_to_eq <- function(model, yname = "Percent", xname = "P_2500") {
  b <- coef(model)
  
  b0 <- unname(b["(Intercept)"])
  bt <- b[names(b) != "(Intercept)"]
  
  # order polynomial terms by the trailing number ... )1 )2 )3
  term_id <- str_match(names(bt), "\\)(\\d+)$")[,2]
  ord <- order(as.integer(term_id))
  bt <- bt[ord]
  
  eq <- paste0(yname, " = ", format(b0, digits = 4, scientific = TRUE))
  
  for (i in seq_along(bt)) {
    bi <- unname(bt[i])
    sign_str <- if (bi >= 0) " + " else " - "
    term <- if (i == 1) xname else paste0(xname, "^", i)
    eq <- paste0(eq, sign_str, format(abs(bi), digits = 4, scientific = TRUE), " * ", term)
  }
  eq
}

safe_lm <- function(dat, meta_var, deg) {
  # checks
  if (!("Percent_subset" %in% names(dat))) return(NULL)
  if (!(meta_var %in% names(dat))) return(NULL)
  
  y <- dat$Percent_subset
  x <- dat[[meta_var]]
  
  if (nrow(dat) < (deg + 2)) return(NULL)
  if (all(is.na(y))) return(NULL)
  if (sd(y, na.rm = TRUE) == 0) return(NULL)
  if (all(is.na(x))) return(NULL)
  if (sd(x, na.rm = TRUE) == 0) return(NULL)
  
  # build formula safely: Percent_subset ~ poly(P_2500, deg, raw=TRUE)
  fml <- as.formula(paste0("Percent_subset ~ poly(", meta_var, ", ", deg, ", raw = TRUE)"))
  lm(fml, data = dat)
}

# ------------------------------------------------------------
# FIT MODELS per (Set × NB_Category)
# ------------------------------------------------------------
mods <- df_nb %>%
  group_by(Set, NB_Category) %>%
  nest() %>%
  mutate(
    n = map_int(data, nrow),
    model = map(data, safe_lm, meta_var = meta_var, deg = deg),
    
    equation = map_chr(
      model,
      ~ if (is.null(.x)) NA_character_ else coef_to_eq(.x, yname = "Percent", xname = meta_var)
    ),
    
    glance = map(
      model,
      ~ if (is.null(.x)) {
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
    )
  ) %>%
  unnest(glance)

mods_eq <- mods %>%
  transmute(
    Set, NB_Category, n,
    equation,
    r2 = r.squared,
    adj_r2 = adj.r.squared,
    F = statistic,
    df1 = df,
    df2 = df.residual,
    p_model = p.value,
    AIC = AIC
  )

mods_eq
mods_eq$equation

2# ------------------------------------------------------------
# (Optional) COEFFICIENT TABLES per (Set × NB_Category)
# ------------------------------------------------------------
empty_tidy <- tibble(
  term = character(),
  estimate = numeric(),
  std.error = numeric(),
  statistic = numeric(),
  p.value = numeric()
)

coef_tables <- mods %>%
  transmute(
    Set, NB_Category,
    coef = map(model, ~ if (is.null(.x)) empty_tidy else broom::tidy(.x))
  ) %>%
  unnest(coef)

coef_tables









