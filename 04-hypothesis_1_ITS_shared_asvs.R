# =============================================================================
# Script: 04_hypothesis_1_ITS_shared_asvs.R
# Project: From forest to fields: The role of soil microbiome spillover in
#          agroecosystem sustainability
# Author: Pedro Mondaca
# Contact: pedromondaca@outlook.com
# Submitted to: Science Advances
#
# Description:
# This script evaluates Hypothesis 1 for fungal communities:
# orchards with greater surrounding forest cover share more soil microbial taxa
# with adjacent forests. It identifies ASVs shared between orchard and forest
# soils, quantifies the proportion of shared ASVs and reads per sample, fits
# frequentist and Bayesian models against forest cover and edaphic/spatial
# variables, produces exploratory plots, and exports summary tables and model
# comparisons.
#
# Required files in the working directory:
#   - results_ITS/phyloseq_ITS_asv_filtered.rds
#   - results_ITS/mapping_ITS_ordered.csv
#
# Output files:
#   - written to the folder "results_ITS_hyp1"
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(phyloseq)
  library(ggplot2)
  library(ggpmisc)
  library(lme4)
  library(lmerTest)
  library(report)
  library(performance)
  library(brms)
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(loo)
})

# =============================================================================
# Paths and setup
# =============================================================================

setwd("c:\\Users\\pedro\\Research\\microbiome_orchard_forest\\Sci Adv - v2\\GitHub")

input_phyloseq <- "results_ITS/phyloseq_ITS_asv_filtered.rds"
input_mapping  <- "results_ITS/mapping_ITS_ordered.csv"
outdir         <- "results_ITS_hyp1"

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

required_files <- c(input_phyloseq, input_mapping)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required files: ", paste(missing_files, collapse = ", "))
}

# =============================================================================
# Load inputs
# =============================================================================

data_phylo_filt <- readRDS(input_phyloseq)
mapping <- read.csv(input_mapping, row.names = 1, check.names = FALSE)

mapping <- mapping[sample_names(data_phylo_filt), , drop = FALSE]

if (!identical(rownames(mapping), sample_names(data_phylo_filt))) {
  stop("Mapping file row names do not match phyloseq sample names.")
}

sample_data(data_phylo_filt) <- sample_data(mapping)

# =============================================================================
# Identify orchard and forest ASVs
# =============================================================================

orchard_phy <- subset_samples(data_phylo_filt, Group == "Orchard")
forest_phy  <- subset_samples(data_phylo_filt, Group == "Forest")

orchard_otu <- as(otu_table(orchard_phy), "matrix")
forest_otu  <- as(otu_table(forest_phy), "matrix")

if (!taxa_are_rows(orchard_phy)) orchard_otu <- t(orchard_otu)
if (!taxa_are_rows(forest_phy))  forest_otu  <- t(forest_otu)

orchard_asvs <- rownames(orchard_otu)[rowSums(orchard_otu) > 0]
forest_asvs  <- rownames(forest_otu)[rowSums(forest_otu) > 0]

shared_asvs         <- intersect(orchard_asvs, forest_asvs)
orchard_unique_asvs <- setdiff(orchard_asvs, forest_asvs)
forest_unique_asvs  <- setdiff(forest_asvs, orchard_asvs)

shared_counts <- data.frame(
  Metric = c("Shared_ASVs", "Orchard_unique_ASVs", "Forest_unique_ASVs"),
  Value  = c(length(shared_asvs), length(orchard_unique_asvs), length(forest_unique_asvs))
)

cat("Shared ASVs:", length(shared_asvs), "\n")
cat("Orchard unique ASVs:", length(orchard_unique_asvs), "\n")
cat("Forest unique ASVs:", length(forest_unique_asvs), "\n")

write.csv(shared_counts, file.path(outdir, "shared_asv_counts.csv"), row.names = FALSE)

# =============================================================================
# Taxonomy table for shared ASVs
# =============================================================================

tax_data <- as.data.frame(tax_table(data_phylo_filt))
tax_data$ASV <- rownames(tax_data)

shared_taxonomy <- tax_data %>%
  filter(ASV %in% shared_asvs)

write.csv(shared_taxonomy, file.path(outdir, "shared_asvs_taxonomy.csv"), row.names = FALSE)

# =============================================================================
# Quantify shared ASVs and shared reads per sample
# =============================================================================

shared_phy <- prune_taxa(shared_asvs, data_phylo_filt)

shared_otu_df <- as(otu_table(shared_phy), "matrix")
total_otu_df  <- as(otu_table(data_phylo_filt), "matrix")

if (taxa_are_rows(shared_phy)) {
  shared_otu_df <- t(shared_otu_df)
}
if (taxa_are_rows(data_phylo_filt)) {
  total_otu_df <- t(total_otu_df)
}

sample_names_all <- rownames(total_otu_df)
sample_names_shared <- rownames(shared_otu_df)

if (!identical(sample_names_all, sample_names_shared)) {
  stop("Sample names in shared and total OTU matrices do not match.")
}

num_shared_asvs   <- rowSums(shared_otu_df > 0)
num_total_asvs    <- rowSums(total_otu_df > 0)
reads_shared_asvs <- rowSums(shared_otu_df)
reads_total_asvs  <- rowSums(total_otu_df)

percentage_shared_asvs  <- (num_shared_asvs / num_total_asvs) * 100
percentage_shared_reads <- (reads_shared_asvs / reads_total_asvs) * 100

sample_metrics <- data.frame(
  SampleID = sample_names_all,
  Num_Shared_ASVs = num_shared_asvs,
  Num_Total_ASVs = num_total_asvs,
  Reads_Shared_ASVs = reads_shared_asvs,
  Reads_Total_ASVs = reads_total_asvs,
  Percentage_Shared_ASVs = percentage_shared_asvs,
  Percentage_Shared_Reads = percentage_shared_reads,
  stringsAsFactors = FALSE
)

write.csv(sample_metrics, file.path(outdir, "sample_metrics_shared_asvs.csv"), row.names = FALSE)

# =============================================================================
# Build orchard and forest summary tables
# =============================================================================

metadata_all <- data.frame(sample_data(data_phylo_filt), check.names = FALSE)
metadata_all$SampleID <- rownames(metadata_all)

summary_df <- dplyr::left_join(metadata_all, sample_metrics, by = "SampleID")
summary_df <- data.frame(summary_df, check.names = FALSE)

orchard_summary_df <- summary_df %>%
  dplyr::filter(Group == "Orchard")

forest_summary_df <- summary_df %>%
  dplyr::filter(Group == "Forest")

cat("Median orchard shared reads (%):",
    median(orchard_summary_df$Percentage_Shared_Reads, na.rm = TRUE), "\n")
cat("Median forest shared reads (%):",
    median(forest_summary_df$Percentage_Shared_Reads, na.rm = TRUE), "\n")

write.csv(orchard_summary_df, file.path(outdir, "orchard_summary_ITS.csv"), row.names = FALSE)
write.csv(forest_summary_df,  file.path(outdir, "forest_summary_ITS.csv"), row.names = FALSE)

# =============================================================================
# Exploratory plot: shared ASVs vs forest cover
# =============================================================================

pH_colors <- c("red", "yellow", "blue")

p_shared_asv_eq <- ggplot(orchard_summary_df, aes(x = P_2500, y = Percentage_Shared_ASVs)) +
  geom_point(aes(color = pH, size = SOMTXT), alpha = 0.9) +
  scale_color_gradientn(
    colors = pH_colors,
    name = "pH",
    limits = c(
      min(orchard_summary_df$pH, na.rm = TRUE),
      max(orchard_summary_df$pH, na.rm = TRUE)
    )
  ) +
  geom_smooth(
    formula = y ~ poly(x, 2, raw = TRUE),
    method = "lm",
    se = TRUE,
    color = "darkblue",
    linewidth = 1.2,
    linetype = "solid"
  ) +
  stat_poly_eq(
    formula = y ~ poly(x, 2, raw = TRUE),
    aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
    parse = TRUE,
    size = 5,
    label.x = 0.8,
    label.y = 30
  ) +
  labs(
    title = "",
    x = expression(italic(P[2500])),
    y = "Shared ASVs abundance (%)",
    size = "SOMTXT"
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.position = "right"
  )

ggsave(
  filename = file.path(outdir, "shared_asvs_vs_P2500_with_eq.png"),
  plot = p_shared_asv_eq,
  width = 7,
  height = 4.5,
  dpi = 300
)

p_shared_asv <- ggplot(orchard_summary_df, aes(x = P_2500, y = Percentage_Shared_ASVs)) +
  geom_point(aes(color = pH, size = SOMTXT), alpha = 0.9) +
  scale_color_gradientn(
    colors = pH_colors,
    name = "pH",
    limits = c(
      min(orchard_summary_df$pH, na.rm = TRUE),
      max(orchard_summary_df$pH, na.rm = TRUE)
    )
  ) +
  geom_smooth(
    formula = y ~ poly(x, 2, raw = TRUE),
    method = "lm",
    se = TRUE,
    color = "darkblue",
    linewidth = 1.2,
    linetype = "solid"
  ) +
  labs(
    title = "",
    x = expression(italic(P[2500])),
    y = "Shared ASVs abundance (%)",
    size = "SOMTXT"
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.position = "right"
  )

ggsave(
  filename = file.path(outdir, "shared_asvs_vs_P2500.png"),
  plot = p_shared_asv,
  width = 7,
  height = 4.5,
  dpi = 300
)

# =============================================================================
# Frequentist models
# =============================================================================

model_lm <- lm(
  Percentage_Shared_ASVs ~ poly(P_2500, 2, raw = TRUE),
  data = orchard_summary_df
)

model_lmer <- lmer(
  Percentage_Shared_ASVs ~ poly(P_2500, 2, raw = TRUE) + (1 | Block),
  data = orchard_summary_df
)

model_glmer <- glmer(
  cbind(Num_Shared_ASVs, Num_Total_ASVs - Num_Shared_ASVs) ~
    poly(P_2500, 2, raw = TRUE) + (1 | Block),
  family = binomial(link = "logit"),
  data = orchard_summary_df
)

writeLines(capture.output(summary(model_lm)),
           file.path(outdir, "model_lm_summary.txt"))
writeLines(capture.output(report(model_lm)),
           file.path(outdir, "model_lm_report.txt"))

writeLines(capture.output(summary(model_lmer)),
           file.path(outdir, "model_lmer_summary.txt"))
writeLines(capture.output(report(model_lmer)),
           file.path(outdir, "model_lmer_report.txt"))
writeLines(capture.output(ranova(model_lmer)),
           file.path(outdir, "model_lmer_ranova.txt"))

writeLines(capture.output(summary(model_glmer)),
           file.path(outdir, "model_glmer_summary.txt"))
writeLines(capture.output(performance::check_model(model_glmer)),
           file.path(outdir, "model_glmer_diagnostics.txt"))

# =============================================================================
# Prepare bounded response for beta models
# =============================================================================

orchard_summary_df <- orchard_summary_df %>%
  mutate(
    Shared_prop = Percentage_Shared_Reads / 100
  )

n_obs <- nrow(orchard_summary_df)

orchard_summary_df <- orchard_summary_df %>%
  mutate(
    Shared_prop_beta = (Shared_prop * (n_obs - 1) + 0.5) / n_obs
  )

# =============================================================================
# Spatial coordinates
# =============================================================================

geo_df <- orchard_summary_df %>%
  dplyr::select(SampleID, Coo_X, Coo_Y)

utm_sf <- st_as_sf(geo_df, coords = c("Coo_X", "Coo_Y"), crs = 32719)
latlon <- st_transform(utm_sf, crs = 4326)
coord_matrix <- as.data.frame(st_coordinates(latlon))

orchard_summary_df$lon <- coord_matrix$X
orchard_summary_df$lat <- coord_matrix$Y

orchard_sf <- st_as_sf(
  orchard_summary_df,
  coords = c("lon", "lat"),
  crs = 4326,
  remove = FALSE
)

# =============================================================================
# Map of orchard samples
# =============================================================================

chile_map <- ne_countries(scale = "medium", returnclass = "sf") %>%
  dplyr::filter(admin == "Chile")

p_map <- ggplot(data = chile_map) +
  geom_sf(fill = "antiquewhite") +
  geom_sf(data = orchard_sf, aes(color = Shared_prop), size = 2, alpha = 0.8) +
  scale_color_viridis_c(option = "D", name = "Shared reads") +
  coord_sf(xlim = c(-71.5, -70.5), ylim = c(-35.5, -34), expand = FALSE) +
  theme_minimal() +
  labs(title = "Location of orchard samples", x = "Longitude", y = "Latitude")

ggsave(
  filename = file.path(outdir, "orchard_map_shared_reads.png"),
  plot = p_map,
  width = 7,
  height = 5,
  dpi = 300
)

# =============================================================================
# Bayesian models
# =============================================================================

set.seed(99)

model_block <- brm(
  Shared_prop_beta ~ t2(as.numeric(Block)),
  data = orchard_summary_df,
  family = Beta(),
  chains = 4,
  cores = 4,
  iter = 4000,
  seed = 99,
  control = list(adapt_delta = 0.99)
)
saveRDS(model_block, file.path(outdir, "brms_block.rds"))
writeLines(capture.output(summary(model_block)),
           file.path(outdir, "brms_block_summary.txt"))

model_P_100 <- brm(
  Shared_prop_beta ~ t2(P_100),
  data = orchard_summary_df,
  family = Beta(),
  chains = 4,
  cores = 4,
  iter = 4000,
  seed = 99,
  control = list(adapt_delta = 0.99)
)
saveRDS(model_P_100, file.path(outdir, "brms_P_100.rds"))
writeLines(capture.output(summary(model_P_100)),
           file.path(outdir, "brms_P_100_summary.txt"))

model_P_250 <- brm(
  Shared_prop_beta ~ t2(P_250),
  data = orchard_summary_df,
  family = Beta(),
  chains = 4,
  cores = 4,
  iter = 4000,
  seed = 99,
  control = list(adapt_delta = 0.99)
)
saveRDS(model_P_250, file.path(outdir, "brms_P_250.rds"))
writeLines(capture.output(summary(model_P_250)),
           file.path(outdir, "brms_P_250_summary.txt"))

model_P_500 <- brm(
  Shared_prop_beta ~ t2(P_500),
  data = orchard_summary_df,
  family = Beta(),
  chains = 4,
  cores = 4,
  iter = 4000,
  seed = 99,
  control = list(adapt_delta = 0.99)
)
saveRDS(model_P_500, file.path(outdir, "brms_P_500.rds"))
writeLines(capture.output(summary(model_P_500)),
           file.path(outdir, "brms_P_500_summary.txt"))

model_P_1000 <- brm(
  Shared_prop_beta ~ t2(P_1000),
  data = orchard_summary_df,
  family = Beta(),
  chains = 4,
  cores = 4,
  iter = 4000,
  seed = 99,
  control = list(adapt_delta = 0.99)
)
saveRDS(model_P_1000, file.path(outdir, "brms_P_1000.rds"))
writeLines(capture.output(summary(model_P_1000)),
           file.path(outdir, "brms_P_1000_summary.txt"))

model_P_2500 <- brm(
  Shared_prop_beta ~ t2(P_2500),
  data = orchard_summary_df,
  family = Beta(),
  chains = 4,
  cores = 4,
  iter = 4000,
  seed = 99,
  control = list(adapt_delta = 0.99)
)
saveRDS(model_P_2500, file.path(outdir, "brms_P_2500.rds"))
writeLines(capture.output(summary(model_P_2500)),
           file.path(outdir, "brms_P_2500_summary.txt"))

model_P_100_re <- brm(
  Shared_prop_beta ~ t2(P_100) + (1 | Block),
  data = orchard_summary_df,
  family = Beta(),
  chains = 4,
  cores = 4,
  iter = 4000,
  seed = 99,
  control = list(adapt_delta = 0.99)
)
saveRDS(model_P_100_re, file.path(outdir, "brms_P_100_re.rds"))
writeLines(capture.output(summary(model_P_100_re)),
           file.path(outdir, "brms_P_100_re_summary.txt"))

model_P_250_re <- brm(
  Shared_prop_beta ~ t2(P_250) + (1 | Block),
  data = orchard_summary_df,
  family = Beta(),
  chains = 4,
  cores = 4,
  iter = 4000,
  seed = 99,
  control = list(adapt_delta = 0.99)
)
saveRDS(model_P_250_re, file.path(outdir, "brms_P_250_re.rds"))
writeLines(capture.output(summary(model_P_250_re)),
           file.path(outdir, "brms_P_250_re_summary.txt"))

model_P_500_re <- brm(
  Shared_prop_beta ~ t2(P_500) + (1 | Block),
  data = orchard_summary_df,
  family = Beta(),
  chains = 4,
  cores = 4,
  iter = 4000,
  seed = 99,
  control = list(adapt_delta = 0.99)
)
saveRDS(model_P_500_re, file.path(outdir, "brms_P_500_re.rds"))
writeLines(capture.output(summary(model_P_500_re)),
           file.path(outdir, "brms_P_500_re_summary.txt"))

model_P_1000_re <- brm(
  Shared_prop_beta ~ t2(P_1000) + (1 | Block),
  data = orchard_summary_df,
  family = Beta(),
  chains = 4,
  cores = 4,
  iter = 4000,
  seed = 99,
  control = list(adapt_delta = 0.99)
)
saveRDS(model_P_1000_re, file.path(outdir, "brms_P_1000_re.rds"))
writeLines(capture.output(summary(model_P_1000_re)),
           file.path(outdir, "brms_P_1000_re_summary.txt"))

model_P_2500_re <- brm(
  Shared_prop_beta ~ t2(P_2500) + (1 | Block),
  data = orchard_summary_df,
  family = Beta(),
  chains = 4,
  cores = 4,
  iter = 4000,
  seed = 99,
  control = list(adapt_delta = 0.99)
)
saveRDS(model_P_2500_re, file.path(outdir, "brms_P_2500_re.rds"))
writeLines(capture.output(summary(model_P_2500_re)),
           file.path(outdir, "brms_P_2500_re_summary.txt"))

model_lat <- brm(
  Shared_prop_beta ~ t2(lat) + (1 | Block),
  data = orchard_summary_df,
  family = Beta(),
  chains = 4,
  cores = 4,
  iter = 4000,
  seed = 99,
  control = list(adapt_delta = 0.99)
)
saveRDS(model_lat, file.path(outdir, "brms_lat.rds"))
writeLines(capture.output(summary(model_lat)),
           file.path(outdir, "brms_lat_summary.txt"))

model_lon <- brm(
  Shared_prop_beta ~ t2(lon),
  data = orchard_summary_df,
  family = Beta(),
  chains = 4,
  cores = 4,
  iter = 4000,
  seed = 99,
  control = list(adapt_delta = 0.99)
)
saveRDS(model_lon, file.path(outdir, "brms_lon.rds"))
writeLines(capture.output(summary(model_lon)),
           file.path(outdir, "brms_lon_summary.txt"))

model_XY <- brm(
  Shared_prop_beta ~ t2(Coo_X, Coo_Y),
  data = orchard_summary_df,
  family = Beta(),
  chains = 4,
  cores = 4,
  iter = 4000,
  seed = 99,
  control = list(adapt_delta = 0.99)
)
saveRDS(model_XY, file.path(outdir, "brms_XY.rds"))
writeLines(capture.output(summary(model_XY)),
           file.path(outdir, "brms_XY_summary.txt"))

model_pH <- brm(
  Shared_prop_beta ~ t2(pH) + (1 | Block),
  data = orchard_summary_df,
  family = Beta(),
  chains = 4,
  cores = 4,
  iter = 4000,
  seed = 99,
  control = list(adapt_delta = 0.99)
)
saveRDS(model_pH, file.path(outdir, "brms_pH.rds"))
writeLines(capture.output(summary(model_pH)),
           file.path(outdir, "brms_pH_summary.txt"))

model_SOM <- brm(
  Shared_prop_beta ~ t2(SOM) + (1 | Block),
  data = orchard_summary_df,
  family = Beta(),
  chains = 4,
  cores = 4,
  iter = 4000,
  seed = 99,
  control = list(adapt_delta = 0.99)
)
saveRDS(model_SOM, file.path(outdir, "brms_SOM.rds"))
writeLines(capture.output(summary(model_SOM)),
           file.path(outdir, "brms_SOM_summary.txt"))

model_N <- brm(
  Shared_prop_beta ~ t2(N) + (1 | Block),
  data = orchard_summary_df,
  family = Beta(),
  chains = 4,
  cores = 4,
  iter = 4000,
  seed = 99,
  control = list(adapt_delta = 0.99)
)
saveRDS(model_N, file.path(outdir, "brms_N.rds"))
writeLines(capture.output(summary(model_N)),
           file.path(outdir, "brms_N_summary.txt"))

model_claysilt <- brm(
  Shared_prop_beta ~ t2(claysilt) + (1 | Block),
  data = orchard_summary_df,
  family = Beta(),
  chains = 4,
  cores = 4,
  iter = 4000,
  seed = 99,
  control = list(adapt_delta = 0.99)
)
saveRDS(model_claysilt, file.path(outdir, "brms_claysilt.rds"))
writeLines(capture.output(summary(model_claysilt)),
           file.path(outdir, "brms_claysilt_summary.txt"))

model_P250_pH <- brm(
  Shared_prop_beta ~ t2(P_250, pH) + (1 | Block),
  data = orchard_summary_df,
  family = Beta(),
  chains = 4,
  cores = 4,
  iter = 4000,
  seed = 99,
  control = list(adapt_delta = 0.99)
)
saveRDS(model_P250_pH, file.path(outdir, "brms_P250_pH.rds"))
writeLines(capture.output(summary(model_P250_pH)),
           file.path(outdir, "brms_P250_pH_summary.txt"))

model_P250_SOM <- brm(
  Shared_prop_beta ~ t2(P_250, SOM) + (1 | Block),
  data = orchard_summary_df,
  family = Beta(),
  chains = 4,
  cores = 4,
  iter = 4000,
  seed = 99,
  control = list(adapt_delta = 0.99)
)
saveRDS(model_P250_SOM, file.path(outdir, "brms_P250_SOM.rds"))
writeLines(capture.output(summary(model_P250_SOM)),
           file.path(outdir, "brms_P250_SOM_summary.txt"))

model_P250_N <- brm(
  Shared_prop_beta ~ t2(P_250, N) + (1 | Block),
  data = orchard_summary_df,
  family = Beta(),
  chains = 4,
  cores = 4,
  iter = 4000,
  seed = 99,
  control = list(adapt_delta = 0.99)
)
saveRDS(model_P250_N, file.path(outdir, "brms_P250_N.rds"))
writeLines(capture.output(summary(model_P250_N)),
           file.path(outdir, "brms_P250_N_summary.txt"))

model_P250_clay <- brm(
  Shared_prop_beta ~ t2(P_250, claysilt) + (1 | Block),
  data = orchard_summary_df,
  family = Beta(),
  chains = 4,
  cores = 4,
  iter = 4000,
  seed = 99,
  control = list(adapt_delta = 0.99)
)
saveRDS(model_P250_clay, file.path(outdir, "brms_P250_clay.rds"))
writeLines(capture.output(summary(model_P250_clay)),
           file.path(outdir, "brms_P250_clay_summary.txt"))

model_P2500_pH <- brm(
  Shared_prop_beta ~ t2(P_2500, pH) + (1 | Block),
  data = orchard_summary_df,
  family = Beta(),
  chains = 4,
  cores = 4,
  iter = 4000,
  seed = 99,
  control = list(adapt_delta = 0.99)
)
saveRDS(model_P2500_pH, file.path(outdir, "brms_P2500_pH.rds"))
writeLines(capture.output(summary(model_P2500_pH)),
           file.path(outdir, "brms_P2500_pH_summary.txt"))

model_P2500_SOM <- brm(
  Shared_prop_beta ~ t2(P_2500, SOM) + (1 | Block),
  data = orchard_summary_df,
  family = Beta(),
  chains = 4,
  cores = 4,
  iter = 4000,
  seed = 99,
  control = list(adapt_delta = 0.99)
)
saveRDS(model_P2500_SOM, file.path(outdir, "brms_P2500_SOM.rds"))
writeLines(capture.output(summary(model_P2500_SOM)),
           file.path(outdir, "brms_P2500_SOM_summary.txt"))

model_P2500_N <- brm(
  Shared_prop_beta ~ t2(P_2500, N) + (1 | Block),
  data = orchard_summary_df,
  family = Beta(),
  chains = 4,
  cores = 4,
  iter = 4000,
  seed = 99,
  control = list(adapt_delta = 0.99)
)
saveRDS(model_P2500_N, file.path(outdir, "brms_P2500_N.rds"))
writeLines(capture.output(summary(model_P2500_N)),
           file.path(outdir, "brms_P2500_N_summary.txt"))

model_P2500_clay <- brm(
  Shared_prop_beta ~ t2(P_2500, claysilt) + (1 | Block),
  data = orchard_summary_df,
  family = Beta(),
  chains = 4,
  cores = 4,
  iter = 4000,
  seed = 99,
  control = list(adapt_delta = 0.99)
)
saveRDS(model_P2500_clay, file.path(outdir, "brms_P2500_clay.rds"))
writeLines(capture.output(summary(model_P2500_clay)),
           file.path(outdir, "brms_P2500_clay_summary.txt"))

# =============================================================================
# LOO comparison
# =============================================================================

loo_block <- loo(model_block)
loo_P_100 <- loo(model_P_100)
loo_P_250 <- loo(model_P_250)
loo_P_500 <- loo(model_P_500)
loo_P_1000 <- loo(model_P_1000)
loo_P_2500 <- loo(model_P_2500)
loo_P_100_re <- loo(model_P_100_re)
loo_P_250_re <- loo(model_P_250_re)
loo_P_500_re <- loo(model_P_500_re)
loo_P_1000_re <- loo(model_P_1000_re)
loo_P_2500_re <- loo(model_P_2500_re)
loo_lat <- loo(model_lat)
loo_lon <- loo(model_lon)
loo_XY <- loo(model_XY)
loo_pH <- loo(model_pH)
loo_SOM <- loo(model_SOM)
loo_N <- loo(model_N)
loo_claysilt <- loo(model_claysilt)
loo_P250_pH <- loo(model_P250_pH)
loo_P250_SOM <- loo(model_P250_SOM)
loo_P250_N <- loo(model_P250_N)
loo_P250_clay <- loo(model_P250_clay)
loo_P2500_pH <- loo(model_P2500_pH)
loo_P2500_SOM <- loo(model_P2500_SOM)
loo_P2500_N <- loo(model_P2500_N)
loo_P2500_clay <- loo(model_P2500_clay)

loo_summary_table <- data.frame(
  model = c(
    "block", "P_100", "P_250", "P_500", "P_1000", "P_2500",
    "P_100_re", "P_250_re", "P_500_re", "P_1000_re", "P_2500_re",
    "lat", "lon", "XY", "pH", "SOM", "N", "claysilt",
    "P250_pH", "P250_SOM", "P250_N", "P250_clay",
    "P2500_pH", "P2500_SOM", "P2500_N", "P2500_clay"
  ),
  elpd_loo = c(
    loo_block$estimates["elpd_loo", "Estimate"],
    loo_P_100$estimates["elpd_loo", "Estimate"],
    loo_P_250$estimates["elpd_loo", "Estimate"],
    loo_P_500$estimates["elpd_loo", "Estimate"],
    loo_P_1000$estimates["elpd_loo", "Estimate"],
    loo_P_2500$estimates["elpd_loo", "Estimate"],
    loo_P_100_re$estimates["elpd_loo", "Estimate"],
    loo_P_250_re$estimates["elpd_loo", "Estimate"],
    loo_P_500_re$estimates["elpd_loo", "Estimate"],
    loo_P_1000_re$estimates["elpd_loo", "Estimate"],
    loo_P_2500_re$estimates["elpd_loo", "Estimate"],
    loo_lat$estimates["elpd_loo", "Estimate"],
    loo_lon$estimates["elpd_loo", "Estimate"],
    loo_XY$estimates["elpd_loo", "Estimate"],
    loo_pH$estimates["elpd_loo", "Estimate"],
    loo_SOM$estimates["elpd_loo", "Estimate"],
    loo_N$estimates["elpd_loo", "Estimate"],
    loo_claysilt$estimates["elpd_loo", "Estimate"],
    loo_P250_pH$estimates["elpd_loo", "Estimate"],
    loo_P250_SOM$estimates["elpd_loo", "Estimate"],
    loo_P250_N$estimates["elpd_loo", "Estimate"],
    loo_P250_clay$estimates["elpd_loo", "Estimate"],
    loo_P2500_pH$estimates["elpd_loo", "Estimate"],
    loo_P2500_SOM$estimates["elpd_loo", "Estimate"],
    loo_P2500_N$estimates["elpd_loo", "Estimate"],
    loo_P2500_clay$estimates["elpd_loo", "Estimate"]
  ),
  se_elpd_loo = c(
    loo_block$estimates["elpd_loo", "SE"],
    loo_P_100$estimates["elpd_loo", "SE"],
    loo_P_250$estimates["elpd_loo", "SE"],
    loo_P_500$estimates["elpd_loo", "SE"],
    loo_P_1000$estimates["elpd_loo", "SE"],
    loo_P_2500$estimates["elpd_loo", "SE"],
    loo_P_100_re$estimates["elpd_loo", "SE"],
    loo_P_250_re$estimates["elpd_loo", "SE"],
    loo_P_500_re$estimates["elpd_loo", "SE"],
    loo_P_1000_re$estimates["elpd_loo", "SE"],
    loo_P_2500_re$estimates["elpd_loo", "SE"],
    loo_lat$estimates["elpd_loo", "SE"],
    loo_lon$estimates["elpd_loo", "SE"],
    loo_XY$estimates["elpd_loo", "SE"],
    loo_pH$estimates["elpd_loo", "SE"],
    loo_SOM$estimates["elpd_loo", "SE"],
    loo_N$estimates["elpd_loo", "SE"],
    loo_claysilt$estimates["elpd_loo", "SE"],
    loo_P250_pH$estimates["elpd_loo", "SE"],
    loo_P250_SOM$estimates["elpd_loo", "SE"],
    loo_P250_N$estimates["elpd_loo", "SE"],
    loo_P250_clay$estimates["elpd_loo", "SE"],
    loo_P2500_pH$estimates["elpd_loo", "SE"],
    loo_P2500_SOM$estimates["elpd_loo", "SE"],
    loo_P2500_N$estimates["elpd_loo", "SE"],
    loo_P2500_clay$estimates["elpd_loo", "SE"]
  ),
  looic = c(
    loo_block$estimates["looic", "Estimate"],
    loo_P_100$estimates["looic", "Estimate"],
    loo_P_250$estimates["looic", "Estimate"],
    loo_P_500$estimates["looic", "Estimate"],
    loo_P_1000$estimates["looic", "Estimate"],
    loo_P_2500$estimates["looic", "Estimate"],
    loo_P_100_re$estimates["looic", "Estimate"],
    loo_P_250_re$estimates["looic", "Estimate"],
    loo_P_500_re$estimates["looic", "Estimate"],
    loo_P_1000_re$estimates["looic", "Estimate"],
    loo_P_2500_re$estimates["looic", "Estimate"],
    loo_lat$estimates["looic", "Estimate"],
    loo_lon$estimates["looic", "Estimate"],
    loo_XY$estimates["looic", "Estimate"],
    loo_pH$estimates["looic", "Estimate"],
    loo_SOM$estimates["looic", "Estimate"],
    loo_N$estimates["looic", "Estimate"],
    loo_claysilt$estimates["looic", "Estimate"],
    loo_P250_pH$estimates["looic", "Estimate"],
    loo_P250_SOM$estimates["looic", "Estimate"],
    loo_P250_N$estimates["looic", "Estimate"],
    loo_P250_clay$estimates["looic", "Estimate"],
    loo_P2500_pH$estimates["looic", "Estimate"],
    loo_P2500_SOM$estimates["looic", "Estimate"],
    loo_P2500_N$estimates["looic", "Estimate"],
    loo_P2500_clay$estimates["looic", "Estimate"]
  ),
  se_looic = c(
    loo_block$estimates["looic", "SE"],
    loo_P_100$estimates["looic", "SE"],
    loo_P_250$estimates["looic", "SE"],
    loo_P_500$estimates["looic", "SE"],
    loo_P_1000$estimates["looic", "SE"],
    loo_P_2500$estimates["looic", "SE"],
    loo_P_100_re$estimates["looic", "SE"],
    loo_P_250_re$estimates["looic", "SE"],
    loo_P_500_re$estimates["looic", "SE"],
    loo_P_1000_re$estimates["looic", "SE"],
    loo_P_2500_re$estimates["looic", "SE"],
    loo_lat$estimates["looic", "SE"],
    loo_lon$estimates["looic", "SE"],
    loo_XY$estimates["looic", "SE"],
    loo_pH$estimates["looic", "SE"],
    loo_SOM$estimates["looic", "SE"],
    loo_N$estimates["looic", "SE"],
    loo_claysilt$estimates["looic", "SE"],
    loo_P250_pH$estimates["looic", "SE"],
    loo_P250_SOM$estimates["looic", "SE"],
    loo_P250_N$estimates["looic", "SE"],
    loo_P250_clay$estimates["looic", "SE"],
    loo_P2500_pH$estimates["looic", "SE"],
    loo_P2500_SOM$estimates["looic", "SE"],
    loo_P2500_N$estimates["looic", "SE"],
    loo_P2500_clay$estimates["looic", "SE"]
  )
)

write.csv(loo_summary_table, file.path(outdir, "loo_summary_table.csv"), row.names = FALSE)

loo_comp <- loo_compare(
  loo_block,
  loo_P_100, loo_P_250, loo_P_500, loo_P_1000, loo_P_2500,
  loo_P_100_re, loo_P_250_re, loo_P_500_re, loo_P_1000_re, loo_P_2500_re,
  loo_lat, loo_lon, loo_XY, loo_pH, loo_SOM, loo_N, loo_claysilt,
  loo_P250_pH, loo_P250_SOM, loo_P250_N, loo_P250_clay,
  loo_P2500_pH, loo_P2500_SOM, loo_P2500_N, loo_P2500_clay
)

write.csv(as.data.frame(unclass(loo_comp)), file.path(outdir, "loo_compare_table.csv"))

# =============================================================================
# Save workspace objects
# =============================================================================

saveRDS(shared_asvs,         file.path(outdir, "shared_asvs_vector.rds"))
saveRDS(orchard_unique_asvs, file.path(outdir, "orchard_unique_asvs_vector.rds"))
saveRDS(forest_unique_asvs,  file.path(outdir, "forest_unique_asvs_vector.rds"))
saveRDS(orchard_summary_df,  file.path(outdir, "orchard_summary_ITS.rds"))
saveRDS(forest_summary_df,   file.path(outdir, "forest_summary_ITS.rds"))

writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo_ITS_hyp1.txt"))
