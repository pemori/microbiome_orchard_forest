# =============================================================================
# Script: 10_hypothesis_4_respiration_assays.R
# Project: From forest to fields: The role of soil microbiome spillover in
#          agroecosystem sustainability
# Author: Pedro Mondaca
# Contact: pedromondaca@outlook.com
# Submitted to: Science Advances
#
# Description:
# This script evaluates hypothesis 4 using soil respiration assays. It:
#   (i) classifies orchard samples into low, medium and high forest cover (LFC,
#       MFC, HFC) based on terciles of P_2500, while keeping forest samples as F,
#   (ii) prepares respiration variables,
#   (iii) summarizes rapid respiration and accumulated CO2 at 20, 30 and 40°C,
#   (iv) generates bar plots with standard errors,
#   (v) compares alternative statistical models for each respiration response,
#   (vi) performs post hoc contrasts among forest-cover categories.
#
# Required objects in the workspace:
#   - BGQ
#
# Main input variables expected in BGQ:
#   - Group
#   - P_2500
#   - pH, N
#   - resp_rapida, rrC, r20, r30, r40
#   - 20C, 30C, 40C
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(car)
  library(emmeans)
  library(performance)
  library(patchwork)
})

# =============================================================================
# 0. Output directories
# =============================================================================

outdir_h4 <- file.path(getwd(), "results_hyp4_respiration")
figdir_h4 <- file.path(outdir_h4, "figures")
tabdir_h4 <- file.path(outdir_h4, "tables")
rdsdir_h4 <- file.path(outdir_h4, "rds")

dir.create(outdir_h4, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir_h4, recursive = TRUE, showWarnings = FALSE)
dir.create(tabdir_h4, recursive = TRUE, showWarnings = FALSE)
dir.create(rdsdir_h4, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. Base data and forest-cover classes
# =============================================================================

map_crop <- BGQ

quantiles <- quantile(
  map_crop$P_2500[map_crop$Group == "Orchard"],
  probs = c(0, 1/3, 2/3, 1),
  na.rm = TRUE
)

quantiles

map_crop <- map_crop %>%
  mutate(
    FC = case_when(
      Group == "Forest" ~ "F",
      P_2500 <= quantiles[2] & Group == "Orchard" ~ "LFC",
      P_2500 <= quantiles[3] & Group == "Orchard" ~ "MFC",
      P_2500 >  quantiles[3] & Group == "Orchard" ~ "HFC",
      TRUE ~ NA_character_
    )
  )

map_crop$FC <- factor(map_crop$FC, levels = c("LFC", "MFC", "HFC", "F"))

write.csv(map_crop, file.path(tabdir_h4, "BGQ_with_FC.csv"), row.names = FALSE)

# =============================================================================
# 2. General summaries
# =============================================================================

summary_by_group <- map_crop %>%
  group_by(FC) %>%
  summarize(
    across(where(is.numeric), list(
      mediana = ~ median(.x, na.rm = TRUE),
      min = ~ min(.x, na.rm = TRUE),
      max = ~ max(.x, na.rm = TRUE)
    )),
    .groups = "drop"
  )

print(summary_by_group)
dplyr::glimpse(summary_by_group)

write.csv(summary_by_group,
          file.path(tabdir_h4, "summary_by_FC_numeric_variables.csv"),
          row.names = FALSE)

# =============================================================================
# 3. Respiration dataset
# =============================================================================

resp <- map_crop[, c("pH", "N", "P_2500", "resp_rapida", "rrC",
                     "r20", "r30", "r40", "20C", "30C", "40C", "FC")]

resp_filt <- na.omit(resp)
resp_filt$FC <- factor(resp_filt$FC, levels = c("LFC", "MFC", "HFC", "F"))

write.csv(resp_filt, file.path(tabdir_h4, "respiration_dataset_complete_cases.csv"),
          row.names = FALSE)

# =============================================================================
# 4. Color palette
# =============================================================================

resp_colors <- c(
  "LFC" = "#1f78b4",
  "MFC" = "gold",
  "HFC" = "#33a02c",
  "F"   = "#006400"
)

# =============================================================================
# 5. Summary statistics for plotting
# =============================================================================

summary_stats <- map_crop %>%
  group_by(FC) %>%
  summarise(
    rrC_mean = mean(rrC, na.rm = TRUE),
    rrC_se   = sd(rrC, na.rm = TRUE) / sqrt(sum(!is.na(rrC))),
    C20_mean = mean(`20C`, na.rm = TRUE),
    C20_se   = sd(`20C`, na.rm = TRUE) / sqrt(sum(!is.na(`20C`))),
    C30_mean = mean(`30C`, na.rm = TRUE),
    C30_se   = sd(`30C`, na.rm = TRUE) / sqrt(sum(!is.na(`30C`))),
    C40_mean = mean(`40C`, na.rm = TRUE),
    C40_se   = sd(`40C`, na.rm = TRUE) / sqrt(sum(!is.na(`40C`))),
    .groups = "drop"
  )

write.csv(summary_stats,
          file.path(tabdir_h4, "respiration_summary_stats.csv"),
          row.names = FALSE)

# =============================================================================
# 6. Respiration plots
# =============================================================================

# -----------------------------------------------------------------------------
# 6.1 Rapid respiration
# -----------------------------------------------------------------------------

plot_rrC <- ggplot(summary_stats, aes(x = FC, y = rrC_mean, fill = FC)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.8, width = 0.6) +
  geom_errorbar(
    aes(ymin = rrC_mean - rrC_se, ymax = rrC_mean + rrC_se),
    width = 0.2
  ) +
  scale_fill_manual(values = resp_colors) +
  labs(x = "", y = expression("Rapid Respiration (" ~ mu ~ "gCO"[2] ~ "/kg soil C)")) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12),
    legend.position = "none"
  )

plot_rrC
# plot_rrC # 350*500

ggsave(file.path(figdir_h4, "plot_rrC.tiff"),
       plot = plot_rrC, width = 3.5, height = 5, dpi = 600,
       compression = "lzw", device = "tiff", bg = "white")

# -----------------------------------------------------------------------------
# 6.2 Accumulated CO2 at 20°C
# -----------------------------------------------------------------------------

plot_20C <- ggplot(summary_stats, aes(x = FC, y = C20_mean, fill = FC)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.8, width = 0.6) +
  geom_errorbar(
    aes(ymin = C20_mean - C20_se, ymax = C20_mean + C20_se),
    width = 0.2
  ) +
  scale_fill_manual(values = resp_colors) +
  labs(x = "", y = expression("Accumulated CO"[2] ~ "(mg CO"[2] ~ " at 20°C)")) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12),
    legend.position = "none"
  )

plot_20C

ggsave(file.path(figdir_h4, "plot_20C.tiff"),
       plot = plot_20C, width = 3.5, height = 5, dpi = 600,
       compression = "lzw", device = "tiff", bg = "white")

# -----------------------------------------------------------------------------
# 6.3 Accumulated CO2 at 30°C
# -----------------------------------------------------------------------------

plot_30C <- ggplot(summary_stats, aes(x = FC, y = C30_mean, fill = FC)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.8, width = 0.6) +
  geom_errorbar(
    aes(ymin = C30_mean - C30_se, ymax = C30_mean + C30_se),
    width = 0.2
  ) +
  scale_fill_manual(values = resp_colors) +
  labs(x = "", y = expression("Accumulated CO"[2] ~ "(mg CO"[2] ~ " at 30°C)")) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12),
    legend.position = "none"
  )

plot_30C

ggsave(file.path(figdir_h4, "plot_30C.tiff"),
       plot = plot_30C, width = 3.5, height = 5, dpi = 600,
       compression = "lzw", device = "tiff", bg = "white")

# -----------------------------------------------------------------------------
# 6.4 Accumulated CO2 at 40°C
# -----------------------------------------------------------------------------

plot_40C <- ggplot(summary_stats, aes(x = FC, y = C40_mean, fill = FC)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.8, width = 0.6) +
  geom_errorbar(
    aes(ymin = C40_mean - C40_se, ymax = C40_mean + C40_se),
    width = 0.2
  ) +
  scale_fill_manual(values = resp_colors) +
  labs(x = "", y = expression("Accumulated CO"[2] ~ "(mg CO"[2] ~ " at 40°C)")) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12),
    legend.position = "none"
  )

plot_40C

ggsave(file.path(figdir_h4, "plot_40C.tiff"),
       plot = plot_40C, width = 3.5, height = 5, dpi = 600,
       compression = "lzw", device = "tiff", bg = "white")

# -----------------------------------------------------------------------------
# 6.5 Combined figure
# -----------------------------------------------------------------------------

final_plot <- plot_rrC | plot_20C | plot_30C | plot_40C
final_plot
# final_plot # 1200*500

ggsave(file.path(figdir_h4, "respiration_combined.tiff"),
       plot = final_plot, width = 12, height = 5, dpi = 600,
       compression = "lzw", device = "tiff", bg = "white")

# =============================================================================
# 7. Statistical analysis: rapid respiration (rrC)
# =============================================================================

rrC_mod_lm <- lm(rrC ~ FC, data = resp_filt)
rrC_mod_invgauss <- glm(rrC ~ FC, data = resp_filt,
                        family = inverse.gaussian(link = "log"))
rrC_mod_lognorm <- lm(log(rrC) ~ FC, data = resp_filt)

rrC_model_comparison <- compare_performance(
  rrC_mod_lm, rrC_mod_invgauss, rrC_mod_lognorm,
  rank = TRUE, verbose = TRUE
)

print(rrC_model_comparison)
check_model(rrC_mod_lognorm)

rrC_emm <- emmeans(rrC_mod_lognorm, pairwise ~ FC, type = "response")

capture.output(rrC_model_comparison,
               file = file.path(tabdir_h4, "rrC_model_comparison.txt"))
capture.output(summary(rrC_emm),
               file = file.path(tabdir_h4, "rrC_emmeans.txt"))

saveRDS(rrC_mod_lm,       file.path(rdsdir_h4, "rrC_mod_lm.rds"))
saveRDS(rrC_mod_invgauss, file.path(rdsdir_h4, "rrC_mod_invgauss.rds"))
saveRDS(rrC_mod_lognorm,  file.path(rdsdir_h4, "rrC_mod_lognorm.rds"))

# =============================================================================
# 8. Statistical analysis: accumulated CO2 at 20°C
# =============================================================================

C20_mod_lm <- lm(`20C` ~ FC, data = resp_filt)
C20_mod_invgauss <- glm(`20C` ~ FC, data = resp_filt,
                        family = inverse.gaussian(link = "log"))
C20_mod_lognorm <- lm(log(`20C`) ~ FC, data = resp_filt)

C20_model_comparison <- compare_performance(
  C20_mod_lm, C20_mod_invgauss, C20_mod_lognorm,
  rank = TRUE, verbose = TRUE
)

print(C20_model_comparison)
check_model(C20_mod_lognorm)

C20_emm <- emmeans(C20_mod_lognorm, pairwise ~ FC, type = "response")

capture.output(C20_model_comparison,
               file = file.path(tabdir_h4, "20C_model_comparison.txt"))
capture.output(summary(C20_emm),
               file = file.path(tabdir_h4, "20C_emmeans.txt"))

saveRDS(C20_mod_lm,       file.path(rdsdir_h4, "20C_mod_lm.rds"))
saveRDS(C20_mod_invgauss, file.path(rdsdir_h4, "20C_mod_invgauss.rds"))
saveRDS(C20_mod_lognorm,  file.path(rdsdir_h4, "20C_mod_lognorm.rds"))

# =============================================================================
# 9. Statistical analysis: accumulated CO2 at 30°C
# =============================================================================

C30_mod_lm <- lm(`30C` ~ FC, data = resp_filt)
C30_mod_invgauss <- glm(`30C` ~ FC, data = resp_filt,
                        family = inverse.gaussian(link = "log"))
C30_mod_lognorm <- lm(log(`30C`) ~ FC, data = resp_filt)

C30_model_comparison <- compare_performance(
  C30_mod_lm, C30_mod_invgauss, C30_mod_lognorm,
  rank = TRUE, verbose = TRUE
)

print(C30_model_comparison)
check_model(C30_mod_lognorm)

C30_emm <- emmeans(C30_mod_lognorm, pairwise ~ FC, type = "response")

capture.output(C30_model_comparison,
               file = file.path(tabdir_h4, "30C_model_comparison.txt"))
capture.output(summary(C30_emm),
               file = file.path(tabdir_h4, "30C_emmeans.txt"))

saveRDS(C30_mod_lm,       file.path(rdsdir_h4, "30C_mod_lm.rds"))
saveRDS(C30_mod_invgauss, file.path(rdsdir_h4, "30C_mod_invgauss.rds"))
saveRDS(C30_mod_lognorm,  file.path(rdsdir_h4, "30C_mod_lognorm.rds"))

# =============================================================================
# 10. Statistical analysis: accumulated CO2 at 40°C
# =============================================================================

C40_mod_lm <- lm(`40C` ~ FC, data = resp_filt)
C40_mod_invgauss <- glm(`40C` ~ FC, data = resp_filt,
                        family = inverse.gaussian(link = "log"))
C40_mod_lognorm <- lm(log(`40C`) ~ FC, data = resp_filt)

C40_model_comparison <- compare_performance(
  C40_mod_lm, C40_mod_invgauss, C40_mod_lognorm,
  rank = TRUE, verbose = TRUE
)

print(C40_model_comparison)
check_model(C40_mod_lognorm)

C40_emm <- emmeans(C40_mod_lognorm, pairwise ~ FC, type = "response")

capture.output(C40_model_comparison,
               file = file.path(tabdir_h4, "40C_model_comparison.txt"))
capture.output(summary(C40_emm),
               file = file.path(tabdir_h4, "40C_emmeans.txt"))

saveRDS(C40_mod_lm,       file.path(rdsdir_h4, "40C_mod_lm.rds"))
saveRDS(C40_mod_invgauss, file.path(rdsdir_h4, "40C_mod_invgauss.rds"))
saveRDS(C40_mod_lognorm,  file.path(rdsdir_h4, "40C_mod_lognorm.rds"))

# =============================================================================
# 11. Session info
# =============================================================================

writeLines(capture.output(sessionInfo()),
           con = file.path(outdir_h4, "sessionInfo_hypothesis_4.txt"))
