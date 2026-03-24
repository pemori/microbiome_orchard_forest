# =============================================================================
# Script: 11_hypothesis_4_qpcr_sem.R
# Project: From forest to fields: The role of soil microbiome spillover in
#          agroecosystem sustainability
# Author: Pedro Mondaca
# Contact: pedromondaca@outlook.com
# Submitted to: Science Advances
#
# Description:
# This script evaluates qPCR-based microbial abundance patterns and their
# relationship with respiration and niche-breadth categories. It:
#   (i) prepares bacterial and fungal niche-breadth matrices,
#   (ii) removes qPCR_ITS outliers using the IQR rule,
#   (iii) classifies samples by forest cover category,
#   (iv) plots qPCR 16S, qPCR ITS and fungal:bacterial biomass ratio,
#   (v) tests differences among forest-cover categories,
#   (vi) builds a combined data frame for SEM analyses,
#   (vii) fits a final piecewise SEM linking respiration, qPCR and generalist taxa.
#
# Required objects in the workspace:
#   - nb_bac_df
#   - nb_fun_df
#   - mapping
#
# Expected variables in `mapping`:
#   - Group, P_2500
#   - qPCR_16s, qPCR_ITS
#   - pH, SOM, N, C, CN, P, clay, claysilt
#   - resp_rapida, rrC, r20, r30, r40, 20C, 30C, 40C
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(emmeans)
  library(performance)
  library(piecewiseSEM)
  library(patchwork)
})

# =============================================================================
# 0. Output directories
# =============================================================================

outdir_h4_sem <- file.path(getwd(), "results_hyp4_qpcr_sem")
figdir_h4_sem <- file.path(outdir_h4_sem, "figures")
tabdir_h4_sem <- file.path(outdir_h4_sem, "tables")
rdsdir_h4_sem <- file.path(outdir_h4_sem, "rds")

dir.create(outdir_h4_sem, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir_h4_sem, recursive = TRUE, showWarnings = FALSE)
dir.create(tabdir_h4_sem, recursive = TRUE, showWarnings = FALSE)
dir.create(rdsdir_h4_sem, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. Prepare niche-breadth matrices
# =============================================================================
# =============================================================================
# 1. Load niche-breadth summary tables from hypothesis 3
# =============================================================================

nb_bac_df <- read.csv(
  "results_16S_hyp3/nb_bacterial_model_dataframe.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

nb_fun_df <- read.csv(
  "results_ITS_hyp3/nb_fungal_model_dataframe.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

stopifnot("SampleID" %in% colnames(nb_bac_df))
stopifnot("SampleID" %in% colnames(nb_fun_df))

rownames(nb_bac_df) <- nb_bac_df$SampleID
rownames(nb_fun_df) <- nb_fun_df$SampleID

bac_rel_nb <- nb_bac_df %>%
  dplyr::select(SampleID, Generalist, Opportunist, Specialist)

rownames(bac_rel_nb) <- bac_rel_nb$SampleID
bac_rel_nb$SampleID <- NULL

colnames(bac_rel_nb)[colnames(bac_rel_nb) == "Generalist"]  <- "Generalist_b"
colnames(bac_rel_nb)[colnames(bac_rel_nb) == "Opportunist"] <- "Opportunist_b"
colnames(bac_rel_nb)[colnames(bac_rel_nb) == "Specialist"]  <- "Specialist_b"

fun_rel_nb <- nb_fun_df %>%
  dplyr::select(SampleID, Generalist, Opportunist, Specialist)

rownames(fun_rel_nb) <- fun_rel_nb$SampleID
fun_rel_nb$SampleID <- NULL

colnames(fun_rel_nb)[colnames(fun_rel_nb) == "Generalist"]  <- "Generalist_f"
colnames(fun_rel_nb)[colnames(fun_rel_nb) == "Opportunist"] <- "Opportunist_f"
colnames(fun_rel_nb)[colnames(fun_rel_nb) == "Specialist"]  <- "Specialist_f"


# =============================================================================
# 2. Base mapping object and outlier filtering
# =============================================================================

plot(mapping$qPCR_16s)
plot(mapping$qPCR_ITS)

q1 <- quantile(mapping$qPCR_ITS, 0.25, na.rm = TRUE)
q3 <- quantile(mapping$qPCR_ITS, 0.75, na.rm = TRUE)
iqr <- q3 - q1

lower_bound <- q1 - 1.5 * iqr
upper_bound <- q3 + 1.5 * iqr

outliers_qPCR_ITS <- mapping$qPCR_ITS[
  mapping$qPCR_ITS < lower_bound | mapping$qPCR_ITS > upper_bound
]

print(outliers_qPCR_ITS)

mapping_clean <- mapping[
  mapping$qPCR_ITS >= lower_bound & mapping$qPCR_ITS <= upper_bound, 
]

plot(mapping_clean$qPCR_ITS)

write.csv(mapping_clean,
          file.path(tabdir_h4_sem, "mapping_clean_qPCRITS_outliers_removed.csv"),
          row.names = TRUE)

# =============================================================================
# 3. Forest-cover categories and fungal:bacterial ratio
# =============================================================================

mapping_clean <- mapping_clean %>%
  mutate(
    FC = case_when(
      Group == "Forest" ~ "F",
      P_2500 <= 0.3584376 ~ "LFC",
      P_2500 <= 0.5224324 ~ "MFC",
      P_2500 >  0.5224324 ~ "HFC"
    )
  )

mapping_clean$FC <- factor(mapping_clean$FC, levels = c("LFC", "MFC", "HFC", "F"))

mapping_clean <- mapping_clean %>%
  mutate(fun_bac = log(qPCR_ITS) / log(qPCR_16s))

write.csv(mapping_clean,
          file.path(tabdir_h4_sem, "mapping_clean_with_FC_fun_bac.csv"),
          row.names = TRUE)

# =============================================================================
# 4. Plot settings
# =============================================================================

resp_colors <- c(
  "LFC" = "#1f78b4",
  "MFC" = "gold",
  "HFC" = "#33a02c",
  "F"   = "#006400"
)

# =============================================================================
# 5. qPCR 16S by forest-cover category
# =============================================================================

df_summary_16S <- mapping_clean %>%
  group_by(FC) %>%
  summarise(
    mean_qPCR = mean(qPCR_16s, na.rm = TRUE),
    se_qPCR   = sd(qPCR_16s, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

df_summary_16S$FC <- factor(df_summary_16S$FC, levels = c("LFC", "MFC", "HFC", "F"))

plot_qPCR_16S <- ggplot(df_summary_16S, aes(x = FC, y = mean_qPCR, fill = FC)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.1),
    width = 0.9,
    color = "black",
    alpha = 0.8
  ) +
  geom_errorbar(
    aes(ymin = mean_qPCR - se_qPCR, ymax = mean_qPCR + se_qPCR),
    position = position_dodge(width = 0.1),
    width = 0.2,
    color = "black"
  ) +
  scale_fill_manual(values = resp_colors) +
  labs(x = "", y = expression("qPCR 16S (copies/g)")) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )

plot_qPCR_16S
# plot_qPCR_16S # 250*450

ggsave(file.path(figdir_h4_sem, "plot_qPCR_16S.tiff"),
       plot = plot_qPCR_16S, width = 2.5, height = 4.5, dpi = 600,
       compression = "lzw", device = "tiff", bg = "white")

# =============================================================================
# 6. qPCR ITS by forest-cover category
# =============================================================================

df_summary_ITS <- mapping_clean %>%
  group_by(FC) %>%
  summarise(
    mean_qPCR = mean(qPCR_ITS, na.rm = TRUE),
    se_qPCR   = sd(qPCR_ITS, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

df_summary_ITS$FC <- factor(df_summary_ITS$FC, levels = c("LFC", "MFC", "HFC", "F"))

plot_qPCR_ITS <- ggplot(df_summary_ITS, aes(x = FC, y = mean_qPCR, fill = FC)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.1),
    width = 0.9,
    color = "black",
    alpha = 0.8
  ) +
  geom_errorbar(
    aes(ymin = mean_qPCR - se_qPCR, ymax = mean_qPCR + se_qPCR),
    position = position_dodge(width = 0.1),
    width = 0.2,
    color = "black"
  ) +
  scale_fill_manual(values = resp_colors) +
  labs(x = "", y = expression("qPCR ITS (copies/g)")) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )

plot_qPCR_ITS
# plot_qPCR_ITS # 250*450

ggsave(file.path(figdir_h4_sem, "plot_qPCR_ITS.tiff"),
       plot = plot_qPCR_ITS, width = 2.5, height = 4.5, dpi = 600,
       compression = "lzw", device = "tiff", bg = "white")

# =============================================================================
# 7. Fungal:bacterial ratio by forest-cover category
# =============================================================================

df_summary_funbac <- mapping_clean %>%
  group_by(FC) %>%
  summarise(
    mean_qPCR = mean(fun_bac, na.rm = TRUE),
    se_qPCR   = sd(fun_bac, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

df_summary_funbac$FC <- factor(df_summary_funbac$FC, levels = c("LFC", "MFC", "HFC", "F"))

plot_fun_bac <- ggplot(df_summary_funbac, aes(x = FC, y = mean_qPCR, fill = FC)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.1),
    width = 0.9,
    color = "black",
    alpha = 0.8
  ) +
  geom_errorbar(
    aes(ymin = mean_qPCR - se_qPCR, ymax = mean_qPCR + se_qPCR),
    position = position_dodge(width = 0.1),
    width = 0.2,
    color = "black"
  ) +
  scale_fill_manual(values = resp_colors) +
  labs(x = "", y = expression("log(qPCR ITS)/log(qPCR16S)")) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )

plot_fun_bac
# plot_fun_bac # 583*505

ggsave(file.path(figdir_h4_sem, "plot_fun_bac_ratio.tiff"),
       plot = plot_fun_bac, width = 5.83, height = 5.05, dpi = 600,
       compression = "lzw", device = "tiff", bg = "white")

# =============================================================================
# 8. Group comparisons for qPCR variables
# =============================================================================

mod_16 <- lm(log(qPCR_16s) ~ FC, data = mapping_clean)
emm_16 <- emmeans(mod_16, pairwise ~ FC, type = "response")

mod_IT <- lm(log(qPCR_ITS) ~ FC, data = mapping_clean)
emm_IT <- emmeans(mod_IT, pairwise ~ FC, type = "response")

mod_fb <- lm(fun_bac ~ FC, data = mapping_clean)
emm_fb <- emmeans(mod_fb, pairwise ~ FC, type = "response")

capture.output(summary(emm_16),
               file = file.path(tabdir_h4_sem, "emmeans_qPCR_16S.txt"))
capture.output(summary(emm_IT),
               file = file.path(tabdir_h4_sem, "emmeans_qPCR_ITS.txt"))
capture.output(summary(emm_fb),
               file = file.path(tabdir_h4_sem, "emmeans_fun_bac_ratio.txt"))

saveRDS(mod_16, file.path(rdsdir_h4_sem, "mod_qPCR_16S.rds"))
saveRDS(mod_IT, file.path(rdsdir_h4_sem, "mod_qPCR_ITS.rds"))
saveRDS(mod_fb, file.path(rdsdir_h4_sem, "mod_fun_bac_ratio.rds"))

# =============================================================================
# 9. Build SEM dataset
# =============================================================================

samples_bac <- rownames(bac_rel_nb)
samples_fun <- rownames(fun_rel_nb)
samples_map <- rownames(mapping_clean)

samples_common <- Reduce(intersect, list(samples_bac, samples_fun, samples_map))

bac_filt <- bac_rel_nb[samples_common, ]
fun_filt <- fun_rel_nb[samples_common, ]
map_filt <- mapping_clean[samples_common, ]

df_sem <- cbind(map_filt, bac_filt, fun_filt)

df_sem <- df_sem %>%
  filter(
    !is.na(rrC),
    !is.na(`20C`),
    !is.na(`30C`),
    !is.na(`40C`)
  )

df_sem <- df_sem %>%
  mutate(
    log_qPCR_16s = log(qPCR_16s),
    log_qPCR_ITS = log(qPCR_ITS)
  )

df_sem <- df_sem %>%
  mutate(
    OppGen_f   = Opportunist_f + Generalist_f,
    bact_gen   = log(qPCR_16s) * Generalist_b,
    fung_gen   = log(qPCR_ITS) * (Generalist_f + Opportunist_f),
    fun_bac    = log(qPCR_ITS) / log(qPCR_16s),
    rr         = resp_rapida,
    LT         = r20,
    MT         = r30,
    HT         = r40,
    LTC        = `20C`,
    MTC        = `30C`,
    HTC        = `40C`,
    gen_esp_b  = Generalist_b / (Opportunist_b + Specialist_b),
    gen_esp_f  = Generalist_f / Specialist_f,
    qgen_esp_b = (Generalist_b / (Opportunist_b + Specialist_b)) * log(qPCR_16s),
    qgen_esp_f = (Generalist_f / Specialist_f) * log(qPCR_ITS)
  )

write.csv(df_sem,
          file.path(tabdir_h4_sem, "df_sem.csv"),
          row.names = TRUE)

# =============================================================================
# 10. Exploratory plots for SEM variables
# =============================================================================

png(file.path(figdir_h4_sem, "exploratory_rrC_vs_N.png"), width = 1200, height = 1200, res = 200)
plot(df_sem$rrC, df_sem$N)
dev.off()

png(file.path(figdir_h4_sem, "exploratory_rrC_vs_Generalist_b.png"), width = 1200, height = 1200, res = 200)
plot(df_sem$rrC, df_sem$Generalist_b)
dev.off()

png(file.path(figdir_h4_sem, "exploratory_rrC_vs_gen_esp_b.png"), width = 1200, height = 1200, res = 200)
plot(df_sem$rrC, df_sem$gen_esp_b)
dev.off()

png(file.path(figdir_h4_sem, "exploratory_rrC_vs_Generalist_f.png"), width = 1200, height = 1200, res = 200)
plot(df_sem$rrC, df_sem$Generalist_f)
dev.off()

png(file.path(figdir_h4_sem, "exploratory_logqPCRITS_vs_N.png"), width = 1200, height = 1200, res = 200)
plot(df_sem$log_qPCR_ITS, df_sem$N)
dev.off()

# =============================================================================
# 11. Final SEM
# =============================================================================

# Based on the latest explicit SEM block in the uploaded workflow:
# HTC ~ Generalist_b + Generalist_f + log_qPCR_16s + log_qPCR_ITS + clay
# Generalist_b ~ pH + P_2500
# Generalist_f ~ pH + P_2500
# log_qPCR_16s ~ C + N + P
# log_qPCR_ITS ~ C + N + P
# plus correlated residuals among microbial predictors
# This final structure follows the last model block in the provided script.
# :contentReference[oaicite:3]{index=3}

#rapid respiration
m1 <- lm(rrC ~ Generalist_b + Generalist_f + log_qPCR_16s + log_qPCR_ITS + clay,
         data = df_sem)
m2 <- lm(Generalist_b ~ pH + P_2500, data = df_sem)
m3 <- lm(Generalist_f ~ pH + P_2500, data = df_sem)
m4 <- lm(log_qPCR_16s ~ C + N + P, data = df_sem)
m5 <- lm(log_qPCR_ITS ~ C + N + P, data = df_sem)

sem_mod <- psem(
  m1, m2, m3, m4, m5,
  Generalist_b %~~% Generalist_f,
  log_qPCR_16s %~~% log_qPCR_ITS,
  Generalist_b %~~% log_qPCR_16s,
  Generalist_f %~~% log_qPCR_16s,
  Generalist_b %~~% log_qPCR_ITS,
  Generalist_f %~~% log_qPCR_ITS,
  df_sem
)

sem_summary <- summary(sem_mod)
sem_summary_conserve <- summary(sem_mod, conserve = TRUE)

print(sem_summary)
print(sem_summary_conserve)

capture.output(sem_summary,
               file = file.path(tabdir_h4_sem, "sem_summary_rrC.txt"))
capture.output(sem_summary_conserve,
               file = file.path(tabdir_h4_sem, "sem_summary_rrC_conservative.txt"))

saveRDS(m1, file.path(rdsdir_h4_sem, "sem_m1_rrC.rds"))
saveRDS(m2, file.path(rdsdir_h4_sem, "sem_m2_Generalist_b_rrC.rds"))
saveRDS(m3, file.path(rdsdir_h4_sem, "sem_m3_Generalist_f_rrC.rds"))
saveRDS(m4, file.path(rdsdir_h4_sem, "sem_m4_log_qPCR_16s_rrC.rds"))
saveRDS(m5, file.path(rdsdir_h4_sem, "sem_m5_log_qPCR_ITS_rrC.rds"))
saveRDS(sem_mod, file.path(rdsdir_h4_sem, "sem_mod_rrC.rds"))


#incubation at 20C
m1 <- lm(LTC ~ Generalist_b + Generalist_f + log_qPCR_16s + log_qPCR_ITS + clay,
         data = df_sem)
m2 <- lm(Generalist_b ~ pH + P_2500, data = df_sem)
m3 <- lm(Generalist_f ~ pH + P_2500, data = df_sem)
m4 <- lm(log_qPCR_16s ~ C + N + P, data = df_sem)
m5 <- lm(log_qPCR_ITS ~ C + N + P, data = df_sem)

sem_mod <- psem(
  m1, m2, m3, m4, m5,
  Generalist_b %~~% Generalist_f,
  log_qPCR_16s %~~% log_qPCR_ITS,
  Generalist_b %~~% log_qPCR_16s,
  Generalist_f %~~% log_qPCR_16s,
  Generalist_b %~~% log_qPCR_ITS,
  Generalist_f %~~% log_qPCR_ITS,
  df_sem
)

sem_summary <- summary(sem_mod)
sem_summary_conserve <- summary(sem_mod, conserve = TRUE)

print(sem_summary)
print(sem_summary_conserve)

capture.output(sem_summary,
               file = file.path(tabdir_h4_sem, "sem_summary_LTC.txt"))
capture.output(sem_summary_conserve,
               file = file.path(tabdir_h4_sem, "sem_summary_LTC_conservative.txt"))

saveRDS(m1, file.path(rdsdir_h4_sem, "sem_m1_LTC.rds"))
saveRDS(m2, file.path(rdsdir_h4_sem, "sem_m2_Generalist_b_LTC.rds"))
saveRDS(m3, file.path(rdsdir_h4_sem, "sem_m3_Generalist_f_LTC.rds"))
saveRDS(m4, file.path(rdsdir_h4_sem, "sem_m4_log_qPCR_16s_LTC.rds"))
saveRDS(m5, file.path(rdsdir_h4_sem, "sem_m5_log_qPCR_ITS_LTC.rds"))
saveRDS(sem_mod, file.path(rdsdir_h4_sem, "sem_mod_LTC.rds"))

#incubation at 30C
m1 <- lm(MTC ~ Generalist_b + Generalist_f + log_qPCR_16s + log_qPCR_ITS + clay,
         data = df_sem)
m2 <- lm(Generalist_b ~ pH + P_2500, data = df_sem)
m3 <- lm(Generalist_f ~ pH + P_2500, data = df_sem)
m4 <- lm(log_qPCR_16s ~ C + N + P, data = df_sem)
m5 <- lm(log_qPCR_ITS ~ C + N + P, data = df_sem)

sem_mod <- psem(
  m1, m2, m3, m4, m5,
  Generalist_b %~~% Generalist_f,
  log_qPCR_16s %~~% log_qPCR_ITS,
  Generalist_b %~~% log_qPCR_16s,
  Generalist_f %~~% log_qPCR_16s,
  Generalist_b %~~% log_qPCR_ITS,
  Generalist_f %~~% log_qPCR_ITS,
  df_sem
)

sem_summary <- summary(sem_mod)
sem_summary_conserve <- summary(sem_mod, conserve = TRUE)

print(sem_summary)
print(sem_summary_conserve)

capture.output(sem_summary,
               file = file.path(tabdir_h4_sem, "sem_summary_MTC.txt"))
capture.output(sem_summary_conserve,
               file = file.path(tabdir_h4_sem, "sem_summary_MTC_conservative.txt"))

saveRDS(m1, file.path(rdsdir_h4_sem, "sem_m1_MTC.rds"))
saveRDS(m2, file.path(rdsdir_h4_sem, "sem_m2_Generalist_b_MTC.rds"))
saveRDS(m3, file.path(rdsdir_h4_sem, "sem_m3_Generalist_f_MTC.rds"))
saveRDS(m4, file.path(rdsdir_h4_sem, "sem_m4_log_qPCR_16s_MTC.rds"))
saveRDS(m5, file.path(rdsdir_h4_sem, "sem_m5_log_qPCR_ITS_MTC.rds"))
saveRDS(sem_mod, file.path(rdsdir_h4_sem, "sem_mod_MTC.rds"))

#incubation at 40C
m1 <- lm(HTC ~ Generalist_b + Generalist_f + log_qPCR_16s + log_qPCR_ITS + clay,
         data = df_sem)
m2 <- lm(Generalist_b ~ pH + P_2500, data = df_sem)
m3 <- lm(Generalist_f ~ pH + P_2500, data = df_sem)
m4 <- lm(log_qPCR_16s ~ C + N + P, data = df_sem)
m5 <- lm(log_qPCR_ITS ~ C + N + P, data = df_sem)

sem_mod <- psem(
  m1, m2, m3, m4, m5,
  Generalist_b %~~% Generalist_f,
  log_qPCR_16s %~~% log_qPCR_ITS,
  Generalist_b %~~% log_qPCR_16s,
  Generalist_f %~~% log_qPCR_16s,
  Generalist_b %~~% log_qPCR_ITS,
  Generalist_f %~~% log_qPCR_ITS,
  df_sem
)

sem_summary <- summary(sem_mod)
sem_summary_conserve <- summary(sem_mod, conserve = TRUE)

print(sem_summary)
print(sem_summary_conserve)

capture.output(sem_summary,
               file = file.path(tabdir_h4_sem, "sem_summary_HTC.txt"))
capture.output(sem_summary_conserve,
               file = file.path(tabdir_h4_sem, "sem_summary_HTC_conservative.txt"))

saveRDS(m1, file.path(rdsdir_h4_sem, "sem_m1_HTC.rds"))
saveRDS(m2, file.path(rdsdir_h4_sem, "sem_m2_Generalist_b_HTC.rds"))
saveRDS(m3, file.path(rdsdir_h4_sem, "sem_m3_Generalist_f_HTC.rds"))
saveRDS(m4, file.path(rdsdir_h4_sem, "sem_m4_log_qPCR_16s_HTC.rds"))
saveRDS(m5, file.path(rdsdir_h4_sem, "sem_m5_log_qPCR_ITS_HTC.rds"))
saveRDS(sem_mod, file.path(rdsdir_h4_sem, "sem_mod_HTC.rds"))


# =============================================================================
# 12. Optional alternative exploratory SEMs retained as commented reference
# =============================================================================

# m1 <- lm(LT ~ bact_gen + fung_gen + N, data = df_sem)
# m2 <- lm(bact_gen ~ pH + SOM + N + claysilt + P_2500, data = df_sem)
# m3 <- lm(fung_gen ~ pH + SOM + N + claysilt + P_2500, data = df_sem)
# sem_mod_LT <- psem(m1, m2, m3, bact_gen %~~% fung_gen, df_sem)

# m1 <- lm(MT ~ Generalist_b + Generalist_f + CN, data = df_sem)
# m2 <- glm(Generalist_b ~ pH + C + P_2500, data = df_sem, family = binomial(link = "logit"))
# m3 <- glm(Generalist_f ~ pH + C + P_2500, data = df_sem, family = binomial(link = "logit"))

# m1 <- lm(LTC ~ Generalist_b + Generalist_f + claysilt + log(qPCR_16s) + log(qPCR_ITS), data = df_sem)
# m2 <- glm(Generalist_b ~ pH + P_2500, data = df_sem, family = binomial(link = "logit"))
# m3 <- glm(Generalist_f ~ pH + P_2500, data = df_sem, family = binomial(link = "logit"))

# =============================================================================
# 13. Session info
# =============================================================================

writeLines(capture.output(sessionInfo()),
           con = file.path(outdir_h4_sem, "sessionInfo_hypothesis_4_qpcr_sem.txt"))

