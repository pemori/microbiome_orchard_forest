# =============================================================================
# Script: 12_biogeochemical_differences.R
# Project: From forest to fields: The role of soil microbiome spillover in
#          agroecosystem sustainability
# Author: Pedro Mondaca
# Contact: pedromondaca@outlook.com
# Submitted to: Science Advances
#
# Description:
# This script evaluates biogeochemical differences:
#   (i) between Forest and Orchard soils, and
#   (ii) among Orchard soils classified by surrounding forest cover (LFC, MFC, HFC).
#
# It:
#   - loads and prepares the BGQ dataset,
#   - computes summary statistics for Tables S1 and S2,
#   - performs Welch ANOVA tests,
#   - performs Games-Howell post hoc comparisons,
#   - exports summary and statistical results.
#
# Input file:
#   - BGQ5.xlsx, sheet = "All"
#
# Expected variables:
#   - Group
#   - P_2500
#   - pH, SOM, N, C, CN, P, clay, claysilt, SOMTXT
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(openxlsx)
  library(userfriendlyscience)
})

# =============================================================================
# 0. Output directories
# =============================================================================

outdir_bgq <- file.path(getwd(), "results_biogeochemical_differences")
tabdir_bgq <- file.path(outdir_bgq, "tables")
rdsdir_bgq <- file.path(outdir_bgq, "rds")

dir.create(outdir_bgq, recursive = TRUE, showWarnings = FALSE)
dir.create(tabdir_bgq, recursive = TRUE, showWarnings = FALSE)
dir.create(rdsdir_bgq, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. Load data
# =============================================================================

BGQ <- read_excel("BGQ5.xlsx", sheet = "All")
BGQ <- as.data.frame(BGQ)

BGQ <- BGQ %>%
  mutate(across(where(is.character), as.factor))

bgq_vars <- c("pH", "SOM", "N", "C", "CN", "P", "clay", "claysilt", "SOMTXT")
BGQ[, bgq_vars] <- lapply(BGQ[, bgq_vars], as.numeric)

write.csv(BGQ, file.path(tabdir_bgq, "BGQ_clean.csv"), row.names = FALSE)

# =============================================================================
# 2. Utility functions
# =============================================================================

calculate_SE <- function(x) {
  n <- sum(!is.na(x))
  if (n > 1) {
    sd(x, na.rm = TRUE) / sqrt(n)
  } else {
    NA_real_
  }
}

run_games_howell <- function(data, response, group_var) {
  x <- data[[response]]
  g <- data[[group_var]]
  out <- userfriendlyscience::posthocTGH(x, g, method = "games-howell")
  capture.output(out)
}

# =============================================================================
# 3. Table S1: Forest vs Orchard
# =============================================================================

mapping_group <- BGQ %>%
  filter(Group %in% c("Forest", "Orchard")) %>%
  droplevels()

# -----------------------------------------------------------------------------
# 3.1 Summary table (Table S1)
# -----------------------------------------------------------------------------

summary_S1 <- mapping_group %>%
  group_by(Group) %>%
  summarise(
    n = n(),
    missing_pH = sum(is.na(pH)),
    mean_pH = mean(pH, na.rm = TRUE),
    median_pH = median(pH, na.rm = TRUE),
    min_pH = min(pH, na.rm = TRUE),
    max_pH = max(pH, na.rm = TRUE),
    
    missing_SOM = sum(is.na(SOM)),
    mean_SOM = mean(SOM, na.rm = TRUE),
    median_SOM = median(SOM, na.rm = TRUE),
    min_SOM = min(SOM, na.rm = TRUE),
    max_SOM = max(SOM, na.rm = TRUE),
    
    missing_C = sum(is.na(C)),
    mean_C = mean(C, na.rm = TRUE),
    median_C = median(C, na.rm = TRUE),
    min_C = min(C, na.rm = TRUE),
    max_C = max(C, na.rm = TRUE),
    
    missing_N = sum(is.na(N)),
    mean_N = mean(N, na.rm = TRUE),
    median_N = median(N, na.rm = TRUE),
    min_N = min(N, na.rm = TRUE),
    max_N = max(N, na.rm = TRUE),
    
    missing_CN = sum(is.na(CN)),
    mean_CN = mean(CN, na.rm = TRUE),
    median_CN = median(CN, na.rm = TRUE),
    min_CN = min(CN, na.rm = TRUE),
    max_CN = max(CN, na.rm = TRUE),
    
    missing_P = sum(is.na(P)),
    mean_P = mean(P, na.rm = TRUE),
    median_P = median(P, na.rm = TRUE),
    min_P = min(P, na.rm = TRUE),
    max_P = max(P, na.rm = TRUE),
    
    missing_clay = sum(is.na(clay)),
    mean_clay = mean(clay, na.rm = TRUE),
    median_clay = median(clay, na.rm = TRUE),
    min_clay = min(clay, na.rm = TRUE),
    max_clay = max(clay, na.rm = TRUE),
    
    missing_claysilt = sum(is.na(claysilt)),
    mean_claysilt = mean(claysilt, na.rm = TRUE),
    median_claysilt = median(claysilt, na.rm = TRUE),
    min_claysilt = min(claysilt, na.rm = TRUE),
    max_claysilt = max(claysilt, na.rm = TRUE),
    
    missing_SOMTXT = sum(is.na(SOMTXT)),
    mean_SOMTXT = mean(SOMTXT, na.rm = TRUE),
    median_SOMTXT = median(SOMTXT, na.rm = TRUE),
    min_SOMTXT = min(SOMTXT, na.rm = TRUE),
    max_SOMTXT = max(SOMTXT, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(summary_S1, file.path(tabdir_bgq, "Table_S1_summary_forest_vs_orchard.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 3.2 Mean ± SE table
# -----------------------------------------------------------------------------

summary_S1_mean_se <- mapping_group %>%
  group_by(Group) %>%
  summarise(
    across(
      all_of(bgq_vars),
      list(mean = ~ mean(.x, na.rm = TRUE),
           SE   = ~ calculate_SE(.x)),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

write.csv(summary_S1_mean_se, file.path(tabdir_bgq, "Table_S1_mean_se_forest_vs_orchard.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 3.3 Welch ANOVA + Games-Howell
# -----------------------------------------------------------------------------

welch_S1 <- lapply(bgq_vars, function(v) {
  fit <- oneway.test(reformulate("Group", response = v), data = mapping_group, var.equal = FALSE)
  data.frame(
    Variable = v,
    Statistic = unname(fit$statistic),
    df_num = unname(fit$parameter[1]),
    df_den = unname(fit$parameter[2]),
    p_value = fit$p.value
  )
}) %>% bind_rows()

write.csv(welch_S1, file.path(tabdir_bgq, "Table_S1_welch_tests.csv"), row.names = FALSE)

games_howell_S1 <- lapply(bgq_vars, function(v) {
  out <- run_games_howell(mapping_group, v, "Group")
  data.frame(
    Variable = v,
    Output = out,
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

write.csv(games_howell_S1, file.path(tabdir_bgq, "Table_S1_games_howell.csv"), row.names = FALSE)

# =============================================================================
# 4. Table S11: Orchard soils by forest-cover category
# =============================================================================

map_crop <- BGQ %>%
  filter(Group == "Orchard")

quantiles <- quantile(
  map_crop$P_2500,
  probs = c(0, 1/3, 2/3, 1),
  na.rm = TRUE
)

map_crop <- map_crop %>%
  mutate(
    FC = case_when(
      P_2500 <= quantiles[2] ~ "LFC",
      P_2500 <= quantiles[3] ~ "MFC",
      P_2500 >  quantiles[3] ~ "HFC",
      TRUE ~ NA_character_
    )
  )

map_crop$FC <- factor(map_crop$FC, levels = c("LFC", "MFC", "HFC"))

write.csv(map_crop, file.path(tabdir_bgq, "orchard_BGQ_with_FC.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 4.1 Summary table (Table S11)
# -----------------------------------------------------------------------------

summary_S11 <- map_crop %>%
  group_by(FC) %>%
  summarise(
    n = n(),
    missing_pH = sum(is.na(pH)),
    mean_pH = mean(pH, na.rm = TRUE),
    median_pH = median(pH, na.rm = TRUE),
    min_pH = min(pH, na.rm = TRUE),
    max_pH = max(pH, na.rm = TRUE),
    
    missing_SOM = sum(is.na(SOM)),
    mean_SOM = mean(SOM, na.rm = TRUE),
    median_SOM = median(SOM, na.rm = TRUE),
    min_SOM = min(SOM, na.rm = TRUE),
    max_SOM = max(SOM, na.rm = TRUE),
    
    missing_C = sum(is.na(C)),
    mean_C = mean(C, na.rm = TRUE),
    median_C = median(C, na.rm = TRUE),
    min_C = min(C, na.rm = TRUE),
    max_C = max(C, na.rm = TRUE),
    
    missing_N = sum(is.na(N)),
    mean_N = mean(N, na.rm = TRUE),
    median_N = median(N, na.rm = TRUE),
    min_N = min(N, na.rm = TRUE),
    max_N = max(N, na.rm = TRUE),
    
    missing_CN = sum(is.na(CN)),
    mean_CN = mean(CN, na.rm = TRUE),
    median_CN = median(CN, na.rm = TRUE),
    min_CN = min(CN, na.rm = TRUE),
    max_CN = max(CN, na.rm = TRUE),
    
    missing_P = sum(is.na(P)),
    mean_P = mean(P, na.rm = TRUE),
    median_P = median(P, na.rm = TRUE),
    min_P = min(P, na.rm = TRUE),
    max_P = max(P, na.rm = TRUE),
    
    missing_clay = sum(is.na(clay)),
    mean_clay = mean(clay, na.rm = TRUE),
    median_clay = median(clay, na.rm = TRUE),
    min_clay = min(clay, na.rm = TRUE),
    max_clay = max(clay, na.rm = TRUE),
    
    missing_claysilt = sum(is.na(claysilt)),
    mean_claysilt = mean(claysilt, na.rm = TRUE),
    median_claysilt = median(claysilt, na.rm = TRUE),
    min_claysilt = min(claysilt, na.rm = TRUE),
    max_claysilt = max(claysilt, na.rm = TRUE),
    
    missing_SOMTXT = sum(is.na(SOMTXT)),
    mean_SOMTXT = mean(SOMTXT, na.rm = TRUE),
    median_SOMTXT = median(SOMTXT, na.rm = TRUE),
    min_SOMTXT = min(SOMTXT, na.rm = TRUE),
    max_SOMTXT = max(SOMTXT, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(summary_S11, file.path(tabdir_bgq, "Table_S11_summary_orchards_by_FC.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 4.2 Mean ± SE table
# -----------------------------------------------------------------------------

summary_S11_mean_se <- map_crop %>%
  group_by(FC) %>%
  summarise(
    across(
      all_of(bgq_vars),
      list(mean = ~ mean(.x, na.rm = TRUE),
           SE   = ~ calculate_SE(.x)),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

write.csv(summary_S11_mean_se, file.path(tabdir_bgq, "Table_S11_mean_se_orchards_by_FC.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 4.3 Welch ANOVA + Games-Howell
# -----------------------------------------------------------------------------

welch_S11 <- lapply(bgq_vars, function(v) {
  fit <- oneway.test(reformulate("FC", response = v), data = map_crop, var.equal = FALSE)
  data.frame(
    Variable = v,
    Statistic = unname(fit$statistic),
    df_num = unname(fit$parameter[1]),
    df_den = unname(fit$parameter[2]),
    p_value = fit$p.value
  )
}) %>% bind_rows()

write.csv(welch_S11, file.path(tabdir_bgq, "Table_S11_welch_tests.csv"), row.names = FALSE)

games_howell_S11 <- lapply(bgq_vars, function(v) {
  out <- run_games_howell(map_crop, v, "FC")
  data.frame(
    Variable = v,
    Output = out,
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

write.csv(games_howell_S11, file.path(tabdir_bgq, "Table_S11_games_howell.csv"), row.names = FALSE)

# =============================================================================
# 5. Export Excel workbook
# =============================================================================

wb <- createWorkbook()

addWorksheet(wb, "S1_summary")
writeData(wb, "S1_summary", summary_S1)

addWorksheet(wb, "S1_mean_se")
writeData(wb, "S1_mean_se", summary_S1_mean_se)

addWorksheet(wb, "S1_welch")
writeData(wb, "S1_welch", welch_S1)

addWorksheet(wb, "S1_games_howell")
writeData(wb, "S1_games_howell", games_howell_S1)

addWorksheet(wb, "S11_summary")
writeData(wb, "S11_summary", summary_S11)

addWorksheet(wb, "S11_mean_se")
writeData(wb, "S11_mean_se", summary_S11_mean_se)

addWorksheet(wb, "S11_welch")
writeData(wb, "S11_welch", welch_S11)

addWorksheet(wb, "S11_games_howell")
writeData(wb, "S11_games_howell", games_howell_S11)

saveWorkbook(wb,
             file = file.path(outdir_bgq, "Tables_S1_S11_biogeochemistry.xlsx"),
             overwrite = TRUE)

# =============================================================================
# 7. Save objects
# =============================================================================

saveRDS(BGQ, file.path(rdsdir_bgq, "BGQ_clean.rds"))
saveRDS(mapping_group, file.path(rdsdir_bgq, "mapping_forest_orchard.rds"))
saveRDS(map_crop, file.path(rdsdir_bgq, "mapping_orchard_FC.rds"))

# =============================================================================
# 8. Session info
# =============================================================================

writeLines(capture.output(sessionInfo()),
           con = file.path(outdir_bgq, "sessionInfo_biogeochemical_differences.txt"))
