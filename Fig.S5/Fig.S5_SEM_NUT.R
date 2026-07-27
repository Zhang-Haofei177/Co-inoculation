# ============================================================================
# Fig.S5 SEM: EPS synthesis → Macroaggregates/Soil nutrients → Crop performance (HY vs YC drought gating)
# ============================================================================
# Model: CPI ~ AGG + NUT; AGG ~ EPS; NUT ~ EPS; AGG ~~ NUT
# Latent variable definitions:
#   EPS = mean(GT4, GT2)              — EPS synthesis capacity
#   AGG = X250um                       — Macroaggregate proportion
#   NUT = mean(TN,SOC,AP,TP,NH4,NO3)  — Soil nutrients
#   CPI = mean(Yield, WUE, NUE)       — Crop performance index
# ============================================================================
library(lavaan)
library(dplyr)
library(tidyr)

# ---------- 0. Paths ----------
out_root <- "C:/Users/ADMIN/Desktop/Co-inoculation/Fig.S5"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
setwd(out_root)

# ---------- 1. Load data ----------
soil <- read.csv(file.path(out_root, "data/soil_properties_topsoil.csv"),
                 stringsAsFactors = FALSE)
soil <- soil[soil$Treatment %in% c("AMF", "Str", "AMF+Str"), ]
soil <- soil[, c("Sample", "Site", "Treatment", "TN", "SOC", "AP", "TP", "NH4", "NO3",
                 "BG", "AKP", "Urease", "LAP", "NAG", "X250um")]

pheno <- read.csv(file.path(out_root, "data/maize_phenotype_summary.csv"), stringsAsFactors = FALSE)
pheno <- pheno[pheno$Treatment %in% c("AMF", "Str", "AMF+Str"), ]
pheno <- pheno[, c("Sample", "Yield_t_ha", "WUE..kg.mm.ha.", "NUE.", "PUE_kg_mm")]
colnames(pheno)[3:5] <- c("WUE", "NUE", "PUE")

cazy_raw <- read.table(file.path(out_root, "data/Unigenes.absolute.level2.xls"),
                        header = TRUE, sep = "\t", row.names = 1, quote = "", comment.char = "")
desc_col <- which(tolower(colnames(cazy_raw)) == "description")
if (length(desc_col) > 0) cazy_raw <- cazy_raw[, -desc_col]
cazy_t <- as.data.frame(t(cazy_raw))
cazy_t$Sample <- rownames(cazy_t)

sig14 <- c("GH94", "GH130", "GT4", "GT2", "GT1", "GH32", "GT31",
           "GT14", "GT7", "AA3", "AA5", "GH13", "GT61", "GH18")
cazy_avail <- intersect(sig14, colnames(cazy_t))
cazy <- cazy_t[, c("Sample", cazy_avail)]

# ---------- 2. Data integration ----------
make_key <- function(x) {
  p <- strsplit(x, "\\.")
  sapply(p, function(v) paste0(v[1], ".", v[3], ".", v[4]))
}

df_ts <- merge(soil, pheno, by = "Sample", all = FALSE)
df_ts$PlotKey <- make_key(df_ts$Sample)

cazy$PlotKey <- make_key(cazy$Sample)
cazy_plot <- cazy %>%
  group_by(PlotKey) %>%
  summarise(across(all_of(cazy_avail), mean), .groups = "drop")

df <- merge(df_ts, cazy_plot, by = "PlotKey", all = FALSE)

# ---------- 3. Standardization + construct latent variables ----------
num_cols <- sapply(df, is.numeric)
df_s <- df
df_s[, num_cols] <- as.data.frame(scale(df[, num_cols]))

df_s$EPS <- rowMeans(df_s[, c("GT4", "GT2")])
df_s$CPI <- rowMeans(df_s[, c("Yield_t_ha", "WUE", "NUE")])
df_s$AGG <- df_s$X250um
df_s$NUT <- rowMeans(df_s[,c("TN","SOC","AP","TP","NH4","NO3")])

# ---------- 4. SEM model definition ----------
syn_model <- '
  CPI ~ b1*AGG + b2*NUT
  AGG ~ a1*EPS
  NUT ~ a2*EPS
  AGG ~~ NUT
  indir_agg  := a1*b1
  indir_nut  := a2*b2
  total_indir := indir_agg + indir_nut
'

# ---------- 5. Fit models (HY + YC) ----------
fit_models <- function(data, site_name) {
  d <- data[data$Site == site_name, ]
  cat(sprintf("\n========== %s (N=%d) ==========\n", site_name, nrow(d)))

  fit <- sem(syn_model, data = d, meanstructure = TRUE)

  # Model fit indices
  fit_measures <- fitMeasures(fit)
  cat(sprintf("  CFI = %.3f\n",   fit_measures["cfi"]))
  cat(sprintf("  RMSEA = %.3f\n", fit_measures["rmsea"]))
  cat(sprintf("  SRMR = %.3f\n",  fit_measures["srmr"]))
  cat(sprintf("  χ²(%d) = %.2f, p = %.3f\n",
              fit_measures["df"], fit_measures["chisq"], fit_measures["pvalue"]))

  # R²
  rsq <- lavInspect(fit, "rsquare")
  cat(sprintf("  R²: AGG=%.3f  NUT=%.3f  CPI=%.3f\n",
              rsq["AGG"], rsq["NUT"], rsq["CPI"]))

  return(fit)
}

fit_hy <- fit_models(df_s, "HY")
fit_yc <- fit_models(df_s, "YC")

# ---------- 6. Extract effects ----------
extract_effects <- function(fit, site_name) {
  sm <- standardizedSolution(fit)

  effects <- data.frame()

  # Direct effects on CPI
  paths_cpi <- sm[sm$op == "~" & sm$lhs == "CPI", ]
  for (i in seq_len(nrow(paths_cpi))) {
    effects <- rbind(effects, data.frame(
      predictor = paths_cpi$rhs[i], outcome = "CPI",
      effect_type = "Direct", effect = paths_cpi$est.std[i],
      p = paths_cpi$pvalue[i], stringsAsFactors = FALSE
    ))
  }

  # Indirect via AGG
  agg_row <- sm[sm$lhs == "indir_agg", ]
  if (nrow(agg_row) > 0) {
    effects <- rbind(effects, data.frame(
      predictor = "EPS", outcome = "CPI",
      effect_type = "Indirect\n(via AGG)", effect = agg_row$est.std,
      p = agg_row$pvalue, stringsAsFactors = FALSE
    ))
  }

  # Indirect via NUT
  nut_row <- sm[sm$lhs == "indir_nut", ]
  if (nrow(nut_row) > 0) {
    effects <- rbind(effects, data.frame(
      predictor = "EPS", outcome = "CPI",
      effect_type = "Indirect\n(via NUT)", effect = nut_row$est.std,
      p = nut_row$pvalue, stringsAsFactors = FALSE
    ))
  }

  # Total indirect
  tot_row <- sm[sm$lhs == "total_indir", ]
  if (nrow(tot_row) > 0) {
    effects <- rbind(effects, data.frame(
      predictor = "EPS", outcome = "CPI",
      effect_type = "Total\nIndirect", effect = tot_row$est.std,
      p = tot_row$pvalue, stringsAsFactors = FALSE
    ))
  }

  # EPS → AGG, EPS → NUT
  agg_path <- sm[sm$lhs == "AGG" & sm$op == "~", ]
  nut_path <- sm[sm$lhs == "NUT" & sm$op == "~", ]
  effects <- rbind(effects, data.frame(
    predictor = "EPS", outcome = "AGG",
    effect_type = "Direct", effect = agg_path$est.std[1],
    p = agg_path$pvalue[1], stringsAsFactors = FALSE
  ))
  effects <- rbind(effects, data.frame(
    predictor = "EPS", outcome = "NUT",
    effect_type = "Direct", effect = nut_path$est.std[1],
    p = nut_path$pvalue[1], stringsAsFactors = FALSE
  ))

  # Significance stars
  effects$stars <- ifelse(effects$p < 0.001, "***",
                   ifelse(effects$p < 0.01, "**",
                   ifelse(effects$p < 0.05, "*", "")))
  effects$stars[effects$stars == ""] <- "ns"

  # Factor ordering
  effects$outcome <- factor(effects$outcome, levels = c("AGG", "NUT", "CPI"))
  effects$effect_type <- factor(effects$effect_type,
    levels = c("Direct", "Indirect\n(via AGG)", "Indirect\n(via NUT)", "Total\nIndirect"))

  effects$Site <- site_name
  return(effects)
}

effects_hy <- extract_effects(fit_hy, "HY")
effects_yc <- extract_effects(fit_yc, "YC")
effects_all <- rbind(effects_hy, effects_yc)

# ---------- 7. Save effects table ----------
write.csv(effects_all, file.path(out_root, "Fig.S5_SEM_effects_HY_YC.csv"), row.names = FALSE)

# ---------- 8. Model fit summary table ----------
fit_summary <- data.frame(
  Site = c("HY", "YC"),
  N = c(12, 12),
  CFI = c(fitMeasures(fit_hy)["cfi"], fitMeasures(fit_yc)["cfi"]),
  RMSEA = c(fitMeasures(fit_hy)["rmsea"], fitMeasures(fit_yc)["rmsea"]),
  SRMR = c(fitMeasures(fit_hy)["srmr"], fitMeasures(fit_yc)["srmr"]),
  Chi2 = c(fitMeasures(fit_hy)["chisq"], fitMeasures(fit_yc)["chisq"]),
  df = c(fitMeasures(fit_hy)["df"], fitMeasures(fit_yc)["df"]),
  Chi2_p = c(fitMeasures(fit_hy)["pvalue"], fitMeasures(fit_yc)["pvalue"]),
  R2_AGG = c(lavInspect(fit_hy, "rsquare")["AGG"], lavInspect(fit_yc, "rsquare")["AGG"]),
  R2_NUT = c(lavInspect(fit_hy, "rsquare")["NUT"], lavInspect(fit_yc, "rsquare")["NUT"]),
  R2_CPI = c(lavInspect(fit_hy, "rsquare")["CPI"], lavInspect(fit_yc, "rsquare")["CPI"]),
  stringsAsFactors = FALSE
)
write.csv(fit_summary, file.path(out_root, "Fig.S5_SEM_model_fit_summary.csv"), row.names = FALSE)

cat("\n========== SEM Model Fit Summary ==========\n")
print(fit_summary)

cat("\n========== SEM Effects ==========\n")
print(effects_all[, c("Site", "predictor", "outcome", "effect_type", "effect", "stars")])

cat("\n========== Fig.S5 SEM (NUT) Analysis Complete ==========\n")
cat("Output files saved to:", out_root, "\n")
cat("  Fig.S5_SEM_effects_HY_YC.csv           — Effects table\n")
cat("  Fig.S5_SEM_model_fit_summary.csv       — Model fit indices summary (CFI/RMSEA/SRMR/χ²/R²)\n")
