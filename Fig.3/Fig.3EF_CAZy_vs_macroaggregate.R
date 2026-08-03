# ============================================================================
# Fig.3E-F: CAZy family (GT2, GT4) vs macroaggregate scatter plot
# ============================================================================
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ---------- 0. Paths ----------
out_root <- "C:/Users/ADMIN/Desktop/Co-inoculation/Fig.3"
script_dir <- out_root
dir.create(script_dir, showWarnings = FALSE, recursive = TRUE)
setwd(script_dir)

# ---------- 1. Load data ----------
family_abund <- readRDS(file.path(out_root, "data/family_abund.rds"))

soil <- read.csv(file.path(out_root, "data/soil_properties_topsoil.csv"), row.names = 1)

group <- read.csv(file.path(out_root, "data/group.csv"), row.names = 1)
group$Treatment <- factor(group$Treatment, levels = c("Ctrl", "AMF", "Str", "AMF+Str"))
group$Site <- factor(group$Site, levels = c("HY", "YC"))
group$Sample <- rownames(group)
group$Rep <- gsub(".*\\.", "", rownames(group))

# ---------- 2. Match rhizosphere metagenome with corresponding plot Topsoil ----------
rhizo_meta <- group[group$Niche == "Rhizosphere", ]
rhizo_meta <- rhizo_meta[rhizo_meta$Sample %in% rownames(family_abund), ]

soil_meta <- group[group$Niche == "Topsoil", ]
rhizo_meta$MatchKey <- paste(rhizo_meta$Site, rhizo_meta$Treatment, rhizo_meta$Rep, sep = "_")
soil_meta$MatchKey <- paste(soil_meta$Site, soil_meta$Treatment, soil_meta$Rep, sep = "_")

soil_data <- soil[rownames(soil_meta), ]
soil_match <- data.frame(MatchKey = soil_meta$MatchKey, soil_data, stringsAsFactors = FALSE)

matched <- inner_join(rhizo_meta %>% select(Sample, MatchKey), soil_match, by = "MatchKey")
rownames(matched) <- matched$Sample

y_vec <- as.numeric(as.character(matched$X250um))
names(y_vec) <- matched$Sample

fam_mat <- family_abund[matched$Sample, ]

# ---------- Keep only GT2 and GT4 ----------
target_families <- c("GT2", "GT4")
available_targets <- intersect(target_families, colnames(fam_mat))
if (length(available_targets) == 0) stop("Neither GT2 nor GT4 found in family_abund!")
fam_mat <- fam_mat[, available_targets, drop = FALSE]
cat("CAZy families analyzed:", paste(available_targets, collapse = ", "), "\n")

common_samples <- intersect(names(y_vec), rownames(fam_mat))
y_vec <- y_vec[common_samples]
fam_mat <- fam_mat[common_samples, ]

plot_df <- cbind(data.frame(Sample = common_samples, X250um = y_vec), fam_mat)
plot_df$Site <- rhizo_meta[common_samples, "Site"]

cazy_families <- colnames(fam_mat)

# ---------- 3. Compute statistics and plot ----------
stats_df <- data.frame(Family = character(), Site = character(), R2 = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

plot_single <- function(family_name) {
  sub <- plot_df[, c("X250um", family_name, "Site")]
  colnames(sub) <- c("X250um", "Abundance", "Site")
  sub <- na.omit(sub)
  if (nrow(sub) < 3) return(NULL)

  fit_stats <- data.frame(Family = family_name, Site = c("HY", "YC"), R2 = NA, p_value = NA, stringsAsFactors = FALSE)
  for (s in c("HY", "YC")) {
    sdata <- sub[sub$Site == s, ]
    if (nrow(sdata) >= 3) {
      model <- lm(X250um ~ Abundance, data = sdata)
      ssum <- summary(model)
      fstat <- ssum$fstatistic
      p_val <- pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
      fit_stats[fit_stats$Site == s, c("R2", "p_value")] <- c(ssum$r.squared, p_val)
    }
  }
  stats_df <<- rbind(stats_df, fit_stats)

  hy_text <- ifelse(!is.na(fit_stats$R2[1]),
                    sprintf("HY: R²=%.2f, p=%.3f", fit_stats$R2[1], fit_stats$p_value[1]),
                    "HY: insufficient data")
  yc_text <- ifelse(!is.na(fit_stats$R2[2]),
                    sprintf("YC: R²=%.2f, p=%.3f", fit_stats$R2[2], fit_stats$p_value[2]),
                    "YC: insufficient data")

  ggplot(sub, aes(x = Abundance, y = X250um, color = Site)) +
    geom_point(size = 2.5, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
    scale_color_manual(values = c(HY = "#1F78B4", YC = "#E31A23")) +
    labs(x = paste(family_name, "abundance"), y = "Macroaggregate (%)", title = family_name) +
    annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.5, label = hy_text, size = 3, color = "#1F78B4") +
    annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 3.0, label = yc_text, size = 3, color = "#E31A23") +
    theme_bw(base_size = 12) +
    theme(legend.position = "top")
}

# Generate single plots (populating stats_df)
plots <- list()
for (fam in cazy_families) {
  p <- plot_single(fam)
  if (!is.null(p)) {
    plots[[length(plots) + 1]] <- p
  }
}

# Combined overview (Fig.3E-F)
if (length(plots) > 0) {
  combined <- wrap_plots(plots, ncol = 2, guides = "collect") +
    plot_annotation(title = "CAZy family abundance vs. macroaggregate proportion (HY vs YC)")
  ggsave(file.path(script_dir, "Fig.3EF_CAZy_vs_macroaggregate_scatter.pdf"), combined, width = 10, height = 5)
  ggsave(file.path(script_dir, "Fig.3EF_CAZy_vs_macroaggregate_scatter.png"), combined, width = 10, height = 5, dpi = 300)
  cat("Fig.3EF scatter saved\n")
}

# Save statistics table
stats_df$R2 <- as.numeric(stats_df$R2)
stats_df$p_value <- as.numeric(stats_df$p_value)
write.csv(stats_df, file.path(script_dir, "Fig.3EF_CAZy_linear_fit_stats.csv"), row.names = FALSE)

cat("Analysis complete! Results saved to", script_dir, "\n")
