# ============================================================================
# CAZy-microbe attribution analysis -- Fig.3B specific
# ============================================================================
# Retained:
#   p1 (Step 5): Log2FC dotplot → Fig.3B_Log2FC_dotplot
#   Fig2 (Step 6): Log2FC heatmap → Fig.3B_Log2FC_heatmap
# Removed: correlation heatmap, regression scatter, species boxplot, phylum-level heatmap, overview plot, data table export
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggpubr)
  library(ggpmisc)
  library(ggsignif)
  library(ggrepel)
  library(patchwork)
  library(pheatmap)
  library(reshape2)
  library(scales)
  library(rstatix)
  library(multcompView)
  library(showtext)
  library(scico)
  library(MetBrewer)
  library(corrplot)
})

# ============================================================================
# 0. Global settings
# ============================================================================
in_dir  <- "C:/Users/ADMIN/Desktop/Co-inoculation/Fig.3/data"
out_dir <- "C:/Users/ADMIN/Desktop/Co-inoculation/Fig.3"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Font settings (sans-serif, similar to Arial/Helvetica)
showtext_auto()
if (file.exists("C:/Windows/Fonts/arial.ttf")) {
  font_add("Arial", regular = "C:/Windows/Fonts/arial.ttf",
           bold = "C:/Windows/Fonts/arialbd.ttf",
           italic = "C:/Windows/Fonts/ariali.ttf")
  base_font <- "Arial"
} else {
  base_font <- "sans"
}
if (!file.exists("C:/Windows/Fonts/arial.ttf") && file.exists("C:/Windows/Fonts/arial.ttf")) {
  base_font <- "sans"
}
cat(sprintf("Using font: %s\n", base_font))

# Color scheme
trt_cols   <- c("AMF"="#FCB255", "Str"="#42B382", "AMF+Str"="#107BA4")
site_cols  <- c("HY"="#E41A1C", "YC"="#377EB8")
cat5_cols  <- c("EPS synthesis & osmoprotectant"="#74C476",
                "Carbon source utilization"="#6BAED6",
                "ROS scavenging"="#FB6A4A",
                "Phytohormone homeostasis"="#9E9AC8",
                "Cell wall remodeling & signaling"="#FDBF6F")

# 14 significant CAZy families
sig_families <- c("GH94","GH130","GT4","GT2","GT1","GH32","GT31",
                  "GT14","GT7","AA3","AA5","GH13","GT61","GH18")

family_category <- data.frame(
  Family = sig_families,
  Category = c("Carbon source utilization","Carbon source utilization",
    "EPS synthesis & osmoprotectant","EPS synthesis & osmoprotectant",
    "Phytohormone homeostasis","Carbon source utilization",
    "Cell wall remodeling & signaling","Cell wall remodeling & signaling",
    "Cell wall remodeling & signaling","ROS scavenging","ROS scavenging",
    "Carbon source utilization","Cell wall remodeling & signaling",
    "Carbon source utilization"),
  Direction = c("Up_HY","Up_HY","Up_YC","Up_YC","Up_YC","Down_YC",
    "Up_YC","Down_YC","Down_YC","Up_YC","Up_YC","Down_YC","Up_YC","Up_YC"),
  stringsAsFactors = FALSE
)

# ggplot theme: Nature/Science style clean theme
theme_journal <- theme_bw(base_size = 9) +
  theme(
    text = element_text(family = base_font, color = "black"),
    axis.text = element_text(family = base_font, size = 8, color = "black"),
    axis.title = element_text(family = base_font, size = 9, face = "bold"),
    plot.title = element_text(family = base_font, size = 11, face = "bold"),
    plot.subtitle = element_text(family = base_font, size = 9),
    legend.text = element_text(family = base_font, size = 7),
    legend.title = element_text(family = base_font, size = 8, face = "bold"),
    strip.text = element_text(family = base_font, size = 8, face = "bold"),
    strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.3),
    panel.grid.major = element_line(linewidth = 0.15, color = "grey90"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 0.5, color = "black"),
    legend.key.size = unit(0.4, "cm"),
    legend.margin = margin(2, 2, 2, 2),
    plot.margin = margin(8, 8, 8, 8)
  )

# ============================================================================
# Helper: read abundance matrix
# ============================================================================
read_abund_matrix <- function(filepath) {
  df <- fread(filepath, sep = "\t", header = TRUE, data.table = FALSE)
  rn <- df[, 1]; df <- df[, -1, drop = FALSE]
  df[] <- lapply(df, function(x) as.numeric(as.character(x)))
  rownames(df) <- rn; return(df)
}

# ============================================================================
# Helper: Log2FC with exhaustive pairwise ratios
# ============================================================================
compute_log2fc_exhaustive <- function(dual_vals, single_vals) {
  n_dual <- length(dual_vals); n_single <- length(single_vals)
  ratios <- numeric(n_dual * n_single); k <- 1
  for (i in seq_len(n_dual)) {
    for (j in seq_len(n_single)) {
      r <- (dual_vals[i] + 0.5) / (single_vals[j] + 0.5)
      ratios[k] <- log2(r); k <- k + 1
    }
  }
  list(mean_log2fc = mean(ratios), sd_log2fc = sd(ratios),
       se_log2fc = sd(ratios)/sqrt(length(ratios)),
       median_log2fc = median(ratios),
       ci95_low = quantile(ratios, 0.025), ci95_high = quantile(ratios, 0.975),
       ratios = ratios, n_pairs = length(ratios))
}

# ============================================================================
# 1. Read basic data
# ============================================================================
cat("\n========== Step 1: Reading metadata ==========\n")
group <- read.csv(file.path(in_dir, "group.csv"), row.names = 1)
group_rhizo <- group[group$Niche == "Rhizosphere" &
                       group$Treatment %in% c("AMF","Str","AMF+Str"), ]
group_rhizo$Treatment <- factor(group_rhizo$Treatment, levels = c("AMF","Str","AMF+Str"))
group_rhizo$Site <- factor(group_rhizo$Site, levels = c("HY","YC"))
rhizo_samples <- rownames(group_rhizo)
cat(sprintf("Rhizosphere samples: %d\n", nrow(group_rhizo)))

# ============================================================================
# 2. Read CAZy family abundance + Log2FC calculation
# ============================================================================
cat("\n========== Step 2: CAZy family abundance & Log2FC ==========\n")

cazy_level2 <- read_abund_matrix(file.path(in_dir, "Unigenes.absolute.level2.xls"))
cazy_abund <- as.data.frame(t(cazy_level2))
cazy_abund$Sample <- rownames(cazy_abund)

available_fams <- intersect(sig_families, colnames(cazy_abund))
cat(sprintf("Families in level2: %d / 14\n", length(available_fams)))

# Fill missing families from gene counts
tmp_filt <- file.path(in_dir, "cazy_tax_sig_families.txt")
if (file.exists(tmp_filt)) {
  cazy_tax_full <- fread(tmp_filt, sep = "\t", header = FALSE, quote = "",
                         col.names = c("Description","CAZy_Family","Gene_ID","Taxonomy"),
                         stringsAsFactors = FALSE, showProgress = FALSE)
  cazy_tax_full[, Sample := sub("_[0-9]+$", "", Gene_ID)]
  cazy_tax_full <- cazy_tax_full[Sample %in% rhizo_samples]
} else {
  cazy_tax_full <- data.table()
}

missing_fams <- setdiff(sig_families, colnames(cazy_abund))
sample_cols <- setdiff(colnames(cazy_abund), "Sample")
for (mf in missing_fams) {
  gene_counts <- cazy_tax_full[CAZy_Family == mf, .(Count = .N), by = Sample]
  gene_counts <- merge(data.frame(Sample = sample_cols, stringsAsFactors = FALSE),
                       gene_counts, by = "Sample", all.x = TRUE)
  gene_counts$Count[is.na(gene_counts$Count)] <- 0
  cazy_abund[[mf]] <- gene_counts$Count[match(cazy_abund$Sample, gene_counts$Sample)]
}

# ---- Calculate Log2FC (exhaustive pairing) ----
log2fc_results <- data.frame()

for (site in c("HY","YC")) {
  site_samples <- group_rhizo[group_rhizo$Site == site, ]
  dual_samples <- rownames(site_samples[site_samples$Treatment == "AMF+Str", ])
  amf_samples  <- rownames(site_samples[site_samples$Treatment == "AMF", ])
  str_samples  <- rownames(site_samples[site_samples$Treatment == "Str", ])

  if (length(dual_samples) < 3 || length(amf_samples) < 3 || length(str_samples) < 3) next

  for (fam in sig_families) {
    if (!fam %in% colnames(cazy_abund)) next

    dual_vals <- cazy_abund[dual_samples, fam]
    amf_vals  <- cazy_abund[amf_samples, fam]
    str_vals  <- cazy_abund[str_samples, fam]

    # Log2FC vs AMF
    res_amf <- compute_log2fc_exhaustive(dual_vals, amf_vals)
    log2fc_results <- rbind(log2fc_results, data.frame(
      Site = site, CAZy_Family = fam, Reference = "AMF",
      mean_log2fc = res_amf$mean_log2fc, se_log2fc = res_amf$se_log2fc,
      median_log2fc = res_amf$median_log2fc,
      ci95_low = res_amf$ci95_low, ci95_high = res_amf$ci95_high,
      n_pairs = res_amf$n_pairs, stringsAsFactors = FALSE
    ))

    # Log2FC vs Str
    res_str <- compute_log2fc_exhaustive(dual_vals, str_vals)
    log2fc_results <- rbind(log2fc_results, data.frame(
      Site = site, CAZy_Family = fam, Reference = "Str",
      mean_log2fc = res_str$mean_log2fc, se_log2fc = res_str$se_log2fc,
      median_log2fc = res_str$median_log2fc,
      ci95_low = res_str$ci95_low, ci95_high = res_str$ci95_high,
      n_pairs = res_str$n_pairs, stringsAsFactors = FALSE
    ))
  }
}

# Add function category
log2fc_results <- merge(log2fc_results, family_category, by.x = "CAZy_Family", by.y = "Family")

# Combine Log2FC of vs_AMF and vs_Str (mean synergy effect)
log2fc_combined <- log2fc_results %>%
  group_by(Site, CAZy_Family, Category, Direction) %>%
  summarise(
    mean_log2fc = mean(mean_log2fc),
    se_log2fc   = sqrt(mean(se_log2fc^2)),
    ci95_low    = mean(ci95_low),
    ci95_high   = mean(ci95_high),
    .groups = 'drop'
  )

# ---- Significance test ----
cat("Computing exhaustive pairwise Log2FC + significance...\n")

log2fc_stats <- data.frame()

for (site in c("HY","YC")) {
  site_samples <- group_rhizo[group_rhizo$Site == site, ]
  dual_samples <- rownames(site_samples[site_samples$Treatment == "AMF+Str", ])
  amf_samples  <- rownames(site_samples[site_samples$Treatment == "AMF", ])
  str_samples  <- rownames(site_samples[site_samples$Treatment == "Str", ])

  for (fam in sig_families) {
    if (!fam %in% colnames(cazy_abund)) next

    dual_vals <- cazy_abund[dual_samples, fam]
    amf_vals  <- cazy_abund[amf_samples, fam]
    str_vals  <- cazy_abund[str_samples, fam]

    all_ratios <- numeric(32); k <- 1
    for (i in 1:4) {
      for (j in 1:4) {
        all_ratios[k] <- log2((dual_vals[i] + 0.5) / (amf_vals[j] + 0.5)); k <- k + 1
        all_ratios[k] <- log2((dual_vals[i] + 0.5) / (str_vals[j] + 0.5)); k <- k + 1
      }
    }

    wt <- wilcox.test(all_ratios, mu = 0, exact = FALSE)
    tt <- t.test(all_ratios, mu = 0)

    log2fc_stats <- rbind(log2fc_stats, data.frame(
      Site = site, CAZy_Family = fam,
      mean_log2fc = mean(all_ratios), se_log2fc = sd(all_ratios)/sqrt(32),
      wilcoxon_P = wt$p.value, ttest_P = tt$p.value,
      n_pairs = 32, stringsAsFactors = FALSE
    ))
  }
}

# FDR per site
log2fc_stats <- log2fc_stats %>%
  group_by(Site) %>%
  mutate(wilcoxon_P_adj = p.adjust(wilcoxon_P, method = "BH"),
         ttest_P_adj = p.adjust(ttest_P, method = "BH")) %>%
  ungroup()

# Between-site comparison
between_site_stats <- data.frame()
for (fam in sig_families) {
  if (!fam %in% colnames(cazy_abund)) next

  ratios_hy <- c(); ratios_yc <- c()
  for (site in c("HY","YC")) {
    site_samples <- group_rhizo[group_rhizo$Site == site, ]
    dual_samples <- rownames(site_samples[site_samples$Treatment == "AMF+Str", ])
    amf_samples  <- rownames(site_samples[site_samples$Treatment == "AMF", ])
    str_samples  <- rownames(site_samples[site_samples$Treatment == "Str", ])
    dual_vals <- cazy_abund[dual_samples, fam]
    amf_vals  <- cazy_abund[amf_samples, fam]
    str_vals  <- cazy_abund[str_samples, fam]

    r <- numeric(32); k <- 1
    for (i in 1:4) for (j in 1:4) {
      r[k] <- log2((dual_vals[i]+0.5)/(amf_vals[j]+0.5)); k <- k+1
      r[k] <- log2((dual_vals[i]+0.5)/(str_vals[j]+0.5)); k <- k+1
    }
    if (site == "HY") ratios_hy <- r else ratios_yc <- r
  }

  mw <- wilcox.test(ratios_hy, ratios_yc, exact = FALSE)
  between_site_stats <- rbind(between_site_stats, data.frame(
    CAZy_Family = fam,
    mean_log2fc_HY = mean(ratios_hy), mean_log2fc_YC = mean(ratios_yc),
    delta_log2fc = mean(ratios_yc) - mean(ratios_hy),
    MW_P = mw$p.value, stringsAsFactors = FALSE
  ))
}
between_site_stats$MW_P_adj <- p.adjust(between_site_stats$MW_P, method = "BH")

cat(sprintf("Within-site significant (Wilcoxon P_adj<0.05): HY=%d, YC=%d\n",
  sum(log2fc_stats$wilcoxon_P_adj[log2fc_stats$Site=="HY"] < 0.05),
  sum(log2fc_stats$wilcoxon_P_adj[log2fc_stats$Site=="YC"] < 0.05)))
cat(sprintf("Between-site significant (MW P_adj<0.05): %d / 14\n",
  sum(between_site_stats$MW_P_adj < 0.05)))

# ============================================================================
# 3. Fig.3B dotplot: Log2FC lollipop chart -- Core result
# ============================================================================
cat("\n========== Step 3: Fig.3B - Log2FC dot plot ==========\n")

plot_data <- merge(log2fc_combined, between_site_stats[, c("CAZy_Family","MW_P_adj")],
                   by = "CAZy_Family")
plot_data$Stars <- ifelse(plot_data$MW_P_adj < 0.001, "***",
                          ifelse(plot_data$MW_P_adj < 0.01, "**",
                                 ifelse(plot_data$MW_P_adj < 0.05, "*", "")))
plot_data$Stars[is.na(plot_data$Stars)] <- ""

fam_order_yc <- plot_data %>%
  filter(Site == "YC") %>%
  arrange(mean_log2fc) %>%
  pull(CAZy_Family)
plot_data$CAZy_Family <- factor(plot_data$CAZy_Family, levels = fam_order_yc)
plot_data$Category <- factor(plot_data$Category,
                             levels = c("EPS synthesis & osmoprotectant","Carbon source utilization",
                                        "ROS scavenging","Phytohormone homeostasis","Cell wall remodeling & signaling"))
plot_data$CAZy_Family <- factor(plot_data$CAZy_Family,levels = rev(c("GT7","GT14","GT61","GT31","GT1","AA5","AA3","GH32","GH13","GH94","GH18","GH130","GT2","GT4")))

plot_data <- plot_data %>%
  mutate(
    y_num = as.numeric(CAZy_Family),
    y_num = ifelse(Site == "HY", y_num - 0.2, y_num + 0.2)
  )

bg_data <- data.frame(
  y_pos = seq_along(levels(plot_data$CAZy_Family)),
  fill = ifelse(seq_along(levels(plot_data$CAZy_Family)) %% 2 == 1, "grey95", "white")
)

# ---------- Plot ----------
p1 <- ggplot(plot_data, aes(x = mean_log2fc, y = y_num, color = Site)) +
  geom_rect(data = bg_data,
            aes(xmin = -Inf, xmax = Inf, ymin = y_pos - 0.5, ymax = y_pos + 0.5,
                fill = fill),
            inherit.aes = FALSE, alpha = 0.8) +
  scale_fill_identity() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_errorbar(aes(xmin = ci95_low, xmax = ci95_high),
                height = 0.25, linewidth = 0.2, alpha = 1,
                position = "identity") +
  geom_point(aes(shape = Site), size = 1, position = "identity") +
  scale_color_manual(values = site_cols, name = "Site") +
  scale_shape_manual(values = c(19, 17), name = "Site") +
  scale_y_continuous(
    breaks = seq_along(levels(plot_data$CAZy_Family)),
    labels = levels(plot_data$CAZy_Family),
    limits = c(0.5, length(levels(plot_data$CAZy_Family)) + 0.5),
    expand = c(0.02, 0.02)
  ) +
  labs(
    title = expression("Synergy effect: Log"[2]*"(AMF+Str / Single inoculant)"),
    subtitle = "14 CAZy families significantly enriched by GSEA; error bars = 95% CI (32 pairs)",
    x = expression("Log"[2]*" fold change (dual vs single)"),
    y = NULL
  ) +
  theme_bw(base_size = 6) +
  theme(
    axis.text = element_text(size = 6, color = "black"),
    axis.title = element_text(size = 6, color = "black"),
    axis.text.y = element_text(size = 6, face = "plain", color = "black"),
    legend.text = element_text(size = 6, color = "black"),
    legend.title = element_text(size = 6, color = "black"),
    plot.title = element_text(size = 6, hjust = 0.5, color = "black"),
    plot.subtitle = element_text(size = 6, hjust = 0.5, color = "black"),
    strip.text = element_text(size = 6, color = "black"),
    text = element_text(family = "sans", color = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.2),
    axis.line = element_line(linewidth = 0.2),
    legend.key.size = unit(0.4, "lines"),
    plot.margin = margin(5,5,5,5, unit = "pt")
  )

# ---------- Significance stars ----------
sig_fams <- unique(plot_data$CAZy_Family[plot_data$Stars != "" & !is.na(plot_data$Stars)])
if (length(sig_fams) > 0) {
  sig_ann <- plot_data %>%
    filter(CAZy_Family %in% sig_fams) %>%
    group_by(CAZy_Family) %>%
    summarise(
      y_pos = mean(y_num),
      x_pos = max(ci95_high) + 0.15,
      Stars = first(Stars),
      .groups = 'drop'
    )
  p1 <- p1 + geom_text(data = sig_ann, aes(x = x_pos, y = y_pos, label = Stars),
                       inherit.aes = FALSE, size = 3, color = "black", hjust = 0)
}

print(p1)
ggsave(file.path(out_dir, "Fig.3B_Log2FC_dotplot.pdf"), p1,
       width = 60, height = 60, units = "mm", device = "pdf")
ggsave(file.path(out_dir, "Fig.3B_Log2FC_dotplot.png"), p1,
       width = 60, height = 60, units = "mm", dpi = 300, device = "png")
cat("  Fig.3B dotplot saved\n")

# ============================================================================
# 4. Fig.3B heatmap: Log2FC heatmap (with significance stars)
# ============================================================================
cat("  Fig.3B heatmap: Log2FC heatmap with significance stars...\n")

heat_mat <- log2fc_combined %>%
  select(CAZy_Family, Site, mean_log2fc) %>%
  pivot_wider(names_from = Site, values_from = mean_log2fc) %>%
  column_to_rownames("CAZy_Family") %>%
  as.matrix()

star_mat <- log2fc_stats %>%
  mutate(Stars = ifelse(wilcoxon_P_adj < 0.001, "***",
                        ifelse(wilcoxon_P_adj < 0.01, "**",
                               ifelse(wilcoxon_P_adj < 0.05, "*", "")))) %>%
  select(CAZy_Family, Site, Stars) %>%
  pivot_wider(names_from = Site, values_from = Stars) %>%
  column_to_rownames("CAZy_Family") %>%
  as.matrix()

row_ann <- family_category %>%
  column_to_rownames("Family") %>%
  select(Category)
row_ann$Category <- factor(row_ann$Category,
                           levels = c("EPS synthesis & osmoprotectant","Carbon source utilization",
                                      "ROS scavenging","Phytohormone homeostasis","Cell wall remodeling & signaling"))

ann_colors <- list(Category = cat5_cols)

fam_order <- c("GT7","GT14","GT61","GT31","GT1","AA5","AA3","GH32","GH13","GH94","GH18","GH130","GT2","GT4")
heat_mat <- heat_mat[fam_order, , drop = FALSE]
row_ann <- row_ann[fam_order, , drop = FALSE]
star_mat <- star_mat[fam_order, , drop = FALSE]

if (nrow(heat_mat) > 1) {
  pheatmap(heat_mat,
           color = colorRampPalette(c("#2166AC","#F5F5F5","#B2182B"))(100),
           cluster_rows = FALSE, cluster_cols = FALSE,
           annotation_row = row_ann, annotation_colors = ann_colors,
           display_numbers = star_mat,
           number_format = "%s",
           fontsize_number = 6, number_color = "red4",
           main = "Log2FC (Dual/Single) of CAZy families",
           fontsize = 6, fontfamily = "sans",
           cellwidth = 10, cellheight = 10,
           filename = file.path(out_dir, "Fig.3B_Log2FC_heatmap.pdf"))

  pheatmap(heat_mat,
           color = colorRampPalette(c("#2166AC","#F5F5F5","#B2182B"))(100),
           cluster_rows = FALSE, cluster_cols = FALSE,
           annotation_row = row_ann, annotation_colors = ann_colors,
           display_numbers = star_mat,
           number_format = "%s",
           fontsize_number = 6, number_color = "red4",
           main = "Log2FC (Dual/Single) of CAZy families",
           fontsize = 6, fontfamily = "sans",
           cellwidth = 10, cellheight = 10,
           filename = file.path(out_dir, "Fig.3B_Log2FC_heatmap.png"))

  cat("  Fig.3B heatmap saved\n")
}

# ============================================================================
# 5. Analysis summary
# ============================================================================
cat("\n" ); cat(paste(rep("=",70), collapse="")); cat("\n")
cat("ANALYSIS SUMMARY\n")
cat(paste(rep("=",70), collapse="")); cat("\n\n")

cat("--- Synergy significance (within-site Wilcoxon P_adj<0.05) ---\n")
for (site in c("HY","YC")) {
  sig <- log2fc_stats[log2fc_stats$Site == site & log2fc_stats$wilcoxon_P_adj < 0.05, ]
  cat(sprintf("  %s: %d / 14 families significantly different from zero\n",
              site, nrow(sig)))
  if (nrow(sig) > 0)
    cat(sprintf("    %s\n", paste(sig$CAZy_Family, collapse=", ")))
}

cat("\n--- Synergy shift HY vs YC (Mann-Whitney P_adj<0.05) ---\n")
sig_bs <- between_site_stats[between_site_stats$MW_P_adj < 0.05, ]
cat(sprintf("  %d / 14 families show significant site difference\n", nrow(sig_bs)))
if (nrow(sig_bs) > 0) {
  for (i in 1:nrow(sig_bs))
    cat(sprintf("  %s: delta_Log2FC=%.3f, P_adj=%.4f\n",
                sig_bs$CAZy_Family[i], sig_bs$delta_log2fc[i], sig_bs$MW_P_adj[i]))
}

cat(paste(rep("=",70), collapse="")); cat("\n")
cat("Analysis files saved to: ", out_dir, "\n")
