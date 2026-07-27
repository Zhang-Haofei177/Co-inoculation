# ============================================================================
# Fig.4A-C: YAS strategy analysis (PNAS method)
#   Fig.4A: YAS pathway log2(YC/HY) bar chart
#   Fig.4B: YAS strategy independent normalization boxplot
#   Fig.4C: Soil enzyme activity synergy bar chart
# ============================================================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(patchwork)

out_root <- "C:/Users/ADMIN/Desktop/Co-inoculation/Fig.4"
script_dir <- out_root
dir.create(script_dir, showWarnings = FALSE, recursive = TRUE)
setwd(script_dir)

# ============================================================================
# PART 1: YAS strategy independent normalization + boxplot (Fig.4B)
# ============================================================================
cat("\n========== Fig.4B: YAS independent scores boxplot ==========\n")

dat <- read.csv(file.path(out_root, "data/YAS_pathway_detail.csv"), stringsAsFactors = FALSE)
dat <- dat %>% filter(Strategy %in% c("Y", "A", "S"))
dat$Site <- factor(dat$Site, levels = c("HY", "YC"))

# Calculate total absolute abundance per strategy per sample
strategy_abund <- dat %>%
  group_by(Sample, Site, Treatment, Strategy) %>%
  summarise(Total = sum(Abundance), .groups = "drop")

# Apply independent min-max normalization (0-1) per strategy
strategy_norm <- strategy_abund %>%
  group_by(Strategy) %>%
  mutate(Score = (Total - min(Total)) / (max(Total) - min(Total))) %>%
  ungroup()

write.csv(strategy_norm, file.path(script_dir, "Fig.4B_YAS_independent_scores.csv"), row.names = FALSE)

# Visualization: boxplot
site_cols <- c("HY" = "#1F78B4", "YC" = "#E31A23")

plot_data <- strategy_norm %>%
  select(Sample, Site, Treatment, Strategy, Score) %>%
  mutate(Strategy = factor(Strategy, levels = c("Y", "A", "S")))

p_box <- ggplot(plot_data, aes(x = Strategy, y = Score, fill = Site)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6, linewidth = 0.2) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1),
              size = 0.8, alpha = 0.6) +
  scale_fill_manual(values = site_cols) +
  labs(x = "Strategy", y = "Normalized score (0-1)",
       title = "Microbial Y-A-S strategies: independent normalization",
       subtitle = "Each strategy independently min-max normalized across all samples") +
  theme_bw(base_size = 6) +
  theme(
    axis.line = element_blank(),
    legend.position = "top",
    axis.text = element_text(size = 6, color = "black"),
    axis.title = element_text(size = 6, color = "black"),
    plot.title = element_text(size = 6, hjust = 0.5),
    plot.subtitle = element_text(size = 6, hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    text = element_text(family = "sans", color = "black")
  )

strategies <- levels(plot_data$Strategy)
pvals <- sapply(strategies, function(s) {
  test <- wilcox.test(Score ~ Site, data = plot_data[plot_data$Strategy == s, ])
  test$p.value
})
stars <- sapply(pvals, function(p) {
  if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
})

star_data <- data.frame(
  Strategy = factor(strategies, levels = strategies),
  star = stars,
  y = sapply(strategies, function(s) max(plot_data$Score[plot_data$Strategy == s], na.rm = TRUE) * 1.1)
)
p_box <- p_box + geom_text(data = star_data, aes(x = Strategy, y = y, label = star),
                           inherit.aes = FALSE, size = 6, color = "black")

ggsave(file.path(script_dir, "Fig.4B_YAS_scores_boxplot.pdf"), p_box, width = 50, height = 50, units = "mm")
ggsave(file.path(script_dir, "Fig.4B_YAS_scores_boxplot.png"), p_box, width = 50, height = 50, units = "mm", dpi = 300)

cat("Wilcoxon test p-values (independent scores):\n")
for (i in seq_along(strategies)) {
  cat(strategies[i], ": p =", signif(pvals[i], 4), "\n")
}

# ============================================================================
# PART 2: YAS pathway log2(YC/HY) bar chart (Fig.4A)
# ============================================================================
cat("\n========== Fig.4A: YAS pathway log2 ratio (YC vs HY) ==========\n")

pval_list <- list()
for (pid in unique(dat$Pathway_ID)) {
  sub <- dat %>% filter(Pathway_ID == pid)
  if (length(unique(sub$Site)) < 2) next
  pval <- wilcox.test(Abundance ~ Site, data = sub)$p.value
  pval_list[[pid]] <- pval
}
pval_df <- data.frame(Pathway_ID = names(pval_list), p_value = unlist(pval_list), stringsAsFactors = FALSE)
pval_df$star <- ifelse(pval_df$p_value < 0.001, "***",
                       ifelse(pval_df$p_value < 0.01, "**",
                              ifelse(pval_df$p_value < 0.05, "*", "")))

mean_abund <- dat %>%
  group_by(Pathway_ID, Description, Strategy, Site) %>%
  summarise(Mean = mean(Abundance, na.rm = TRUE), .groups = "drop")

wide_abund <- mean_abund %>%
  pivot_wider(names_from = Site, values_from = Mean, values_fill = 0) %>%
  mutate(Ratio = ifelse(HY > 0, YC / HY, NA),
         Ratio = ifelse(is.infinite(Ratio), NA, Ratio)) %>%
  filter(!is.na(Ratio))

wide_abund <- left_join(wide_abund, pval_df, by = "Pathway_ID")
wide_abund$Strategy <- factor(wide_abund$Strategy, levels = c("Y", "A", "S"))
wide_abund <- wide_abund %>% arrange(Strategy, desc(Ratio))

wide_abund$Pathway_label <- sapply(strsplit(wide_abund$Description, ";\\s*"), function(x) tail(x, 1))
wide_abund$Pathway_label <- factor(wide_abund$Pathway_label, levels = unique(wide_abund$Pathway_label))

wide_abund <- wide_abund %>%
  mutate(Log2Ratio = log2(Ratio)) %>%
  filter(is.finite(Log2Ratio))

wide_abund <- wide_abund %>% arrange(Strategy, desc(Log2Ratio))
wide_abund$Pathway_label <- factor(wide_abund$Pathway_label,
                                   levels = unique(wide_abund$Pathway_label))

strategy_colors <- c("Y" = "#F7D092", "A" = "#96C6E0", "S" = "#EDA8A6")

p_left <- ggplot(wide_abund, aes(x = 1, y = Pathway_label, fill = Strategy)) +
  geom_tile(width = 0.8, color = NA) +
  scale_fill_manual(values = strategy_colors, guide = "none") +
  theme_void()

ratio_range <- max(wide_abund$Log2Ratio, na.rm = TRUE) - min(wide_abund$Log2Ratio, na.rm = TRUE)
offset <- ratio_range * 0.02

star_data2 <- wide_abund %>%
  filter(star != "") %>%
  mutate(star_x = ifelse(Log2Ratio >= 0, Log2Ratio + offset, Log2Ratio - offset))

p_right <- ggplot(wide_abund, aes(x = Log2Ratio, y = Pathway_label, fill = Strategy)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
  scale_fill_manual(values = strategy_colors, guide = "none") +
  labs(x = expression(log[2]~"(YC / HY)"), y = NULL,
       title = "KEGG level-3 pathway abundance: YC vs HY\n(* p<0.05, ** p<0.01, *** p<0.001, Wilcoxon test)") +
  theme_bw(base_size = 10) +
  theme(axis.text.y = element_text(size = 8),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_text(data = star_data2, aes(x = star_x, label = star),
            hjust = 0.5, size = 4, color = "black") +
  expand_limits(x = c(min(wide_abund$Log2Ratio) * 1.1, max(wide_abund$Log2Ratio) * 1.1))

p_combined <- p_left + p_right + plot_layout(widths = c(0.05, 1))

ggsave(file.path(script_dir, "Fig.4A_YAS_pathway_log2ratio.pdf"), p_combined, width = 8, height = 6)
ggsave(file.path(script_dir, "Fig.4A_YAS_pathway_log2ratio.png"), p_combined, width = 8, height = 6, dpi = 300)

cat("Pathway label order:\n")
print(levels(wide_abund$Pathway_label))

# ============================================================================
# PART 3: Soil enzyme activity synergy bar chart (Fig.4C)
# ============================================================================
cat("\n========== Fig.4C: Enzyme synergy barplot ==========\n")

soil <- read.csv(file.path(out_root, "data/soil_properties_topsoil.csv"),
                 row.names = 1, stringsAsFactors = FALSE)
group <- read.csv(file.path(out_root, "data/group.csv"), row.names = 1)
group$Site <- factor(group$Site, levels = c("HY", "YC"))
group$Treatment <- factor(group$Treatment, levels = c("Ctrl", "AMF", "Str", "AMF+Str"))

group$Niche <- trimws(group$Niche)
group$Niche[group$Niche == "Bottom Soil"] <- "Bottom soil"
group$Niche[group$Niche == "Top Soil"] <- "Topsoil"

depth <- "Topsoil"

samples_use <- rownames(group)[group$Niche == depth & group$Treatment != "Ctrl"]
common <- intersect(samples_use, rownames(soil))
soil_sub <- soil[common, ]
meta_sub <- group[common, ]

enzymes <- c("NAG", "BG", "LAP", "AKP", "Urease")
enzymes <- intersect(enzymes, colnames(soil_sub))

enzyme_data <- soil_sub[, enzymes, drop = FALSE]
enzyme_data[] <- lapply(enzyme_data, as.numeric)

# Extract replicate ID
meta_sub$Rep <- gsub(".*\\.", "", rownames(meta_sub))

synergy_individual <- data.frame()

for (site in c("HY", "YC")) {
  for (rep_id in unique(meta_sub$Rep[meta_sub$Site == site])) {
    sample_amf <- rownames(meta_sub)[meta_sub$Site == site & meta_sub$Treatment == "AMF" & meta_sub$Rep == rep_id]
    sample_str <- rownames(meta_sub)[meta_sub$Site == site & meta_sub$Treatment == "Str" & meta_sub$Rep == rep_id]
    sample_comb <- rownames(meta_sub)[meta_sub$Site == site & meta_sub$Treatment == "AMF+Str" & meta_sub$Rep == rep_id]

    if (length(sample_comb) == 0) next

    for (enz in enzymes) {
      if (length(sample_amf) > 0) {
        ratio_amf <- enzyme_data[sample_comb, enz] / enzyme_data[sample_amf, enz]
        synergy_individual <- rbind(synergy_individual, data.frame(
          Site = site, Rep = rep_id, Enzyme = enz,
          Contrast = "vs_AMF", Ratio = ratio_amf, stringsAsFactors = FALSE
        ))
      }
      if (length(sample_str) > 0) {
        ratio_str <- enzyme_data[sample_comb, enz] / enzyme_data[sample_str, enz]
        synergy_individual <- rbind(synergy_individual, data.frame(
          Site = site, Rep = rep_id, Enzyme = enz,
          Contrast = "vs_Str", Ratio = ratio_str, stringsAsFactors = FALSE
        ))
      }
    }
  }
}

synergy_individual$log2FC <- log2(synergy_individual$Ratio)
synergy_individual <- synergy_individual[is.finite(synergy_individual$log2FC), ]

# Wilcoxon test
results <- data.frame()
for (enz in enzymes) {
  sub <- synergy_individual %>% filter(Enzyme == enz)
  if (nrow(sub) < 6) next
  mean_hy <- mean(sub$log2FC[sub$Site == "HY"], na.rm = TRUE)
  mean_yc <- mean(sub$log2FC[sub$Site == "YC"], na.rm = TRUE)
  test <- wilcox.test(log2FC ~ Site, data = sub, exact = FALSE)
  results <- rbind(results, data.frame(
    Enzyme = enz, HY_mean = mean_hy, YC_mean = mean_yc,
    pvalue = test$p.value, stringsAsFactors = FALSE
  ))
}

results$star <- ifelse(results$pvalue < 0.001, "***",
                       ifelse(results$pvalue < 0.01, "**",
                              ifelse(results$pvalue < 0.05, "*", "ns")))

plot_stats <- synergy_individual %>%
  group_by(Enzyme, Site) %>%
  summarise(
    Mean = mean(log2FC, na.rm = TRUE),
    SE   = sd(log2FC, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  left_join(results[, c("Enzyme", "pvalue", "star")], by = "Enzyme")

plot_stats$Enzyme <- factor(plot_stats$Enzyme, levels = c("NAG", "BG", "LAP", "AKP", "Urease"))
plot_stats$Site <- factor(plot_stats$Site, levels = c("HY", "YC"))

plot_stats <- plot_stats %>%
  group_by(Enzyme) %>%
  mutate(y_star = max(Mean + SE) * 1.15) %>%
  ungroup()

site_colors <- c("HY" = "#1F78B4", "YC" = "#E31A23")

p_enzyme <- ggplot(plot_stats, aes(x = Enzyme, y = Mean, fill = Site)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7,
           color = "black", linewidth = 0.4) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                position = position_dodge(width = 0.8), width = 0.2, linewidth = 0.6) +
  geom_text(aes(x = Enzyme, y = y_star, label = star),
            position = position_dodge(width = 0.8),
            size = 6, vjust = 0, color = "black", family = "sans") +
  scale_fill_manual(values = site_colors, name = "Site") +
  labs(x = NULL, y = expression(log[2]~"(Fold Change: AMF+Str / Single)"),
       title = "Soil enzyme activity synergy: HY vs YC") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  theme_classic(base_size = 14) +
  theme(axis.line = NULL,
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.title = element_text(size = 16, hjust = 0.5, face = "plain", color = "black"),
        axis.text.x = element_text(angle = 0, vjust = 0.5, face = "plain", color = "black", size = 13),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, face = "plain", color = "black"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        text = element_text(family = "sans"),
        plot.margin = margin(10, 10, 10, 10, unit = "pt")
  )

ggsave(file.path(script_dir, "Fig.4C_Enzyme_synergy_barplot.pdf"), p_enzyme, width = 5, height = 4.5)
ggsave(file.path(script_dir, "Fig.4C_Enzyme_synergy_barplot.png"), p_enzyme, width = 5, height = 4.5, dpi = 300)

cat("\nEnzyme synergy results:\n")
print(results[, c("Enzyme", "HY_mean", "YC_mean", "pvalue", "star")], row.names = FALSE)

cat("\n========== Fig.4A-C Analysis Complete ==========\n")
cat("Output files saved to:", script_dir, "\n")
