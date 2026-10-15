# ============================================================================
# Fig.S6: Soil moisture by growth stage (VT & R6) — HY vs YC + treatment effects
# ============================================================================
library(ggplot2)
library(dplyr)
library(tidyr)
library(agricolae)
library(patchwork)

out_root <- "C:/Users/ADMIN/Desktop/Co-inoculation/Fig.S6"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
setwd(out_root)

# ---------- 1. Load data ----------
df <- read.csv("data/soil_Moisture.csv", header = TRUE, stringsAsFactors = FALSE)
df <- df[df$Niche == "Topsoil" & df$Treatment %in% c("AMF", "Str", "AMF+Str"), ]
df$Site <- factor(df$Site, levels = c("HY", "YC"))
df$Treatment <- factor(df$Treatment, levels = c("AMF", "Str", "AMF+Str"))

# Reshape to long format
df_long <- df %>%
  pivot_longer(cols = c(Moisture_VT, Moisture_R6),
               names_to = "Stage", values_to = "Moisture") %>%
  mutate(Stage = factor(Stage,
                        levels = c("Moisture_VT", "Moisture_R6"),
                        labels = c("VT stage", "R6 stage")))

site_colors <- c("HY" = "#1F78B4", "YC" = "#E31A23")
trt_colors  <- c("AMF" = "#F0C986", "Str" = "#3FBDA7", "AMF+Str" = "#F08080")

# ============================================================================
# Part A: HY vs YC comparison (per stage) — boxplot, unified theme
# ============================================================================
cat("\n========== HY vs YC comparison (Wilcoxon rank-sum test) ==========\n")

for (st in levels(df_long$Stage)) {
  sub <- df_long[df_long$Stage == st, ]
  n_hy <- sum(sub$Site == "HY")
  n_yc <- sum(sub$Site == "YC")
  med_hy <- median(sub$Moisture[sub$Site == "HY"], na.rm = TRUE)
  med_yc <- median(sub$Moisture[sub$Site == "YC"], na.rm = TRUE)
  test <- wilcox.test(Moisture ~ Site, data = sub, exact = FALSE, conf.int = TRUE)
  star <- ifelse(test$p.value < 0.001, "***",
                 ifelse(test$p.value < 0.01, "**",
                        ifelse(test$p.value < 0.05, "*", "ns")))
  cat(sprintf("\n%s:\n", st))
  cat(sprintf("  n(HY) = %d, n(YC) = %d\n", n_hy, n_yc))
  cat(sprintf("  Median(HY) = %.2f %%, Median(YC) = %.2f %%\n", med_hy, med_yc))
  cat(sprintf("  Median difference (YC - HY) = %.2f %%\n", test$estimate))
  cat(sprintf("  95%% confidence interval: [%.2f, %.2f]\n",
              test$conf.int[1], test$conf.int[2]))
  cat(sprintf("  Wilcoxon W = %.2f, p = %.4f %s\n", test$statistic, test$p.value, star))
}

star_data <- df_long %>%
  group_by(Stage) %>%
  summarise(y_max = max(Moisture, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    y_star = y_max * 1.08,
    y_line = y_max * 1.04,
    star = sapply(levels(Stage), function(st) {
      test <- wilcox.test(Moisture ~ Site,
                          data = df_long[df_long$Stage == st, ], exact = FALSE)
      ifelse(test$p.value < 0.001, "***",
             ifelse(test$p.value < 0.01, "**",
                    ifelse(test$p.value < 0.05, "*", "ns")))
    })
  )

p_hy_vs_yc <- ggplot(df_long, aes(x = Site, y = Moisture, fill = Site)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.5, linewidth = 0.3) +
  geom_jitter(width = 0.12, size = 1.2, alpha = 0.5, color = "grey30") +
  facet_wrap(~ Stage) +
  geom_text(data = star_data, aes(x = 1.5, y = y_star, label = star),
            inherit.aes = FALSE, size = 3, color = "black") +
  geom_segment(data = star_data,
               aes(x = 1, xend = 2, y = y_line, yend = y_line),
               inherit.aes = FALSE, linewidth = 0.3) +
  scale_fill_manual(values = site_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Soil moisture (%)") +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.15),
    axis.line = element_blank(),
    axis.ticks = element_line(linewidth = 0.15),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5,
                               face = "plain", color = "black", size = 6),
    axis.text.y = element_text(color = "black", size = 6),
    axis.title.y = element_text(size = 6, color = "black"),
    axis.title.x = element_blank(),
    strip.background = element_rect(color = "black", fill = NA, linewidth = 0.15),
    strip.text = element_text(size = 6, face = "bold", color = "black"),
    plot.title = element_blank(),
    legend.position = "none",
    text = element_text(family = "sans", size = 6),
    plot.margin = margin(2, 2, 2, 2, unit = "pt")
  )

# ============================================================================
# Part B: Treatment comparison within each Site x Stage — four panels in one row
# ============================================================================
cat("\n========== Treatment comparisons (ANOVA + Tukey HSD) ==========\n")

compute_letters <- function(data) {
  model <- aov(Moisture ~ Treatment, data = data)
  hsd <- HSD.test(model, "Treatment", group = TRUE, alpha = 0.05)
  groups <- hsd$groups
  letter_col <- colnames(groups)[2]
  means <- hsd$means
  res <- data.frame(
    Treatment = rownames(groups),
    mean = groups[[1]],
    groups = groups[[letter_col]],
    se = means[rownames(groups), "se"],
    stringsAsFactors = FALSE
  )
  rownames(res) <- NULL
  return(res)
}

letter_list <- list()
for (st in levels(df_long$Stage)) {
  for (site in c("HY", "YC")) {
    sub <- df_long[df_long$Stage == st & df_long$Site == site, ]
    if (nrow(sub) < 6) next
    res <- compute_letters(sub)
    res$Stage <- st
    res$Site <- site
    letter_list[[length(letter_list) + 1]] <- res
  }
}
letters_all <- bind_rows(letter_list)
letters_all$Treatment <- factor(letters_all$Treatment, levels = c("AMF", "Str", "AMF+Str"))

letters_all <- letters_all %>%
  group_by(Stage, Site) %>%
  mutate(label_y = mean + se * 1.5) %>%
  ungroup()

cat("Tukey HSD results:\n")
print(letters_all)

letters_all <- letters_all %>%
  mutate(
    Site_Stage = factor(
      paste(Site, Stage, sep = "_"),
      levels = c("HY_VT stage", "HY_R6 stage", "YC_VT stage", "YC_R6 stage")
    )
  )

p_trt <- ggplot(letters_all,
                aes(x = Treatment, y = mean, fill = Treatment)) +
  geom_col(color = "black", width = 0.6, linewidth = 0.15) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.15, linewidth = 0.15) +
  geom_text(aes(y = label_y, label = groups),
            size = 1.8, fontface = "plain", color = "black", family = "sans") +
  facet_wrap(~ Site_Stage, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = trt_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = NULL, y = "Soil moisture (%)") +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.15),
    axis.line = element_blank(),
    axis.ticks = element_line(linewidth = 0.15),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1,
                               face = "plain", color = "black", size = 6),
    axis.text.y = element_text(color = "black", size = 6),
    axis.title.y = element_text(size = 6, color = "black"),
    axis.title.x = element_blank(),
    strip.background = element_rect(color = "black", fill = NA, linewidth = 0.15),
    strip.text = element_text(size = 6, face = "bold", color = "black"),
    plot.title = element_blank(),
    legend.position = "none",
    text = element_text(family = "sans", size = 6),
    plot.margin = margin(2, 2, 2, 2, unit = "pt")
  )

# ============================================================================
# Combine Part B (left) and Part A (right) into one figure
# ============================================================================
combined <- p_trt + p_hy_vs_yc + plot_layout(widths = c(2, 1))

ggsave(file.path(out_root, "Fig.S6_moisture_combined.pdf"),
       combined, width = 140, height = 40, units = "mm")
ggsave(file.path(out_root, "Fig.S6_moisture_combined.png"),
       combined, width = 140, height = 40, units = "mm", dpi = 300)

cat("\n========== Fig.S6 combined plot saved ==========\n")
# Save Tukey results table (optional)
write.csv(letters_all, file.path(out_root, "Fig.S6_moisture_tukey_results.csv"), row.names = FALSE)

