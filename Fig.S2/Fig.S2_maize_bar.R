# ============================================================================
# Fig.S2: Maize phenotype bar chart + 差异分析结果导出
# ============================================================================
library(ggplot2)
library(dplyr)
library(tidyr)
library(agricolae)
library(patchwork)

out_root <- "C:/Users/ADMIN/Desktop/Co-inoculation/Fig.S2"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
setwd(out_root)

# ---------- 1. Data preparation ----------
df <- read.csv("data/maize_phenotype_summary.csv", header = TRUE)
df <- df[grepl("Topsoil", df$Niche), ]
df <- df[!grepl("Ctrl", df$Treatment), ]
df$Treatment <- factor(df$Treatment, levels = c("AMF", "Str", "AMF+Str"))

var_units <- c(
  "Yield_t_ha"            = "Yield (t/ha)",
  "Biomass_t_ha"          = "Biomass (t/ha)",
  "WUE..kg.mm.ha."        = "Water use efficiency (kg/mm/ha)",
  "NUE."                  = "Nitrogen use efficiency (%)",
  "PUE_kg_mm"             = "Precipitation use efficiency (kg/mm/ha)",
  "Grain_N_perc"          = "Grain nitrogen content (%)",
  "Straw_N_perc"          = "Straw nitrogen content (%)",
  "Crude.fat.content"     = "Crude fat content (%)",
  "Starch.content"        = "Starch content (%)",
  "Crude.protein.content" = "Crude protein content (%)"
)

# 简洁变量标签（用于导出，便于转述给 AI）
var_labels <- c(
  "Yield_t_ha"            = "Yield",
  "Biomass_t_ha"          = "Biomass",
  "WUE..kg.mm.ha."        = "WUE",
  "NUE."                  = "NUE",
  "PUE_kg_mm"             = "PUE",
  "Grain_N_perc"          = "Grain_N",
  "Straw_N_perc"          = "Straw_N",
  "Crude.fat.content"     = "Crude_fat",
  "Starch.content"        = "Starch",
  "Crude.protein.content" = "Crude_protein"
)

vars <- names(var_units)

# ---------- 2. Multiple comparisons + ANOVA p-value ----------
compute_letters <- function(data, variable) {
  formula <- as.formula(paste(variable, "~ Treatment"))
  model <- aov(formula, data = data)
  hsd <- HSD.test(model, "Treatment", group = TRUE, alpha = 0.05)
  groups <- hsd$groups
  letter_col <- colnames(groups)[2]
  means <- hsd$means
  
  # 提取 ANOVA 整体 p 值
  anova_p <- summary(model)[[1]][["Pr(>F)"]][1]
  
  res <- data.frame(
    Treatment = rownames(groups),
    mean = groups[[1]],
    groups = groups[[letter_col]],
    se = means[rownames(groups), "se"],
    anova_p = anova_p,
    stringsAsFactors = FALSE
  )
  rownames(res) <- NULL
  return(res)
}

result_list <- list()
for (site in c("HY", "YC")) {
  site_data <- filter(df, Site == site)
  for (var in vars) {
    stats <- compute_letters(site_data, var)
    stats$Variable <- var
    stats$Site <- site
    result_list[[length(result_list) + 1]] <- stats
  }
}
plot_data <- bind_rows(result_list)
plot_data$Variable <- factor(plot_data$Variable, levels = vars)

plot_data <- plot_data %>%
  group_by(Site, Variable) %>%
  mutate(label_y = mean + se * 1.2) %>%
  ungroup()

plot_data <- plot_data %>%
  mutate(Site_Treatment = paste(Site, Treatment, sep = "_")) %>%
  mutate(Site_Treatment = factor(Site_Treatment,
                                 levels = c(paste("HY", levels(df$Treatment), sep = "_"),
                                            paste("YC", levels(df$Treatment), sep = "_"))))

# ============================================================================
# 新增：差异分析结果导出
# ============================================================================
dir.create(file.path(out_root, "stats"), showWarnings = FALSE, recursive = TRUE)

# --- 导出 1：长表（每行 = 位点 × 变量 × 处理）---
stats_long <- plot_data %>%
  mutate(Variable_label = var_labels[as.character(Variable)]) %>%
  dplyr::select(Site, Variable, Variable_label, Treatment,
                mean, se, groups, anova_p) %>%
  mutate(
    mean    = round(mean, 4),
    se      = round(se, 4),
    anova_p = signif(anova_p, 4)
  ) %>%
  arrange(Site, Variable, Treatment)

write.csv(stats_long,
          file.path(out_root, "stats", "FigS2_maize_stats_long.csv"),
          row.names = FALSE)

# --- 导出 2：宽表（行 = 位点 × 变量，列 = 每个处理的 mean±se 和字母 + ANOVA p）---
stats_wide <- stats_long %>%
  mutate(
    Val = sprintf("%.3f ± %.3f (%s)", mean, se, groups)
  ) %>%
  dplyr::select(Site, Variable, Variable_label, Treatment, Val, anova_p) %>%
  pivot_wider(names_from = Treatment, values_from = Val) %>%
  mutate(ANOVA_p = signif(anova_p, 4)) %>%
  dplyr::select(Site, Variable, Variable_label, AMF, Str, `AMF+Str`, ANOVA_p) %>%
  arrange(Site, Variable)

write.csv(stats_wide,
          file.path(out_root, "stats", "FigS2_maize_stats_wide.csv"),
          row.names = FALSE)

# --- 终端打印宽表
cat("\n\n")
cat(strrep("=", 100), "\n", sep = "")
cat("Fig.S2 玉米表型差异分析结果（每个位点内 Treatment 间 ANOVA + Tukey HSD）\n")
cat("格式：均值 ± 标准误 (显著性字母)；ANOVA_p 为整体 F 检验 p 值\n")
cat(strrep("=", 100), "\n", sep = "")

for (site in c("HY", "YC")) {
  cat("\n\n########## Site:", site, "##########\n\n")
  sub_wide <- stats_wide %>% filter(Site == site) %>%
    dplyr::select(Variable_label, AMF, Str, `AMF+Str`, ANOVA_p)
  print(as.data.frame(sub_wide), row.names = FALSE)
}

cat("\n\nCSV 文件已保存至：", file.path(out_root, "stats"), "\n", sep = "")

# ---------- 3. Plot function ----------
plot_bar_combined <- function(data, var_name) {
  sub_data <- data %>% filter(Variable == var_name)
  sub_data$Site_Treatment <- factor(sub_data$Site_Treatment, levels = levels(data$Site_Treatment))
  
  colors <- c("#F0C986", "#3FBDA7", "#F08080")
  names(colors) <- levels(df$Treatment)
  
  ggplot(sub_data, aes(x = Site_Treatment, y = mean, fill = Treatment)) +
    annotate("rect", xmin = 3.5, xmax = 6.5, ymin = -Inf, ymax = Inf,
             fill = "grey90", alpha = 0.3) +
    geom_vline(xintercept = 3.5, linetype = "solid", color = "black", linewidth = 0.15) +
    geom_col(color = "black", width = 0.6, linewidth = 0.15) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                  width = 0.15, linewidth = 0.15) +
    geom_text(aes(y = label_y, label = groups),
              vjust = -0.3, size = 6, fontface = "plain",
              color = "black", family = "sans", size.unit = "pt") +
    scale_fill_manual(values = colors) +
    scale_x_discrete(labels = c("AMF","Str","AMF+Str", "AMF","Str","AMF+Str")) +
    labs(x = NULL, y = var_units[var_name]) +
    theme_classic() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.15),
      axis.line = element_blank(),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1,
                                 face = "plain", color = "black", size = 6),
      axis.text.y = element_text(color = "black", size = 6),
      axis.title.y = element_text(size = 6, color = "black"),
      axis.title.x = element_blank(),
      axis.ticks = element_line(linewidth = 0.1),
      plot.title = element_blank(),
      legend.position = "none",
      text = element_text(family = "sans", size = 6),
      plot.margin = margin(2, 2, 2, 2, unit = "pt")
    )
}

# ---------- 4. Combine and save ----------
plot_list <- lapply(vars, function(var) plot_bar_combined(plot_data, var))

combined <- wrap_plots(plot_list, ncol = 5, nrow = 2) +
  plot_annotation(title = "HY and YC Combined",
                  theme = theme(plot.title = element_text(size = 8, hjust = 0.5,
                                                          face = "plain", family = "sans")))

ggsave(file.path(out_root, "Fig.S2_maize_bar_test.pdf"), combined, width = 180, height = 90, units = "mm", dpi = 300)
ggsave(file.path(out_root, "Fig.S2_maize_bar.png"), combined, width = 180, height = 90, units = "mm", dpi = 300)

cat("\nFig.S2 saved\n")