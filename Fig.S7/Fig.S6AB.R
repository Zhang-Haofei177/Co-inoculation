# ============================================================================
# 宏基因组物种丰度箱线图：基于 Unigenes.absolute.g.xls（相对丰度修正版）
# 真菌：Glomeromycota（门/纲）+ Rhizophagus（属）
# 细菌：Actinomycetota（门）+ Streptomyces（属）
# 稳健版 get_letters + 终端显著性表 + 检出率表
# ============================================================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(agricolae)

# ---------- 0. 路径 ----------
setwd("C:/Users/ADMIN/Desktop/Co-inoculation/Fig.S7/")
out_dir <- "C:/Users/ADMIN/Desktop/Co-inoculation/Fig.S7"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- 1. 读取分组信息 ----------
group <- read.csv("data/group.csv", row.names = 1)
group$Treatment <- factor(group$Treatment, levels = c("Ctrl", "AMF", "Str", "AMF+Str"))
group$Site      <- factor(group$Site,      levels = c("HY", "YC"))

target_trt   <- c("AMF", "Str", "AMF+Str")
target_niche <- c("Rhizosphere")
group_sel <- group[group$Treatment %in% target_trt & group$Niche %in% target_niche, ]
group_sel$Treatment <- droplevels(group_sel$Treatment)

# ---------- 2. 读取属绝对丰度表 ----------
raw <- read.table("data/Unigenes.absolute.g.xls",
                  header = TRUE, sep = "\t", row.names = 1,
                  quote = "", comment.char = "", stringsAsFactors = FALSE)

if ("Others" %in% rownames(raw)) {
  raw <- raw[rownames(raw) != "Others", ]
}

desc_col <- which(tolower(colnames(raw)) == "description")
if (length(desc_col) != 1) stop("未找到唯一的 Description 列")
sample_cols <- setdiff(colnames(raw), "Description")

desc   <- raw$Description
abund  <- raw[, sample_cols, drop = FALSE]
abund[] <- lapply(abund, as.numeric)

# 转置为 样本×属，然后每行（样本）除以该样本总丰度，再乘以 100
# 注意：宏基因组的相对丰度分母是「全部微生物总丰度」，
#       而扩增子是真菌、细菌各自群落内计算，二者不能直接对比
otu_t     <- as.data.frame(t(abund))
otu_t_rel <- otu_t / rowSums(otu_t) * 100
otu_t_rel[is.na(otu_t_rel)] <- 0

# ---------- 3. 识别目标类群 ----------
# 3.1 Glomeromycota（门/纲，作背景参考）
glom_idx <- grepl("Glomeromycota|Glomeromycetes|Glomerales",
                  desc, ignore.case = TRUE)

# 3.2 解析属名（兼容 g__ 与 g: 两种写法）
parse_genus <- function(tax) {
  if (is.na(tax) || tax == "") return("")
  m <- regmatches(tax, regexec("g__([^;]+)", tax))[[1]]
  if (length(m) > 1) return(trimws(m[2]))
  m <- regmatches(tax, regexec("g:([^;]+)", tax))[[1]]
  if (length(m) > 1) return(trimws(m[2]))
  ""
}
genera <- sapply(desc, parse_genus)

# 3.3 Rhizophagus（属，AMF 指示类群）
rhizo_idx <- grepl("^Rhizophagus$", genera, ignore.case = TRUE)

# 3.4 Streptomyces（属）
str_idx <- grepl("^Streptomyces$", genera, ignore.case = TRUE)

# 3.5 Actinomycetota（门，作背景参考）
act_idx <- grepl("Actinomycetota|Actinobacteria", desc, ignore.case = TRUE)

cat("Glomeromycota（门/纲） 匹配到属数：", sum(glom_idx), "\n")
cat("Rhizophagus     （属） 匹配到属数：", sum(rhizo_idx), "\n")
cat("Streptomyces    （属） 匹配到属数：", sum(str_idx), "\n")
cat("Actinomycetota  （门） 匹配到属数：", sum(act_idx), "\n")

glom_abund  <- rowSums(otu_t_rel[, glom_idx,  drop = FALSE], na.rm = TRUE)
rhizo_abund <- rowSums(otu_t_rel[, rhizo_idx, drop = FALSE], na.rm = TRUE)
str_abund   <- rowSums(otu_t_rel[, str_idx,   drop = FALSE], na.rm = TRUE)
act_abund   <- rowSums(otu_t_rel[, act_idx,   drop = FALSE], na.rm = TRUE)

# ---------- 4. 合并数据 ----------
common_samples <- intersect(rownames(group_sel), rownames(otu_t_rel))
meta <- group_sel[common_samples, c("Site", "Niche", "Treatment")]
meta$Sample <- rownames(meta)

plot_data <- meta %>%
  mutate(
    Glomeromycota  = glom_abund[Sample],
    Rhizophagus    = rhizo_abund[Sample],
    Streptomyces   = str_abund[Sample],
    Actinomycetota = act_abund[Sample]
  ) %>%
  replace_na(list(Glomeromycota = 0, Rhizophagus = 0,
                  Streptomyces = 0, Actinomycetota = 0)) %>%
  mutate(
    Site      = factor(Site,      levels = c("HY", "YC")),
    Treatment = factor(Treatment, levels = c("AMF", "Str", "AMF+Str"))
  )

cat("验证相对丰度总和（样本1）：", sum(otu_t_rel[1, ]), "\n")

# ---------- 5. 绘图设置 + 稳健字母函数 ----------
colors_no <- c("AMF" = "#3FBDA7", "Str" = "#F08080", "AMF+Str" = "#F0C986")
expected_levels <- c("AMF", "Str", "AMF+Str")

# 稳健版 get_letters（带三层保护）
#   情形 A：响应值全相同（含全 0）——直接全部标 "a"
#   情形 B：正常 ANOVA + HSD
#   情形 C：HSD 因某组方差为 0 失败——退回全部标 "a"
# 想让退化格子里留空字母，把 letters = "a" 改成 letters = ""
get_letters_robust <- function(data, response, all_levels = expected_levels) {
  data$Treatment <- droplevels(data$Treatment)
  present <- levels(data$Treatment)
  
  # 空 / 单组
  if (nrow(data) < 2 || length(present) < 2) {
    res <- data.frame(Treatment = all_levels, letters = "", stringsAsFactors = FALSE)
    means <- tapply(data[[response]], data$Treatment, mean, na.rm = TRUE)
    res$Mean <- as.numeric(means[all_levels])
    res$letters[res$Treatment %in% present] <- "a"
    return(res)
  }
  
  y <- data[[response]]
  
  # 情形 A：整组响应值全相同（含全 0），方差为 0
  if (length(unique(y[!is.na(y)])) <= 1) {
    means <- tapply(y, data$Treatment, mean, na.rm = TRUE)
    return(data.frame(Treatment = all_levels,
                      Mean      = as.numeric(means[all_levels]),
                      letters   = "a",
                      stringsAsFactors = FALSE))
  }
  
  # 情形 B：正常 ANOVA + HSD，tryCatch 兜底
  formula <- as.formula(paste(response, "~ Treatment"))
  model   <- aov(formula, data = data)
  hsd <- tryCatch(
    agricolae::HSD.test(model, "Treatment", group = TRUE),
    error   = function(e) NULL,
    warning = function(w) NULL
  )
  
  # 情形 C：HSD 失败
  if (is.null(hsd) || is.null(hsd$groups)) {
    means <- tapply(y, data$Treatment, mean, na.rm = TRUE)
    return(data.frame(Treatment = all_levels,
                      Mean      = as.numeric(means[all_levels]),
                      letters   = "a",
                      stringsAsFactors = FALSE))
  }
  
  letters_df <- hsd$groups
  mean_col   <- setdiff(names(letters_df), "groups")[1]
  letters_df <- data.frame(
    Treatment = rownames(letters_df),
    Mean      = letters_df[[mean_col]],
    letters   = letters_df$groups,
    stringsAsFactors = FALSE
  )
  rownames(letters_df) <- NULL
  
  missing <- setdiff(all_levels, letters_df$Treatment)
  if (length(missing) > 0) {
    letters_df <- bind_rows(letters_df,
                            data.frame(Treatment = missing, Mean = NA,
                                       letters = "", stringsAsFactors = FALSE))
  }
  letters_df <- letters_df[match(all_levels, letters_df$Treatment), ]
  rownames(letters_df) <- NULL
  letters_df
}

# ---------- 6. 终端输出函数 ----------
print_sig_table <- function(data, var_name) {
  combo_list <- split(data, data$Site)
  results <- bind_rows(lapply(names(combo_list), function(nm) {
    sub <- combo_list[[nm]]
    if (nrow(sub) < 2) return(NULL)
    res <- get_letters_robust(sub, var_name)
    res$Site <- nm
    res
  }))
  results <- results %>%
    dplyr::select(Site, Treatment, Mean, letters) %>%
    mutate(Mean = round(Mean, 4)) %>%
    arrange(Site, Treatment)
  
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("指标：", var_name, "  （Site 内 Treatment 间 ANOVA + Tukey HSD）\n", sep = "")
  cat(strrep("=", 70), "\n", sep = "")
  print(as.data.frame(results), row.names = FALSE)
  
  wide <- results %>%
    mutate(Val = paste0(Mean, " (", letters, ")")) %>%
    dplyr::select(Site, Treatment, Val) %>%
    tidyr::pivot_wider(names_from = Treatment, values_from = Val)
  cat("\n宽表（Mean (letter)）：\n")
  print(as.data.frame(wide), row.names = FALSE)
  cat("\n")
  
  invisible(results)
}

print_detect_table <- function(data, var_name) {
  tab <- data %>%
    group_by(Site, Treatment) %>%
    summarise(
      n          = n(),
      n_detected = sum(.data[[var_name]] > 0, na.rm = TRUE),
      mean       = round(mean(.data[[var_name]], na.rm = TRUE), 4),
      median     = round(median(.data[[var_name]], na.rm = TRUE), 4),
      .groups    = "drop"
    ) %>%
    arrange(Site, Treatment)
  
  cat("\n", strrep("-", 70), "\n", sep = "")
  cat("检出率：", var_name, "  （n_detected / n）\n", sep = "")
  cat(strrep("-", 70), "\n", sep = "")
  print(as.data.frame(tab), row.names = FALSE)
  cat("\n")
  
  invisible(tab)
}

# ---------- 7. 绘图函数 ----------
plot_one_var <- function(data, var_name, y_label, file_name) {
  combo_list <- split(data, data$Site)
  letters_all <- bind_rows(lapply(names(combo_list), function(nm) {
    sub <- combo_list[[nm]]
    res <- get_letters_robust(sub, var_name)
    res$Site <- nm
    res
  }))
  
  ymax_data <- data %>%
    group_by(Site, Treatment) %>%
    summarise(y_max = max(.data[[var_name]], na.rm = TRUE), .groups = "drop")
  
  sig_data <- ymax_data %>%
    left_join(letters_all, by = c("Site", "Treatment")) %>%
    mutate(y_label = y_max * 1.05)
  
  plot_df <- data %>%
    dplyr::select(Site, Treatment, Value = all_of(var_name))
  
  p <- ggplot(plot_df, aes(x = Treatment, y = Value)) +
    geom_boxplot(aes(color = Treatment), outlier.shape = NA, alpha = 0.7) +
    geom_jitter(aes(color = Treatment), width = 0.15, size = 1.5, alpha = 0.6) +
    facet_wrap(~ Site, scales = "free_y", ncol = 2) +
    geom_text(data = sig_data,
              aes(x = Treatment, y = y_label, label = letters),
              size = 3.5, vjust = 0, inherit.aes = FALSE) +
    scale_color_manual(values = colors_no) +
    labs(title = var_name, x = NULL, y = y_label) +
    theme_classic() +
    theme(
      text             = element_text(family = "sans"),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.title       = element_text(size = 14, hjust = 0.5, color = "#333333",
                                      face = "bold",
                                      margin = ggplot2::margin(b = 10)),
      axis.text.x      = element_text(angle = 0, vjust = 0.5,
                                      face = "bold", color = "#333333"),
      plot.margin      = ggplot2::margin(10, 10, 10, 10, unit = "pt"),
      legend.position  = "none"
    )
  
  ggsave(file.path(out_dir, paste0(file_name, ".pdf")), p, width = 4, height = 2.5)
  ggsave(file.path(out_dir, paste0(file_name, ".png")), p, width = 4, height = 2.5, dpi = 300)
  return(p)
}

# ---------- 8. 绘制四张图 ----------
p_glom  <- plot_one_var(plot_data, "Glomeromycota",
                        "Relative abundance of Glomeromycota (%)",
                        "Glomeromycota_metagenome_boxplot")
p_rhizo <- plot_one_var(plot_data, "Rhizophagus",
                        "Relative abundance of Rhizophagus (%)",
                        "Rhizophagus_metagenome_boxplot")
p_str   <- plot_one_var(plot_data, "Streptomyces",
                        "Relative abundance of Streptomyces (%)",
                        "Streptomyces_metagenome_boxplot")
p_act   <- plot_one_var(plot_data, "Actinomycetota",
                        "Relative abundance of Actinomycetota (%)",
                        "Actinomycetota_metagenome_boxplot")

# ---------- 9. 终端输出显著性表 + 检出率表 ----------
cat("\n########## 检出率概览 ##########\n")
det_glom  <- print_detect_table(plot_data, "Glomeromycota")
det_rhizo <- print_detect_table(plot_data, "Rhizophagus")
det_str   <- print_detect_table(plot_data, "Streptomyces")
det_act   <- print_detect_table(plot_data, "Actinomycetota")

cat("\n########## 显著性表（ANOVA + Tukey HSD）##########\n")
sig_glom  <- print_sig_table(plot_data, "Glomeromycota")
sig_rhizo <- print_sig_table(plot_data, "Rhizophagus")
sig_str   <- print_sig_table(plot_data, "Streptomyces")
sig_act   <- print_sig_table(plot_data, "Actinomycetota")

# ---------- 10. 保存 CSV ----------
write.csv(sig_glom,  file.path(out_dir, "sig_Glomeromycota_metagenome.csv"),  row.names = FALSE)
write.csv(sig_rhizo, file.path(out_dir, "sig_Rhizophagus_metagenome.csv"),    row.names = FALSE)
write.csv(sig_str,   file.path(out_dir, "sig_Streptomyces_metagenome.csv"),   row.names = FALSE)
write.csv(sig_act,   file.path(out_dir, "sig_Actinomycetota_metagenome.csv"), row.names = FALSE)

write.csv(det_glom,  file.path(out_dir, "detect_Glomeromycota_metagenome.csv"),  row.names = FALSE)
write.csv(det_rhizo, file.path(out_dir, "detect_Rhizophagus_metagenome.csv"),    row.names = FALSE)
write.csv(det_str,   file.path(out_dir, "detect_Streptomyces_metagenome.csv"),   row.names = FALSE)
write.csv(det_act,   file.path(out_dir, "detect_Actinomycetota_metagenome.csv"), row.names = FALSE)

cat("\n四张箱线图 + 显著性/检出率表格已保存至：", out_dir, "\n")
