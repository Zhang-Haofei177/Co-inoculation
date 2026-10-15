# ============================================================================
# 突出AMF和链霉菌：门/属指标分别绘制箱线图 + 终端显著性表格
# 真菌：Glomeromycota（门）+ Rhizophagus（属）
# 细菌：Actinomycetota（门）+ Streptomyces（属）
# 修复：Rhizophagus 全 0 导致 HSD.test 内部 pvalue=NA 报错
# ============================================================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(agricolae)

# ---------- 0. 路径 ----------
setwd("C:/Users/ADMIN/Desktop/Co-inoculation/Fig.S7/")
out_dir <- "C:/Users/ADMIN/Desktop/Co-inoculation/Fig.S7"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- 1. 分组 ----------
group <- read.csv("data/group.csv", row.names = 1)
group$Treatment <- factor(group$Treatment, levels = c("Ctrl", "AMF", "Str", "AMF+Str"))
group$Site <- factor(group$Site, levels = c("HY", "YC"))

target_trt   <- c("AMF", "Str", "AMF+Str")
target_niche <- c("Topsoil", "Rhizosphere", "Endosphere")

group_sel <- group[group$Treatment %in% target_trt & group$Niche %in% target_niche, ]
group_sel$Treatment <- droplevels(group_sel$Treatment)

# ---------- 2. 真菌：Glomeromycota（门）+ Rhizophagus（属） ----------
fun_otu <- read.csv("data/FUN_rarefied.csv", row.names = 1)
rownames(fun_otu) <- paste0("Fun_", rownames(fun_otu))
fun_otu_t <- as.data.frame(t(fun_otu))

fun_tax <- read.csv("data/FUN_Tax.csv", row.names = 1)
rownames(fun_tax) <- paste0("Fun_", rownames(fun_tax))

# 2.1 Glomeromycota 门
glom_asv <- rownames(fun_tax)[grepl("Glomeromycota", fun_tax$Phylum, ignore.case = TRUE)]
if (length(glom_asv) == 0) {
  warning("Phylum 中未匹配到 Glomeromycota，门水平结果将为空。")
}

# 2.2 Rhizophagus 属
rhizo_asv <- rownames(fun_tax)[grepl("Rhizophagus", fun_tax$Genus, ignore.case = TRUE)]
cat("Glomeromycota 门匹配到 ASV 数：", length(glom_asv), "\n")
cat("Rhizophagus   属匹配到 ASV 数：", length(rhizo_asv), "\n")

fun_total  <- rowSums(fun_otu_t)
common_fun <- intersect(rownames(group_sel), rownames(fun_otu_t))

fun_glom  <- fun_otu_t[common_fun, colnames(fun_otu_t) %in% glom_asv,  drop = FALSE]
fun_rhizo <- fun_otu_t[common_fun, colnames(fun_otu_t) %in% rhizo_asv, drop = FALSE]

glom_rel  <- rowSums(fun_glom)  / fun_total[common_fun] * 100
rhizo_rel <- rowSums(fun_rhizo) / fun_total[common_fun] * 100

# ---------- 3. 细菌：Actinomycetota（门）+ Streptomyces（属） ----------
bac_otu <- read.csv("data/BAC_rarefied.csv", row.names = 1)
rownames(bac_otu) <- paste0("Bac_", rownames(bac_otu))
bac_otu_t <- as.data.frame(t(bac_otu))

bac_tax <- read.csv("data/BAC_Tax.csv", row.names = 1)
rownames(bac_tax) <- paste0("Bac_", rownames(bac_tax))

str_asv <- rownames(bac_tax)[grepl("Streptomyces", bac_tax$Genus, ignore.case = TRUE)]
act_asv <- rownames(bac_tax)[grepl("Actinomycetota|Actinobacteria",
                                   bac_tax$Phylum, ignore.case = TRUE)]
cat("Streptomyces   属匹配到 ASV 数：", length(str_asv), "\n")
cat("Actinomycetota 门匹配到 ASV 数：", length(act_asv), "\n")

bac_total  <- rowSums(bac_otu_t)
common_bac <- intersect(rownames(group_sel), rownames(bac_otu_t))

bac_str <- bac_otu_t[common_bac, colnames(bac_otu_t) %in% str_asv, drop = FALSE]
bac_act <- bac_otu_t[common_bac, colnames(bac_otu_t) %in% act_asv, drop = FALSE]

str_rel <- rowSums(bac_str) / bac_total[common_bac] * 100
act_rel <- rowSums(bac_act) / bac_total[common_bac] * 100

# ---------- 4. 合并数据 ----------
meta <- group_sel[, c("Site", "Niche", "Treatment")]
meta$Sample <- rownames(meta)

plot_data <- meta %>%
  mutate(
    Glomeromycota  = glom_rel[Sample],
    Rhizophagus    = rhizo_rel[Sample],
    Streptomyces   = str_rel[Sample],
    Actinomycetota = act_rel[Sample]
  ) %>%
  replace_na(list(Glomeromycota = 0, Rhizophagus = 0,
                  Streptomyces = 0, Actinomycetota = 0)) %>%
  mutate(
    Site      = factor(Site, levels = c("HY", "YC")),
    Niche     = factor(Niche, levels = c("Topsoil", "Rhizosphere", "Endosphere")),
    Treatment = factor(Treatment, levels = c("AMF", "Str", "AMF+Str"))
  )

colors_no <- c("AMF" = "#3FBDA7", "Str" = "#F08080", "AMF+Str" = "#F0C986")

# ---------- 5. 字母计算函数（稳健版） ----------
# 说明：
#   - 若整组响应值为常数（含全 0）→ 直接全部标 "a"，避免 HSD.test 内部 pvalue=NA 报错；
#   - 若 HSD 因某组方差为 0 失败 → 退回全部标 "a"；
#   - 想在这些退化格子里留空，把 letters = "a" 改成 letters = "" 即可。
get_letters <- function(data, response) {
  data$Treatment <- droplevels(data$Treatment)
  all_trts <- levels(data$Treatment)
  
  # 空 / 单组情形
  if (nrow(data) < 2 || length(all_trts) < 2) {
    means <- tapply(data[[response]], data$Treatment, mean, na.rm = TRUE)
    return(data.frame(Treatment = all_trts,
                      Mean      = as.numeric(means[all_trts]),
                      letters   = "",
                      stringsAsFactors = FALSE))
  }
  
  y <- data[[response]]
  
  # 情形 A：响应值全相同（含全 0）——方差为 0
  if (length(unique(y[!is.na(y)])) <= 1) {
    means <- tapply(y, data$Treatment, mean, na.rm = TRUE)
    return(data.frame(Treatment = all_trts,
                      Mean      = as.numeric(means[all_trts]),
                      letters   = "a",
                      stringsAsFactors = FALSE))
  }
  
  # 情形 B：正常 ANOVA + HSD，用 tryCatch 兜底
  formula <- as.formula(paste(response, "~ Treatment"))
  model   <- aov(formula, data = data)
  
  hsd <- tryCatch(
    agricolae::HSD.test(model, "Treatment", group = TRUE),
    error   = function(e) NULL,
    warning = function(w) NULL
  )
  
  # 情形 C：HSD 失败 → 退回全部标 "a"
  if (is.null(hsd) || is.null(hsd$groups)) {
    means <- tapply(y, data$Treatment, mean, na.rm = TRUE)
    return(data.frame(Treatment = all_trts,
                      Mean      = as.numeric(means[all_trts]),
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
  
  missing <- setdiff(all_trts, letters_df$Treatment)
  if (length(missing) > 0) {
    letters_df <- bind_rows(
      letters_df,
      data.frame(Treatment = missing, Mean = NA, letters = "",
                 stringsAsFactors = FALSE)
    )
  }
  letters_df <- letters_df[match(all_trts, letters_df$Treatment), ]
  rownames(letters_df) <- NULL
  letters_df
}

# ---------- 6. 终端打印显著性表格 ----------
print_sig_table <- function(data, var_name) {
  combo_list <- split(data, list(data$Site, data$Niche))
  results <- lapply(names(combo_list), function(combo_name) {
    sub <- combo_list[[combo_name]]
    if (nrow(sub) < 2) return(NULL)
    res <- get_letters(sub, var_name)
    res$Site  <- sub$Site[1]
    res$Niche <- sub$Niche[1]
    res
  })
  results <- bind_rows(results)
  results <- results %>%
    dplyr::select(Site, Niche, Treatment, Mean, letters) %>%
    mutate(Mean = round(Mean, 4)) %>%
    arrange(Site, Niche, Treatment)
  
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("指标：", var_name, "  （Site × Niche 内 Treatment 间 ANOVA + Tukey HSD）\n", sep = "")
  cat(strrep("=", 70), "\n", sep = "")
  print(as.data.frame(results), row.names = FALSE)
  cat("\n")
  
  # 宽表（每个 Site × Niche 一行）
  wide <- results %>%
    mutate(Cell = paste0(Site, "_", Niche),
           Val  = paste0(Mean, " (", letters, ")")) %>%
    dplyr::select(Cell, Treatment, Val) %>%
    tidyr::pivot_wider(names_from = Treatment, values_from = Val)
  cat("宽表（Mean (letter)）：\n")
  print(as.data.frame(wide), row.names = FALSE)
  cat("\n")
  
  invisible(results)
}

# ---------- 6b. 检出率表（判断哪些格子不宜统计） ----------
print_detect_table <- function(data, var_name) {
  tab <- data %>%
    group_by(Site, Niche, Treatment) %>%
    summarise(
      n          = n(),
      n_detected = sum(.data[[var_name]] > 0, na.rm = TRUE),
      mean       = round(mean(.data[[var_name]], na.rm = TRUE), 4),
      median     = round(median(.data[[var_name]], na.rm = TRUE), 4),
      .groups    = "drop"
    ) %>%
    arrange(Site, Niche, Treatment)
  
  cat("\n", strrep("-", 70), "\n", sep = "")
  cat("检出率：", var_name, "  （n_detected / n）\n", sep = "")
  cat(strrep("-", 70), "\n", sep = "")
  print(as.data.frame(tab), row.names = FALSE)
  cat("\n")
  
  invisible(tab)
}

# ---------- 7. 绘图函数 ----------
plot_one_variable <- function(data, var_name, y_label, file_name) {
  combo_list <- split(data, list(data$Site, data$Niche))
  letters_all <- lapply(names(combo_list), function(combo_name) {
    sub <- combo_list[[combo_name]]
    if (nrow(sub) < 2) return(NULL)
    res <- get_letters(sub, var_name)
    res$Site  <- sub$Site[1]
    res$Niche <- sub$Niche[1]
    res
  })
  letters_all <- bind_rows(letters_all)
  
  ymax_data <- data %>%
    group_by(Site, Niche, Treatment) %>%
    summarise(y_max = max(.data[[var_name]], na.rm = TRUE), .groups = "drop")
  
  sig_data <- ymax_data %>%
    left_join(letters_all, by = c("Site", "Niche", "Treatment")) %>%
    mutate(y_label = y_max * 1.05)
  
  plot_df <- data %>%
    dplyr::select(Site, Niche, Treatment, Value = all_of(var_name)) %>%
    mutate(Site_Niche = paste(Site, Niche, sep = "_"))
  
  plot_df$Site_Niche <- factor(
    plot_df$Site_Niche,
    levels = c("HY_Topsoil", "HY_Rhizosphere", "HY_Endosphere",
               "YC_Topsoil", "YC_Rhizosphere", "YC_Endosphere")
  )
  sig_data$Site_Niche <- factor(
    paste(sig_data$Site, sig_data$Niche, sep = "_"),
    levels = levels(plot_df$Site_Niche)
  )
  
  p <- ggplot(plot_df, aes(x = Treatment, y = Value)) +
    geom_boxplot(aes(color = Treatment), outlier.shape = NA, alpha = 0.7) +
    geom_jitter(aes(color = Treatment), width = 0.15, size = 1.5, alpha = 0.6) +
    facet_wrap(~ Site_Niche, scales = "free_y", ncol = 6) +
    geom_text(data = sig_data,
              aes(x = Treatment, y = y_label, label = letters),
              size = 3.5, vjust = 0, inherit.aes = FALSE) +
    scale_color_manual(values = colors_no) +
    labs(title = var_name, x = NULL, y = y_label) +
    theme_classic() +
    theme(
      axis.line        = NULL,
      text             = element_text(family = "sans"),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.title       = element_text(size = 14, hjust = 0.5, color = "black",
                                      face = "plain",
                                      margin = ggplot2::margin(b = 10)),
      axis.text.x      = element_text(angle = 0, vjust = 0.5,
                                      face = "plain", color = "black"),
      plot.margin      = ggplot2::margin(10, 10, 10, 10, unit = "pt"),
      legend.position  = "none"
    )
  
  ggsave(file.path(out_dir, file_name), p, width = 9, height = 2.2, dpi = 300)
  return(p)
}

# ---------- 8. 绘图 ----------
p_glom  <- plot_one_variable(plot_data, "Glomeromycota",
                             "Relative abundance of Glomeromycota (%)",
                             "Glomeromycota_abundance_boxplot.pdf")

p_rhizo <- plot_one_variable(plot_data, "Rhizophagus",
                             "Relative abundance of Rhizophagus (%)",
                             "Rhizophagus_abundance_boxplot.pdf")

p_str   <- plot_one_variable(plot_data, "Streptomyces",
                             "Relative abundance of Streptomyces (%)",
                             "Streptomyces_abundance_boxplot.pdf")

p_act   <- plot_one_variable(plot_data, "Actinomycetota",
                             "Relative abundance of Actinomycetota (%)",
                             "Actinomycetota_abundance_boxplot.pdf")

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

# ---------- 10. 保存 CSV（方便直接进补充材料） ----------
write.csv(sig_glom,  file.path(out_dir, "sig_Glomeromycota.csv"),  row.names = FALSE)
write.csv(sig_rhizo, file.path(out_dir, "sig_Rhizophagus.csv"),    row.names = FALSE)
write.csv(sig_str,   file.path(out_dir, "sig_Streptomyces.csv"),   row.names = FALSE)
write.csv(sig_act,   file.path(out_dir, "sig_Actinomycetota.csv"), row.names = FALSE)

write.csv(det_glom,  file.path(out_dir, "detect_Glomeromycota.csv"),  row.names = FALSE)
write.csv(det_rhizo, file.path(out_dir, "detect_Rhizophagus.csv"),    row.names = FALSE)
write.csv(det_str,   file.path(out_dir, "detect_Streptomyces.csv"),   row.names = FALSE)
write.csv(det_act,   file.path(out_dir, "detect_Actinomycetota.csv"), row.names = FALSE)

cat("\n四张箱线图 + 显著性/检出率表格已保存至：", out_dir, "\n")


# ---------- 11. 输出 Table S6（含 SE） ----------
make_table_S6 <- function(plot_data, vars) {
  out <- list()
  for (v in vars) {
    combo_list <- split(plot_data, list(plot_data$Site, plot_data$Niche))
    res <- bind_rows(lapply(names(combo_list), function(nm) {
      sub <- combo_list[[nm]]
      if (nrow(sub) < 2) return(NULL)
      let <- get_letters(sub, v)
      stats <- sub %>%
        group_by(Treatment) %>%
        summarise(
          Mean = mean(.data[[v]], na.rm = TRUE),
          SE   = sd(.data[[v]], na.rm = TRUE) / sqrt(n()),
          n    = n(),
          n_detected = sum(.data[[v]] > 0, na.rm = TRUE),
          .groups = "drop"
        )
      stats <- stats %>%
        left_join(let %>% select(Treatment, letters),
                  by = "Treatment") %>%
        mutate(Site  = sub$Site[1],
               Niche = sub$Niche[1],
               Taxon = v) %>%
        select(Site, Niche, Taxon, Treatment, Mean, SE,
               n, n_detected, letters)
      stats
    }))
    out[[v]] <- res
  }
  bind_rows(out)
}

table_S6 <- make_table_S6(plot_data,
                          c("Glomeromycota", "Rhizophagus", "Streptomyces", "Actinomycetota"))

table_S6 <- table_S6 %>%
  mutate(across(c(Mean, SE), ~ round(.x, 5))) %>%
  rename(`n_detected / n` = n_detected, `Tukey letter` = letters) %>%
  mutate(`n_detected / n` = paste0(`n_detected / n`, " / ", n)) %>%
  select(-n)

write.csv(table_S6,
          file.path(out_dir, "Table_S6_amplicon_colonization.csv"),
          row.names = FALSE)
