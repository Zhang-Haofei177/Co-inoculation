# ============================================================================
# β多样性分析：PERMANOVA + 6组Treatment效应宽表
# 仅保留处理：AMF, Str, AMF+Str
# 仅保留生态位：Topsoil, Rhizosphere, Endosphere
# 输出：每个 Site×Niche 组合中 Treatment 的 PERMANOVA R² 和 p 值宽表
# ============================================================================
rm(list = ls())

# 1. 加载包
library(vegan)
library(dplyr)
library(tidyr)

# 2. 设置路径和输出目录
out_root <- "C:/Users/ADMIN/Desktop/Co-inoculation/Fig.2"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
# 3. 定义分析函数（单个数据集）
analyze_beta <- function(otu_file, prefix) {
  cat("\n========== 正在分析", prefix, "==========\n")
  
  # 读取抽平后的绝对丰度表（行为ASV，列为样本）
  otu <- read.csv(otu_file, row.names = 1)
  
  # 过滤ASV：至少在2个样本中出现且总丰度≥10
  asv_presence <- rowSums(otu > 0)
  asv_total <- rowSums(otu)
  keep_asv <- asv_presence >= 2 & asv_total >= 10
  otu_filtered <- otu[keep_asv, ]
  
  # 转置为样本×ASV
  otu_t <- t(otu_filtered)
  
  # 读取分组信息
  group <- read.csv("data/group.csv", row.names = 1)
  
  # ---------- 关键过滤：移除 Ctrl，移除 Subsoil/Bottom soil ----------
  group <- group[group$Treatment != "Ctrl", ]
  group <- group[group$Niche %in% c("Topsoil", "Rhizosphere", "Endosphere"), ]
  
  # 重新设置因子水平
  treatment_order <- c("AMF", "Str", "AMF+Str")
  group$Treatment <- factor(group$Treatment, levels = treatment_order)
  
  niche_order <- c("Topsoil", "Rhizosphere", "Endosphere")
  group$Niche <- factor(group$Niche, levels = niche_order)
  
  group$Site <- factor(group$Site, levels = c("HY", "YC"))
  
  # 确保样本顺序一致
  common_samples <- intersect(rownames(group), rownames(otu_t))
  otu_t <- otu_t[common_samples, ]
  group <- group[common_samples, ]
  
  # 计算 Bray-Curtis 距离
  bray_dist <- vegdist(otu_t, method = "bray")
  
  # -------------------- 整体 PERMANOVA（可选，保留） --------------------
  set.seed(123)
  perm_full <- adonis2(bray_dist ~ Site * Niche * Treatment, data = group, permutations = 999, by = "terms")
  perm_full <- as.data.frame(perm_full)
  perm_full$Factor <- rownames(perm_full)
  write.csv(perm_full, file.path(out_root, paste0(prefix, "_PERMANOVA_full.csv")), row.names = FALSE)
  
  # -------------------- 核心：2地点 × 3生态位 = 6 组，分别计算 Treatment 的 PERMANOVA --------------------
  comb <- expand.grid(Site = levels(group$Site), Niche = levels(group$Niche), stringsAsFactors = FALSE)
  results_long <- data.frame()
  
  for (i in 1:nrow(comb)) {
    sit <- comb$Site[i]
    nic <- comb$Niche[i]
    idx <- group$Site == sit & group$Niche == nic
    if (sum(idx) < 4) next   # 样本太少则跳过
    
    sub_otu <- otu_t[idx, , drop = FALSE]
    sub_group <- group[idx, , drop = FALSE]
    
    # 至少需要两个处理水平
    if (length(unique(sub_group$Treatment)) < 2) next
    
    sub_dist <- vegdist(sub_otu, method = "bray")
    
    set.seed(123)
    res <- adonis2(sub_dist ~ Treatment, data = sub_group, permutations = 999)
    r2 <- res$R2[1]
    pval <- res$`Pr(>F)`[1]
    
    # 显著性星号
    signif <- ""
    if (!is.na(pval)) {
      if (pval < 0.001) signif <- "***"
      else if (pval < 0.01) signif <- "**"
      else if (pval < 0.05) signif <- "*"
    }
    
    results_long <- rbind(results_long, data.frame(
      Site = sit, Niche = nic, R2 = round(r2, 4), P = round(pval, 4), Signif = signif,
      stringsAsFactors = FALSE
    ))
  }
  
  # 保存长表
  write.csv(results_long, file.path(out_root, paste0(prefix, "_Treatment_PERMANOVA_long.csv")), row.names = FALSE)
  
  # 整理成宽表：行=Site，列=Niche，单元格为 "R²=..., p=...*"
  if (nrow(results_long) > 0) {
    results_long$Value <- sprintf("R²=%.4f, p=%.4f%s", results_long$R2, results_long$P, results_long$Signif)
    results_wide <- results_long %>%
      dplyr::select(Site, Niche, Value) %>%
      tidyr::pivot_wider(names_from = Niche, values_from = Value)
    
    write.csv(results_wide, file.path(out_root, paste0(prefix, "_Treatment_PERMANOVA_wide.csv")), row.names = FALSE)
    
    cat("\n", prefix, " Treatment PERMANOVA 宽表：\n", sep = "")
    print(as.data.frame(results_wide), row.names = FALSE)
  } else {
    cat("\n", prefix, " 没有足够的样本进行 Treatment PERMANOVA 分析。\n", sep = "")
  }
  
  # -------------------- 其他分析（可选，按需保留） --------------------
  # 分层PERMANOVA（按生态位）
  niches <- levels(group$Niche)
  perm_niche <- list()
  for (nic in niches) {
    idx <- group$Niche == nic
    sub_otu <- otu_t[idx, , drop = FALSE]
    sub_group <- group[idx, , drop = FALSE]
    if (nrow(sub_otu) >= 3 & length(unique(sub_group$Treatment)) >= 2) {
      sub_dist <- vegdist(sub_otu, method = "bray")
      res <- adonis2(sub_dist ~ Site * Treatment, data = sub_group, permutations = 999, by = "terms")
      res <- as.data.frame(res)
      res$Niche <- nic
      perm_niche[[nic]] <- res
    }
  }
  perm_niche_df <- bind_rows(perm_niche)
  write.csv(perm_niche_df, file.path(out_root, paste0(prefix, "_PERMANOVA_by_niche.csv")), row.names = FALSE)
  
  # 分层PERMANOVA（按地点）
  sites <- levels(group$Site)
  perm_site <- list()
  for (st in sites) {
    idx <- group$Site == st
    sub_otu <- otu_t[idx, , drop = FALSE]
    sub_group <- group[idx, , drop = FALSE]
    if (nrow(sub_otu) >= 3 & length(unique(sub_group$Treatment)) >= 2) {
      sub_dist <- vegdist(sub_otu, method = "bray")
      res <- adonis2(sub_dist ~ Niche * Treatment, data = sub_group, permutations = 999, by = "terms")
      res <- as.data.frame(res)
      res$Site <- st
      perm_site[[st]] <- res
    }
  }
  perm_site_df <- bind_rows(perm_site)
  write.csv(perm_site_df, file.path(out_root, paste0(prefix, "_PERMANOVA_by_site.csv")), row.names = FALSE)
  
  # PERMDISP（按处理、按生态位）
  beta_disp_treat <- betadisper(bray_dist, group$Treatment)
  sink(file.path(out_root, paste0(prefix, "_PERMDISP_treatment.txt")))
  print(permutest(beta_disp_treat, permutations = 999))
  sink()
  
  beta_disp_niche <- betadisper(bray_dist, group$Niche)
  sink(file.path(out_root, paste0(prefix, "_PERMDISP_niche.txt")))
  print(permutest(beta_disp_niche, permutations = 999))
  sink()
  
  cat(prefix, "分析完成！结果保存至：", out_root, "\n")
  return(results_wide)
}

p_bac <- analyze_beta("data/BAC_rarefied.csv", "Bacteria")
p_fun <- analyze_beta("data/FUN_rarefied.csv", "Fungi")
cat("\n所有分析完成！Table S4 Fungi_Treatment_PERMANOVA_long.csv AND Bacteria_Treatment_PERMANOVA_long.csv\n")

