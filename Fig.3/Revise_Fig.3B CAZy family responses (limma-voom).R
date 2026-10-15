# ============================================================================
# CAZyme limma-voom: export raw tables + rounded supplementary tables
#   Raw outputs (GitHub):
#     CAZyme_limma_all_contrasts.csv
#     CAZyme_limma_focal_families.csv
#     CAZyme_limma_separate_contrasts.csv
#   Supplementary tables (manuscript):
#     Table_S2_CAZyme_focal_families.csv       (AMF+Str vs pooled single)
#     Table_S3_CAZyme_separate_contrasts.csv   (AMF+Str vs AMF; AMF+Str vs Str)
#   Rounding rule:
#     P.Value, adj.P.Val: signif(x, 5)
#     all other numeric:   round(x, 5)
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(limma)
  library(edgeR)
})

# ---------------------------------------------------------------------------
# 0. Paths
# ---------------------------------------------------------------------------
in_dir  <- "C:/Users/ADMIN/Desktop/Co-inoculation/Fig.3/data"
out_dir <- "C:/Users/ADMIN/Desktop/Co-inoculation/Fig.3"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# 1. Metadata
# ---------------------------------------------------------------------------
group <- read.csv(file.path(in_dir, "group.csv"), row.names = 1)
group_rhizo <- group[group$Niche == "Rhizosphere" &
                       group$Treatment %in% c("AMF", "Str", "AMF+Str"), ]
group_rhizo$Treatment <- factor(group_rhizo$Treatment,
                                levels = c("AMF", "Str", "AMF+Str"))
group_rhizo$Site <- factor(group_rhizo$Site, levels = c("HY", "YC"))

group_rhizo$Treatment_safe <- ifelse(group_rhizo$Treatment == "AMF+Str",
                                     "AMF_Str",
                                     as.character(group_rhizo$Treatment))
group_rhizo$Group <- factor(
  paste(group_rhizo$Site, group_rhizo$Treatment_safe, sep = "_"),
  levels = c("HY_AMF", "HY_Str", "HY_AMF_Str",
             "YC_AMF", "YC_Str", "YC_AMF_Str"))

rhizo_samples <- rownames(group_rhizo)

# ---------------------------------------------------------------------------
# 2. CAZy family abundance (level2)
# ---------------------------------------------------------------------------
cazy_level2 <- fread(file.path(in_dir, "Unigenes.absolute.level2.xls"),
                     sep = "\t", header = TRUE, data.table = FALSE)
rn <- cazy_level2[, 1]
cazy_level2 <- cazy_level2[, -1, drop = FALSE]
rownames(cazy_level2) <- rn

cazy_abund <- as.data.frame(t(cazy_level2))
cazy_abund[] <- lapply(cazy_abund, function(x) as.numeric(as.character(x)))
cazy_abund$Sample <- rownames(cazy_abund)

common_samples <- intersect(rhizo_samples, cazy_abund$Sample)
cazy_abund <- cazy_abund[cazy_abund$Sample %in% common_samples, ]
cazy_abund <- cazy_abund[match(common_samples, cazy_abund$Sample), ]
group_rhizo <- group_rhizo[common_samples, ]

count_mat <- t(as.matrix(cazy_abund[, setdiff(colnames(cazy_abund), "Sample")]))
count_mat <- round(count_mat)
mode(count_mat) <- "integer"

keep <- rowSums(count_mat) >= 10 & rowSums(count_mat > 0) >= 2
count_mat <- count_mat[keep, , drop = FALSE]

cat(sprintf("CAZy families retained after filtering: %d\n", nrow(count_mat)))

# ---------------------------------------------------------------------------
# 3. limma-voom
# ---------------------------------------------------------------------------
dge <- DGEList(counts = count_mat, group = group_rhizo$Group)
dge <- calcNormFactors(dge)

design <- model.matrix(~ 0 + Group, data = group_rhizo)
colnames(design) <- make.names(colnames(design))

v <- voom(dge, design, plot = FALSE)
fit <- lmFit(v, design)

cont <- makeContrasts(
  HY_Dual_vs_Single = GroupHY_AMF_Str - (GroupHY_AMF + GroupHY_Str)/2,
  YC_Dual_vs_Single = GroupYC_AMF_Str - (GroupYC_AMF + GroupYC_Str)/2,
  HY_Dual_vs_AMF    = GroupHY_AMF_Str - GroupHY_AMF,
  HY_Dual_vs_Str    = GroupHY_AMF_Str - GroupHY_Str,
  YC_Dual_vs_AMF    = GroupYC_AMF_Str - GroupYC_AMF,
  YC_Dual_vs_Str    = GroupYC_AMF_Str - GroupYC_Str,
  Diff_YC_vs_HY     = (GroupYC_AMF_Str - (GroupYC_AMF + GroupYC_Str)/2) -
    (GroupHY_AMF_Str - (GroupHY_AMF + GroupHY_Str)/2),
  levels = design
)

fit2 <- contrasts.fit(fit, cont)
fit2 <- eBayes(fit2)

# ---------------------------------------------------------------------------
# 4. Extract all results
# ---------------------------------------------------------------------------
all_res <- data.frame()
for (coef in colnames(cont)) {
  res <- topTable(fit2, coef = coef, number = Inf, sort.by = "none")
  res$Contrast <- coef
  res$CAZy_Family <- rownames(res)
  all_res <- rbind(all_res, res)
}

all_res <- all_res %>%
  select(CAZy_Family, Contrast, logFC, AveExpr, t, P.Value, adj.P.Val, B) %>%
  mutate(
    SE      = logFC / t,
    CI_low  = logFC - 1.96 * SE,
    CI_high = logFC + 1.96 * SE
  ) %>%
  select(CAZy_Family, Contrast, logFC, SE, CI_low, CI_high,
         AveExpr, t, P.Value, adj.P.Val, B)

# ---------------------------------------------------------------------------
# 5. Focal families (14) and separate contrasts subsets
# ---------------------------------------------------------------------------
sig_families <- c("GH94", "GH130", "GT4", "GT2", "GT1", "GH32", "GT31",
                  "GT14", "GT7", "AA3", "AA5", "GH13", "GT61", "GH18")

focal_res <- all_res %>% filter(CAZy_Family %in% sig_families)

sep_res <- focal_res %>%
  filter(Contrast %in% c("HY_Dual_vs_AMF", "HY_Dual_vs_Str",
                         "YC_Dual_vs_AMF", "YC_Dual_vs_Str")) %>%
  mutate(
    Site      = ifelse(grepl("^HY", Contrast), "HY", "YC"),
    Reference = ifelse(grepl("vs_AMF$", Contrast), "AMF", "Str")
  ) %>%
  select(CAZy_Family, Site, Reference, logFC, SE, CI_low, CI_high,
         P.Value, adj.P.Val) %>%
  arrange(CAZy_Family, Site, Reference)

# ---------------------------------------------------------------------------
# 6. Export RAW tables (full precision) -> GitHub
# ---------------------------------------------------------------------------
write.csv(all_res,
          file.path(out_dir, "CAZyme_limma_all_contrasts.csv"),
          row.names = FALSE)

write.csv(focal_res,
          file.path(out_dir, "CAZyme_limma_focal_families.csv"),
          row.names = FALSE)

write.csv(sep_res,
          file.path(out_dir, "CAZyme_limma_separate_contrasts.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# 7. Rounding helper
#   P.Value / adj.P.Val: signif(x, 5)
#   all other numeric:   round(x, 5)
# ---------------------------------------------------------------------------
round_table <- function(df) {
  df %>%
    mutate(across(where(is.numeric), ~ round(.x, 5))) %>%
    mutate(
      P.Value   = signif(P.Value, 5),
      adj.P.Val = signif(adj.P.Val, 5)
    )
}

# ---------------------------------------------------------------------------
# 8. Table S2: focal families (main contrast) - rounded
#    Keep only the two main contrasts for the supplementary table
# ---------------------------------------------------------------------------
table_S2 <- focal_res %>%
  filter(Contrast %in% c("HY_Dual_vs_Single", "YC_Dual_vs_Single")) %>%
  mutate(
    Site = ifelse(grepl("^HY", Contrast), "HY", "YC"),
    Contrast = "AMF+Str_vs_pooled_single"
  ) %>%
  select(CAZy_Family, Site, Contrast, logFC, SE, CI_low, CI_high,
         P.Value, adj.P.Val) %>%
  arrange(CAZy_Family, Site)

table_S2 <- round_table(table_S2)

write.csv(table_S2,
          file.path(out_dir, "Table_S2_CAZyme_focal_families.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# 9. Table S3: separate contrasts - rounded
# ---------------------------------------------------------------------------
table_S3 <- round_table(sep_res)

write.csv(table_S3,
          file.path(out_dir, "Table_S3_CAZyme_separate_contrasts.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# 10. Summary
# ---------------------------------------------------------------------------
cat("\n===== Export summary =====\n")
cat("RAW tables (GitHub):\n")
cat("  CAZyme_limma_all_contrasts.csv\n")
cat("  CAZyme_limma_focal_families.csv\n")
cat("  CAZyme_limma_separate_contrasts.csv\n\n")
cat("Supplementary tables (rounded):\n")
cat("  Table_S2_CAZyme_focal_families.csv\n")
cat("  Table_S3_CAZyme_separate_contrasts.csv\n\n")

cat("Table S2 preview:\n"); print(head(table_S2, 10))
cat("\nTable S3 preview:\n"); print(head(table_S3, 10))
# ============================================================================
# 11. Fig.3B visualization: heatmap (left) + dotplot (right)
#     基于 limma-voom 主对比 AMF+Str vs pooled single-inoculant
#     星号用未校正 P 值；FDR 不显著写入图注
#     固定 family 顺序，与原始脚本一致
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(pheatmap)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

# ---------- 11.1 取主对比，准备绘图数据 ----------
plot_dat <- focal_res %>%
  filter(Contrast %in% c("HY_Dual_vs_Single", "YC_Dual_vs_Single")) %>%
  mutate(
    Site  = ifelse(grepl("^HY", Contrast), "HY", "YC"),
    Stars = ifelse(P.Value < 0.001, "***",
                   ifelse(P.Value < 0.01, "**",
                          ifelse(P.Value < 0.05, "*", "")))
  )

# ---------- 11.2 热图矩阵 ----------
heat_mat <- plot_dat %>%
  select(CAZy_Family, Site, logFC) %>%
  pivot_wider(names_from = Site, values_from = logFC) %>%
  column_to_rownames("CAZy_Family") %>%
  as.matrix()

star_mat <- plot_dat %>%
  select(CAZy_Family, Site, Stars) %>%
  pivot_wider(names_from = Site, values_from = Stars) %>%
  column_to_rownames("CAZy_Family") %>%
  as.matrix()

# 固定 family 顺序（与原始脚本一致）
fam_fixed <- c("GT7", "GT14", "GT61", "GT31", "GT1", "AA5", "AA3",
               "GH32", "GH13", "GH94", "GH18", "GH130", "GT2", "GT4")

# 取 14 个 focal family，缺失的 family 自动跳过
fam_keep <- intersect(fam_fixed, rownames(heat_mat))
heat_mat <- heat_mat[fam_keep, , drop = FALSE]
star_mat <- star_mat[fam_keep, , drop = FALSE]

# ---------- 11.3 热图（左） ----------
pheatmap(
  heat_mat,
  color = colorRampPalette(c("#2166AC", "#F5F5F5", "#B2182B"))(100),
  cluster_rows = FALSE, cluster_cols = FALSE,
  display_numbers = star_mat,
  number_format = "%s",
  fontsize_number = 6,
  number_color = "red4",
  fontsize = 7,
  main = "log2FC (AMF+Str vs pooled single inoculant)",
  filename = file.path(out_dir, "Revised_Fig.3B_heatmap_limma.pdf"),
  width = 1.3, height = 2
)

pheatmap(
  heat_mat,
  color = colorRampPalette(c("#2166AC", "#F5F5F5", "#B2182B"))(100),
  cluster_rows = FALSE, cluster_cols = FALSE,
  display_numbers = star_mat,
  number_format = "%s",
  fontsize_number = 6,
  number_color = "red4",
  fontsize = 7,
  main = "log2FC (AMF+Str vs pooled single inoculant)",
  filename = file.path(out_dir, "Revised_Fig.3B_heatmap_limma.png"),
  width = 1.3, height = 2,
  res = 300
)

# ---------- 11.4 点图（右） ----------
# 固定 family 顺序：点图 y 轴从下到上 = fam_fixed
plot_dat$CAZy_Family <- factor(plot_dat$CAZy_Family, levels = rev(fam_fixed))

p_dot <- ggplot(plot_dat, aes(x = logFC, y = CAZy_Family, color = Site)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high),
                orientation = "y",
                height = 0.2, linewidth = 0.3,
                position = position_dodge(width = 0.4)) +
  geom_point(aes(shape = Site), size = 1.8,
             position = position_dodge(width = 0.4)) +
  scale_color_manual(values = c("HY" = "#E41A1C", "YC" = "#377EB8")) +
  scale_shape_manual(values = c("HY" = 19, "YC" = 17)) +
  labs(
    x = expression("log"[2]*" FC (co-inoc / single)"),
    y = NULL
  ) +
  theme_bw(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    axis.text.y = element_text(size = 7),
    legend.position = "right"
  )

# 未校正 P 星号
star_ann <- plot_dat %>%
  filter(Stars != "") %>%
  group_by(CAZy_Family) %>%
  reframe(
    x_pos = max(CI_high, na.rm = TRUE) + 0.02,
    y_pos = mean(as.numeric(CAZy_Family)),
    Stars = first(Stars)
  )

if (nrow(star_ann) > 0) {
  p_dot <- p_dot + geom_text(
    data = star_ann,
    aes(x = x_pos, y = y_pos, label = Stars),
    inherit.aes = FALSE,
    size = 3, color = "black", hjust = 0
  )
}

print(p_dot)

ggsave(file.path(out_dir, "Revised_Fig.3B_dotplot_limma.pdf"),
       p_dot, width = 50, height = 60, units = "mm", device = "pdf")
ggsave(file.path(out_dir, "Revised_Fig.3B_dotplot_limma.png"),
       p_dot, width = 50, height = 60, units = "mm",
       dpi = 300, device = "png")

cat("\nFig.3B visualization finished.\n")
cat("  Fig.3B_heatmap_limma.pdf / .png\n")
cat("  Fig.3B_dotplot_limma.pdf / .png\n")

