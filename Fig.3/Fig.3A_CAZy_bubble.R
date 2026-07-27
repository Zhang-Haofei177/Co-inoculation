# ============================================================================
# Fig.3A: CAZy GSEA bubble plot (HY vs YC)
# ============================================================================

library(limma)
library(edgeR)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(tidyr)

out_root <- "C:/Users/ADMIN/Desktop/Co-inoculation/Fig.3"
script_dir <- out_root
dir.create(script_dir, showWarnings = FALSE, recursive = TRUE)
setwd(script_dir)

# ---------- 1. Read group info ----------
group <- read.csv(file.path(out_root, "data/group.csv"), row.names = 1)
group <- group[group$Niche == "Rhizosphere" & group$Treatment %in% c("Ctrl","AMF","Str","AMF+Str"), ]
group$Treatment <- factor(
  ifelse(group$Treatment == "AMF+Str", "AMF.Str", as.character(group$Treatment)),
  levels = c("AMF", "Str", "AMF.Str")
)
group$Site <- factor(group$Site, levels = c("HY", "YC"))

# ---------- 2. Read function (remove Others) ----------
read_absolute_table <- function(file_path) {
  df <- read.table(file_path, header = TRUE, sep = "\t", row.names = 1,
                   check.names = FALSE, quote = "", comment.char = "",
                   stringsAsFactors = FALSE)
  last_col <- colnames(df)[ncol(df)]
  if (tolower(last_col) == "description") {
    df <- df[, -ncol(df), drop = FALSE]
  }
  df <- df[rownames(df) != "Others", , drop = FALSE]
  abund <- as.data.frame(t(df))
  abund[] <- lapply(abund, as.numeric)
  return(abund)
}

# ---------- 3. Read and filter ----------
ko_abund <- read_absolute_table(file.path(out_root, "data/Unigenes.absolute.ko.xls"))
og_abund <- read_absolute_table(file.path(out_root, "data/Unigenes.absolute.og.xls"))
ec_abund <- read_absolute_table(file.path(out_root, "data/Unigenes.absolute.ec.xls"))

common_samples <- Reduce(intersect, list(rownames(group), rownames(ko_abund),
                                         rownames(og_abund), rownames(ec_abund)))
group <- group[common_samples, ]
ko_abund <- ko_abund[common_samples, ]
og_abund <- og_abund[common_samples, ]
ec_abund <- ec_abund[common_samples, ]

filter_features <- function(mat) {
  keep <- colSums(mat) >= 10 & colSums(mat > 0) >= 2
  mat[, keep, drop = FALSE]
}
ko_abund <- filter_features(ko_abund)
og_abund <- filter_features(og_abund)
ec_abund <- filter_features(ec_abund)

# ---------- 4. Build term2gene ----------
# KEGG
ko_anno <- read.table(file.path(out_root, "data/Unigenes.KEGG.tax.xls"),
                      header=TRUE, sep="\t", quote="", comment.char="")
ko_term2gene <- ko_anno %>%
  select(KO_ID, Map_ID) %>%
  distinct() %>%
  filter(KO_ID %in% colnames(ko_abund)) %>%
  rename(gene = KO_ID, term = Map_ID) %>%
  as.data.frame()

# CAZy: match families via Description
cazy_anno <- read.table(file.path(out_root, "data/Unigenes.CAZy.tax.xls"),
                        header=TRUE, sep="\t", quote="", comment.char="")
cazy_desc_to_family <- cazy_anno %>%
  select(Description, CAZy_Family) %>%
  distinct() %>%
  group_by(Description) %>%
  slice(1) %>%
  ungroup()

cazy_term2gene <- data.frame(gene = colnames(ec_abund), stringsAsFactors = FALSE) %>%
  left_join(cazy_desc_to_family, by = c("gene" = "Description")) %>%
  filter(!is.na(CAZy_Family)) %>%
  rename(term = CAZy_Family) %>%
  select(term, gene)

cat("CAZy term2gene rows:", nrow(cazy_term2gene), "\n")

# ---------- 5. Differential analysis function (combined contrast) ----------
get_logFC_combined <- function(abund, meta, site) {
  idx <- meta$Site == site & meta$Treatment %in% c("AMF", "Str", "AMF.Str")
  sub_meta <- meta[idx, , drop = FALSE]
  sub_abund <- abund[rownames(sub_meta), , drop = FALSE]
  if (!all(c("AMF", "Str", "AMF.Str") %in% sub_meta$Treatment)) return(NULL)

  count_mat <- t(sub_abund)
  mode(count_mat) <- "integer"

  dge <- DGEList(counts = count_mat, group = sub_meta$Treatment)
  keep <- filterByExpr(dge, min.count = 2, min.total.count = 10)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)

  design <- model.matrix(~ 0 + Treatment, data = sub_meta)
  colnames(design) <- make.names(colnames(design))
  if (!all(c("TreatmentAMF.Str", "TreatmentAMF", "TreatmentStr") %in% colnames(design))) return(NULL)

  v <- voom(dge, design, plot = FALSE)
  fit <- lmFit(v, design)

  # Combined contrast: AMF+Str vs mean of two single inoculants
  cont_matrix <- makeContrasts(
    Comb_vs_AvgSingle = TreatmentAMF.Str - (TreatmentAMF + TreatmentStr)/2,
    levels = design
  )
  fit2 <- contrasts.fit(fit, cont_matrix)
  fit2 <- eBayes(fit2)

  res <- topTable(fit2, coef = "Comb_vs_AvgSingle", number = Inf, adjust.method = "none", sort.by = "none")
  logFC <- res$logFC
  names(logFC) <- rownames(res)
  return(sort(logFC, decreasing = TRUE))
}

# ---------- 6. Calculate logFC for all contrasts (one list per site) ----------
# Replace original contrast_plans with simple site names
sites <- c("HY", "YC")

compute_logFC_list_combined <- function(abund, db_name) {
  logFC_list <- list()
  for (site in sites) {
    cat(db_name, "-", site, "\n")
    lfc <- get_logFC_combined(abund, group, site)
    if (!is.null(lfc)) logFC_list[[site]] <- lfc
  }
  return(logFC_list)
}

ko_logFCs <- compute_logFC_list_combined(ko_abund, "KEGG")
ec_logFCs <- compute_logFC_list_combined(ec_abund, "CAZy")

# ---------- 7. GSEA ----------
run_gsea <- function(logFC, term2gene, db_name, contrast_name) {
  if (is.null(logFC) || length(logFC) == 0 || nrow(term2gene) == 0) return(NULL)
  set.seed(123)
  gsea_res <- tryCatch({
    GSEA(
      geneList = logFC,
      TERM2GENE = term2gene,
      minGSSize = 5,
      maxGSSize = 500,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      eps = 0
    )
  }, error = function(e) NULL)

  if (!is.null(gsea_res) && nrow(gsea_res@result) > 0) {
    res_df <- as.data.frame(gsea_res@result)
    res_df$Database <- db_name
    res_df$Contrast <- contrast_name
    return(res_df)
  }
  return(NULL)
}

gsea_tables <- list()
for (contr in names(ko_logFCs)) {
  res <- run_gsea(ko_logFCs[[contr]], ko_term2gene, "KEGG", contr)
  if (!is.null(res)) gsea_tables[[length(gsea_tables)+1]] <- res
}
if (nrow(cazy_term2gene) > 0) {
  for (contr in names(ec_logFCs)) {
    res <- run_gsea(ec_logFCs[[contr]], cazy_term2gene, "CAZy", contr)
    if (!is.null(res)) gsea_tables[[length(gsea_tables)+1]] <- res
  }
}

gsea_all <- bind_rows(gsea_tables)
cat("Total significant enriched pathways:", nrow(gsea_all), "\n")
write.csv(gsea_all, file.path(script_dir, "GSEA_all_significant_pathways.csv"), row.names = FALSE)
############# Plotting ##########################
# ============================================================================
# CAZy GSEA visualization (adapted for combined contrast HY vs YC)
# Fig.3A
# ============================================================================
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Ensure gsea_all exists (CAZy results only)
if (!exists("gsea_all") || nrow(gsea_all) == 0) {
  stop("gsea_all is empty, please run GSEA analysis first.")
}

cazy_res <- gsea_all[gsea_all$Database == "CAZy", ]
if (nrow(cazy_res) == 0) stop("No CAZy results.")

# ---------- 1. Set contrast order ----------
cazy_res$Contrast <- factor(cazy_res$Contrast, levels = c("HY", "YC"))

# ---------- 2. Dynamically build function annotation table ----------
# ---------- Concise function annotation table (14 families) ----------
family_func <- data.frame(
  ID = c("GH94", "GH130", "GT4", "GT2", "GT1", "GH32", "GT31", "GT14", "GT7", "AA3", "AA5", "GH13", "GT61", "GH18"),
  Function = c(
    "Cellobiose/cellodextrin phosphorylase (cellulose utilization)",
    "β-mannoside phosphorylase (fungal cell wall utilization)",
    "Sucrose/trehalose/EPS synthase (osmoprotectant & EPS synthesis)",
    "Cellulose/EPS/β-glucan synthase (EPS & biofilm formation)",
    "UDP-glycosyltransferase (phytohormone glycosylation)",
    "Fructan/inulin hydrolase (fructan degradation, downregulated)",
    "β1,3-galactosyltransferase (AGP glycosylation, cell wall remodeling)",
    "β-1,6-GlcNAc transferase (glycoprotein branching, stress response)",
    "Lactose/LacNAc synthase (β-1,4-galactosyltransferase, downregulated)",
    "Sugar/alcohol oxidase (H₂O₂ production, redox signaling)",
    "Glyoxal/alcohol oxidase (ROS scavenging)",
    "α-amylase (starch/glycogen degradation, downregulated)",
    "Xylan/arabinoxylan sidechain synthase (cell wall reinforcement)",
    "Chitinase (chitin degradation, fungal wall remodeling)"
  ),
  stringsAsFactors = FALSE
)

# Merge into cazy_res (remove Function column first if exists)
cazy_res <- cazy_res[, !colnames(cazy_res) %in% "Function"]
cazy_res <- left_join(cazy_res, family_func, by = "ID")

# ---------- Family function category mapping table ----------
# 14 significant families, each assigned a drought-related function category
# ---------- Family function category mapping table ----------
family_category <- data.frame(
  ID = c("GH94", "GH130", "GT4", "GT2", "GT1", "GH32", "GT31", "GT14", "GT7", "AA3", "AA5", "GH13", "GT61", "GH18"),
  Category = c(
    "Carbon source utilization",        # GH94 cellodextrin phosphorylase
    "Carbon source utilization",        # GH130 beta-mannoside phosphorylase
    "EPS synthesis & osmoprotectant",   # GT4 sucrose/trehalose/EPS synthase
    "EPS synthesis & osmoprotectant",   # GT2 cellulose/EPS/beta-glucan synthase
    "Phytohormone homeostasis",         # GT1 UDP-glycosyltransferase
    "Carbon source utilization",        # GH32 fructan/inulin hydrolase (downregulated)
    "Cell wall remodeling & signaling", # GT31 beta-1,3-galactosyltransferase
    "Cell wall remodeling & signaling", # GT14 beta-1,6-GlcNAc transferase
    "Cell wall remodeling & signaling", # GT7 lactose/LacNAc synthase
    "ROS scavenging",                   # AA3 sugar/alcohol oxidase
    "ROS scavenging",                   # AA5 glyoxal/alcohol oxidase
    "Carbon source utilization",        # GH13 alpha-amylase (downregulated)
    "Cell wall remodeling & signaling", # GT61 xylan/arabinoxylan sidechain synthase
    "Carbon source utilization"         # GH18 chitinase
  ),
  stringsAsFactors = FALSE
)

# Merge into cazy_res
cazy_res <- cazy_res %>%
  left_join(family_category, by = "ID")

# ---------- Category color definitions ----------
cat_colors <- c(
  "EPS synthesis & osmoprotectant"      = "#74C476",
  "Carbon source utilization"           = "#6BAED6",
  "ROS scavenging"                      = "#FB6A4A",
  "Phytohormone homeostasis"            = "#9E9AC8",
  "Cell wall remodeling & signaling"    = "#FDBF6F"
)


# ---------- 4. Sort rows by function category and YC NES ----------
# Specify category order
category_order <- c("EPS synthesis & osmoprotectant",
                    "Carbon source utilization",
                    "ROS scavenging",
                    "Phytohormone homeostasis",
                    "Cell wall remodeling & signaling")

# Calculate NES in YC for each ID (0 if missing)
yc_nes <- cazy_res %>%
  filter(Contrast == "YC") %>%
  select(ID, NES_YC = NES)

cazy_res <- cazy_res %>%
  left_join(yc_nes, by = "ID") %>%
  mutate(NES_YC = ifelse(is.na(NES_YC), 0, NES_YC))

# Order by category, then by YC NES descending
id_order <- cazy_res %>%
  distinct(ID, Category, NES_YC) %>%
  mutate(Category = factor(Category, levels = category_order)) %>%
  arrange(Category, desc(NES_YC)) %>%
  pull(ID)

cazy_res$ID <- factor(cazy_res$ID, levels = id_order)

# ---------- 5. Left color tile ----------
p_color <- ggplot(cazy_res, aes(x = 1, y = ID, fill = Category)) +
  geom_tile(width = 0.9, color = NA) +
  scale_fill_manual(values = cat_colors, name = "Category") +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())

# ---------- 6. Middle function description ----------
p_text <- ggplot(cazy_res, aes(x = 0, y = ID, label = Function)) +
  geom_text(hjust = 0, size = 2.8, color = "black") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())

# ---------- 7. Bubble plot ----------
p_dot <- ggplot(cazy_res, aes(x = Contrast, y = ID)) +
  geom_point(aes(size = -log10(p.adjust), color = NES)) +
  scale_color_gradient2(low = "#2166AC", mid = "#F5F5F5", high = "#B2182B", midpoint = 0, name = "NES") +
  scale_size_continuous(name = "-log10(p.adjust)", range = c(2, 7)) +
  labs(title = "CAZy families enriched by dual vs. single inoculum (HY vs YC)",
       x = NULL, y = NULL) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0, face = "plain", size = 12),
        axis.text.y = element_text(size = 10),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(0, 0, 0, 0))

# ---------- 8. Combine ----------
combined <- p_color + p_text + p_dot +
  plot_layout(widths = c(0.08, 1.4,0.38), guides = "collect") &
  theme(plot.margin = margin(2, 2, 2, 2))
combined
ggsave(file.path(script_dir, "Fig.3A_CAZy_dotplot.pdf"), combined, width = 7.5, height = max(4, length(unique(cazy_res$ID)) * 0.25))
ggsave(file.path(script_dir, "Fig.3A_CAZy_dotplot.png"), combined, width = 7.5, height = max(4, length(unique(cazy_res$ID)) * 0.25), dpi = 300)
