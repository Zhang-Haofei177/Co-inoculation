# ============================================================================
# Fig.4A: YAS pathway synergy differential heatmap (Dual vs Single, HY vs YC)
# ============================================================================
library(dplyr)
library(tidyr)
library(pheatmap)

# ---------- 0. Paths ----------
out_root <- "C:/Users/ADMIN/Desktop/Co-inoculation/Fig.4"
script_dir <- out_root
dir.create(script_dir, showWarnings = FALSE, recursive = TRUE)
setwd(script_dir)

# ---------- 1. Read data ----------
kegg3_file <- file.path(out_root, "data/Unigenes.absolute.level3.xls")
kegg3_abund <- read.table(kegg3_file, header = TRUE, sep = "\t", row.names = 1,
                          check.names = FALSE, quote = "", comment.char = "")
desc_col <- which(tolower(colnames(kegg3_abund)) == "description")
sample_cols <- setdiff(colnames(kegg3_abund), "Description")
abund_mat <- as.data.frame(t(kegg3_abund[, sample_cols, drop = FALSE]))
abund_mat[] <- lapply(abund_mat, as.numeric)

group <- read.csv(file.path(out_root, "data/group.csv"), row.names = 1)
group$Treatment <- factor(group$Treatment, levels = c("Ctrl", "AMF", "Str", "AMF+Str"))
group$Site <- factor(group$Site, levels = c("HY", "YC"))

# Filter rhizosphere + non-Ctrl
rhizo_s <- rownames(group)[group$Niche == "Rhizosphere" & group$Treatment != "Ctrl"]
abund_rhizo <- abund_mat[intersect(rownames(abund_mat), rhizo_s), ]
meta <- group[rownames(abund_rhizo), ]

# Filter low-abundance pathways (total abundance >= 10 and present in >= 2 samples)
abund_rhizo <- abund_rhizo[, colSums(abund_rhizo) >= 10 & colSums(abund_rhizo > 0) >= 2]

# ---------- 2. Define YAS pathways and match ----------
Y_paths <- c("Fatty acid biosynthesis", "Carbon fixation pathways in prokaryotes",
             "Citrate cycle (TCA cycle)", "Ribosome", "Pentose phosphate pathway",
             "Glycolysis / Gluconeogenesis", "Valine, leucine and isoleucine biosynthesis",
             "Lysine biosynthesis", "Phenylalanine, tyrosine and tryptophan biosynthesis",
             "One carbon pool by folate")
A_paths <- c("ABC transporters", "Phosphotransferase system (PTS)",
             "Biosynthesis of siderophore group nonribosomal peptides", "Flagellar assembly",
             "Regulation of actin cytoskeleton", "Bacterial chemotaxis",
             "RNA degradation", "Biotin metabolism", "cAMP signaling pathway",
             "Bacterial secretion system")
S_paths <- c("Homologous recombination", "Base excision repair",
             "Nucleotide excision repair", "Lipopolysaccharide biosynthesis",
             "Peptidoglycan biosynthesis", "Ascorbate and aldarate metabolism",
             "Mismatch repair", "Porphyrin and chlorophyll metabolism",
             "Porphyrin metabolism", "Thiamine metabolism", "Proteasome")

path_desc <- data.frame(Pathway_ID = rownames(kegg3_abund),
                        Description = kegg3_abund$Description, stringsAsFactors = FALSE)

Y_lower <- tolower(trimws(Y_paths))
A_lower <- tolower(trimws(A_paths))
S_lower <- tolower(trimws(S_paths))

path_hits <- path_desc %>%
  rowwise() %>%
  mutate(
    hit_Y = any(sapply(Y_lower, function(pat) grepl(pat, tolower(Description), fixed = TRUE))),
    hit_A = any(sapply(A_lower, function(pat) grepl(pat, tolower(Description), fixed = TRUE))),
    hit_S = any(sapply(S_lower, function(pat) grepl(pat, tolower(Description), fixed = TRUE)))
  ) %>%
  ungroup() %>%
  mutate(Strategy = case_when(hit_Y ~ "Y", hit_A ~ "A", hit_S ~ "S", TRUE ~ "None"))

yas_ids <- path_hits %>% filter(Strategy != "None") %>% pull(Pathway_ID) %>% unique()

path_labels <- path_hits %>%
  filter(Pathway_ID %in% yas_ids) %>%
  mutate(Label = sapply(strsplit(Description, ";\\s*"), function(x) tail(x, 1))) %>%
  select(Pathway_ID, Strategy, Label)

# ---------- 3. Calculate log2FC and p-values (Dual vs Single, per site) ----------
sites <- c("HY", "YC")
log2fc_mat <- matrix(NA, nrow = length(yas_ids), ncol = length(sites),
                     dimnames = list(yas_ids, sites))
pval_mat <- matrix(NA, nrow = length(yas_ids), ncol = length(sites),
                   dimnames = list(yas_ids, sites))

for (pid in yas_ids) {
  vals <- abund_rhizo[, pid]
  names(vals) <- rownames(abund_rhizo)
  for (site in sites) {
    dual_idx <- meta$Site == site & meta$Treatment == "AMF+Str"
    single_idx <- meta$Site == site & meta$Treatment %in% c("AMF", "Str")
    x_dual <- vals[dual_idx]
    x_single <- vals[single_idx]
    if (length(x_dual) < 3 | length(x_single) < 3) next
    fc <- mean(x_dual) / mean(x_single)
    log2fc_mat[pid, site] <- log2(fc)
    pval_mat[pid, site] <- wilcox.test(x_dual, x_single, exact = FALSE)$p.value
  }
}

# ---------- 4. Prepare heatmap data ----------
keep <- rowSums(is.na(log2fc_mat)) < ncol(log2fc_mat)
log2fc_mat <- log2fc_mat[keep, , drop = FALSE]
pval_mat <- pval_mat[keep, , drop = FALSE]

info <- path_labels %>% filter(Pathway_ID %in% rownames(log2fc_mat))
info <- info[match(rownames(log2fc_mat), info$Pathway_ID), ]
rownames(log2fc_mat) <- info$Label
rownames(pval_mat) <- info$Label

# Strategy factor
strategies <- factor(info$Strategy, levels = c("Y", "A", "S"))

# Sort by strategy, then by mean log2FC within each strategy
ord <- order(strategies, rowMeans(log2fc_mat, na.rm = TRUE))
log2fc_mat <- log2fc_mat[ord, , drop = FALSE]
pval_mat <- pval_mat[ord, , drop = FALSE]
strategies <- strategies[ord]

# Significance matrix
sig_mat <- matrix("", nrow = nrow(pval_mat), ncol = ncol(pval_mat),
                  dimnames = dimnames(pval_mat))
sig_mat[pval_mat < 0.05] <- "*"
sig_mat[pval_mat < 0.01] <- "**"
sig_mat[pval_mat < 0.001] <- "***"

# Row annotation
annotation_row <- data.frame(Strategy = strategies)
rownames(annotation_row) <- rownames(log2fc_mat)
ann_colors <- list(Strategy = c("Y" = "#E41A1C", "A" = "#377EB8", "S" = "#4DAF4A"))

# ---------- 5. Draw heatmap (Fig.4A) ----------
max_abs <- max(abs(log2fc_mat), na.rm = TRUE) * 1.05
breaks <- seq(-max_abs, max_abs, length.out = 101)
colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)

pheatmap(log2fc_mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         cellwidth = 20, cellheight = 10,
         display_numbers = sig_mat,
         number_color = "black",
         fontsize_number = 14,
         fontface_number = "bold",
         color = colors,
         breaks = breaks,
         annotation_row = annotation_row,
         annotation_colors = ann_colors,
         main = "Synergistic regulation of YAS pathways (Dual vs Single, log2FC, * p<0.05)",
         filename = file.path(script_dir, "Fig.4A_heatmap.pdf"),
         width = 8, height = 8)

pheatmap(log2fc_mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         cellwidth = 20, cellheight = 10,
         display_numbers = sig_mat,
         number_color = "black",
         fontsize_number = 14,
         fontface_number = "bold",
         color = colors,
         breaks = breaks,
         annotation_row = annotation_row,
         annotation_colors = ann_colors,
         main = "Synergistic regulation of YAS pathways (Dual vs Single, log2FC, * p<0.05)",
         filename = file.path(script_dir, "Fig.4A_heatmap.png"),
         width = 8, height = 8, dpi = 300)

cat("Fig.4A heatmap saved to", script_dir, "\n")
