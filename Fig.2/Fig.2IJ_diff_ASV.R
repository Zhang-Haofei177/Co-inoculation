# ============================================================================
# Fig.2I+J: Rhizosphere differential ASV scatter plot (DESeq2 Dual vs Single)
# Bacteria (Fig.2I) + Fungi (Fig.2J) combined
# ============================================================================
library(DESeq2)
library(dplyr)
library(ggplot2)
library(tibble)
library(patchwork)

out_root <- "C:/Users/ADMIN/Desktop/Co-inoculation/Fig.2"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
setwd(out_root)

# ---------- 1. Read group metadata ----------
group <- read.csv("data/group.csv", row.names = 1)
group$Treatment <- factor(group$Treatment, levels = c("Ctrl", "AMF", "Str", "AMF+Str"))
group$Niche     <- factor(group$Niche)
group$Site      <- factor(group$Site, levels = c("HY", "YC"))

# ---------- 2. DESeq2 analysis function ----------
analyze_dual_vs_single <- function(otu_file, tax_file, microbe_name) {
  cat("\n========== Analyzing", microbe_name, "==========\n")
  otu <- read.csv(otu_file, row.names = 1)
  rownames(otu) <- paste0(microbe_name, "_", rownames(otu))
  tax <- read.csv(tax_file, row.names = 1)
  rownames(tax) <- paste0(microbe_name, "_", rownames(tax))
  tax$Phylum <- as.character(tax$Phylum)
  tax$Phylum[is.na(tax$Phylum) | tax$Phylum == ""] <- "Other"
  tax$Phylum[grepl("subgroup|unclassified|uncultured", tax$Phylum, ignore.case = TRUE)] <- "Other"

  keep <- rowSums(otu) >= 10 & rowSums(otu > 0) >= 2
  otu <- otu[keep, ]
  cat("Filtered ASV count:", nrow(otu), "\n")
  otu_t <- t(otu)

  idx_rhizo <- group$Niche == "Rhizosphere"
  meta_rhizo <- group[idx_rhizo, , drop = FALSE]
  common <- intersect(rownames(meta_rhizo), rownames(otu_t))
  otu_rhizo <- otu_t[common, ]
  meta_rhizo <- meta_rhizo[common, ]

  meta_rhizo$Treatment_bin <- ifelse(meta_rhizo$Treatment == "AMF+Str", "Dual", "Single")
  meta_rhizo$Treatment_bin <- factor(meta_rhizo$Treatment_bin, levels = c("Single", "Dual"))

  sites <- levels(meta_rhizo$Site)
  res_list <- list()
  for (st in sites) {
    cat("  Processing site:", st, "\n")
    idx <- meta_rhizo$Site == st
    sub_meta <- meta_rhizo[idx, , drop = FALSE]
    sub_otu <- otu_rhizo[rownames(sub_meta), , drop = FALSE]
    if (nrow(sub_meta) < 6 || length(unique(sub_meta$Treatment_bin)) < 2) next
    count_mat <- t(sub_otu)
    mode(count_mat) <- "integer"
    dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = sub_meta, design = ~ Treatment_bin)
    dds <- DESeq(dds, quiet = TRUE)
    res <- results(dds, contrast = c("Treatment_bin", "Dual", "Single"))
    sig <- subset(res, pvalue < 0.05)
    if (nrow(sig) > 0) {
      df <- as.data.frame(sig) %>%
        rownames_to_column("ASV") %>%
        mutate(Site = st, Niche = "Rhizosphere", Contrast = "Dual_vs_Single", Microbe = microbe_name)
      res_list[[st]] <- df
    }
  }
  if (length(res_list) == 0) return(NULL)
  res_all <- bind_rows(res_list)
  res_all$Phylum <- tax[res_all$ASV, "Phylum"]
  res_all$Phylum[is.na(res_all$Phylum) | res_all$Phylum == ""] <- "Other"
  write.csv(res_all, file.path(out_root, paste0(microbe_name, "_Rhizosphere_diff_ASV.csv")), row.names = FALSE)
  return(res_all)
}

# ---------- 3. Run analysis for each kingdom ----------
bac_res <- analyze_dual_vs_single("data/BAC_rarefied.csv", "data/BAC_Tax.csv", "Bacteria")
fun_res <- analyze_dual_vs_single("data/FUN_rarefied.csv", "data/FUN_Tax.csv", "Fungi")

# ---------- 4. Phylum color mapping ----------
bac_phylum_colors <- c(
  "Pseudomonadota"="#377EB8","Actinomycetota"="#CCEBC5","Gemmatimonadota"="#FDB462",
  "Bacillota"="#FFED6F","Acidobacteriota"="#FB8072","Bacteroidota"="#FED9A6",
  "Myxococcota"="#B3DE69","Chloroflexota"="#BEBADA","Verrucomicrobiota"="#80B1D3",
  "Nitrospirota"="#A6CEE3","Entotheonellaeota"="#8DD3C7","Bdellovibrionota"="#FCCDE5",
  "Planctomycetota"="#BC80BD","Others"="#D3D3D3"
)
fun_phylum_colors <- c(
  "Ascomycota"="#8dd3c7","Basidiomycota"="#ffffb3","Mortierellomycota"="#bebada",
  "Glomeromycota"="#fb8072","Chytridiomycota"="#80b1d3","Rozellomycota"="#fdb462",
  "Blastocladiomycota"="#b3de69","Others"="gray"
)

# ---------- 5. Plotting function ----------
plot_microbe <- function(res_df, microbe, phylum_colors) {
  if (is.null(res_df) || nrow(res_df) == 0) return(NULL)
  res_df$Phylum_plot <- ifelse(res_df$Phylum %in% names(phylum_colors), res_df$Phylum, "Others")
  phyla_present <- unique(res_df$Phylum_plot)
  colors_sub <- phylum_colors[phyla_present]
  colors_sub[is.na(colors_sub)] <- "gray"

  ggplot(res_df, aes(x = Site, y = log2FoldChange, color = Phylum_plot, size = -log10(pvalue))) +
    geom_vline(xintercept = 1.5, linetype = "solid", color = "gray", linewidth = 0.2) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.2) +
    geom_jitter(alpha = 0.8, width = 0.2, stroke = 0.1) +
    scale_color_manual(values = colors_sub) +
    scale_size_continuous(range = c(0.1, 2)) +
    labs(x = NULL, y = "log2 Fold Change (Dual vs Single)", title = microbe,
         subtitle = "p < 0.05 (unadjusted)") +
    theme_bw(base_size = 6) +
    theme(
      axis.text = element_text(size = 6, color = "black"),
      axis.title = element_text(size = 6, color = "black"),
      legend.text = element_text(size = 6, color = "black"),
      legend.title = element_text(size = 6, color = "black"),
      legend.position = "right",
      panel.grid.major.x = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 6, color = "black"),
      plot.subtitle = element_text(hjust = 0.5, size = 6, color = "black"),
      text = element_text(family = "sans", color = "black")
    ) +
    guides(color = guide_legend(title = "Phylum", override.aes = list(size = 4)),
           size = guide_legend(title = "-log10(p)"))
}

p_bac <- plot_microbe(bac_res, "Bacteria", bac_phylum_colors)
p_fun <- plot_microbe(fun_res, "Fungi", fun_phylum_colors)

# ---------- 6. Combine panels ----------
if (!is.null(p_bac) && !is.null(p_fun)) {
  combined <- p_bac + p_fun +
    plot_annotation(
      title = "Rhizosphere: Significantly altered ASVs (Dual vs Single inoculum)",
      theme = theme(plot.title = element_text(hjust = 0.5, size = 6, color = "black", family = "sans"))
    )
} else if (!is.null(p_bac)) {
  combined <- p_bac
} else {
  combined <- p_fun
}

ggsave(file.path(out_root, "Fig.2IJ_diff_ASV_scatter.pdf"), combined, width = 100, height = 60, units = "mm")
ggsave(file.path(out_root, "Fig.2IJ_diff_ASV_scatter.png"), combined, width = 100, height = 60, units = "mm", dpi = 300)
cat("\nFig.2I+J saved\n")
