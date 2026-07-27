# ============================================================================
# Fig.2C (Bacteria) + Fig.2F (Fungi): Combined beta-diversity PCoA with PERMANOVA annotation
# Faceted PCoA portion removed; single combined PCoA only
# ============================================================================
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)

out_root <- "C:/Users/ADMIN/Desktop/Co-inoculation/Fig.2"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
setwd(out_root)

analyze_beta <- function(otu_file, prefix) {
  cat("\n========== Analyzing", prefix, "==========\n")

  otu <- read.csv(otu_file, row.names = 1)

  asv_presence <- rowSums(otu > 0)
  asv_total <- rowSums(otu)
  keep_asv <- asv_presence >= 2 & asv_total >= 10
  otu_filtered <- otu[keep_asv, ]
  cat("Retained ASV count:", sum(keep_asv), "\n")

  otu_t <- t(otu_filtered)

  group <- read.csv("data/group.csv", row.names = 1)
  group <- group[group$Treatment != "Ctrl", ]
  group <- group[!group$Niche %in% c("Subsoil", "Bottom soil"), ]
  group$Treatment <- factor(group$Treatment, levels = c("AMF", "Str", "AMF+Str"))
  group$Niche <- factor(group$Niche, levels = c("Topsoil", "Rhizosphere", "Endosphere"))
  group$Site <- factor(group$Site, levels = c("HY", "YC"))

  common <- intersect(rownames(group), rownames(otu_t))
  otu_t <- otu_t[common, ]
  group <- group[common, ]

  bray_dist <- vegdist(otu_t, method = "bray")

  # PERMANOVA
  set.seed(123)
  perm_full <- adonis2(bray_dist ~ Site * Niche * Treatment, data = group, permutations = 999, by = "terms")
  perm_full <- as.data.frame(perm_full)
  perm_full$Factor <- rownames(perm_full)
  write.csv(perm_full, file.path(out_root, paste0(prefix, "_PERMANOVA_full.csv")), row.names = FALSE)

  # Stratify by niche
  niches <- levels(group$Niche)
  perm_niche <- list()
  for (nic in niches) {
    idx <- group$Niche == nic
    sub_otu <- otu_t[idx, ]
    sub_group <- group[idx, ]
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

  # Stratify by site
  sites <- levels(group$Site)
  perm_site <- list()
  for (st in sites) {
    idx <- group$Site == st
    sub_otu <- otu_t[idx, ]
    sub_group <- group[idx, ]
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

  # PERMDISP
  beta_disp_treat <- betadisper(bray_dist, group$Treatment)
  sink(file.path(out_root, paste0(prefix, "_PERMDISP_treatment.txt")))
  print(permutest(beta_disp_treat, permutations = 999))
  sink()
  beta_disp_niche <- betadisper(bray_dist, group$Niche)
  sink(file.path(out_root, paste0(prefix, "_PERMDISP_niche.txt")))
  print(permutest(beta_disp_niche, permutations = 999))
  sink()

  # ========== Single combined PCoA plot ==========
  treat_colors <- c("AMF" = "#F0C986", "Str" = "#3FBDA7", "AMF+Str" = "#F08080")

  pcoa_res <- cmdscale(bray_dist, k = 2, eig = TRUE)
  eig <- pcoa_res$eig
  var1 <- round(eig[1]/sum(eig)*100, 2)
  var2 <- round(eig[2]/sum(eig)*100, 2)

  pcoa_df <- data.frame(Sample = rownames(pcoa_res$points),
                        Axis1 = pcoa_res$points[,1],
                        Axis2 = pcoa_res$points[,2])
  pcoa_df <- merge(pcoa_df, group, by.x = "Sample", by.y = "row.names")

  pcoa_df <- pcoa_df %>%
    mutate(
      Site_Treatment = paste(Site, Treatment, sep = "_"),
      Shape = case_when(
        Site_Treatment == "HY_AMF"      ~ 0,
        Site_Treatment == "HY_Str"      ~ 1,
        Site_Treatment == "HY_AMF+Str"  ~ 2,
        Site_Treatment == "YC_AMF"      ~ 15,
        Site_Treatment == "YC_Str"      ~ 16,
        Site_Treatment == "YC_AMF+Str"  ~ 17
      )
    )

  # Extract annotation text from PERMANOVA
  perm_full_terms <- perm_full[!grepl("Residual|Total", perm_full$Factor), ]
  label_lines <- c()
  for (i in 1:nrow(perm_full_terms)) {
    factor <- perm_full_terms$Factor[i]
    r2 <- perm_full_terms$R2[i]
    pv <- perm_full_terms$`Pr(>F)`[i]
    if (pv < 0.001) sig <- "***"
    else if (pv < 0.01) sig <- "**"
    else if (pv < 0.05) sig <- "*"
    else sig <- "ns"
    label_lines <- c(label_lines, sprintf("%s: R-squared=%.3f, p=%.3f%s", factor, r2, pv, sig))
  }
  label_text <- paste(label_lines, collapse = "\n")

  pcoa_df$Shape <- factor(pcoa_df$Shape, levels = c(0, 1, 2, 15, 16, 17))
  shape_labels <- c("HY_Topsoil", "HY_Rhizosphere", "HY_Endosphere",
                    "YC_Topsoil", "YC_Rhizosphere", "YC_Endosphere")

  p <- ggplot(pcoa_df, aes(x = Axis1, y = Axis2, color = Treatment, shape = Shape)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = treat_colors) +
    scale_shape_manual(values = c(0, 1, 2, 15, 16, 17), labels = shape_labels) +
    labs(x = paste0("PCoA1 (", var1, "%)"),
         y = paste0("PCoA2 (", var2, "%)"),
         title = paste(prefix, "beta-diversity (Bray-Curtis)")) +
    annotate("text", x = Inf, y = Inf, label = label_text,
             hjust = 1.1, vjust = 1.1, size = 6/.pt,
             color = "black", family = "sans") +
    theme_bw() +
    theme(axis.text = element_text(size = 6, color = "black"),
          axis.title = element_text(size = 6, color = "black"),
          plot.title = element_text(size = 6, hjust = 0.5, color = "black"),
          legend.text = element_text(size = 6, color = "black"),
          legend.title = element_text(size = 6, color = "black"),
          legend.position = "bottom",
          panel.grid = element_blank(),
          text = element_text(family = "sans", size = 6)) +
    guides(shape = guide_legend(title = "Niche_Site"),
           color = guide_legend(title = "Treatment"))

  panel_label <- ifelse(prefix == "Bacteria", "Fig.2C", "Fig.2F")
  ggsave(file.path(out_root, paste0(panel_label, "_PCoA.pdf")), p, width = 50, height = 70, units = "mm")
  ggsave(file.path(out_root, paste0(panel_label, "_PCoA.png")), p, width = 50, height = 70, units = "mm", dpi = 300)

  cat(prefix, "PCoA done\n")
  return(p)
}

p_bac <- analyze_beta("data/BAC_rarefied.csv", "Bacteria")
p_fun <- analyze_beta("data/FUN_rarefied.csv", "Fungi")
cat("\nFig.2C + Fig.2F done\n")
