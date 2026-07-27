# ============================================================================
# Fig.S4: Rarefaction curves (Bacteria + Fungi, Niche + Treatment, 2×2)
# ============================================================================
library(vegan)
library(ggplot2)
library(patchwork)

out_root <- "C:/Users/ADMIN/Desktop/Co-inoculation/Fig.S4"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
setwd(out_root)

# ---------- 1. Group info ----------
group <- read.csv("data/group.csv", stringsAsFactors = FALSE)
colnames(group)[1] <- "names"
group <- group[group$Treatment %in% c("AMF","Str","AMF+Str") &
               group$Niche %in% c("Topsoil","Rhizosphere","Endosphere"), ]
group$Niche <- factor(group$Niche, levels = c("Topsoil","Rhizosphere","Endosphere"))
group$Site  <- factor(group$Site, levels = c("HY","YC"), labels = c("Heyang","Yuci"))

# ---------- 2. Rarefaction function ----------
ly.rare <- function(a, n) {
  b <- rowSums(a)
  i <- 1
  n.i <- as.vector(c(seq(0, b[i], n), b[i]))
  rare.i <- as.integer(rarefy(a[i, ], n.i))
  rare <- cbind(rep(names(b)[i], length(n.i)), n.i, rare.i)
  for (i in 2:length(b)) {
    n.i <- as.vector(c(seq(0, b[i], n), b[i]))
    rare.i <- as.integer(rarefy(a[i, ], n.i))
    rare.ni <- cbind(rep(names(b)[i], length(n.i)), n.i, rare.i)
    rare <- rbind(rare, rare.ni)
  }
  colnames(rare)[1:3] <- c("names","reads","value")
  rare <- as.data.frame(rare)
  rare$reads <- as.numeric(as.character(rare$reads))
  rare$value <- as.numeric(as.character(rare$value))
  rare
}

# ---------- 3. Processing function ----------
process <- function(asv_file, microbe_name, step = 2000) {
  cat(sprintf("\n=== %s ===\n", microbe_name))
  a1 <- read.csv(asv_file, row.names = 1, header = TRUE, check.names = FALSE)
  a1 <- a1[, group$names]
  a1 <- t(a1)
  cat(sprintf("Samples: %d, ASVs: %d\n", nrow(a1), ncol(a1)))

  rare <- ly.rare(a1, step)
  rare.new <- merge(rare, group, by = "names")

  # By Niche
  p1 <- ggplot(rare.new, aes(reads, value, group = names, colour = Niche)) +
    geom_line(linewidth = 0.3, alpha = 0.7) +
    facet_wrap(~ Site, ncol = 2) +
    scale_color_manual(values = c("Topsoil"="#FCB255","Rhizosphere"="#42B382","Endosphere"="#107BA4")) +
    labs(x = "Number of sequences", y = "Number of ASVs", title = NULL) +
    theme_bw(base_size = 6) +
    theme(panel.grid = element_blank(),
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(size = 6, color = "black"),
          legend.key.size = unit(0.25, "cm"),
          legend.text = element_text(size = 6),
          legend.title = element_blank(),
          axis.text = element_text(size = 6, color = "black"),
          axis.title = element_text(size = 6, color = "black"),
          plot.title = element_blank(),
          text = element_text(family = "sans"))

  # By Treatment
  p2 <- ggplot(rare.new, aes(reads, value, group = names, colour = Treatment)) +
    geom_line(linewidth = 0.3, alpha = 0.7) +
    facet_wrap(~ Site, ncol = 2) +
    scale_color_manual(values = c("AMF"="#FCB255","Str"="#42B382","AMF+Str"="#107BA4")) +
    labs(x = "Number of sequences", y = "Number of ASVs", title = NULL) +
    theme_bw(base_size = 6) +
    theme(panel.grid = element_blank(),
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(size = 6, color = "black"),
          legend.key.size = unit(0.25, "cm"),
          legend.text = element_text(size = 6),
          legend.title = element_blank(),
          axis.text = element_text(size = 6, color = "black"),
          axis.title = element_text(size = 6, color = "black"),
          plot.title = element_blank(),
          text = element_text(family = "sans"))

  return(list(p1 = p1, p2 = p2))
}

# ---------- 4. Generate and combine ----------
bac_plots <- process("data/BAC_Absolute.csv", "Bacteria")
fun_plots <- process("data/FUN_Absolute.csv", "Fungi")

combined <- (bac_plots$p1 + fun_plots$p1) / (bac_plots$p2 + fun_plots$p2) +
  plot_annotation(
    title = "Rarefaction curves",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 6, hjust = 0.5, face = "bold", family = "sans"),
      plot.tag = element_text(size = 6, face = "bold", family = "sans")
    )
  )

ggsave(file.path(out_root, "Fig.S4_rarefaction.pdf"), combined, width = 200, height = 80, units = "mm")
ggsave(file.path(out_root, "Fig.S4_rarefaction.png"), combined, width = 200, height = 80, units = "mm", dpi = 300)
cat("Fig.S4 saved\n")
