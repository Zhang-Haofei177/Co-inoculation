# ============================================================================
# Fig.3C-D: Genus x CAZy family co-occurrence network (HY vs YC) -- Full circle bipartite network
# ============================================================================
# Modifications:
#   1. Remove two-row layout section, keep full circle plot
#   2. New: output topology parameters (degree, positive edge ratio) per enzyme family
#   3. Fix: duplicate function definitions, variable name typos, CSV not written bugs
# ============================================================================
library(dplyr)
library(tidyr)
library(igraph)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(cowplot)

# ---------- 0. Paths ----------
out_root <- "C:/Users/ADMIN/Desktop/Co-inoculation/Fig.3"
script_dir <- out_root
dir.create(script_dir, showWarnings = FALSE, recursive = TRUE)
setwd(script_dir)

# ---------- 1. Load genus abundance ----------
genus_file <- file.path(out_root, "data/genus_abund.rds")
if (file.exists(genus_file)) {
  genus_abund <- readRDS(genus_file)
} else {
  genus_raw <- read.table(file.path(out_root, "data/Unigenes.absolute.g.xls"),
                          header = TRUE, sep = "\t", row.names = 1, quote = "", comment.char = "")
  genus_raw <- genus_raw[rownames(genus_raw) != "Others", ]
  if ("Description" %in% colnames(genus_raw)) genus_raw <- genus_raw[, -which(colnames(genus_raw) == "Description")]
  genus_abund <- as.data.frame(t(genus_raw))
  genus_abund[] <- lapply(genus_abund, as.numeric)
  saveRDS(genus_abund, genus_file)
}

# Load CAZy family abundance
family_abund <- readRDS(file.path(out_root, "data/family_abund.rds"))
family_abund <- family_abund[,c("GH94","GH130","GT4","GT2","GT1","GH32","GT31","GT14","GT7","AA3","AA5","GH13","GT61","GH18")]

# Group info
group <- read.csv(file.path(out_root, "data/group.csv"), row.names = 1)

# ---------- 2. Filter rhizosphere samples and split by site ----------
rhizo_samples <- rownames(group)[group$Niche == "Rhizosphere"]

genus_rhizo <- genus_abund[intersect(rownames(genus_abund), rhizo_samples), ]
genus_rhizo <- genus_rhizo[, colSums(genus_rhizo) > 0]

fam_rhizo <- family_abund[intersect(rownames(family_abund), rhizo_samples), ]
fam_rhizo <- fam_rhizo[, colSums(fam_rhizo) > 0]

common_s <- intersect(rownames(genus_rhizo), rownames(fam_rhizo))
genus_rhizo <- genus_rhizo[common_s, ]
fam_rhizo <- fam_rhizo[common_s, ]

site_vec <- group[common_s, "Site"]

hy_s <- common_s[site_vec == "HY"]
yc_s <- common_s[site_vec == "YC"]

# ---------- 3. Network construction function ----------
build_bipartite_network <- function(genus_mat, family_mat,
                                    cor_threshold = 0.6, p_adj_threshold = 0.05) {
  cor_mat <- matrix(NA, ncol(family_mat), ncol(genus_mat))
  rownames(cor_mat) <- colnames(family_mat)
  colnames(cor_mat) <- colnames(genus_mat)
  p_mat <- cor_mat

  for (i in seq_len(ncol(family_mat))) {
    for (j in seq_len(ncol(genus_mat))) {
      res <- cor.test(family_mat[, i], genus_mat[, j], method = "spearman")
      cor_mat[i, j] <- res$estimate
      p_mat[i, j] <- res$p.value
    }
  }

  p_vals <- p_mat[!is.na(p_mat)]
  p_adj <- p.adjust(p_vals, method = "BH")
  p_mat_adj <- p_mat
  p_mat_adj[!is.na(p_mat)] <- p_adj

  sig_edges <- which(p_mat_adj < p_adj_threshold & abs(cor_mat) > cor_threshold, arr.ind = TRUE)
  if (nrow(sig_edges) == 0) return(NULL)

  edges <- data.frame(
    from = rownames(cor_mat)[sig_edges[, 1]],
    to = colnames(cor_mat)[sig_edges[, 2]],
    weight = cor_mat[sig_edges],
    p_adj = p_mat_adj[sig_edges],
    stringsAsFactors = FALSE
  )

  g <- graph_from_data_frame(edges, directed = FALSE)
  V(g)$type <- ifelse(V(g)$name %in% colnames(family_mat), "Family", "Genus")

  return(g)
}

# Build HY and YC networks
g_hy <- build_bipartite_network(genus_rhizo[hy_s, ], fam_rhizo[hy_s, ])
g_yc <- build_bipartite_network(genus_rhizo[yc_s, ], fam_rhizo[yc_s, ])

# ---------- 4. Network topology parameters ----------
compute_network_params <- function(g, site_name) {
  if (is.null(g)) {
    return(data.frame(Site = site_name, Edges = 0, Nodes_Genus = 0, Nodes_Family = 0,
                      Positive_Edges = 0, Negative_Edges = 0, Positive_Ratio = NA,
                      Avg_Degree = NA, Density = NA, stringsAsFactors = FALSE))
  }
  pos_edges <- sum(E(g)$weight > 0)
  neg_edges <- sum(E(g)$weight < 0)
  data.frame(
    Site = site_name,
    Edges = ecount(g),
    Nodes_Genus = length(unique(V(g)$name[V(g)$type == "Genus"])),
    Nodes_Family = length(unique(V(g)$name[V(g)$type == "Family"])),
    Positive_Edges = pos_edges,
    Negative_Edges = neg_edges,
    Positive_Ratio = ifelse(ecount(g) > 0, pos_edges / ecount(g), NA),
    Avg_Degree = mean(degree(g)),
    Density = edge_density(g),
    stringsAsFactors = FALSE
  )
}

net_params <- rbind(
  compute_network_params(g_hy, "HY"),
  compute_network_params(g_yc, "YC")
)
write.csv(net_params, file.path(script_dir, "Fig.3CD_network_params.csv"), row.names = FALSE)
cat("\n========== Network Topology Parameters ==========\n")
print(net_params)

# ---------- 4b. Topology parameters per enzyme family ----------
compute_family_topology <- function(g, site_name) {
  if (is.null(g)) return(data.frame())
  fam_nodes <- V(g)$name[V(g)$type == "Family"]
  if (length(fam_nodes) == 0) return(data.frame())

  result <- data.frame()
  for (fam in fam_nodes) {
    nbrs <- neighbors(g, fam, mode = "all")
    nbr_names <- V(g)$name[nbrs]
    # Get all edges connected to this family
    incident_edges <- incident(g, fam, mode = "all")
    edge_weights <- E(g)$weight[incident_edges]
    pos_count <- sum(edge_weights > 0)
    neg_count <- sum(edge_weights < 0)
    total_edges <- length(edge_weights)

    result <- rbind(result, data.frame(
      Site = site_name,
      CAZy_Family = fam,
      Degree = degree(g)[fam],
      Positive_Edges = pos_count,
      Negative_Edges = neg_count,
      Total_Edges = total_edges,
      Positive_Ratio = ifelse(total_edges > 0, pos_count / total_edges, NA),
      stringsAsFactors = FALSE
    ))
  }
  return(result)
}

family_topology <- rbind(
  compute_family_topology(g_hy, "HY"),
  compute_family_topology(g_yc, "YC")
)
write.csv(family_topology, file.path(script_dir, "Fig.3CD_family_topology.csv"), row.names = FALSE)
cat("\n========== Topology Parameters per Enzyme Family ==========\n")
print(family_topology)

# ---------- 5. Basic igraph visualization ----------
plot_bipartite <- function(g, site) {
  if (is.null(g)) {
    png(file.path(script_dir, paste0("network_", site, ".png")), width = 10, height = 8, units = "in", res = 300)
    plot.new()
    text(0.5, 0.5, paste(site, "- no significant edges"))
    dev.off()
    return()
  }

  E(g)$weight_abs <- abs(E(g)$weight)
  E(g)$sign <- ifelse(E(g)$weight > 0, "positive", "negative")
  if (max(E(g)$weight_abs) == 0) E(g)$weight_abs <- 1

  layout <- layout_with_fr(g, weights = E(g)$weight_abs)

  V(g)$color <- ifelse(V(g)$type == "Family", "#E41A1C", "#377EB8")
  V(g)$shape <- ifelse(V(g)$type == "Family", "square", "circle")
  V(g)$size <- ifelse(V(g)$type == "Family", 5, 2.5)
  V(g)$label.cex <- 0.5
  V(g)$label <- ifelse(V(g)$type == "Family", V(g)$name, NA)

  edge_colors <- ifelse(E(g)$sign == "positive", "#E41A1C60", "#377EB860")
  edge_widths <- E(g)$weight_abs * 3

  pdf(file.path(script_dir, paste0("network_", site, ".pdf")), width = 10, height = 8)
  plot(g, layout = layout,
       vertex.label.color = "black",
       edge.color = edge_colors,
       edge.width = edge_widths,
       main = paste("Bipartite network (red = positive, blue = negative):", site))
  legend("bottomright", legend = c("Family", "Genus", "Positive", "Negative"),
         col = c("#E41A1C", "#377EB8", "#E41A1C", "#377EB8"),
         pch = c(22, 21, NA, NA), lty = c(NA, NA, 1, 1), pt.bg = c("#E41A1C", "#377EB8", NA, NA))
  dev.off()

  png(file.path(script_dir, paste0("network_", site, ".png")), width = 10, height = 8, units = "in", res = 300)
  plot(g, layout = layout,
       vertex.label.color = "black",
       edge.color = edge_colors,
       edge.width = edge_widths,
       main = paste("Bipartite network (red = positive, blue = negative):", site))
  legend("bottomright", legend = c("Family", "Genus", "Positive", "Negative"),
         col = c("#E41A1C", "#377EB8", "#E41A1C", "#377EB8"),
         pch = c(22, 21, NA, NA), lty = c(NA, NA, 1, 1), pt.bg = c("#E41A1C", "#377EB8", NA, NA))
  dev.off()
}

plot_bipartite(g_hy, "HY")
plot_bipartite(g_yc, "YC")

# ---------- 6. Extract list of genera connected to families ----------
extract_family_genera <- function(g) {
  if (is.null(g)) return(data.frame())
  edges <- as_data_frame(g, what = "edges")
  edges %>%
    filter(from %in% V(g)$name[V(g)$type == "Family"] |
             to %in% V(g)$name[V(g)$type == "Family"]) %>%
    mutate(Family = ifelse(from %in% V(g)$name[V(g)$type == "Family"], from, to),
           Genus = ifelse(from %in% V(g)$name[V(g)$type == "Family"], to, from)) %>%
    dplyr::select(Family, Genus, weight, p_adj)
}

hy_edges <- extract_family_genera(g_hy)
yc_edges <- extract_family_genera(g_yc)

# ============================================================================
# Build genus-to-phylum mapping table (from original absolute abundance file)
# ============================================================================
genus_raw <- read.table(
  file.path(out_root, "data/Unigenes.absolute.g.xls"),
  header = TRUE, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE
)

colnames(genus_raw)[1] <- "Genus_node"
colnames(genus_raw)[ncol(genus_raw)] <- "Taxonomy"

genus_to_phylum <- genus_raw %>%
  dplyr::select(Genus_node, Taxonomy) %>%
  distinct() %>%
  mutate(Phylum = sapply(Taxonomy, function(tax) {
    if (is.na(tax) || tax == "") return("Unclassified")
    m <- regmatches(tax, regexec("p__([^;]+)", tax))
    if (length(m[[1]]) > 1) m[[1]][2] else "Unclassified"
  })) %>%
  dplyr::select(Genus_node, Phylum)

genus_nodes_hy <- V(g_hy)$name[V(g_hy)$type == "Genus"]
genus_nodes_yc <- V(g_yc)$name[V(g_yc)$type == "Genus"]
all_genera <- unique(c(genus_nodes_hy, genus_nodes_yc))

key_genera_phyla <- data.frame(Genus_node = all_genera) %>%
  left_join(genus_to_phylum, by = "Genus_node")

unique_phyla <- unique(key_genera_phyla$Phylum)
cat("Phyla for key genera (total", length(unique_phyla), "):\n")
print(unique_phyla)

# Save genus-to-phylum mapping table
write.csv(key_genera_phyla, file.path(script_dir, "genus_phylum_map.csv"), row.names = FALSE)

# ============================================================================
# Fig.3C-D: Full circle bipartite network (360 degrees)
# ============================================================================
# Color definitions
cat_colors <- c(
  "Carbon source utilization"          = "#6BAED6",
  "EPS synthesis & osmoprotectant"    = "#74C476",
  "ROS scavenging"          = "#FB6A4A",
  "Phytohormone homeostasis"            = "#9E9AC8",
  "Cell wall remodeling & signaling" = "#FDBF6F"
)

family_func <- data.frame(
  ID = c("GH94","GH130","GT4","GT2","GT1","GH32","GT31","GT14","GT7","AA3","AA5","GH13","GT61","GH18"),
  Category = c("Carbon source utilization","Carbon source utilization",
               "EPS synthesis & osmoprotectant","EPS synthesis & osmoprotectant",
               "Phytohormone homeostasis","Carbon source utilization",
               "Cell wall remodeling & signaling","Cell wall remodeling & signaling",
               "Cell wall remodeling & signaling","ROS scavenging","ROS scavenging",
               "Carbon source utilization","Cell wall remodeling & signaling","Carbon source utilization"),
  stringsAsFactors = FALSE
)

phy_colors <- c(
  "Pseudomonadota"        = "#377EB8",
  "Myxococcota"           = "#B3DE69",
  "Planctomycetota"       = "#BC80BD",
  "Verrucomicrobiota"     = "#80B1D3",
  "Ascomycota"            = "#8dd3c7",
  "Uroviricota"           = "#FBB4AE",
  "Rhodothermota"         = "#B2DF8A",
  "Armatimonadota"        = "#9400D3",
  "Acidobacteriota"       = "#FB8072",
  "Actinomycetota"        = "#CCEBC5",
  "Chloroflexota"         = "#BEBADA",
  "Candidatus Tectomicrobia" = "#33A02C",
  "Cyanobacteriota"       = "#E41A1C",
  "Thermodesulfobacteriota" = "#00FFFF",
  "Gemmatimonadota"       = "#FDB462",
  "Bacillota"             = "#FFED6F",
  "Euryarchaeota"         = "#984EA3",
  "Spirochaetota"         = "#FF7F00",
  "Abditibacteriota"      = "#00FF7F",
  "Candidatus Moduliflexota" = "#A65628",
  "Unclassified"          = "grey85",
  "Other"                 = "grey85"
)

# ============================================================================
# Full circle bipartite network plot function
# ============================================================================
plot_bipartite_circle_manual <- function(g, site,
                                         linewidth_range = c(0.1, 0.6),
                                         background_r0 = 0.65,
                                         background_r1 = 1.02,
                                         label_offset = 0.01) {
  if (is.null(g) || vcount(g) == 0) return(NULL)

  # ---------- 1. Family node angle assignment (sorted by degree) ----------
  fam_nodes <- V(g)$name[V(g)$type == "Family"]
  deg_fam <- degree(g)[fam_nodes]
  df_fam <- data.frame(name = fam_nodes, degree = deg_fam, stringsAsFactors = FALSE)
  df_fam <- df_fam[order(df_fam$degree, decreasing = TRUE), ]
  n_fam <- nrow(df_fam)
  if (n_fam > 0) {
    angles_fam <- seq(0, 2*pi*(1 - 1/n_fam), length.out = n_fam)
    mid <- (n_fam + 1) / 2
    order_idx <- order(abs(seq_len(n_fam) - mid))
    df_fam$angle <- angles_fam[order_idx]
  } else {
    df_fam <- data.frame(name = character(), degree = numeric(), angle = numeric())
  }

  # ---------- 2. Microbial node angle: weighted mean of connected family angles ----------
  gen_nodes <- V(g)$name[V(g)$type == "Genus"]
  if (length(gen_nodes) > 0) {
    fam_angle_dict <- setNames(df_fam$angle, df_fam$name)
    mean_angle <- sapply(gen_nodes, function(gn) {
      nbrs <- neighbors(g, gn, mode = "all")
      fam_nbrs <- nbrs[V(g)$type[nbrs] == "Family"]
      if (length(fam_nbrs) == 0) return(NA_real_)
      fam_names <- V(g)$name[fam_nbrs]
      angles <- fam_angle_dict[fam_names]
      angles <- angles[!is.na(angles)]
      if (length(angles) == 0) return(NA_real_)
      edge_weights <- sapply(fam_names, function(fn) {
        e <- get.edge.ids(g, c(fn, gn))
        if (e == 0) e <- get.edge.ids(g, c(gn, fn))
        if (e == 0) return(1)
        abs(E(g)$weight[e])
      })
      weighted.mean(angles, w = edge_weights)
    })
    df_gen <- data.frame(name = gen_nodes, degree = degree(g)[gen_nodes],
                         mean_angle = mean_angle, stringsAsFactors = FALSE)
    df_gen_connected <- df_gen[!is.na(df_gen$mean_angle), ]
    df_gen_isolated <- df_gen[is.na(df_gen$mean_angle), ]
    if (nrow(df_gen_connected) > 0) {
      df_gen_connected <- df_gen_connected[order(df_gen_connected$mean_angle), ]
      n_conn <- nrow(df_gen_connected)
      angles_gen_conn <- seq(0, 2*pi*(1 - 1/n_conn), length.out = n_conn)
      df_gen_connected$angle <- angles_gen_conn
    }
    if (nrow(df_gen_isolated) > 0) {
      df_gen_isolated <- df_gen_isolated[order(df_gen_isolated$degree, decreasing = TRUE), ]
      n_iso <- nrow(df_gen_isolated)
      left_angles <- seq(0, 0.1*pi, length.out = floor(n_iso/2) + 2)
      right_angles <- seq(2*pi - 0.1*pi, 2*pi, length.out = floor(n_iso/2) + 2)
      iso_angles <- c(left_angles[-c(1, length(left_angles))],
                      right_angles[-c(1, length(right_angles))])
      if (length(iso_angles) < n_iso) {
        iso_angles <- c(iso_angles, rep(pi, n_iso - length(iso_angles)))
      }
      df_gen_isolated$angle <- iso_angles[1:n_iso]
      df_gen <- rbind(df_gen_connected, df_gen_isolated)
    } else {
      df_gen <- df_gen_connected
    }
    df_gen$angle <- pmax(0, pmin(2*pi, df_gen$angle))
    df_gen <- df_gen[order(df_gen$angle), ]
  } else {
    df_gen <- data.frame(name = character(), degree = numeric(), angle = numeric())
  }

  # ---------- 3. Merge node coordinates ----------
  r_inner <- 0.7
  r_outer <- 1.0
  node_data <- rbind(
    data.frame(name = df_fam$name, type = "Family", degree = df_fam$degree,
               angle = df_fam$angle, stringsAsFactors = FALSE),
    data.frame(name = df_gen$name, type = "Genus", degree = df_gen$degree,
               angle = df_gen$angle, stringsAsFactors = FALSE)
  )
  node_data$x <- ifelse(node_data$type == "Family",
                        r_inner * cos(node_data$angle),
                        r_outer * cos(node_data$angle))
  node_data$y <- ifelse(node_data$type == "Family",
                        r_inner * sin(node_data$angle),
                        r_outer * sin(node_data$angle))

  # ---------- 4. Node colors ----------
  node_data$color <- NA
  for (i in seq_len(nrow(node_data))) {
    if (node_data$type[i] == "Family") {
      fam <- node_data$name[i]
      cat <- family_func$Category[match(fam, family_func$ID)]
      if (is.na(cat)) cat <- "Other"
      node_data$color[i] <- cat_colors[cat]
    } else {
      raw <- node_data$name[i]
      phy <- genus_to_phylum$Phylum[match(raw, genus_to_phylum$Genus_node)]
      if (is.na(phy) || phy == "") phy <- "Unclassified"
      node_data$color[i] <- ifelse(phy %in% names(phy_colors), phy_colors[phy], "grey70")
    }
  }
  node_data$label <- sapply(node_data$name, function(nm) {
    if (node_data$type[node_data$name == nm] == "Family") nm
    else if (grepl("g__", nm)) sub(".*g__", "", nm) else nm
  })

  # ---------- 5. Edge data ----------
  edge_list <- as_data_frame(g, what = "edges")
  edge_data <- data.frame(from = edge_list$from, to = edge_list$to, weight = edge_list$weight)
  edge_data$x_from <- node_data$x[match(edge_data$from, node_data$name)]
  edge_data$y_from <- node_data$y[match(edge_data$from, node_data$name)]
  edge_data$x_to   <- node_data$x[match(edge_data$to, node_data$name)]
  edge_data$y_to   <- node_data$y[match(edge_data$to, node_data$name)]
  edge_data$sign   <- ifelse(edge_data$weight > 0, "positive", "negative")
  edge_data$abs_weight <- abs(edge_data$weight)
  edge_data <- edge_data[complete.cases(edge_data[, c("x_from", "y_from", "x_to", "y_to")]), ]

  # ---------- 6. Circular background ----------
  theta <- seq(0, 2*pi, length.out = 200)
  x_outer <- background_r1 * cos(theta)
  y_outer <- background_r1 * sin(theta)
  theta_rev <- rev(theta)
  x_inner <- background_r0 * cos(theta_rev)
  y_inner <- background_r0 * sin(theta_rev)
  background_poly <- data.frame(x = c(x_outer, x_inner), y = c(y_outer, y_inner))

  # ---------- 7. Label positions (family labels only) ----------
  fam_labels <- node_data[node_data$type == "Family", ]
  fam_labels$radius_label <- background_r0 - label_offset
  fam_labels$x_label <- fam_labels$radius_label * cos(fam_labels$angle)
  fam_labels$y_label <- fam_labels$radius_label * sin(fam_labels$angle)
  fam_labels$angle_deg <- fam_labels$angle * 180 / pi
  fam_labels$angle_deg <- ifelse(fam_labels$angle_deg > 90 | fam_labels$angle_deg < -90,
                                 fam_labels$angle_deg + 180, fam_labels$angle_deg)
  fam_labels_left <- fam_labels[fam_labels$angle < pi/2 | fam_labels$angle > 3*pi/2, ]
  fam_labels_right <- fam_labels[fam_labels$angle >= pi/2 & fam_labels$angle <= 3*pi/2, ]

  # ---------- 8. Plot ----------
  p <- ggplot() +
    geom_polygon(data = background_poly, aes(x = x, y = y), fill = "#E4EBC1", alpha = 0.6) +
    geom_curve(data = edge_data,
               aes(x = x_from, y = y_from, xend = x_to, yend = y_to,
                   color = sign, linewidth = abs_weight),
               curvature = 0.2, alpha = 1, lineend = "round") +
    geom_point(data = node_data, aes(x = x, y = y, fill = color, size = degree),
               shape = 21, stroke = 0.3, color = "black") +
    scale_fill_identity() +
    scale_color_manual(values = c("positive" = "#FC876D", "negative" = "#91C3E0"),
                       name = "Spearman's r") +
    scale_size_continuous(range = c(2, 30), guide = "none") +
    scale_linewidth_continuous(range = linewidth_range, guide = "none") +
    geom_text(data = fam_labels_left,
              aes(x = x_label, y = y_label, label = label, angle = angle_deg),
              size = 3.5, hjust = 0, vjust = 0.5, color = "black", fontface = "bold") +
    geom_text(data = fam_labels_right,
              aes(x = x_label, y = y_label, label = label, angle = angle_deg),
              size = 3.5, hjust = 1, vjust = 0.5, color = "black", fontface = "bold") +
    coord_fixed(clip = "off", xlim = c(-1.3, 1.3), ylim = c(-1.3, 1.3)) +
    labs(title = paste("Bipartite network (full circle):", site)) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "right")

  # Name by site: Fig.3C (HY) and Fig.3D (YC)
  panel_label <- ifelse(site == "HY", "Fig.3C", "Fig.3D")
  ggsave(file.path(script_dir, paste0(panel_label, "_network_full_", site, ".pdf")),
         p, width = 12, height = 12, device = "pdf")
  ggsave(file.path(script_dir, paste0(panel_label, "_network_full_", site, ".png")),
         p, width = 12, height = 12, dpi = 300)
  return(p)
}

# Call: keep full circle plot only
plot_bipartite_circle_manual(g_hy, "HY")
plot_bipartite_circle_manual(g_yc, "YC")

# ============================================================================
# Fig.3G-H supplement: Hub family degree bar chart + phylum-level stacked bar chart of associated genera
# Using already computed family_topology / hy_edges / yc_edges / genus_to_phylum
# ============================================================================
cat("\n========== Hub analysis: degree bars + phylum composition ==========\n")

theme_nature <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line      = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks     = element_line(linewidth = 0.35, colour = "black"),
      axis.title     = element_text(size = base_size),
      axis.text      = element_text(size = base_size - 0.5, colour = "black"),
      legend.title   = element_text(size = base_size - 0.3),
      legend.text    = element_text(size = base_size - 0.7),
      legend.key.size = unit(0.35, "cm"),
      strip.text     = element_text(size = base_size - 0.3, face = "bold"),
      plot.title     = element_text(size = base_size + 1, face = "bold", hjust = 0.5),
      plot.subtitle  = element_text(size = base_size - 0.5, hjust = 0.5, colour = "grey40"),
      panel.grid     = element_blank()
    )
}

plot_hub_analysis <- function(site_name, edges_df) {
  # --- Panel A: All families degree bar chart ---
  cent_df <- family_topology %>%
    filter(Site == site_name) %>%
    arrange(desc(Degree))

  if (nrow(cent_df) == 0) return(NULL)

  cent_df$Category <- family_func$Category[match(cent_df$CAZy_Family, family_func$ID)]
  cent_df$CAZy_Family <- factor(cent_df$CAZy_Family, levels = cent_df$CAZy_Family)

  # Hub: top 2 enzyme families by degree
  hub_fams <- head(cent_df$CAZy_Family, 2)
  cent_df$is_hub <- cent_df$CAZy_Family %in% hub_fams

  cat(sprintf("  %s: top-2 hub families by degree: %s\n",
              site_name, paste(hub_fams, collapse = ", ")))

  p_a <- ggplot(cent_df, aes(x = CAZy_Family, y = Degree, fill = Category)) +
    geom_col(width = 0.65, color = "black", linewidth = 0.2) +
    geom_text(aes(label = Degree), vjust = -0.5, size = 2.8, family = "sans") +
    scale_fill_manual(values = cat_colors, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(title = site_name, x = NULL, y = "Degree") +
    theme_nature(7) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 6.5))

  # --- Panel B: Phylum-level stacked bar chart of genera associated with hub families ---
  if (length(hub_fams) > 0 && nrow(edges_df) > 0) {
    # Extract genera connected to hub families from edge list
    hub_edges <- edges_df %>%
      filter(Family %in% hub_fams) %>%
      rename(Hub_Family = Family)

    if (nrow(hub_edges) > 0) {
      hub_phy <- hub_edges %>%
        left_join(genus_to_phylum, by = c("Genus" = "Genus_node"))
      hub_phy$Phylum[is.na(hub_phy$Phylum)] <- "Unclassified"

      phy_n <- hub_phy %>% count(Phylum) %>% arrange(desc(n))
      top_p <- phy_n$Phylum[1:min(8, nrow(phy_n))]
      hub_phy$Phylum_plt <- ifelse(hub_phy$Phylum %in% top_p, hub_phy$Phylum, "Other")

      phy_order <- hub_phy %>%
        count(Phylum_plt) %>%
        arrange(desc(n)) %>%
        pull(Phylum_plt)
      phy_order <- c(setdiff(phy_order, "Other"), "Other")
      hub_phy$Phylum_plt <- factor(hub_phy$Phylum_plt, levels = phy_order)

      p_b <- ggplot(hub_phy, aes(x = Hub_Family, fill = Phylum_plt)) +
        geom_bar(position = "fill", color = "black", linewidth = 0.15) +
        scale_fill_manual(values = phy_colors, name = "Phylum") +
        scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
        labs(title = "Phylum of connected genera", x = "Hub family", y = NULL) +
        theme_nature(7) +
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 6.5),
              legend.key.size = unit(0.25, "cm"), legend.text = element_text(size = 5.5))

      p_combined <- p_a + p_b + plot_layout(widths = c(1.6, 1))
    } else {
      p_combined <- p_a
    }
  } else {
    p_combined <- p_a
  }

  # Save
  panel_label <- ifelse(site_name == "HY", "Fig.3G", "Fig.3H")
  ggsave(file.path(script_dir, paste0(panel_label, "_Hub_analysis_", site_name, ".pdf")),
         p_combined, width = 183, height = 80, units = "mm", device = "pdf")
  ggsave(file.path(script_dir, paste0(panel_label, "_Hub_analysis_", site_name, ".png")),
         p_combined, width = 183, height = 80, units = "mm", dpi = 300)

  print(p_combined)
  return(p_combined)
}

plot_hub_analysis("HY", hy_edges)
plot_hub_analysis("YC", yc_edges)

cat("  Hub analysis figures saved\n")

# ============================================================================
# Microbial phylum-level color legend
# ============================================================================
df_legend <- data.frame(Phylum = names(phy_colors), Color = unlist(phy_colors))

p_legend <- ggplot(df_legend, aes(x = 1, y = rev(seq_len(nrow(df_legend))), fill = Phylum)) +
  geom_tile(show.legend = TRUE) +
  scale_fill_manual(values = setNames(df_legend$Color, df_legend$Phylum)) +
  theme_void() +
  theme(legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10, face = "plain"),
        legend.key.size = unit(0.6, "cm")) +
  guides(fill = guide_legend(title = "Phylum", ncol = 1))

legend_only <- get_legend(p_legend)

png(file.path(script_dir, "Fig.3CD_Phylum_legend.png"), width = 800, height = 1000, res = 150)
grid::grid.draw(legend_only)
dev.off()

pdf(file.path(script_dir, "Fig.3CD_Phylum_legend.pdf"), width = 6, height = 8)
grid::grid.draw(legend_only)
dev.off()

cat("\n========== Analysis Complete ==========\n")
cat("Output files saved to:", script_dir, "\n")
