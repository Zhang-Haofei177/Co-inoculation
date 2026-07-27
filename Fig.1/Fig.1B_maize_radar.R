# ============================================================================
# Fig.1B: Maize phenotype overlapping radar plot (AMF / Str / AMF+Str, HY vs YC)
# ============================================================================
library(dplyr)
library(tidyr)
library(agricolae)
library(ggiraphExtra)
library(scales)
library(ggplot2)
library(patchwork)

out_root <- "C:/Users/ADMIN/Desktop/Co-inoculation/Fig.1"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
setwd(out_root)

# ---------- 1. Data preparation ----------
df <- read.csv("data/maize_phenotype_summary.csv", header = TRUE)
df <- df[grepl("Topsoil", df$Niche), ]
df <- df %>% filter(Treatment %in% c("AMF", "Str", "AMF+Str"))
df$Treatment <- factor(df$Treatment, levels = c("AMF", "Str", "AMF+Str"))

var_units <- c(
  "Yield_t_ha"            = "Yield (t/ha)",
  "Biomass_t_ha"          = "Biomass (t/ha)",
  "WUE..kg.mm.ha."        = "WUE (kg/mm/ha)",
  "NUE."                  = "NUE (%)",
  "PUE_kg_mm"             = "PUE (kg/mm/ha)",
  "Grain_N_perc"          = "Grain N (%)",
  "Straw_N_perc"          = "Straw N (%)",
  "Crude.fat.content"     = "Crude fat (%)",
  "Starch.content"        = "Starch (%)",
  "Crude.protein.content" = "Crude protein (%)"
)
vars <- names(var_units)

colors <- c("#F0C986", "#3FBDA7", "#F08080")
names(colors) <- levels(df$Treatment)

# ---------- 2. Multiple comparison ----------
compute_letters <- function(data, variable) {
  formula <- as.formula(paste(variable, "~ Treatment"))
  model <- aov(formula, data = data)
  hsd <- HSD.test(model, "Treatment", group = TRUE, alpha = 0.05)
  groups <- hsd$groups
  letter_col <- colnames(groups)[2]
  means <- hsd$means
  res <- data.frame(
    Treatment = rownames(groups),
    mean = groups[[1]],
    groups = groups[[letter_col]],
    se = means[rownames(groups), "se"],
    stringsAsFactors = FALSE
  )
  rownames(res) <- NULL
  return(res)
}

plot_data <- list()
for (site in c("HY", "YC")) {
  site_data <- filter(df, Site == site)
  for (var in vars) {
    stats <- compute_letters(site_data, var)
    stats$Variable <- var
    stats$Site <- site
    plot_data[[length(plot_data) + 1]] <- stats
  }
}
plot_data <- bind_rows(plot_data)
plot_data$Variable <- factor(plot_data$Variable, levels = vars)

# ---------- 3. Prepare radar plot data ----------
prepare_radar_data <- function(site_name) {
  site_data <- plot_data %>% filter(Site == site_name)
  wide <- site_data %>%
    dplyr::select(Treatment, Variable, mean) %>%
    pivot_wider(names_from = Variable, values_from = mean)

  mat <- as.matrix(wide[, -1])
  rownames(mat) <- wide$Treatment

  scale_to_max <- function(x) {
    if (max(x) == 0) return(rep(0, length(x)))
    x / max(x)
  }
  mat_scaled <- apply(mat, 2, scale_to_max)

  df_radar <- as.data.frame(mat_scaled)
  df_radar$Treatment <- rownames(df_radar)
  rownames(df_radar) <- NULL

  sig_check <- site_data %>%
    group_by(Variable) %>%
    summarise(sig = !(length(unique(groups)) == 1)) %>%
    ungroup() %>%
    mutate(Variable = as.character(Variable))
  sig_vec <- setNames(sig_check$sig, sig_check$Variable)

  new_names <- sapply(colnames(df_radar), function(v) {
    if (v == "Treatment") return("Treatment")
    label <- var_units[v]
    if (!is.na(sig_vec[v]) && sig_vec[v]) label <- paste0(label, " *")
    return(label)
  })
  colnames(df_radar) <- new_names
  return(df_radar)
}

# ---------- 4. Plot overlapping radar ----------
plot_radar_overlay <- function(site_name) {
  df_radar <- prepare_radar_data(site_name)
  df_radar$Treatment <- factor(df_radar$Treatment, levels = c("AMF", "Str", "AMF+Str"))
  treat_levels <- levels(df_radar$Treatment)
  col <- colors[treat_levels]

  p <- ggRadar(data = df_radar,
               aes(color = Treatment, fill = Treatment),
               rescale = FALSE,
               legend.position = "bottom",
               size = 0.4, alpha = 0.4)

  p <- p +
    scale_fill_manual(values = col) +
    scale_color_manual(values = col) +
    scale_y_continuous(limits = c(0, 1), breaks = 1) +
    theme_bw(base_size = 5) +
    theme(
      panel.background = element_blank(),
      panel.border = element_rect(linewidth = 0.2, color = "black"),
      panel.grid.major = element_line(linewidth = 0.5),
      axis.text = element_text(size = 4, color = "grey30"),
      axis.text.y = element_blank(),
      axis.line = element_blank(),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 6),
      plot.margin = margin(0, 0, 0, 0, "cm")
    ) +
    ggtitle(paste("Site:", site_name)) +
    theme(plot.title = element_text(size = 12, hjust = 0.5, face = "plain"))
  return(p)
}

p_hy <- plot_radar_overlay("HY")
p_yc <- plot_radar_overlay("YC")
combined <- p_hy + p_yc + plot_annotation(tag_levels = NULL)

ggsave(file.path(out_root, "Fig.1B_maize_radar.pdf"), combined, width = 5, height = 2.5)
ggsave(file.path(out_root, "Fig.1B_maize_radar.png"), combined, width = 5, height = 2.5, dpi = 300)
cat("Fig.1B saved\n")
