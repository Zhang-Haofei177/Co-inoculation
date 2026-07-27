# ============================================================================
# Fig.S3: Topsoil physicochemical + enzyme activity bar chart (12 indicators x 3 treatments x HY/YC, multiple comparison letters)
# ============================================================================
library(ggplot2)
library(dplyr)
library(agricolae)
library(patchwork)

out_root <- "C:/Users/ADMIN/Desktop/Co-inoculation/Fig.S3"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
setwd(out_root)

# ---------- 1. Data preparation ----------
df <- read.csv("data/soil_properties_topsoil.csv", header = TRUE)
df <- df[grepl("Topsoil", df$Niche), ]
df <- df[grepl("AMF|Str|AMF+Str", df$Treatment), ]
df$Treatment <- factor(df$Treatment, levels = c("AMF", "Str", "AMF+Str"))
df$X52um <- NULL; df$X52.250um <- NULL; df$X250um <- NULL

var_units <- c(
  "TN"       = "Total nitrogen (mg/kg)",
  "SOC"      = "Soil organic carbon (g/kg)",
  "AP"       = "Available phosphorus (mg/kg)",
  "TP"       = "Total phosphorus (mg/kg)",
  "NH4"      = "Ammonium nitrogen (mg/kg)",
  "NO3"      = "Nitrate nitrogen (mg/kg)",
  "Moisture" = "Soil moisture (%)",
  "NAG"      = "NAG activity (nmol/g/h)",
  "BG"       = "BG activity (nmol/g/h)",
  "LAP"      = "LAP activity (nmol/g/h)",
  "AKP"      = "AKP activity (nmol/g/h)",
  "Urease"   = "Urease activity (μg/g/h)"
)
vars <- names(var_units)

# ---------- 2. Multiple comparisons ----------
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

result_list <- list()
for (site in c("HY", "YC")) {
  site_data <- filter(df, Site == site)
  for (var in vars) {
    stats <- compute_letters(site_data, var)
    stats$Variable <- var
    stats$Site <- site
    result_list[[length(result_list) + 1]] <- stats
  }
}
plot_data <- bind_rows(result_list)
plot_data$Variable <- factor(plot_data$Variable, levels = vars)

plot_data <- plot_data %>%
  group_by(Site, Variable) %>%
  mutate(label_y = mean + se * 1.2) %>%
  ungroup()

plot_data <- plot_data %>%
  mutate(Site_Treatment = paste(Site, Treatment, sep = "_")) %>%
  mutate(Site_Treatment = factor(Site_Treatment,
                                 levels = c(paste("HY", levels(df$Treatment), sep = "_"),
                                            paste("YC", levels(df$Treatment), sep = "_"))))

# ---------- 3. Plot function ----------
plot_bar_combined <- function(data, var_name) {
  sub_data <- data %>% filter(Variable == var_name)
  sub_data$Site_Treatment <- factor(sub_data$Site_Treatment, levels = levels(data$Site_Treatment))

  colors <- c("#F0C986", "#3FBDA7", "#F08080")
  names(colors) <- levels(df$Treatment)

  ggplot(sub_data, aes(x = Site_Treatment, y = mean, fill = Treatment)) +
    annotate("rect", xmin = 3.5, xmax = 6.5, ymin = -Inf, ymax = Inf,
             fill = "grey90", alpha = 0.3) +
    geom_vline(xintercept = 3.5, linetype = "solid", color = "black", linewidth = 0.15) +
    geom_col(color = "black", width = 0.6, linewidth = 0.15) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.15, linewidth = 0.15) +
    geom_text(aes(y = label_y, label = groups),
              vjust = -0.3, size = 1.8, fontface = "plain", color = "black", family = "sans") +
    scale_fill_manual(values = colors) +
    scale_x_discrete(labels = c("AMF","Str","AMF+Str", "AMF","Str","AMF+Str")) +
    labs(x = NULL, y = var_units[var_name]) +
    theme_classic() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.15),
      axis.line = element_blank(),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "plain", color = "black", size = 6),
      axis.text.y = element_text(color = "black", size = 6),
      axis.title.y = element_text(size = 6, color = "black"),
      axis.title.x = element_blank(),
      axis.ticks = element_line(linewidth = 0.1),
      plot.title = element_blank(),
      legend.position = "none",
      text = element_text(family = "sans", size = 6),
      plot.margin = margin(2, 2, 2, 2, unit = "pt")
    )
}

# ---------- 4. Combine ----------
plot_list <- lapply(vars, function(var) plot_bar_combined(plot_data, var))
combined <- wrap_plots(plot_list, ncol = 6, nrow = 2) +
  plot_annotation(title = "HY and YC Combined",
                  theme = theme(plot.title = element_text(size = 8, hjust = 0.5, face = "plain", family = "sans")))

ggsave(file.path(out_root, "Fig.S3_soil_bar.pdf"), combined, width = 180, height = 90, units = "mm")
ggsave(file.path(out_root, "Fig.S3_soil_bar.png"), combined, width = 180, height = 90, units = "mm", dpi = 300)

cat("\n========== Data means and significance letters ==========\n")
for (site in c("HY", "YC")) {
  cat(sprintf("\n--- Site: %s ---\n", site))
  for (var in vars) {
    sub_data <- plot_data %>%
      filter(Variable == var, Site == site) %>%
      dplyr::select(Treatment, mean, se, groups) %>%
      arrange(Treatment)
    cat(sprintf("\nVariable: %s (%s)\n", var, var_units[var]))
    print(as.data.frame(sub_data), row.names = FALSE)
  }
}
cat("\nFig.S3 saved\n")
