# ============================================================================
# Fig.1A: Weather data (temperature + precipitation) + growth stage bars
# Fig.1D: Stacked bar chart of total precipitation by growth stage
# ============================================================================

Sys.setlocale("LC_TIME", "C")

library(ggplot2)
library(dplyr)
library(lubridate)
library(patchwork)

out_root <- "C:/Users/ADMIN/Desktop/Co-inoculation/Fig.1"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
setwd(out_root)

# ============================================================================
# Fig.1A: Combined weather plot with growth stage bars
# ============================================================================
cat("\n========== Fig.1A ==========\n")

hy <- read.csv("data/2025_HY.csv", stringsAsFactors = FALSE)
yc <- read.csv("data/2025_YC.csv", stringsAsFactors = FALSE)

hy <- hy %>% mutate(Location = "Heyang", data = ymd(data))
yc <- yc %>% mutate(Location = "Yvci", data = ymd(data))
weather <- bind_rows(hy, yc)

# Precipitation scaling factor
max_temp <- max(weather$MT, na.rm = TRUE)
max_precip <- max(weather$P, na.rm = TRUE)
scale_factor <- max_temp / max_precip * 0.45

colors <- list(
  Heyang_ribbon = "#A6CEE3",
  Heyang_line   = "#1F78B4",
  Heyang_bar    = "#1F78B4",
  Yvci_ribbon   = "#FDBF6F",
  Yvci_line     = "#E31A23",
  Yvci_bar      = "#E31A23"
)

precip_data <- weather %>%
  select(data, Location, P) %>%
  mutate(P_scaled = P * scale_factor)

# Main plot (weather data)
p_main <- ggplot() +
  geom_ribbon(data = filter(weather, Location == "Heyang"),
              aes(x = data, ymin = LT, ymax = MT, fill = "Heyang Range"), alpha = 0.25) +
  geom_ribbon(data = filter(weather, Location == "Yvci"),
              aes(x = data, ymin = LT, ymax = MT, fill = "Yvci Range"), alpha = 0.25) +
  geom_line(data = filter(weather, Location == "Heyang"),
            aes(x = data, y = AT, color = "Heyang Avg. Temp"), linewidth = 0.9) +
  geom_line(data = filter(weather, Location == "Yvci"),
            aes(x = data, y = AT, color = "Yvci Avg. Temp"), linewidth = 0.9) +
  geom_col(data = precip_data,
           aes(x = data, y = P_scaled, fill = Location, color = Location),
           alpha = 0.8, width = 0.8, position = position_dodge(0.9), linewidth = 0.2) +
  scale_y_continuous(
    name = "Temperature (°C)",
    sec.axis = sec_axis(~ . / scale_factor, name = "Precipitation (mm)")
  ) +
  scale_fill_manual(
    name = NULL,
    values = c("Heyang Range" = colors$Heyang_ribbon,
               "Yvci Range"   = colors$Yvci_ribbon,
               "Heyang"       = colors$Heyang_bar,
               "Yvci"         = colors$Yvci_bar),
    labels = c("Heyang Range" = "Heyang temp. range",
               "Yvci Range"   = "Yvci temp. range",
               "Heyang"       = "Heyang precipitation",
               "Yvci"         = "Yvci precipitation")
  ) +
  scale_color_manual(
    name = NULL,
    values = c("Heyang Avg. Temp" = colors$Heyang_line,
               "Yvci Avg. Temp"   = colors$Yvci_line,
               "Heyang"           = colors$Heyang_bar,
               "Yvci"             = colors$Yvci_bar),
    labels = c("Heyang Avg. Temp" = "Heyang avg. temp",
               "Yvci Avg. Temp"   = "Yvci avg. temp",
               "Heyang"           = "Heyang precipitation",
               "Yvci"             = "Yvci precipitation")
  ) +
  labs(x = NULL, y = "Temperature (°C)",
       title = "Daily weather variables during maize growing season") +
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  ) +
  guides(
    fill = guide_legend(order = 1, nrow = 2, override.aes = list(alpha = 0.7)),
    color = guide_legend(order = 2)
  )

# Growth stage intervals
get_stage_intervals <- function(df) {
  df %>% group_by(Stage) %>% summarise(start = min(data), end = max(data), .groups = "drop")
}
hy_stages <- get_stage_intervals(hy)
yc_stages <- get_stage_intervals(yc)

stage_colors <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69")
hy_stages$color <- stage_colors[1:nrow(hy_stages)]
yc_stages$color <- stage_colors[1:nrow(yc_stages)]

create_stage_bar <- function(stage_df, site_name, y_min = 0, y_max = 5) {
  ggplot(stage_df) +
    geom_rect(aes(xmin = start, xmax = end, ymin = y_min, ymax = y_max, fill = Stage),
              alpha = 0.8) +
    geom_text(aes(x = start + (end - start)/2, y = (y_min + y_max)/2, label = Stage),
              size = 3.5, color = "black") +
    scale_fill_manual(values = setNames(stage_df$color, stage_df$Stage)) +
    labs(x = NULL, y = NULL, title = paste(site_name, "Growth stages")) +
    theme_void() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold", size = 10))
}

p_hy_stage <- create_stage_bar(hy_stages, "Heyang", 0, 5)
p_yc_stage <- create_stage_bar(yc_stages, "Yvci", 0, 5)

date_range <- range(weather$data)
p_main <- p_main + xlim(date_range)
p_hy_stage <- p_hy_stage + xlim(date_range) + theme(axis.text.x = element_blank())
p_yc_stage <- p_yc_stage + xlim(date_range) + theme(axis.text.x = element_blank())

combined_plot <- (p_hy_stage / p_yc_stage / p_main) +
  plot_annotation(title = NULL) &
  theme(plot.margin = margin(5, 5, 5, 5))

ggsave(file.path(out_root, "Fig.1A_weather_combined.pdf"), combined_plot, width = 9, height = 8)
ggsave(file.path(out_root, "Fig.1A_weather_combined.png"), combined_plot, width = 9, height = 8, dpi = 300)
cat("Fig.1A saved\n")

# ============================================================================
# Fig.1D: Stacked bar chart of total precipitation by growth stage
# ============================================================================
cat("\n========== Fig.1D ==========\n")

stage_order <- c("VE_V5", "V6_V12", "VT_R5", "R6")

hy_precip <- hy %>%
  group_by(Stage) %>%
  summarise(Precip = sum(P, na.rm = TRUE), .groups = "drop") %>%
  mutate(Location = "HY")

yc_precip <- yc %>%
  group_by(Stage) %>%
  summarise(Precip = sum(P, na.rm = TRUE), .groups = "drop") %>%
  mutate(Location = "YC")

precip_data <- bind_rows(hy_precip, yc_precip) %>%
  mutate(Stage = factor(Stage, levels = stage_order))

total_precip <- precip_data %>%
  group_by(Location) %>%
  summarise(Total = sum(Precip), .groups = "drop")

stage_colors_p <- c("#E8F4B4", "#78CCB0", "#FFDD94", "#F27A49")
names(stage_colors_p) <- stage_order

p <- ggplot(precip_data, aes(x = Location, y = Precip, fill = Stage)) +
  geom_col(position = "stack", width = 0.6, color = "black", linewidth = 0.2) +
  geom_text(data = total_precip,
            aes(x = Location, y = Total + 5, label = paste0(round(Total, 1), " mm")),
            inherit.aes = FALSE, size = 4, vjust = 0) +
  scale_fill_manual(values = stage_colors_p) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Location", y = "Total precipitation (mm)",
       title = "Total precipitation by growth stage") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(out_root, "Fig.1D_precipitation_by_stage.pdf"), p, width = 3.5, height = 4)
ggsave(file.path(out_root, "Fig.1D_precipitation_by_stage.png"), p, width = 3.5, height = 4, dpi = 300)
cat("Fig.1D saved\n")

cat("\n========== Fig.1A + Fig.1D done ==========\n")
