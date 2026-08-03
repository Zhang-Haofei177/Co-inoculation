# ============================================================================
# Fig.1C: Topsoil aggregate size distribution (3 treatments x HY/YC)
# ============================================================================
library(ggplot2)
library(dplyr)
library(tidyr)

out_root <- "C:/Users/ADMIN/Desktop/Co-inoculation/Fig.1"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
setwd(out_root)

# ---------- 1. Read data ----------
df <- read.csv("data/soil_properties_topsoil.csv", header = TRUE, row.names = 1)
df <- df[grepl("T1|T2|T6", rownames(df)), c(2:5, 18:20)]
df <- df[grepl("Topsoil", df$Niche), ]
df$Site <- factor(df$Site, levels = c("HY", "YC"), labels = c("Heyang", "Yvci"))

# ---------- 2. Calculate mean of each component ----------
df_summary <- df %>%
  group_by(Group, Niche, Site, Treatment) %>%
  summarise(
    `Clay_Silt` = mean(X52um, na.rm = TRUE),
    `Microaggregates` = mean(X52.250um, na.rm = TRUE),
    `Macroaggregates` = mean(X250um, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(`Clay_Silt`, `Microaggregates`, `Macroaggregates`),
    names_to = "component",
    values_to = "percentage"
  )

df_summary$component <- factor(df_summary$component,
                               levels = c("Clay_Silt", "Microaggregates", "Macroaggregates"))
df_summary$Treatment <- factor(df_summary$Treatment,
                               levels = c("AMF", "Str", "AMF+Str"),
                               labels = c("AMF", "Str", "AMF+Str"))

# ---------- 3. Plot ----------
p <- ggplot(df_summary, aes(x = Treatment, y = percentage, fill = component)) +
  geom_col(position = "stack", width = 0.7) +
  facet_grid(~Site, scales = "free_x", space = "free") +
  labs(x = "", y = "Aggregate size distribution (%)", title = "") +
  scale_fill_manual(values = c("Clay_Silt" = "#FCB255",
                               "Microaggregates" = "#42B382",
                               "Macroaggregates" = "#107BA4"),
                    labels = c("Clay_Silt" = "Clay + silt (≤52 μm)",
                               "Microaggregates" = "Microaggregates (52-250 μm)",
                               "Macroaggregates" = "Macroaggregates (≥250 μm)")) +
  theme_bw(base_size = 6) +
  theme(legend.title = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.text = element_text(color = "black", size = 6),
        axis.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        strip.text = element_text(size = 6),
        strip.background = element_rect(fill = "white"))
p
ggsave(file.path(out_root, "Fig.1C_soil_aggregate.pdf"), p, width = 50, height = 40, units = "mm")
ggsave(file.path(out_root, "Fig.1C_soil_aggregate.png"), p, width = 50, height = 40, units = "mm", dpi = 300)
write.csv(df_summary, file.path(out_root, "Fig.1C_aggregate_data.csv"), row.names = FALSE)
cat("Fig.1C saved\n")
