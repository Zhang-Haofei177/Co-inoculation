# ============================================================================
# Fig.2D: Bacterial iCAMP assembly process stacked bar plot (merged single inocula, 3 niches)
# ============================================================================
library(dplyr)
library(tidyr)
library(ggplot2)

out_root <- "C:/Users/ADMIN/Desktop/Co-inoculation/Fig.2"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
setwd(out_root)

# ---------- 1. Read iCAMP results ----------
df <- read.csv("data/Bacteria_treatment.ProcessImportance_EachGroup.csv", header = TRUE)

# ---------- 2. Extract dual vs single comparisons within the same niche ----------
df_vs <- df %>%
  filter(grepl("_vs_", Group)) %>%
  filter(grepl("T6", Group)) %>%
  filter(grepl("T1", Group) | grepl("T2", Group)) %>%
  mutate(
    left_part  = sub("_vs_.*", "", Group),
    right_part = sub(".*_vs_", "", Group),
    left_site_niche  = sub("\\.T[126]$", "", left_part),
    right_site_niche = sub("\\.T[126]$", "", right_part)
  ) %>%
  filter(left_site_niche == right_site_niche) %>%
  select(-left_part, -right_part, -left_site_niche, -right_site_niche) %>%
  mutate(Contrast = case_when(
    grepl("T1", Group) & grepl("T6", Group) ~ "vs_AMF",
    grepl("T2", Group) & grepl("T6", Group) ~ "vs_Str",
    TRUE ~ "Other"
  )) %>%
  filter(Contrast != "Other")

# ---------- 3. Parse site and niche ----------
df_vs <- df_vs %>%
  mutate(
    first_part = sub("_vs_.*", "", Group),
    Site = case_when(
      substr(first_part, 1, 1) == "H" ~ "HY",
      substr(first_part, 1, 1) == "T" ~ "YC",
      TRUE ~ "Other"
    ),
    Niche = case_when(
      grepl("\\.T\\.", first_part) ~ "Topsoil",
      grepl("\\.R\\.", first_part) ~ "Rhizosphere",
      grepl("\\.E\\.", first_part) ~ "Endosphere",
      TRUE ~ "Other"
    )
  ) %>%
  filter(Site != "Other", Niche != "Other") %>%
  filter(Niche %in% c("Topsoil", "Rhizosphere", "Endosphere")) %>%
  select(-first_part)

# ---------- 4. Convert to long format ----------
process_cols <- c("HeS", "HoS", "DL", "HD", "DR")
df_long <- df_vs %>%
  pivot_longer(cols = all_of(process_cols), names_to = "Process", values_to = "Proportion")

df_long$Process <- factor(df_long$Process, levels = c("HeS", "HoS", "DL", "HD", "DR"))

# ---------- 5. Merge the two single-inoculum contrasts ----------
df_merged <- df_long %>%
  group_by(Site, Niche, Process) %>%
  summarise(Proportion = mean(Proportion, na.rm = TRUE), .groups = "drop")

df_merged <- df_merged %>%
  mutate(
    Niche_short = case_when(
      Niche == "Topsoil"    ~ "Top",
      Niche == "Rhizosphere"~ "Rhi",
      Niche == "Endosphere" ~ "End"
    ),
    Site_Niche = paste(Niche_short, Site, sep = "_")
  ) %>%
  mutate(Site_Niche = factor(Site_Niche,
                             levels = c("Top_HY", "Top_YC", "Rhi_HY", "Rhi_YC", "End_HY", "End_YC")))

# ---------- 6. Plot ----------
process_colors <- c("HeS"="#F9CCCC","HoS"="#FFE4A9","DL"="#D8F3FF","HD"="#8DD3C7","DR"="#E4E8EE")

p <- ggplot(df_merged, aes(x = Site_Niche, y = Proportion, fill = Process)) +
  annotate("rect", xmin = 2.5, xmax = 4.5, ymin = -Inf, ymax = Inf, fill = "#CCD0DD", alpha = 0.5) +
  geom_vline(xintercept = 2.5, linetype = "solid", color = "black", linewidth = 0.2) +
  geom_vline(xintercept = 4.5, linetype = "solid", color = "black", linewidth = 0.2) +
  geom_col(position = "fill", color = "black", linewidth = 0.25, width = 0.7) +
  scale_fill_manual(values = process_colors) +
  scale_x_discrete(labels = c("HY", "YC", "HY", "YC", "HY", "YC")) +
  labs(x = NULL, y = "Proportion of assembly processes",
       title = "Bacterial synergistic assembly processes (merged single inocula)") +
  theme_classic(base_size = 6) +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
    plot.title = element_text(size = 6, hjust = 0.5, face = "plain", color = "black"),
    axis.text = element_text(face = "plain", color = "black", size = 6),
    axis.title = element_text(face = "plain", color = "black", size = 6),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
    strip.text = element_text(face = "plain", size = 6, color = "black"),
    legend.position = "bottom",
    legend.title = element_text(face = "plain", size = 6, color = "black"),
    legend.text = element_text(size = 6, color = "black"),
    text = element_text(family = "sans", size = 6, color = "black"),
    plot.margin = margin(5, 5, 5, 5, unit = "pt")
  )

ggsave(file.path(out_root, "Fig.2D_bacteria_iCAMP.pdf"), p, width = 40, height = 50, units = "mm")
ggsave(file.path(out_root, "Fig.2D_bacteria_iCAMP.png"), p, width = 40, height = 50, units = "mm", dpi = 300)
write.csv(df_merged, file.path(out_root, "Fig.2D_bacteria_iCAMP_data.csv"), row.names = FALSE)
cat("Fig.2D saved\n")
