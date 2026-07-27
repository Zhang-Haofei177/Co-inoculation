# ============================================================================
# Fig.2E: Bacterial phylum-level HoS + DL contribution stacked bar plot (3 niches)
# ============================================================================
library(dplyr)
library(tidyr)
library(ggplot2)

out_root <- "C:/Users/ADMIN/Desktop/Co-inoculation/Fig.2"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
setwd(out_root)

# ---------- 1. Fixed color mapping ----------
phylum_colors <- c(
  "Pseudomonadota"           = "#377EB8",
  "Actinomycetota"           = "#CCEBC5",
  "Gemmatimonadota"          = "#FDB462",
  "Bacillota"                = "#FFED6F",
  "Acidobacteriota"          = "#FB8072",
  "Bacteroidota"             = "#FED9A6",
  "Myxococcota"              = "#B3DE69",
  "Chloroflexota"            = "#BEBADA",
  "Verrucomicrobiota"        = "#80B1D3",
  "Nitrospirota"             = "#A6CEE3",
  "Entotheonellaeota"        = "#8DD3C7",
  "Bdellovibrionota"         = "#FCCDE5",
  "Planctomycetota"          = "#BC80BD",
  "NB1-j"                    = "#FFB6C1",
  "Thermodesulfobacteriota"  = "#00FFFF",
  "Armatimonadota"           = "#9400D3",
  "Candidatus_Kapabacteria"  = "#228B22",
  "Methylomirabilota"        = "#FF6347",
  "Abditibacteriota"         = "#00FF7F",
  "Latescibacterota"         = "#EEE8AA",
  "Others"                   = "#D3D3D3"
)

# ---------- 2. Read data ----------
df_vs <- read.csv("data/Bacteria_bin_filtered_synergy.csv")
df_vs <- df_vs[!grepl("S|B", df_vs$Group), ]
bin_tax <- read.csv("data/Bacteria_treatment.Taxon_Bin.csv", header = TRUE)

# ---------- 3. Preprocessing ----------
bin_cols <- grep("^bin[0-9]+$", colnames(df_vs), value = TRUE)
bin_tax$Bin_num <- gsub("Bin", "", bin_tax$Bin)
bin_tax$Bin_num <- as.character(bin_tax$Bin_num)

bin_size <- bin_tax %>%
  group_by(Bin_num) %>%
  summarise(ASV_count = n(), .groups = "drop")

df_long <- df_vs %>%
  pivot_longer(cols = all_of(bin_cols), names_to = "Bin", values_to = "Contribution") %>%
  mutate(Bin_num = gsub("bin", "", Bin)) %>%
  left_join(bin_size, by = "Bin_num") %>%
  mutate(ASV_contrib = Contribution / ASV_count)

df_asv <- df_long %>%
  left_join(bin_tax %>% select(ID, Bin_num, Phylum), by = "Bin_num", relationship = "many-to-many") %>%
  filter(!is.na(ID))

df_asv$Phylum[is.na(df_asv$Phylum) | df_asv$Phylum == ""] <- "Unclassified"

# ---------- 4. Aggregate by Site x Niche x Process x Phylum ----------
df_plot <- df_asv %>%
  filter(Process %in% c("HoS", "DL")) %>%
  group_by(Site, Niche, Process, Phylum) %>%
  summarise(TotalContrib = sum(ASV_contrib, na.rm = TRUE), .groups = "drop")

df_plot$Niche <- factor(df_plot$Niche, levels = c("Topsoil", "Rhizosphere", "Endosphere"))
df_plot$Process <- factor(df_plot$Process, levels = c("HoS", "DL"))

# ---------- 5. Top 10 phyla ----------
df_plot_no_unclass <- df_plot %>% filter(Phylum != "Unclassified")
total_contrib_by_phylum <- df_plot_no_unclass %>%
  group_by(Phylum) %>%
  summarise(Total = sum(TotalContrib, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Total))

top10_phyla <- total_contrib_by_phylum$Phylum[1:10]

# ---------- 6. Unify to Top 10 + other ----------
df_plot_clean <- df_plot %>%
  mutate(Phylum_plot = case_when(
    Phylum == "Unclassified" ~ "other",
    Phylum %in% top10_phyla ~ Phylum,
    TRUE ~ "other"
  ))

phyla_order <- c(top10_phyla, "other")
df_plot_clean$Phylum_plot <- factor(df_plot_clean$Phylum_plot, levels = phyla_order)

# ---------- 7. Color mapping ----------
phyla_order_col <- phyla_order
phyla_order_col[phyla_order_col == "other"] <- "Others"
plot_colors <- phylum_colors[phyla_order_col]
plot_colors[is.na(plot_colors)] <- "#D3D3D3"
names(plot_colors) <- phyla_order

# ---------- 8. Site_Niche ----------
df_plot_clean <- df_plot_clean %>%
  mutate(
    Niche_short = case_when(
      Niche == "Topsoil"    ~ "Top",
      Niche == "Rhizosphere"~ "Rhi",
      Niche == "Endosphere" ~ "End"
    ),
    Site_Niche = paste(Niche_short, Site, sep = "_")
  ) %>%
  mutate(Site_Niche = factor(Site_Niche,
                             levels = c("Top_HY","Top_YC","Rhi_HY","Rhi_YC","End_HY","End_YC")))

# ---------- 9. Plot ----------
p <- ggplot(df_plot_clean, aes(x = Site_Niche, y = TotalContrib/2, fill = Phylum_plot)) +
  annotate("rect", xmin = 2.5, xmax = 4.5, ymin = -Inf, ymax = Inf, fill = "grey90", alpha = 0.5) +
  geom_vline(xintercept = 2.5, linetype = "solid", color = "black", linewidth = 0.2) +
  geom_vline(xintercept = 4.5, linetype = "solid", color = "black", linewidth = 0.2) +
  geom_col(position = "stack", width = 0.6, linewidth = 0.2) +
  scale_fill_manual(values = plot_colors, name = "Phylum") +
  facet_grid(~Process) +
  scale_x_discrete(labels = c("HY", "YC", "HY", "YC", "HY", "YC")) +
  labs(x = NULL, y = "Sum of ASV-level contribution",
       title = "Bacterial phyla driving HoS and DL (Top 10 + other)") +
  theme_classic(base_size = 6) +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.2),
    strip.text = element_text(face = "plain", size = 6, color = "black"),
    axis.text = element_text(face = "plain", color = "black", size = 6),
    axis.title = element_text(size = 6, color = "black"),
    legend.text = element_text(size = 5, color = "black"),
    legend.title = element_text(size = 5, color = "black"),
    legend.position = "bottom",
    text = element_text(family = "sans", size = 6, color = "black")
  )

ggsave(file.path(out_root, "Fig.2E_bacteria_phylum_HoS_DL.pdf"), p, width = 70, height = 60, units = "mm")
ggsave(file.path(out_root, "Fig.2E_bacteria_phylum_HoS_DL.png"), p, width = 70, height = 60, units = "mm", dpi = 300)
cat("Fig.2E saved\n")
