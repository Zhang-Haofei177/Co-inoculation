# =============================================================================
# Fig.2A: Bacterial phylum relative abundance (3 niches x HY/YC)
# =============================================================================
library(dplyr)
library(reshape2)
library(ggplot2)

set.seed(123)

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

# ---------- 2. Data loading ----------
otu <- read.csv("data/BAC_Relative.csv", header = TRUE)
tax <- read.csv("data/BAC_Tax.csv", header = TRUE)
group <- read.csv("data/group.csv", header = TRUE)

group <- group[!grepl("Ctrl", group$Treatment), ]
group <- group[!grepl("Subsoil|Bottom soil", group$Niche), ]
otu <- otu[c("X", group$Sample)]

df <- otu
df_rel <- df
df_rel[, 2:ncol(df)] <- 100 * df[, 2:ncol(df)] / colSums(df[, 2:ncol(df)])
otu <- df_rel

# ---------- 3. Summarize by phylum ----------
Order <- merge(otu, tax[c("X","Phylum")], by = "X")
rownames(Order) <- Order$X
Order$X <- NULL

sum_Order <- group_by(Order, Phylum) %>%
  summarise_all(sum)
sum_Order$Phylum[is.na(sum_Order$Phylum) | sum_Order$Phylum == ""] <- "Others"

sum_Order$row_sum <- rowSums(sum_Order[, -1], na.rm = TRUE)
sorted_data <- sum_Order[order(-sum_Order$row_sum), ]
sorted_data$row_sum <- NULL

Ord_20 <- sorted_data[c(1:21), ]
Ord_20 <- subset(Ord_20, Ord_20$Phylum != "Others")
name_20 <- Ord_20$Phylum

Ord_other <- subset(sorted_data, !grepl(paste(name_20, collapse = "|"), Phylum))
Ord_other <- colSums(Ord_other[2:ncol(Ord_other)])
Ord_other <- c("Others", Ord_other)

Ord_20_other <- rbind(Ord_20, Ord_other)
Ord_21_name <- Ord_20_other$Phylum

# ---------- 4. Extract colors ----------
CS_cols <- phylum_colors[Ord_21_name]
CS_cols[is.na(CS_cols)] <- "#D3D3D3"
names(CS_cols) <- Ord_21_name

# ---------- 5. Convert to long format ----------
Ord_long <- melt(Ord_20_other, id.vars = "Phylum", variable.name = "Sample", value.name = "Value")
Ord_long_group <- merge(Ord_long, group, by = "Sample")

Ord_long_group$Treatment <- factor(Ord_long_group$Treatment, levels = c("AMF","Str","AMF+Str"))
Ord_long_group$Niche <- factor(Ord_long_group$Niche, levels = c("Topsoil","Rhizosphere","Endosphere"))
Ord_long_group$Site <- factor(Ord_long_group$Site, levels = c("HY","YC"))
Ord_long_group$Phylum <- factor(Ord_long_group$Phylum, levels = Ord_21_name)
Ord_long_group$Value <- as.numeric(Ord_long_group$Value)

Ord_long_group <- Ord_long_group %>%
  mutate(
    Niche_short = case_when(
      Niche == "Topsoil"    ~ "Top",
      Niche == "Rhizosphere"~ "Rhi",
      Niche == "Endosphere" ~ "End"
    ),
    Treatment_Niche = paste(Niche_short, Treatment, sep = "_")
  )

treatment_niche_order <- c("Top_AMF", "Top_Str", "Top_AMF+Str",
                           "Rhi_AMF", "Rhi_Str", "Rhi_AMF+Str",
                           "End_AMF", "End_Str", "End_AMF+Str")
Ord_long_group$Treatment_Niche <- factor(Ord_long_group$Treatment_Niche, levels = treatment_niche_order)
Ord_long_group$Phylum <- factor(Ord_long_group$Phylum, levels = Ord_21_name)

# ---------- 6. Plot ----------
p <- ggplot(Ord_long_group, aes(x = Treatment_Niche, y = Value/4, fill = Phylum)) +
  annotate("rect", xmin = 3.5, xmax = 6.5, ymin = -Inf, ymax = Inf,
           fill = "grey90", alpha = 0.5) +
  geom_vline(xintercept = 3.5, linetype = "solid", color = "black", linewidth = 0.2) +
  geom_vline(xintercept = 6.5, linetype = "solid", color = "black", linewidth = 0.2) +
  geom_col(position = "stack", width = 0.6) +
  facet_grid(~Site, scales = "free_x") +
  scale_fill_manual(values = CS_cols) +
  scale_x_discrete(labels = rep(c("AMF", "Str", "AMF+Str"), times = 3)) +
  labs(x = NULL, y = 'Bacteria Relative Abundance(%)') +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(color = 'black', fill = 'transparent'),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
    axis.ticks = element_line(linewidth = 0.2),
    strip.background = element_rect(fill = "white", color = "black", size = 0.2),
    strip.text = element_text(family = "sans", size = 6, colour = "black"),
    axis.text = element_text(family = "sans", size = 6, colour = "black"),
    axis.title = element_text(family = "sans", size = 6, colour = "black"),
    legend.text = element_text(family = "sans", size = 5, face = "plain", colour = "black"),
    legend.position = "bottom",
    legend.key.size = unit(0.5, "lines"),
    legend.key = element_rect(fill = NA, linewidth = 0.2),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

ggsave(file.path(out_root, "Fig.2A_bacteria_phylum.pdf"), p, width = 90, height = 75, units = "mm")
ggsave(file.path(out_root, "Fig.2A_bacteria_phylum.png"), p, width = 90, height = 75, units = "mm", dpi = 300)
cat("Fig.2A saved\n")
