# ============================================================================
# Fungi bin data filtering: extract synergy comparisons from BinContributeToProcess
# (adapted from Bacteria version)
# ============================================================================
library(dplyr)
library(tidyr)

# ---------- 0. Paths ----------
wd <- "C:/Users/ADMIN/Desktop/Co-inoculation/Fig.2/iCAMP/Output_simple"
setwd(wd)

# ---------- 1. Read raw iCAMP bin contribution ----------
df <- read.csv("Fungi_treatment.BinContributeToProcess_EachGroup.csv", header = TRUE)

# ---------- 2. Extract dual vs. single inoculation within the same niche ----------
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
      grepl("\\.S\\.", first_part) ~ "Subsoil",
      grepl("\\.B\\.", first_part) ~ "Bottom soil",
      grepl("\\.R\\.", first_part) ~ "Rhizosphere",
      grepl("\\.E\\.", first_part) ~ "Endosphere",
      TRUE ~ "Other"
    )
  ) %>%
  filter(Site != "Other", Niche != "Other") %>%
  select(-first_part)

cat("\n===== Row count per niche =====\n")
print(table(df_vs$Niche))

# Save filtered results
write.csv(df_vs, "Fungi_bin_filtered_synergy_process.csv", row.names = FALSE)
cat("Filtering complete:", nrow(df_vs), "rows\n")
