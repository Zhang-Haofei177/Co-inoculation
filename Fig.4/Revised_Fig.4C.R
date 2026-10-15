# ============================================================================
# Fig.4C: Soil enzyme activity — direct linear model contrasts (fixed)
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(emmeans)
})

out_root <- "C:/Users/ADMIN/Desktop/Co-inoculation/Fig.4"
setwd(out_root)

# ---------- 1. Read data ----------
soil  <- read.csv(file.path(out_root, "data/soil_properties_topsoil.csv"),
                  row.names = 1, stringsAsFactors = FALSE)
group <- read.csv(file.path(out_root, "data/group.csv"), row.names = 1)

group$Site <- factor(group$Site, levels = c("HY", "YC"))
group$Treatment <- factor(group$Treatment,
                          levels = c("Ctrl", "AMF", "Str", "AMF+Str"))
group$Niche <- trimws(group$Niche)
group$Niche[group$Niche == "Top Soil"]   <- "Topsoil"
group$Niche[group$Niche == "Bottom Soil"] <- "Bottom soil"

samples_use <- rownames(group)[group$Niche == "Topsoil" &
                                 group$Treatment %in% c("AMF", "Str", "AMF+Str")]
common <- intersect(samples_use, rownames(soil))

soil_sub <- soil[common, ]
meta_sub <- group[common, ]
meta_sub$Treatment <- droplevels(meta_sub$Treatment)

# 确保 Treatment 顺序：AMF, Str, AMF+Str
stopifnot(identical(levels(meta_sub$Treatment),
                    c("AMF", "Str", "AMF+Str")))

enzymes <- c("NAG", "BG", "LAP", "AKP", "Urease")
enzymes <- intersect(enzymes, colnames(soil_sub))
soil_sub[enzymes] <- lapply(soil_sub[enzymes], as.numeric)

# ---------- 2. Per-enzyme analysis ----------
run_enzyme_lm <- function(df, meta, enzyme) {
  y <- df[[enzyme]]
  if (any(y <= 0, na.rm = TRUE)) return(NULL)
  y <- log2(y)
  
  work <- data.frame(Y = y, Site = meta$Site, Treatment = meta$Treatment)
  work$Treatment <- factor(work$Treatment, levels = c("AMF", "Str", "AMF+Str"))
  
  m <- lm(Y ~ Site * Treatment, data = work)
  
  # 2a. 地点内共接种效应（显式按名称指定系数）
  emm_by_site <- emmeans(m, ~ Treatment | Site)
  site_eff <- contrast(emm_by_site,
                       list(Dual_vs_Single = c(AMF = -0.5, Str = -0.5, `AMF+Str` = 1)))
  site_eff_df <- as.data.frame(summary(site_eff, infer = c(TRUE, TRUE)))
  site_eff_df$Enzyme <- enzyme
  
  # 2b. YC vs HY 交互对比（按名称查格子，不受水平顺序影响）
  emm_grid <- emmeans(m, ~ Site * Treatment)
  grid <- emm_grid@grid
  idx <- function(s, t) which(grid$Site == s & grid$Treatment == t)
  
  cont <- rep(0, nrow(grid))
  cont[idx("YC", "AMF+Str")] <-  1
  cont[idx("YC", "AMF")]     <- -0.5
  cont[idx("YC", "Str")]     <- -0.5
  cont[idx("HY", "AMF+Str")] <- -1
  cont[idx("HY", "AMF")]     <-  0.5
  cont[idx("HY", "Str")]     <-  0.5
  
  c_int <- contrast(emm_grid, list(Diff_YC_vs_HY = cont))
  s_int <- as.data.frame(summary(c_int, infer = c(TRUE, TRUE)))
  s_int$Enzyme <- enzyme
  
  list(site = site_eff_df, inter = s_int)
}

res_all <- lapply(enzymes, function(e) run_enzyme_lm(soil_sub, meta_sub, e))
names(res_all) <- enzymes

# ---------- 3. 整理结果 ----------
site_tab <- do.call(rbind, lapply(res_all, `[[`, "site")) %>%
  select(Enzyme, Site, estimate, SE, lower.CL, upper.CL, p.value) %>%
  rename(logFC = estimate, CI_low = lower.CL, CI_high = upper.CL,
         P.Value = p.value)

inter_tab <- do.call(rbind, lapply(res_all, `[[`, "inter")) %>%
  select(Enzyme, estimate, SE, lower.CL, upper.CL, p.value) %>%
  rename(logFC = estimate, CI_low = lower.CL, CI_high = upper.CL,
         P.Value = p.value) %>%
  mutate(Site = "YC vs HY", Contrast = "Diff_YC_vs_HY")

out_tab <- bind_rows(
  site_tab %>% mutate(Contrast = "Dual_vs_Single") %>%
    select(Enzyme, Site, Contrast, logFC, SE, CI_low, CI_high, P.Value),
  inter_tab %>%
    select(Enzyme, Site, Contrast, logFC, SE, CI_low, CI_high, P.Value)
)

write.csv(out_tab, "Fig.4C_enzyme_lm_contrasts_fixed.csv", row.names = FALSE)

cat("\n=== Fig.4C results (fixed) ===\n")
print(out_tab, row.names = FALSE)

# ---------- 4. 绘图：HY vs YC 两个柱子 + 地点差异星号 ----------
plot_stats <- site_tab %>%
  left_join(inter_tab %>% select(Enzyme, P.Value_inter = P.Value),
            by = "Enzyme") %>%
  mutate(
    Stars = ifelse(P.Value_inter < 0.001, "***",
                   ifelse(P.Value_inter < 0.01, "**",
                          ifelse(P.Value_inter < 0.05, "*", "ns"))),
    Enzyme = factor(Enzyme, levels = c("NAG", "BG", "LAP", "AKP", "Urease")),
    Site   = factor(Site, levels = c("HY", "YC"))
  ) %>%
  group_by(Enzyme) %>%
  mutate(y_star = max(logFC + SE, na.rm = TRUE) * 1.15 + 0.05) %>%
  ungroup()

site_colors <- c("HY" = "#1F78B4", "YC" = "#E31A23")

p_enzyme <- ggplot(plot_stats, aes(x = Enzyme, y = logFC, fill = Site)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7,
           color = "black", linewidth = 0.4) +
  geom_errorbar(aes(ymin = logFC - SE, ymax = logFC + SE),
                position = position_dodge(width = 0.8),
                width = 0.2, linewidth = 0.6) +
  geom_text(aes(x = Enzyme, y = y_star, label = Stars),
            position = position_dodge(width = 0.8),
            size = 5, vjust = 0, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = site_colors, name = "Site") +
  labs(x = NULL,
       y = expression(log[2]*" fold change (AMF+Str / pooled single inoculant)")) +
  theme_classic(base_size = 13) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 11),
    axis.title = element_text(color = "black", size = 12),
    legend.position = "top",
    legend.title = element_blank()
  )
p_enzyme
ggsave("Revised_Fig.4C_enzyme_lm_barplot.pdf", p_enzyme, width = 5, height = 4.5)
ggsave("Revised_Fig.4C_enzyme_lm_barplot.png", p_enzyme, width = 5, height = 4.5, dpi = 300)

cat("\nOutputs:\n")
cat("  Fig.4C_enzyme_lm_contrasts.csv\n")
cat("  Fig.4C_enzyme_lm_barplot.pdf / .png\n")
