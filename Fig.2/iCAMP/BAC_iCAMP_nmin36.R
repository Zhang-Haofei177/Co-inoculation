# ============================================================
# Reference script: iCAMP analysis (Bacteria, nmin=36)
# ⚠ This script cannot be rerun locally (requires server: 40 threads / 200 GB memory)
# Provided as analysis workflow documentation; results have been copied to Output_simple/
# ============================================================

# 0. Load packages ------------------------------------------------------------
library(iCAMP)
library(ape)
library(ggplot2)

# 1. Paths (rewritten from server paths to local paths) ------------------------
wd <- "C:/Users/ADMIN/Desktop/Co-inoculation/Fig.2/iCAMP/data"
com.file <- file.path(wd, "BAC_rarefied.txt")
tree.file <- file.path(wd, "BAC_tree.nwk")
treat.file <- file.path(wd, "group.txt")
clas.file <- file.path(wd, "BAC_Tax.txt")
save.wd <- "C:/Users/ADMIN/Desktop/Co-inoculation/Fig.2/iCAMP/Output_simple"
if(!dir.exists(save.wd)) dir.create(save.wd, recursive = TRUE)

prefix <- "Bacteria_treatment"
rand.time <- 999
nworker <- 40
memory.G <- 200
bin.size.limit <- 36

# 2. Read data ----------------------------------------------------------------
setwd(wd)

comm <- t(read.table(com.file, header = TRUE, sep = "\t", row.names = 1,
                     as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                     check.names = FALSE))

tree <- read.tree(tree.file)

treat <- read.table(treat.file, header = TRUE, sep = "\t", row.names = 1,
                    as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                    check.names = FALSE)

if (file.exists(clas.file)) {
  clas <- read.table(clas.file, header = TRUE, sep = "\t", row.names = 1,
                     as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
                     check.names = FALSE)
} else {
  clas <- NULL
  cat("Classification file not found; skipping clas-related output\n")
}

# 3. Data alignment and filtering ---------------------------------------------
samp_match <- match.name(rn.list = list(comm = comm, treat = treat))
comm <- samp_match$comm
treat <- samp_match$treat

comm <- comm[, colSums(comm) > 0, drop = FALSE]

# Abundance filtering: total abundance >= 10
asv_total <- colSums(comm)
comm <- comm[, asv_total >= 10, drop = FALSE]
cat("Species remaining after abundance filtering:", ncol(comm), "\n")

# Occurrence filtering: present in at least 5% of samples
occurrence <- colSums(comm > 0)
min_occ <- ceiling(0.05 * nrow(comm))
comm <- comm[, occurrence >= min_occ, drop = FALSE]
cat("Species remaining after frequency filtering:", ncol(comm), "\n")

# Clean single quotes from tree labels and match with community data
tree$tip.label <- gsub("'", "", tree$tip.label)
common_sp <- intersect(colnames(comm), tree$tip.label)
comm <- comm[, common_sp, drop = FALSE]
tree <- keep.tip(tree, common_sp)

cat("Final data:", nrow(comm), "samples x", ncol(comm), "species\n")
cat("Number of tips in tree:", length(tree$tip.label), "\n")

# 4. Prepare grouping vector --------------------------------------------------
group_vec <- treat[rownames(comm), "treatment"]
group_vec <- as.factor(group_vec)

# 5. Run iCAMP core analysis --------------------------------------------------
# ⚠ The following steps require server execution; local runs will fail due to insufficient memory
cat("\n===== Running icamp.big (nmin=36) =====\n")
icres <- icamp.big(
  comm = comm, tree = tree,
  bin.size.limit = bin.size.limit,
  rand = rand.time, nworker = nworker, memory.G = memory.G,
  detail.save = TRUE, qp.save = FALSE, rtree.save = FALSE
)

saveRDS(icres, file = file.path(save.wd, paste0(prefix, ".iCAMP.detail.rds")))
save(icres, file = file.path(save.wd, paste0(prefix, ".iCAMP.detail.rda")))

# 6. iCAMP bin-level statistics -----------------------------------------------
cat("\n===== Running icamp.bins =====\n")
icbin <- icamp.bins(
  icamp.detail = icres$detail, treat = treat, clas = clas,
  silent = FALSE, boot = TRUE, rand.time = rand.time, between.group = TRUE
)

save(icbin, file = file.path(save.wd, paste0(prefix, ".iCAMP.Summary.rda")))

write.csv(icbin$Pt, file.path(save.wd, paste0(prefix, ".ProcessImportance_EachGroup.csv")), row.names = FALSE)
write.csv(icbin$Ptk, file.path(save.wd, paste0(prefix, ".ProcessImportance_EachBin_EachGroup.csv")), row.names = FALSE)
write.csv(icbin$Ptuv, file.path(save.wd, paste0(prefix, ".ProcessImportance_EachTurnover.csv")), row.names = FALSE)
write.csv(icbin$BPtk, file.path(save.wd, paste0(prefix, ".BinContributeToProcess_EachGroup.csv")), row.names = FALSE)

if (!is.null(icbin$Class.Bin)) {
  write.csv(data.frame(ID = rownames(icbin$Class.Bin), icbin$Class.Bin, stringsAsFactors = FALSE),
            file.path(save.wd, paste0(prefix, ".Taxon_Bin.csv")), row.names = FALSE)
}
if (!is.null(icbin$Bin.TopClass)) {
  write.csv(icbin$Bin.TopClass, file.path(save.wd, paste0(prefix, ".Bin_TopTaxon.csv")), row.names = FALSE)
}
if (!is.null(icres$detail$taxabin$sp.bin)) {
  write.csv(icres$detail$taxabin$sp.bin, file.path(save.wd, paste0(prefix, ".iCamp.Bin.taxon.csv")), row.names = FALSE)
}

# 7. Bootstrap testing --------------------------------------------------------
cat("\n===== Running icamp.boot =====\n")
treat.use <- treat[rownames(comm), "Treatment", drop = FALSE]
icboot <- icamp.boot(
  icamp.result = icres$CbMPDiCBraya, treat = treat.use,
  rand.time = rand.time, compare = TRUE, silent = FALSE,
  between.group = TRUE, ST.estimation = TRUE
)

save(icboot, file = file.path(save.wd, paste0(prefix, ".iCAMP.Boot.treatment.rda")))
write.csv(icboot$summary, file.path(save.wd, paste0(prefix, ".iCAMP.BootSummary.treatment.csv")), row.names = FALSE)
if (!is.null(icboot$compare)) {
  write.csv(icboot$compare, file.path(save.wd, paste0(prefix, ".iCAMP.Compare.treatment.csv")), row.names = FALSE)
}

cat("\n===== Analysis complete =====\n")
cat("All results saved to:", save.wd, "\n")

# ============================================================
# Input files:
#   BAC_rarefied.txt    - Rarefied ASV abundance table (82M)
#   BAC_tree.nwk        - Phylogenetic tree (5.5M)
#   group.txt           - Grouping information
#   BAC_Tax.txt         - Taxonomic annotation (13M)
#
# Output files (Output_simple/):
#   Bacteria_treatment.iCAMP.detail.rds/.rda  - Full results
#   Bacteria_treatment.iCAMP.Summary.rda       - Bin summary
#   Bacteria_treatment.ProcessImportance_*.csv - Process importance
#   Bacteria_treatment.Taxon_Bin.csv           - Bin taxonomy
#   Bacteria_treatment.iCAMP.Boot*.csv/.rda    - Bootstrap results
#   Bacteria_treatment.iCAMP.Compare.treatment.csv - Between-group comparison
# ============================================================
