library(BiocManager)
library(readr)
library(R.utils)
library(WGCNA)
library(limma)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(patchwork)
library(openxlsx)
library(edgeR)
library(ComplexHeatmap)
library(circlize)
# library(clusterProfiler)
library(grid)
library(ggfortify)
library(EnhancedVolcano)
library(tidyHeatmap)
library(fgsea)
# library(enrichplot)
library(ggnewscale)
library(affy)
library(dplyr)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

###############################################################################
# === CONFIGURATION ===
###############################################################################

# ---- File paths ----
setwd("~/OSD152")
getwd()
counts_path <- "~/OSD152/GLDS-152_array_raw_intensities_GLmicroarray.csv"
meta_path   <- "~/OSD152/OSD152_metadata.csv"
out_dir_p2  <- "~/OSD152/OSD152_WGCNA"   # Phase 2 outputs
ws_fig <- file.path(out_dir_p2, "figures")

project_name <- "OSD152"

# # ---- Metadata column names ----
sample_id_col <- "Sample Name"
# 
# # ---- VST parameters ----
# vst_blind <- TRUE

# ---- Filtering threshold ----
counts_filter_threshold <- 10   # rowSums(counts) > threshold

# ---- PCA coloring variables ----
# pca_color_vars <- c("diagnosis", "age", "sex", "pmi", "rin", "hemisphere")

# # ---- WGCNA parameters (USER-SET) ----
# soft_power       <- 14       # USER DECISION from pickSoftThreshold diagnostics
# network_type     <- "signed"
# tom_type         <- "signed"
# cor_type         <- "bicor"
# max_block_size   <- 35000    # USER SET: all 33,373 genes in one block
# deep_split       <- 3
# min_module_size  <- 25
# merge_cut_height <- 0.2



###############################################################################
# END CONFIGURATION
###############################################################################

# dir.create(out_dir_p1, showWarnings = FALSE, recursive = TRUE)
# dir.create(fig_dir_p1, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_p2, showWarnings = FALSE, recursive = TRUE)
# dir.create(fig_dir_p2, showWarnings = FALSE, recursive = TRUE)
dir.create(ws_fig,  showWarnings = FALSE, recursive = TRUE)

log_lines <- c()
logit <- function(msg) { cat(msg, "\n"); log_lines <<- c(log_lines, msg) }


# ###############################################################################
# # PHASE 1: QC + PCA
# ###############################################################################
# out_dir <- out_dir_p1
# fig_dir <- fig_dir_p1

# STEP 1: Load raw counts (xlsx via openxlsx)
###############################################################################
counts_raw <- read.csv(counts_path, header = TRUE, sep = ",", check.names = FALSE)
norm_counts <- read.csv("GLDS-152_array_normalized_expression_GLmicroarray.csv", header = TRUE, sep = ",", check.names = FALSE)
temp.norm <- norm_counts %>%
  filter(!is.na(SYMBOL), SYMBOL != "")
colnames(temp.norm)
temp.norm <- temp.norm[, -c(1, 3:10)]
rowGroup <- as.character(temp.norm$SYMBOL)
rownames(temp.norm) <- make.unique(temp.norm$SYMBOL)
temp.norm <- temp.norm[, -1]


# gene_id_col <- colnames(counts_raw)[1]
# logit("Step 1 - Loading raw counts from xlsx")
# logit(sprintf("  First column: '%s'", gene_id_col))
# 
# rownames(counts_raw) <- counts_raw[[gene_id_col]]
# counts_raw[[gene_id_col]] <- NULL
# counts <- as.matrix(counts_raw)
# mode(counts) <- "numeric"
# 
# logit(sprintf("  Raw counts: %d genes x %d samples", nrow(counts), ncol(counts)))
# logit(sprintf("  Value range: %.1f - %.1f, all integer: %s",
#               min(counts), max(counts), all(counts == floor(counts), na.rm = TRUE)))
# 
# if (any(counts != floor(counts), na.rm = TRUE)) {
#   counts <- round(counts)
#   logit("  Rounded to integers for DESeq2")
# }

# ###############################################################################
# # STEP 2: Duplicate Entrez ID check (before mapping)
# ###############################################################################
# entrez_ids <- rownames(counts)
# n_total <- length(entrez_ids)
# n_unique <- length(unique(entrez_ids))
# n_dup <- n_total - n_unique
# logit("")
# logit("Step 2 - Duplicate Entrez ID check (before mapping)")
# logit(sprintf("  Total rows: %d, Unique: %d, Duplicates: %d", n_total, n_unique, n_dup))
# if (n_dup > 0) {
#   dup_table <- table(entrez_ids)
#   dup_ids <- names(dup_table[dup_table > 1])
#   logit(sprintf("  Duplicated Entrez IDs: %d", length(dup_ids)))
# }
# 
# ###############################################################################
# # STEP 3: Entrez -> SYMBOL mapping (org.Hs.eg.db)
# ###############################################################################
# symbols <- mapIds(org.Hs.eg.db, keys = as.character(entrez_ids),
#                   column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
# mapped_mask <- !is.na(symbols)
# symbols[!mapped_mask] <- entrez_ids[!mapped_mask]  # unmapped retain Entrez ID
# logit("")
# logit("Step 3 - Entrez-to-Symbol mapping")
# logit(sprintf("  Mapped: %d / %d (%.1f%%), Unmapped: %d",
#               sum(mapped_mask), n_total, sum(mapped_mask)/n_total*100,
#               sum(!mapped_mask)))
# 
# n_dup_symbols <- sum(table(symbols) > 1)
# logit(sprintf("  Duplicate symbols after mapping: %d", n_dup_symbols))

###############################################################################
# STEP 4: Check for duplicated gene symbols -> collapseRows if needed
#          STOP if any duplicates remain after collapseRows (no make.unique)
###############################################################################
logit("")
logit("Step 4 - Duplicate symbol handling")

n_dup_symbols <- sum(table(rowGroup) > 1)
logit(sprintf("  Duplicate symbols after mapping: %d", n_dup_symbols))

if (n_dup_symbols > 0) {
  logit(sprintf("  %d duplicate symbols found — running collapseRows (maxRowVariance)", n_dup_symbols))
  temp.norm <- collapseRows(
    datET = temp.norm, rowGroup = rowGroup, rowID = rownames(temp.norm),
    method = "maxRowVariance", connectivityBasedCollapsing = FALSE,
    selectFewestMissing = TRUE
  )
  counts_collapsed <- temp.norm$datETcollapsed
  logit(sprintf("  collapseRows: %d -> %d genes", nrow(counts), nrow(counts_collapsed)))
} else {
  logit("  No duplicate symbols — skipping collapseRows")
  counts_collapsed <- counts
  rownames(counts_collapsed) <- symbols
}

n_dup_remaining <- anyDuplicated(rownames(counts_collapsed))
if (n_dup_remaining > 0) {
  dup_names <- rownames(counts_collapsed)[duplicated(rownames(counts_collapsed))]
  logit("")
  logit(sprintf("  *** STOP: %d duplicate row names remain after collapseRows ***", n_dup_remaining))
  logit(sprintf("  Duplicated symbols: %s", paste(head(dup_names, 20), collapse = ", ")))
  logit("  Investigate these symbols before proceeding.")
  writeLines(log_lines, file.path(out_dir, "workflow_summary.txt"))
  stop("Duplicate row names remain after collapseRows. Stopping per plan.")
}
logit("  No duplicate row names remaining. Row names are SYMBOLS.")
logit(sprintf("  First 10: %s", paste(head(rownames(counts_collapsed), 10), collapse = ", ")))

dim(counts_collapsed)
temp.norm <- counts_collapsed
# counts_collapsed <- as.data.frame(counts_collapsed)
## PCA Check
# temp.count <- log2(1 + temp.count) # temp.count is already in log2 by this point

###############################################################################
# STEP 5: Verify sample ID alignment between metadata and counts
###############################################################################
## Load metadata
OSD152.meta <- read.csv(meta_path, header = TRUE, sep = ",", check.names = FALSE)
temp.meta <- OSD152.meta
temp.meta$'Sample Name' <- sub(" sample$", "_sample", temp.meta$'Sample Name')
colnames(temp.meta)[colnames(temp.meta) == 'Characteristics[Sex]'] <- 'Sex'
colnames(temp.meta)[colnames(temp.meta) == 'Characteristics[age]'] <- 'Age'
rownames(temp.meta) <- temp.meta$'Sample Name'
temp.meta$Dose <- paste0(as.character(temp.meta$Dose), "_Gy")
table(temp.meta$Dose)
temp.meta$Hour <- paste0(as.character(temp.meta$Hour), "Hrs_Gamma")

colnames(temp.meta)
all.equal(colnames(temp.norm), temp.meta$`Sample Name`)

temp.meta$Sex_Dose_Hour <- paste(temp.meta$Sex, temp.meta$Dose, temp.meta$Hour, sep = "_")
temp.meta$Dose_Hour <- paste(temp.meta$Dose, temp.meta$Hour, sep = "_")
table(temp.meta$Dose_Hour)

all.equal(rownames(temp.meta),colnames(temp.norm))
###############################################################################
# STEP 6: Filter genes (rowSums > threshold)
###############################################################################
keep <- rowSums(temp.norm) > counts_filter_threshold
counts_filt <- temp.norm[keep, ]
logit("")
logit(sprintf("Step 6 - Filter (rowSums > %d): %d -> %d genes",
              counts_filter_threshold, nrow(temp.norm), nrow(counts_filt)))

# ###############################################################################
# # STEP 7: Parse metadata covariates
# ###############################################################################
# logit("")
# logit("Step 7 - Parse metadata covariates")
# 
# strip_prefix <- function(x) {
#   x <- as.character(x)
#   sub("^[^:]+:\\s*", "", x)
# }
# 
# meta <- meta_raw
# meta <- meta[match(colnames(counts_filt), meta[[sample_id_col]]), ]
# rownames(meta) <- meta[[sample_id_col]]
# 
# # openxlsx deduplicates column names; use positional indexing for characteristics
# char_cols <- which(colnames(meta) == "Sample_characteristics_ch1")
# logit(sprintf("  Found %d Sample_characteristics_ch1 columns (positional: %s)",
#               length(char_cols), paste(char_cols, collapse = ", ")))
# 
# meta$diagnosis <- strip_prefix(meta[[char_cols[1]]])
# meta$age       <- as.numeric(strip_prefix(meta[[char_cols[2]]]))
# meta$sex       <- strip_prefix(meta[[char_cols[3]]])
# meta$pmi       <- as.numeric(strip_prefix(meta[[char_cols[4]]]))
# meta$rin       <- as.numeric(strip_prefix(meta[[char_cols[5]]]))
# meta$hemisphere <- strip_prefix(meta[[char_cols[6]]])
# meta$symptom_side <- strip_prefix(meta[[char_cols[7]]])
# meta$group     <- strip_prefix(meta[[char_cols[8]]])
# meta$source    <- meta[["Sample_source_name_ch1"]]
# 
# logit(sprintf("  diagnosis: %s", paste(names(table(meta$diagnosis)), "=", as.integer(table(meta$diagnosis)), collapse = ", ")))
# logit(sprintf("  age range: %.0f - %.0f", min(meta$age, na.rm=TRUE), max(meta$age, na.rm=TRUE)))
# logit(sprintf("  sex: %s", paste(names(table(meta$sex)), "=", as.integer(table(meta$sex)), collapse = ", ")))
# logit(sprintf("  pmi range: %.1f - %.1f", min(meta$pmi, na.rm=TRUE), max(meta$pmi, na.rm=TRUE)))
# logit(sprintf("  rin range: %.1f - %.1f", min(meta$rin, na.rm=TRUE), max(meta$rin, na.rm=TRUE)))
# logit(sprintf("  hemisphere: %s", paste(names(table(meta$hemisphere)), "=", as.integer(table(meta$hemisphere)), collapse = ", ")))
# 
# write.csv(meta, file.path(out_dir, paste0(project_name, "_metadata_parsed.csv")), row.names = FALSE)
# logit("  Saved parsed metadata CSV")
# 
# ###############################################################################
# # STEP 8: 3-way PCA (Raw vs log2+1 vs VST)
# ###############################################################################
# logit("")
# logit("Step 8 - 3-way PCA comparison (Raw vs log2+1 vs VST)")
# 
# # Build three expression matrices (genes x samples)
# # (a) RAW counts (truly unnormalized)
# mat_raw <- counts_filt
# logit(sprintf("  Raw:        range %.0f - %.0f", min(mat_raw), max(mat_raw)))
# 
# # (b) log2(raw counts + 1)
# mat_log2 <- log2(counts_filt + 1)
# logit(sprintf("  log2+1:     range %.2f - %.2f", min(mat_log2), max(mat_log2)))
# 
# # (c) VST (DESeq2 variance-stabilizing transformation, blind = TRUE)
# coldata <- data.frame(row.names = colnames(counts_filt))
# coldata$diagnosis <- factor(meta$diagnosis)
# dds <- DESeqDataSetFromMatrix(countData = counts_filt, colData = coldata, design = ~1)
# vsd <- varianceStabilizingTransformation(dds, blind = vst_blind)
# mat_vst <- assay(vsd)
# logit(sprintf("  VST (blind=%s): range %.2f - %.2f", vst_blind, min(mat_vst), max(mat_vst)))
# 
# saveRDS(mat_vst, file.path(out_dir, "vst_normalized_filtered.rds"))
# 
# # Trait vectors for coloring (aligned to count sample order)
# diag_vec <- meta$diagnosis
# age_vec  <- meta$age
# sex_vec  <- meta$sex
# pmi_vec  <- meta$pmi
# rin_vec  <- meta$rin
# hemi_vec <- meta$hemisphere
# gsm_vec  <- meta[[sample_id_col]]
# 
# # PCA on all three matrices (prcomp on transposed: samples x genes)
# run_pca <- function(mat, label) {
#   pc <- prcomp(t(mat), center = TRUE, scale. = TRUE)
#   ve <- pc$sdev^2 / sum(pc$sdev^2)
#   df <- data.frame(
#     PC1 = pc$x[,1], PC2 = pc$x[,2],
#     sample     = rownames(pc$x),
#     diagnosis  = diag_vec,
#     age        = age_vec,
#     sex        = sex_vec,
#     pmi        = pmi_vec,
#     rin        = rin_vec,
#     hemisphere = hemi_vec,
#     stringsAsFactors = FALSE
#   )
#   list(pca = pc, ve = ve, df = df, label = label)
# }
# 
# pca_raw  <- run_pca(mat_raw,  "Raw counts")
# pca_log2 <- run_pca(mat_log2, "log2(raw+1)")
# pca_vst  <- run_pca(mat_vst,  "VST")
# 
# logit(sprintf("  Raw:     PC1=%.1f%%  PC2=%.1f%%  PC3=%.1f%%",
#               100*pca_raw$ve[1], 100*pca_raw$ve[2], 100*pca_raw$ve[3]))
# logit(sprintf("  log2+1:  PC1=%.1f%%  PC2=%.1f%%  PC3=%.1f%%",
#               100*pca_log2$ve[1], 100*pca_log2$ve[2], 100*pca_log2$ve[3]))
# logit(sprintf("  VST:     PC1=%.1f%%  PC2=%.1f%%  PC3=%.1f%%",
#               100*pca_vst$ve[1], 100*pca_vst$ve[2], 100*pca_vst$ve[3]))
# 
# # Trait-PC correlations (numeric covariates correlated directly; categorical as numeric factor)
# diag_num <- as.numeric(factor(diag_vec))
# sex_num  <- as.numeric(factor(sex_vec))
# hemi_num <- as.numeric(factor(hemi_vec))
# 
# logit("")
# logit("  Trait correlations with PC1/PC2:")
# for (nm in c("pca_raw", "pca_log2", "pca_vst")) {
#   obj <- get(nm)
#   logit(sprintf("    %s:", obj$label))
#   logit(sprintf("      cor(PC1, diagnosis)=%.3f  cor(PC1, age)=%.3f  cor(PC1, sex)=%.3f",
#                 cor(obj$pca$x[,1], diag_num), cor(obj$pca$x[,1], age_vec), cor(obj$pca$x[,1], sex_num)))
#   logit(sprintf("      cor(PC1, pmi)=%.3f  cor(PC1, rin)=%.3f  cor(PC1, hemisphere)=%.3f",
#                 cor(obj$pca$x[,1], pmi_vec), cor(obj$pca$x[,1], rin_vec), cor(obj$pca$x[,1], hemi_num)))
#   logit(sprintf("      cor(PC2, diagnosis)=%.3f  cor(PC2, age)=%.3f  cor(PC2, sex)=%.3f",
#                 cor(obj$pca$x[,2], diag_num), cor(obj$pca$x[,2], age_vec), cor(obj$pca$x[,2], sex_num)))
#   logit(sprintf("      cor(PC2, pmi)=%.3f  cor(PC2, rin)=%.3f  cor(PC2, hemisphere)=%.3f",
#                 cor(obj$pca$x[,2], pmi_vec), cor(obj$pca$x[,2], rin_vec), cor(obj$pca$x[,2], hemi_num)))
# }
# 
# # --- Plotting function (supports 6 covariates) ---
# make_pca_plot <- function(pca_obj, color_by, title_suffix) {
#   df <- pca_obj$df
#   ve <- pca_obj$ve
#   xl <- sprintf("PC1 (%.1f%%)", 100*ve[1])
#   yl <- sprintf("PC2 (%.1f%%)", 100*ve[2])
#   base_title <- paste0("PCA - ", pca_obj$label, " - ", title_suffix)
#   
#   p <- if (color_by == "diagnosis") {
#     ggplot(df, aes(PC1, PC2, color = diagnosis)) +
#       geom_point(size = 3.2, alpha = 0.9) +
#       scale_color_manual(values = diag_colors)
#   } else if (color_by == "sex") {
#     ggplot(df, aes(PC1, PC2, color = sex)) +
#       geom_point(size = 3.2, alpha = 0.9) +
#       scale_color_manual(values = sex_colors)
#   } else if (color_by == "hemisphere") {
#     ggplot(df, aes(PC1, PC2, color = hemisphere)) +
#       geom_point(size = 3.2, alpha = 0.9) +
#       scale_color_manual(values = hemi_colors)
#   } else if (color_by == "age") {
#     ggplot(df, aes(PC1, PC2, color = age)) +
#       geom_point(size = 3.2, alpha = 0.9) +
#       scale_color_gradient(low = "#2166AC", high = "#B2182B")
#   } else if (color_by == "pmi") {
#     ggplot(df, aes(PC1, PC2, color = pmi)) +
#       geom_point(size = 3.2, alpha = 0.9) +
#       scale_color_gradient(low = "#75A025", high = "#FF9400")
#   } else if (color_by == "rin") {
#     ggplot(df, aes(PC1, PC2, color = rin)) +
#       geom_point(size = 3.2, alpha = 0.9) +
#       scale_color_gradient(low = "#0279EE", high = "#FD9BED")
#   }
#   
#   p + geom_text_repel(aes(label = sample), size = 2.4, max.overlaps = 18, show.legend = FALSE) +
#     labs(title = base_title, x = xl, y = yl) +
#     theme_minimal(base_family = "sans") +
#     theme(text = element_text(size = 10), plot.title = element_text(face = "bold", size = 11),
#           legend.position = "bottom")
# }
# 
# save_plot <- function(p, name) {
#   ggsave(file.path(workspace_fig, paste0(name, ".png")), p, width = 8, height = 6, dpi = 150)
#   ggsave(file.path(workspace_fig, paste0(name, ".svg")), p, width = 8, height = 6)
#   system2("cp", args = c(file.path(workspace_fig, paste0(name, ".png")), fig_dir))
#   system2("cp", args = c(file.path(workspace_fig, paste0(name, ".svg")), fig_dir))
# }
# 
# # --- Individual plots: 3 methods x 6 covariates = 18 ---
# color_labels <- c(diagnosis="diagnosis", age="age", sex="sex",
#                   pmi="pmi", rin="rin", hemisphere="hemisphere")
# method_keys <- list(raw="pca_raw", log2="pca_log2", vst="pca_vst")
# method_short <- c(raw="raw", log2="log2", vst="vst")
# 
# for (mkey in names(method_keys)) {
#   pca_obj <- get(method_keys[[mkey]])
#   for (cv in names(color_labels)) {
#     plot_name <- paste0("pca_", method_short[mkey], "_", cv)
#     save_plot(make_pca_plot(pca_obj, cv, paste0("colored by ", cv)), plot_name)
#   }
# }
# logit("  Saved 18 individual PCA plots (3 methods x 6 covariates)")
# 
# # --- Side-by-side comparison panels (patchwork): one per covariate ---
# for (cv in names(color_labels)) {
#   cv_label <- color_labels[cv]
#   compare_p <- (
#     (make_pca_plot(pca_raw,  cv, "Raw")     + labs(title = "Raw counts")) |
#       (make_pca_plot(pca_log2, cv, "log2+1")  + labs(title = "log2(raw+1)")) |
#       (make_pca_plot(pca_vst,  cv, "VST")     + labs(title = "VST"))
#   ) + plot_annotation(
#     title = paste0("PCA comparison - colored by ", cv_label, " (Raw vs log2+1 vs VST)"),
#     theme = theme(plot.title = element_text(face = "bold", size = 13))
#   )
#   fname <- paste0("pca_compare_", cv)
#   ggsave(file.path(workspace_fig, paste0(fname, ".png")), compare_p, width = 20, height = 7, dpi = 150)
#   ggsave(file.path(workspace_fig, paste0(fname, ".svg")), compare_p, width = 20, height = 7)
#   system2("cp", args = c(file.path(workspace_fig, paste0(fname, ".png")), fig_dir))
#   system2("cp", args = c(file.path(workspace_fig, paste0(fname, ".svg")), fig_dir))
# }
# logit("  Saved 6 comparison panels (one per covariate, 3 methods side-by-side)")
# 
# # --- Scree comparison ---
# make_scree <- function(ve, label, fill_col) {
#   df <- data.frame(PC = factor(1:10, levels = 1:10), VarExp = ve[1:10] * 100)
#   ggplot(df, aes(PC, VarExp)) +
#     geom_col(fill = fill_col, alpha = 0.85) +
#     geom_text(aes(label = sprintf("%.1f", VarExp)), vjust = -0.3, size = 2.8) +
#     labs(title = label, x = "PC", y = "% variance") +
#     theme_minimal(base_family = "sans") +
#     theme(text = element_text(size = 9), plot.title = element_text(face = "bold", size = 10)) +
#     ylim(0, max(ve[1:10]*100) * 1.25)
# }
# scree_compare <- (
#   make_scree(pca_raw$ve,  "Raw counts",  "#75A025") |
#     make_scree(pca_log2$ve, "log2(raw+1)", "#0279EE") |
#     make_scree(pca_vst$ve,  "VST",         "#FF9400")
# ) + plot_annotation(title = "Variance explained - Raw vs log2+1 vs VST",
#                     theme = theme(plot.title = element_text(face = "bold", size = 13)))
# ggsave(file.path(workspace_fig, "pca_scree_compare.png"), scree_compare, width = 18, height = 6, dpi = 150)
# ggsave(file.path(workspace_fig, "pca_scree_compare.svg"), scree_compare, width = 18, height = 6)
# system2("cp", args = c(file.path(workspace_fig, "pca_scree_compare.png"), fig_dir))
# system2("cp", args = c(file.path(workspace_fig, "pca_scree_compare.svg"), fig_dir))
# logit("  Saved scree comparison plot")
# 
# # Save variance explained tables
# for (mkey in names(method_keys)) {
#   obj <- get(method_keys[[mkey]])
#   ve_df <- data.frame(
#     PC = paste0("PC", 1:length(obj$ve)),
#     Variance_Pct = round(obj$ve * 100, 2),
#     Cumulative_Pct = round(cumsum(obj$ve) * 100, 2)
#   )
#   write.csv(ve_df, file.path(out_dir, paste0("pca_variance_", method_short[mkey], ".csv")),
#             row.names = FALSE)
# }
# logit("  Saved 3 variance explained CSVs")

###############################################################################
# STEP 9: goodSamplesGenes QC
###############################################################################
logit("")
logit("Step 9 - goodSamplesGenes QC")

datExpr <- t(temp.norm)
gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  logit(sprintf("  Removed %d bad genes, %d bad samples",
                sum(!gsg$goodGenes), sum(!gsg$goodSamples)))
} else {
  logit("  goodSamplesGenes: all samples and genes pass QC")
}
logit(sprintf("  Final datExpr: %d samples x %d genes", nrow(datExpr), ncol(datExpr)))

# Re-align meta to datExpr sample order
meta <- temp.meta[match(rownames(datExpr), temp.meta[[sample_id_col]]), ]

###############################################################################
# STEP 10: Sample dendrogram (outlier check)
###############################################################################
logit("")
logit("Step 10 - Sample dendrogram (outlier check)")

sampleDist <- dist(datExpr, method = "euclidean")
sampleTree <- hclust(sampleDist, method = "average")

get_color_diag <- function(val) {
  if (is.na(val)) return("white")
  if (val == "Female") return("#B2182B")
  if (val == "Male") return("#2166AC")
  return("white")
}
# get_color_hemi <- function(val) {
#   if (is.na(val)) return("white")
#   if (val == "L") return("#0279EE")
#   if (val == "R") return("#FF9400")
#   return("white")
# }

sample_order <- sampleTree$labels[sampleTree$order]
colors_diag  <- sapply(meta[sample_order, "Sex"], get_color_diag)
# colors_hemi  <- sapply(meta[sample_order, "hemisphere"], get_color_hemi)

# PNG
png(file.path(ws_fig, "sample_dendrogram.png"), width = 1400, height = 800, res = 150)
layout(matrix(c(1, 2, 3), nrow = 3, byrow = TRUE), heights = c(5, 0.6, 0.6))
par(mar = c(7, 5, 3, 2))
plot(sampleTree, main = "Sample Dendrogram (Microarray-Norm) - OSD152 Outlier Check",
     xlab = "", sub = "", hang = -1, cex = 0.65, las = 2)
par(mar = c(0, 5, 0.5, 2))
plot(0, type = "n", xlim = c(0.5, length(sample_order)), ylim = c(0, 1),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
rect(seq(0.5, length(sample_order) - 0.5), 0, seq(1.5, length(sample_order) + 0.5), 1,
     col = colors_diag, border = NA)
mtext("Sex", side = 2, line = 0.5, cex = 0.7, las = 2)
# par(mar = c(0, 5, 0.5, 2))
# plot(0, type = "n", xlim = c(0.5, length(sample_order)), ylim = c(0, 1),
#      xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
# rect(seq(0.5, length(sample_order) - 0.5), 0, seq(1.5, length(sample_order) + 0.5), 1,
#      col = colors_hemi, border = NA)
# mtext("hemisphere", side = 2, line = 0.5, cex = 0.7, las = 2)
dev.off()


# ###############################################################################
# # Write summary log
# ###############################################################################
# writeLines(log_lines, file.path(out_dir, "workflow_summary.txt"))
# 
# cat("\n=== Phase 1 (QC + PCA) complete ===\n")
# cat("All outputs in:", out_dir, "\n")
# cat("  Figures: 18 individual PCA, 6 comparison panels, 1 scree, 1 dendrogram\n")
# cat("  Tables: metadata_parsed.csv, pca_variance_*.csv\n")
# cat("  RDS: vst_normalized_filtered.rds\n")
# cat("  Log: workflow_summary.txt\n")


###############################################################################
# PHASE 2: WGCNA on LIMMA-normalized data
###############################################################################
out_dir <- out_dir_p2
fig_dir <- ws_fig

# # Point Phase 2 at the VST matrix saved by Phase 1
# vst_path  <- file.path(out_dir_p1, "vst_normalized_filtered.rds")
# 
# # SECTION 4: Load VST + goodSamplesGenes + sample dendrogram
# ###############################################################################
# logit("=== SECTION 4: Load VST + goodSamplesGenes + sample dendrogram ===")
# 
# mat_vst <- readRDS(vst_path)
# logit(sprintf("  Loaded VST matrix: %d genes x %d samples", nrow(mat_vst), ncol(mat_vst)))
# 
# datExpr <- t(mat_vst)
# logit(sprintf("  datExpr: %d samples x %d genes", nrow(datExpr), ncol(datExpr)))
# 
# # --- Parse metadata (positional indexing) ---
# meta_raw <- read.xlsx(meta_path, sheet = 1, check.names = FALSE)
# colnames(meta_raw) <- trimws(colnames(meta_raw))
# 
# strip_prefix <- function(s) { s <- as.character(s); sub("^[^:]+:\\s*", "", s) }
# 
# meta_raw$diagnosis     <- strip_prefix(meta_raw[[10]])
# meta_raw$age           <- as.numeric(strip_prefix(meta_raw[[11]]))
# meta_raw$sex           <- strip_prefix(meta_raw[[12]])
# meta_raw$pmi           <- as.numeric(strip_prefix(meta_raw[[13]]))
# meta_raw$rin           <- as.numeric(strip_prefix(meta_raw[[14]]))
# meta_raw$hemisphere    <- strip_prefix(meta_raw[[15]])
# meta_raw$symptom_side  <- strip_prefix(meta_raw[[16]])
# meta_raw$group         <- strip_prefix(meta_raw[[17]])
# rownames(meta_raw) <- meta_raw$Sample_geo_accession
# meta_aligned <- meta_raw[rownames(datExpr), ]
# logit(sprintf("  Metadata aligned to datExpr: %s",
#               identical(rownames(datExpr), rownames(meta_aligned))))
# logit(sprintf("  Diagnosis: %s",
#               paste(sprintf("%s=%d", names(table(meta_aligned$diagnosis)),
#                             table(meta_aligned$diagnosis)), collapse=", ")))

# # goodSamplesGenes safety check
# gsg <- goodSamplesGenes(datExpr, verbose = 0)
# if (!gsg$allOK) {
#   datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
#   mat_vst <- t(datExpr)
#   logit(sprintf("  Removed %d bad genes, %d bad samples",
#                 sum(!gsg$goodGenes), sum(!gsg$goodSamples)))
# } else {
#   logit("  goodSamplesGenes: all samples and genes pass QC")
# }
# logit(sprintf("  datExpr after GSG: %d samples x %d genes", nrow(datExpr), ncol(datExpr)))
# 
# # Sample dendrogram with trait color bars
# sampleDist <- dist(datExpr, method = "euclidean")
# sampleTree <- hclust(sampleDist, method = "average")
# sample_order <- sampleTree$labels[sampleTree$order]
# colors_dx  <- dx_colors[meta_aligned[sample_order, "diagnosis"]]
# # colors_grp <- grp_colors[meta_aligned[sample_order, "group"]]
# # colors_hem <- hem_colors[meta_aligned[sample_order, "hemisphere"]]
# 
# png(file.path(ws_fig, "sample_dendrogram.png"), width = 1400, height = 900, res = 150)
# layout(matrix(c(1, 2, 3, 4), nrow = 4, byrow = TRUE), heights = c(5, 0.6, 0.6, 0.6))
# par(mar = c(7, 5, 3, 2))
# plot(sampleTree, main = "Sample Dendrogram (VST) - GSE135036 Outlier Check",
#      xlab = "", sub = "", hang = -1, cex = 0.65, las = 2)
# par(mar = c(0, 5, 0.5, 2))
# plot(0, type = "n", xlim = c(0.5, length(sample_order)), ylim = c(0, 1),
#      xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
# rect(seq(0.5, length(sample_order) - 0.5), 0, seq(1.5, length(sample_order) + 0.5), 1,
#      col = colors_dx, border = NA)
# mtext("diagnosis", side = 2, line = 0.5, cex = 0.7, las = 2)
# par(mar = c(0, 5, 0.5, 2))
# plot(0, type = "n", xlim = c(0.5, length(sample_order)), ylim = c(0, 1),
#      xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
# rect(seq(0.5, length(sample_order) - 0.5), 0, seq(1.5, length(sample_order) + 0.5), 1,
#      col = colors_grp, border = NA)
# mtext("group", side = 2, line = 0.5, cex = 0.7, las = 2)
# par(mar = c(0, 5, 0.5, 2))
# plot(0, type = "n", xlim = c(0.5, length(sample_order)), ylim = c(0, 1),
#      xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
# rect(seq(0.5, length(sample_order) - 0.5), 0, seq(1.5, length(sample_order) + 0.5), 1,
#      col = colors_hem, border = NA)
# mtext("hemisphere", side = 2, line = 0.5, cex = 0.7, las = 2)
# dev.off()
# 
# svg(file.path(ws_fig, "sample_dendrogram.svg"), width = 14, height = 9)
# layout(matrix(c(1, 2, 3, 4), nrow = 4, byrow = TRUE), heights = c(5, 0.6, 0.6, 0.6))
# par(mar = c(7, 5, 3, 2))
# plot(sampleTree, main = "Sample Dendrogram (VST) - GSE135036 Outlier Check",
#      xlab = "", sub = "", hang = -1, cex = 0.65, las = 2)
# par(mar = c(0, 5, 0.5, 2))
# plot(0, type = "n", xlim = c(0.5, length(sample_order)), ylim = c(0, 1),
#      xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
# rect(seq(0.5, length(sample_order) - 0.5), 0, seq(1.5, length(sample_order) + 0.5), 1,
#      col = colors_dx, border = NA)
# mtext("diagnosis", side = 2, line = 0.5, cex = 0.7, las = 2)
# par(mar = c(0, 5, 0.5, 2))
# plot(0, type = "n", xlim = c(0.5, length(sample_order)), ylim = c(0, 1),
#      xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
# rect(seq(0.5, length(sample_order) - 0.5), 0, seq(1.5, length(sample_order) + 0.5), 1,
#      col = colors_grp, border = NA)
# mtext("group", side = 2, line = 0.5, cex = 0.7, las = 2)
# par(mar = c(0, 5, 0.5, 2))
# plot(0, type = "n", xlim = c(0.5, length(sample_order)), ylim = c(0, 1),
#      xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
# rect(seq(0.5, length(sample_order) - 0.5), 0, seq(1.5, length(sample_order) + 0.5), 1,
#      col = colors_hem, border = NA)
# mtext("hemisphere", side = 2, line = 0.5, cex = 0.7, las = 2)
# dev.off()
# 
# system2("cp", args = c(file.path(ws_fig, "sample_dendrogram.png"), fig_dir))
# system2("cp", args = c(file.path(ws_fig, "sample_dendrogram.svg"), fig_dir))
# logit(sprintf("  Sample dendrogram saved | max height: %.2f, mean height: %.2f",
#               max(sampleTree$height), mean(sampleTree$height)))

###############################################################################
# SECTION 5: datExpr save
###############################################################################
logit("")
logit("=== SECTION 5: datExpr save ===")
meta_final <- meta[rownames(datExpr), ]
saveRDS(datExpr, file.path(out_dir, "OSD152_datExpr.rds"))
logit(sprintf("  Saved: OSD152_datExpr.rds (%d samples x %d genes)",
              nrow(datExpr), ncol(datExpr)))

###############################################################################
# SECTION 6: 6-trait matrix construction
###############################################################################
logit("")
logit("=== SECTION 6: 6-trait matrix construction ===")
stopifnot(identical(rownames(datExpr), rownames(meta_final)))

# Note: diagnosis uses Unicode right single quote (U+2019)
is_male <- meta_final$Sex == "Male"
is_female <- meta_final$Sex == "Female"
is_sham_0_Gy   <- meta_final$Dose == "0_Gy"
is_gamma_0.5_Gy   <- meta_final$Dose == "0.5_Gy"
is_gamma_2_Gy   <- meta_final$Dose == "2_Gy"
is_gamma_5_Gy   <- meta_final$Dose == "5_Gy"
is_gamma_8_Gy   <- meta_final$Dose == "8_Gy"
is_24hrs <- meta_final$Hour == "24Hrs_Gamma"
is_6hrs <- meta_final$Hour == "6Hrs_Gamma"

logit(sprintf("is_male=%d is_female=%d is_sham_0_Gy=%d is_gamma_0.5_Gy =%d is_gamma_2_Gy =%d is_gamma_5_Gy =%d is_gamma_8_Gy =%d is_24hrs=%d is_6hrs=%d",
              sum(is_male), sum(is_female),
              sum(is_sham_0_Gy), sum(is_gamma_0.5_Gy), sum(is_gamma_2_Gy), sum(is_gamma_5_Gy), sum(is_gamma_8_Gy),
              sum(is_24hrs), sum(is_6hrs)))

traits <- data.frame(
  Meta1_Male_only = ifelse(is_male, 1L, ifelse(is_female, 0L, NA)),
  Meta2_Female_only   = ifelse(is_female, 1L, ifelse(is_male, 0L, NA)),
  Meta3_Sham_only = ifelse(is_sham_0_Gy, 1L, 0L),
  Meta4_0.5_Gy_Only   = ifelse(is_gamma_0.5_Gy, 1L, 0L),
  Meta5_2_Gy_Only  = ifelse(is_gamma_2_Gy, 1L, 0L),
  Meta6_8_Gy_Only  = ifelse(is_gamma_8_Gy, 1L, 0L),
  Meta7_24hrs_Only = ifelse(is_24hrs, 1L, 0L),
  Meta8_6hrs_Only = ifelse(is_6hrs, 1L, 0L),
  stringsAsFactors = FALSE
)
rownames(traits) <- rownames(meta_final)
traits <- traits[rownames(datExpr), ]

logit("  Trait matrix columns:")
for (col in colnames(traits)) {
  logit(sprintf("    %-22s n=%-3d  1s=%-3d  0s=%-3d  NA=%-3d",
                col, sum(!is.na(traits[[col]])),
                sum(traits[[col]]==1, na.rm=TRUE),
                sum(traits[[col]]==0, na.rm=TRUE),
                sum(is.na(traits[[col]]))))
}

for (col in colnames(traits)) meta_final[[paste0("trait_", col)]] <- traits[[col]]
write.csv(meta_final, file.path(out_dir, "OSD152_metadata_aligned.csv"))
saveRDS(traits, file.path(out_dir, "OSD152_traits.rds"))
write.csv(traits, file.path(out_dir, "OSD152_traits.csv"))
logit("  Saved: traits.rds + traits.csv + metadata_aligned.csv")

logit("")
logit("  Alignment checks:")
logit(sprintf("    rownames(datExpr) == rownames(traits): %s",
              identical(rownames(datExpr), rownames(traits))))
logit(sprintf("    rownames(datExpr) == rownames(meta_final): %s",
              identical(rownames(datExpr), rownames(meta_final))))

###############################################################################
# SECTION 7: Soft power diagnostics
###############################################################################
logit("")
logit("=== SECTION 7: Soft power diagnostics ===")

gc()
enableWGCNAThreads(2)
powers <- c(c(1:10), seq(from = 12, to = 30, by = 2))
sft <- pickSoftThreshold(datExpr, dataIsExpr = TRUE, powerVector = powers,
                         verbose = 5, networkType = "signed",
                         corFnc = "bicor", moreNetworkConcepts = TRUE,
                         blockSize = 10000)

logit("")
logit(sprintf("  Auto estimate (for reference only): %s", sft$powerEstimate))
logit("")
logit("  Power  SFT.R.sq   slope    trunc.R.sq  mean.k.   median.k.  max.k.")
logit(paste("  ", paste(rep("-", 68), collapse = "")))
for (i in 1:nrow(sft$fitIndices)) {
  logit(sprintf("  %-6d %-10.4f %-8.4f %-11.4f %-9.1f %-10.1f %-8.1f",
                sft$fitIndices$Power[i], sft$fitIndices$SFT.R.sq[i],
                sft$fitIndices$slope[i], sft$fitIndices$truncated.R.sq[i],
                sft$fitIndices$mean.k.[i], sft$fitIndices$median.k.[i],
                sft$fitIndices$max.k.[i]))
}

# Diagnostic plots
png(file.path(ws_fig, "soft_power_selection.png"), width = 1200, height = 600, res = 150)
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "SFT, signed R^2", type = "n",
     main = "Scale independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = "red")
abline(h = 0.9, col = "red")
abline(h = 0.85, col = "red", lty = 2)
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], type = "n",
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red")
dev.off()

# svg(file.path(ws_fig, "soft_power_selection.svg"), width = 12, height = 6)
# par(mfrow = c(1, 2))
# plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
#      xlab = "Soft Threshold (power)", ylab = "SFT, signed R^2", type = "n",
#      main = "Scale independence")
# text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
#      labels = powers, col = "red")
# abline(h = 0.9, col = "red")
# abline(h = 0.85, col = "red", lty = 2)
# plot(sft$fitIndices[, 1], sft$fitIndices[, 5], type = "n",
#      xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
#      main = "Mean connectivity")
# text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red")
# dev.off()
# 
# system2("cp", args = c(file.path(ws_fig, "soft_power_selection.png"), fig_dir))
# system2("cp", args = c(file.path(ws_fig, "soft_power_selection.svg"), fig_dir))
saveRDS(sft, file.path(out_dir, "soft_threshold_results.rds"))
logit("  Saved: soft_power_selection.png, soft_threshold_results.rds")


# SECTION 8: Network construction (blockwiseModules)
###############################################################################
logit("")
logit("=== SECTION 8: Network construction (blockwiseModules) ===")

cor <- WGCNA::cor   # required for bicor

# ---- WGCNA parameters (USER-SET) ----
soft_power       <- 14       # USER DECISION from pickSoftThreshold diagnostics
network_type     <- "signed"
tom_type         <- "signed"
cor_type         <- "bicor"
max_block_size   <- 20000    # USER SET: all 17259 genes in one block
deep_split       <- 3
min_module_size  <- 25
merge_cut_height <- 0.2

enableWGCNAThreads()
net <- blockwiseModules(
  datExpr,
  power = soft_power,
  networkType = network_type,
  TOMType = tom_type,
  corType = cor_type,
  maxBlockSize = max_block_size,
  deepSplit = deep_split,
  minModuleSize = min_module_size,
  mergeCutHeight = merge_cut_height,
  detectCutHeight = 1,
  pamStage = TRUE,
  pamRespectsDendro = TRUE,
  saveTOMs = TRUE,
  tomFileBase = file.path(out_dir, "wgcna_TOM"),
  numericLabels = TRUE,
  verbose = 3
)

module_colors <- labels2colors(net$colors)
names(module_colors) <- colnames(datExpr)

cor <- stats::cor   # switch back for downstream

logit("")
logit(sprintf("  Network construction complete (power = %d)", soft_power))
tab <- table(module_colors)
logit(sprintf("  Modules (excl grey): %d", sum(tab[names(tab) != "grey"] > 0)))
logit(sprintf("  Grey (unassigned): %d (%.1f%%)", tab["grey"],
              round(tab["grey"]/sum(tab)*100, 1)))
logit("  Module sizes:")
for (m in names(sort(tab, decreasing=TRUE))) {
  logit(sprintf("    %-15s %5d genes", m, tab[m]))
}

saveRDS(net, file.path(out_dir, "wgcna_network.rds"))
saveRDS(module_colors, file.path(out_dir, "wgcna_module_colors.rds"))
logit("  Saved: wgcna_network.rds, wgcna_module_colors.rds")

###############################################################################
# SECTION 9: Module-trait correlations + dendrogram with legend
###############################################################################
logit("")
logit("=== SECTION 9: Module-trait correlations + dendrogram ===")

# Build color-label mapping for ME column lookup
colorLabelMap <- unique(data.frame(color = module_colors, label = net$colors))
colorLabelMap$MEcol <- paste0("ME", colorLabelMap$label)

# Module eigengenes — rename ME columns from numeric labels to color names
MEs <- net$MEs
MEs <- MEs[rownames(datExpr), ]
me_to_color <- unique(data.frame(MEcol = colorLabelMap$MEcol, color = colorLabelMap$color))
new_me_names <- paste0("ME", me_to_color$color[match(colnames(MEs), me_to_color$MEcol)])
colnames(MEs) <- new_me_names
logit(sprintf("  ME columns renamed to colors: %d modules", length(colnames(MEs))))

# Module-trait correlations
module_trait_cor <- cor(MEs, traits, use = "p", method = "pearson")

# P-values with per-trait n
n_per_trait <- colSums(!is.na(traits))
moduleTraitPvalue <- matrix(NA, nrow = ncol(MEs), ncol = ncol(traits))
rownames(moduleTraitPvalue) <- colnames(MEs)
colnames(moduleTraitPvalue) <- colnames(traits)
for (j in 1:ncol(traits)) {
  moduleTraitPvalue[, j] <- corPvalueStudent(module_trait_cor[, j], n_per_trait[j])
}

logit(sprintf("  Non-NA n per trait: %s",
              paste(sprintf("%s=%d", names(n_per_trait), n_per_trait), collapse=", ")))

# Print correlation table
trait_cols <- colnames(traits)
logit("")
logit("  Module-trait correlations (r values):")
logit(sprintf("  %-20s %s", "Module",
              paste(sprintf("%-22s", trait_cols), collapse="")))
logit(paste("  ", paste(rep("-", 20 + 22*length(trait_cols)), collapse = "")))
for (me in rownames(module_trait_cor)) {
  logit(sprintf("  %-20s %s", me,
                paste(sprintf("%-22.4f", module_trait_cor[me, ]), collapse="")))
}

# Save CSV — full r and p for all traits
mt_results <- data.frame(
  Module = rownames(module_trait_cor),
  module_trait_cor,
  stringsAsFactors = FALSE)
pval_df <- as.data.frame(moduleTraitPvalue)
colnames(pval_df) <- paste0("p_", trait_cols)
mt_results <- cbind(mt_results, pval_df)

write.csv(mt_results, file.path(out_dir, "wgcna_module_trait_cor.csv"), row.names = FALSE)
logit("  Saved: wgcna_module_trait_cor.csv")

# Save eigengenes
write.csv(MEs, file.path(out_dir, "OSD152_eigengenes.csv"))
logit("  Saved: OSD152_eigengenes.csv")

# --- Module-trait correlation heatmap (ComplexHeatmap) ---
text_matrix <- matrix(
  paste0(signif(module_trait_cor, 2), "\n(", signif(moduleTraitPvalue, 1), ")"),
  nrow = nrow(module_trait_cor), ncol = ncol(module_trait_cor))
rownames(text_matrix) <- rownames(module_trait_cor)
colnames(text_matrix) <- colnames(module_trait_cor)

col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
ht <- Heatmap(module_trait_cor,
              name = "Correlation", col = col_fun,
              cluster_rows = FALSE, cluster_columns = FALSE,
              show_row_names = TRUE, show_column_names = TRUE,
              row_names_side = "left", column_names_side = "bottom",
              column_title = "Module-Trait Relationships (OSD152- Gamma Irradiation)",
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(text_matrix[i, j], x, y, gp = gpar(fontsize = 5))
              },
              heatmap_legend_param = list(title = "Corr", at = c(-1, -0.5, 0, 0.5, 1)),
              width = unit(12, "cm"), height = unit(max(8, nrow(module_trait_cor) * 0.6), "cm"))

png(file.path(ws_fig, "module_trait_correlation.png"),
    width = 12, height = 32, units = "in", res = 300)
draw(ht)
dev.off()
logit("  Saved: module_trait_correlation.png")

# --- Gene dendrogram with ComplexHeatmap-style legend ---
cor_to_color <- function(r, col_fun) {
  if (is.na(r)) return("grey90")
  return(col_fun(r))
}

get_module_cor <- function(mc, trait_col, cor_mat) {
  me_name <- paste0("ME", mc)
  if (me_name %in% rownames(cor_mat)) return(cor_mat[me_name, trait_col]) else return(NA)
}

# Trait columns ordered
trait_cols_ordered <- c("Meta1_Male_only", "Meta2_Female_only", "Meta3_Sham_only",
                        "Meta4_0.5_Gy_Only", "Meta5_2_Gy_Only", "Meta6_8_Gy_Only", 
                        "Meta7_24hrs_Only", "Meta8_6hrs_Only")
trait_color_list <- lapply(trait_cols_ordered, function(tc) {
  cor_per_gene <- sapply(module_colors, function(mc)
    get_module_cor(mc, tc, module_trait_cor))
  sapply(cor_per_gene, function(r) cor_to_color(r, col_fun))
})
names(trait_color_list) <- trait_cols_ordered

datColors <- data.frame(
  ModuleColors = module_colors,
  do.call(cbind, trait_color_list),
  stringsAsFactors = FALSE)
short_labels <- c("Modules",
                  "Meta1_Male (n=25)", "Meta2_Female (n=25)", "Meta3_Sham (n=10)",
                  "Meta4_0.5_Gy (n=10)", "Meta5_2_Gy (n=10)", "Meta6_8_Gy (n=10)",
                  "Meta7_24hrs (n=25)", "Meta8_6hrs (n=25)")
colnames(datColors) <- short_labels

# Generate dendrogram and legend as separate images, then stack
b <- 1  # single block
block_genes <- net$blockGenes[[b]]
datColors_block <- datColors[block_genes, , drop = FALSE]

# Panel 1: Dendrogram (standalone)
png(filename = file.path(ws_fig, "dendro_panel.png"),
    width = 2400, height = 900, res = 150)
par(mar = c(7, 5, 3, 2))
plotDendroAndColors(
  dendro = net$dendrograms[[b]], colors = datColors_block,
  groupLabels = colnames(datColors), dendroLabels = FALSE,
  hang = 0.03, addGuide = TRUE, guideHang = 0.05,
  main = "Block 1 - Gene Dendrogram and Module Colors")
dev.off()

# Panel 2: Legend (ComplexHeatmap::Legend, matching heatmap exactly)
lgd <- Legend(col_fun = col_fun, title = "Corr",
              at = c(-1, -0.5, 0, 0.5, 1),
              legend_height = unit(3, "cm"),
              legend_width = unit(0.6, "cm"),
              title_gp = gpar(fontsize = 10, fontface = "bold"),
              labels_gp = gpar(fontsize = 9),
              title_position = "topleft")

png(filename = file.path(ws_fig, "legend_panel.png"),
    width = 300, height = 800, res = 150)
draw(lgd)
dev.off()

# Trim excess whitespace from legend panel
system2("convert", args = c(
  file.path(ws_fig, "legend_panel.png"),
  "-trim",
  file.path(ws_fig, "legend_panel_trimmed.png")
))

# Stack: dendrogram on left, trimmed legend on right
system2("convert", args = c(
  file.path(ws_fig, "dendro_panel.png"),
  file.path(ws_fig, "legend_panel_trimmed.png"),
  "+append",
  file.path(ws_fig, "WGCNA_Dendrogram_Block1.png")
))

# # SVG: dendrogram panel
# svg(filename = file.path(ws_fig, "dendro_panel.svg"),
#     width = 20, height = 8)
# par(mar = c(2, 5, 3, 2))
# plotDendroAndColors(
#   dendro = net$dendrograms[[b]], colors = datColors_block,
#   groupLabels = colnames(datColors), dendroLabels = FALSE,
#   hang = 0.03, addGuide = TRUE, guideHang = 0.05,
#   main = "Block 1 - Gene Dendrogram and Module Colors")
# dev.off()
# system2("cp", args = c(file.path(ws_fig, "dendro_panel.svg"),
#                        file.path(ws_fig, "WGCNA_Dendrogram_Block1.svg")))
# 
# # Copy final outputs to results
# system2("cp", args = c(file.path(ws_fig, "WGCNA_Dendrogram_Block1.png"), fig_dir))
# system2("cp", args = c(file.path(ws_fig, "WGCNA_Dendrogram_Block1.svg"), fig_dir))

# Clean up temp panels
unlink(file.path(ws_fig, "dendro_panel.png"))
unlink(file.path(ws_fig, "legend_panel.png"))
unlink(file.path(ws_fig, "legend_panel_trimmed.png"))
# unlink(file.path(ws_fig, "dendro_panel.svg"))
# unlink(file.path(ws_fig, "legend_panel.svg"))
logit("  Saved: WGCNA_Dendrogram_Block1.png (with ComplexHeatmap legend)")

# --- Module eigengene dendrogram ---
MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")

png(file.path(ws_fig, "module_eigengene_dendrogram.png"),
    width = 1000, height = 800, res = 150)
par(mar = c(4, 5, 3, 2))
plot(METree, main = "Module Eigengene Dendrogram", xlab = "", sub = "", hang = -1)
abline(h = 0.2, col = "red", lty = 2)
dev.off()

# svg(file.path(ws_fig, "module_eigengene_dendrogram.svg"), width = 10, height = 5)
# par(mar = c(4, 5, 3, 2))
# plot(METree, main = "Module Eigengene Dendrogram", xlab = "", sub = "", hang = -1)
# abline(h = 0.2, col = "red", lty = 2)
# dev.off()

# system2("cp", args = c(file.path(ws_fig, "module_eigengene_dendrogram.png"), fig_dir))
# system2("cp", args = c(file.path(ws_fig, "module_eigengene_dendrogram.svg"), fig_dir))
logit("  Saved: module_eigengene_dendrogram.png")

###############################################################################
# SECTION 10: KME + hub genes + GO BP enrichment + lollipop plots + Excel + RDS
###############################################################################
logit("")
logit("=== SECTION 10: KME + hub genes + GO BP enrichment + Excel + RDS ===")

geneNames <- colnames(datExpr)
nGenes <- ncol(datExpr)

# KME for all genes against all module eigengenes (Pearson)
datKME <- stats::cor(datExpr, MEs, use = "p")

# KME_own: each gene's KME to its own module eigengene
KME_own <- rep(NA, nGenes)
names(KME_own) <- geneNames
for (i in 1:nGenes) {
  myColor <- module_colors[i]
  myMEcol <- paste0("ME", myColor)
  if (myMEcol %in% colnames(datKME)) {
    KME_own[i] <- datKME[i, myMEcol]
  }
}

# Total connectivity
cor <- WGCNA::cor
kTotal <- softConnectivity(datExpr, power = soft_power, type = network_type, corFnc = cor_type)
cor <- stats::cor

# kWithin: compute per-module adjacency
kWithin <- rep(0, nGenes)
names(kWithin) <- geneNames

unique_modules <- unique(module_colors)
unique_modules <- unique_modules[unique_modules != "grey"]

for (mod in unique_modules) {
  mod_idx <- which(module_colors == mod)
  if (length(mod_idx) < 2) next
  adj_mod <- adjacency(datExpr[, mod_idx], power = soft_power,
                       type = network_type, corFnc = cor_type)
  kWithin[mod_idx] <- rowSums(adj_mod) - 1
}

kOut <- kTotal - kWithin

# Hub score: |KME_own| * kWithin (grey = NA)
hub_score <- abs(KME_own) * kWithin
hub_score[module_colors == "grey"] <- NA

gene_info <- data.frame(
  gene = geneNames,
  module = module_colors,
  KME_own = round(KME_own, 4),
  kWithin = round(kWithin, 2),
  kOut = round(kOut, 2),
  kTotal = round(kTotal, 2),
  hub_score = round(hub_score, 2),
  stringsAsFactors = FALSE)
gene_info <- gene_info[order(gene_info$module, -gene_info$hub_score), ]

logit(sprintf("  Genes with hub_score (non-grey): %d", sum(!is.na(gene_info$hub_score))))
logit(sprintf("  Grey genes (hub_score=NA): %d", sum(is.na(gene_info$hub_score))))

# Top 5 hub genes per module
logit("  Top 5 hub genes per module:")
for (mod in names(sort(table(module_colors), decreasing=TRUE))) {
  if (mod == "grey") next
  mod_genes <- gene_info[gene_info$module == mod, ]
  top5 <- head(mod_genes[order(-mod_genes$hub_score), ], 5)
  logit(sprintf("    %-15s (%4d genes): %s", mod, nrow(mod_genes),
                paste(top5$gene, collapse=", ")))
}

write.csv(gene_info, file.path(out_dir, "wgcna_gene_modules.csv"), row.names = FALSE)
logit("  Saved: wgcna_gene_modules.csv")

# ---- Thresholds for module-trait significance + core genes ----
r_threshold      <- 0.3
p_threshold      <- 0.05
kme_threshold    <- 0.7

# --- Two-level filtering + GO BP enrichment + lollipop plots ---
logit("")
logit("  --- Two-level filtering + GO BP enrichment ---")

# # Helper: convert gene names to Entrez IDs
# genes_to_entrez <- function(gene_names) {
#   entrez_from_symbol <- mapIds(org.Hs.eg.db, keys = gene_names,
#                                column = "ENTREZID", keytype = "SYMBOL",
#                                multiVals = "first")
#   unmapped <- is.na(entrez_from_symbol)
#   numeric_genes <- gene_names[unmapped]
#   is_numeric <- grepl("^[0-9]+$", numeric_genes)
#   entrez_from_symbol[unmapped][is_numeric] <- numeric_genes[is_numeric]
#   return(entrez_from_symbol)
# }
# 
# all_core_genes <- list()
# all_enrichment <- list()
# 
# for (trait_name in colnames(module_trait_cor)) {
#   logit("")
#   logit(sprintf("  Trait: %s (n = %d)", trait_name, n_per_trait[trait_name]))
#   
#   # Level 1: significant modules
#   sig_modules <- which(abs(module_trait_cor[, trait_name]) > r_threshold &
#                          moduleTraitPvalue[, trait_name] < p_threshold)
#   sig_module_names <- rownames(module_trait_cor)[sig_modules]
#   sig_colors <- gsub("^ME", "", sig_module_names)
#   logit(sprintf("  Level 1 - Significant modules (|r|>%g, p<%g): %d",
#                 r_threshold, p_threshold, length(sig_colors)))
#   
#   if (length(sig_colors) == 0) {
#     logit("  No significant modules. Skipping.")
#     next
#   }
#   
#   for (sm in sig_module_names) {
#     logit(sprintf("    %s: r = %.4f, p = %.2e", sm,
#                   module_trait_cor[sm, trait_name], moduleTraitPvalue[sm, trait_name]))
#   }
#   
#   # Level 2: core genes (|KME| > kme_threshold) within significant modules
#   core_mask <- gene_info$module %in% sig_colors & abs(gene_info$KME_own) > kme_threshold
#   core_genes <- gene_info[core_mask, ]
#   logit(sprintf("  Level 2 - Core genes (|KME|>%g): %d", kme_threshold, nrow(core_genes)))
#   
#   if (nrow(core_genes) == 0) {
#     logit("  No core genes found. Skipping GO enrichment.")
#     next
#   }
#   
#   all_core_genes[[trait_name]] <- core_genes
#   write.csv(core_genes, file.path(out_dir, paste0("wgcna_core_genes_", trait_name, ".csv")),
#             row.names = FALSE)
#   
#   # GO BP enrichment per module
#   for (mod_color in unique(core_genes$module)) {
#     mod_genes <- core_genes$gene[core_genes$module == mod_color]
#     if (length(mod_genes) < 5) {
#       logit(sprintf("    %s: <5 core genes, skipping GO", mod_color))
#       next
#     }
#     
#     entrez <- genes_to_entrez(mod_genes)
#     entrez <- na.omit(entrez)
#     if (length(entrez) < 5) {
#       logit(sprintf("    %s: <5 Entrez-mapped genes, skipping GO", mod_color))
#       next
#     }
#     
#     logit(sprintf("    %s: GO BP enrichment with %d genes...", mod_color, length(entrez)))
#     
#     ego <- enrichGO(gene = entrez, OrgDb = org.Hs.eg.db, ont = "BP",
#                     pAdjustMethod = "BH", pvalueCutoff = 0.05,
#                     qvalueCutoff = 0.2, readable = TRUE)
#     
#     if (nrow(as.data.frame(ego)) == 0) {
#       logit(sprintf("    %s: no significant GO terms", mod_color))
#       next
#     }
#     
#     ego_df <- as.data.frame(ego)
#     ego_df$module <- mod_color
#     ego_df$trait <- trait_name
#     ego_df$r_value <- module_trait_cor[paste0("ME", mod_color), trait_name]
#     all_enrichment[[length(all_enrichment) + 1]] <- ego_df
#     
#     # Lollipop plot (top 15 terms)
#     top15 <- head(ego_df[order(ego_df$p.adjust), ], 15)
#     top15$Description <- factor(top15$Description, levels = rev(top15$Description))
#     
#     p <- ggplot(top15, aes(x = Count, y = Description, fill = -log10(p.adjust))) +
#       geom_segment(aes(x = 0, xend = Count, y = Description, yend = Description), color = "grey70") +
#       geom_point(size = 4, shape = 21, color = "black", stroke = 0.5) +
#       scale_fill_gradient(low = "#0279EE", high = "#FF9400", name = "-log10(padj)") +
#       labs(title = paste0("GO BP - ", mod_color, " (", trait_name, ", r=",
#                           round(module_trait_cor[paste0("ME", mod_color), trait_name], 2), ")"),
#            x = "Gene Count", y = "") +
#       theme_minimal(base_family = "sans") +
#       theme(text = element_text(size = 10), plot.title = element_text(face = "bold", size = 11))
#     
#     ggsave(file.path(ws_fig, paste0("lollipop_", trait_name, "_", mod_color, ".png")),
#            p, width = 10, height = 6, dpi = 150, limitsize = FALSE)
#     ggsave(file.path(ws_fig, paste0("lollipop_", trait_name, "_", mod_color, ".svg")),
#            p, width = 10, height = 6, limitsize = FALSE)
#     system2("cp", args = c(file.path(ws_fig, paste0("lollipop_", trait_name, "_", mod_color, ".png")), fig_dir))
#     system2("cp", args = c(file.path(ws_fig, paste0("lollipop_", trait_name, "_", mod_color, ".svg")), fig_dir))
#     
#     logit(sprintf("    %s: %d significant GO terms", mod_color, nrow(ego_df)))
#   }
# }
# 
# # Combine enrichment results
# if (length(all_enrichment) > 0) {
#   enrichment_combined <- do.call(rbind, all_enrichment)
#   write.csv(enrichment_combined, file.path(out_dir, "wgcna_go_enrichment.csv"), row.names = FALSE)
#   
#   # Combined lollipop: top 5 per module-trait
#   top_combined <- do.call(rbind, lapply(
#     split(enrichment_combined, paste(enrichment_combined$module, enrichment_combined$trait, sep = "_")),
#     function(x) head(x[order(x$p.adjust), ], 5)
#   ))
#   top_combined$label <- paste0(top_combined$module, " [", top_combined$trait, "] | ",
#                                substr(top_combined$Description, 1, 50))
#   top_combined$label <- make.unique(top_combined$label)
#   top_combined$label <- factor(top_combined$label, levels = rev(top_combined$label))
#   
#   p_combined <- ggplot(top_combined, aes(x = Count, y = label, fill = -log10(p.adjust))) +
#     geom_segment(aes(x = 0, xend = Count, y = label, yend = label), color = "grey70") +
#     geom_point(size = 3.5, shape = 21, color = "black", stroke = 0.5) +
#     scale_fill_gradient(low = "#0279EE", high = "#FF9400", name = "-log10(padj)") +
#     labs(title = "GO BP Enrichment - Top 5 per Significant Module-Trait",
#          x = "Gene Count", y = "") +
#     theme_minimal(base_family = "sans") +
#     theme(text = element_text(size = 9), plot.title = element_text(face = "bold", size = 11))
#   
#   ggsave(file.path(ws_fig, "lollipop_combined.png"), p_combined,
#          width = 14, height = 16, dpi = 150, limitsize = FALSE)
#   ggsave(file.path(ws_fig, "lollipop_combined.svg"), p_combined,
#          width = 14, height = 16, dpi = 150, limitsize = FALSE)
#   system2("cp", args = c(file.path(ws_fig, "lollipop_combined.png"), fig_dir))
#   system2("cp", args = c(file.path(ws_fig, "lollipop_combined.svg"), fig_dir))
#   
#   logit(sprintf("  Saved: wgcna_go_enrichment.csv (%d terms), lollipop_combined.png/.svg",
#                 nrow(enrichment_combined)))
# } else {
#   logit("  No enrichment results to combine.")
# }
# 
# # Save combined core genes
# if (length(all_core_genes) > 0) {
#   core_all <- do.call(rbind, all_core_genes)
#   write.csv(core_all, file.path(out_dir, "wgcna_core_genes_filtered.csv"), row.names = FALSE)
#   logit(sprintf("  Saved: wgcna_core_genes_filtered.csv (%d genes)", nrow(core_all)))
# }

# --- Excel workbook ---
logit("")
logit("  --- Excel workbook ---")

wb <- createWorkbook()

# Summary sheet — all 6 traits
summary_df <- data.frame(
  Module = rownames(module_trait_cor),
  Module_Color = gsub("^ME", "", rownames(module_trait_cor)),
  Module_Size = as.numeric(table(module_colors)[match(
    gsub("^ME", "", rownames(module_trait_cor)), names(table(module_colors)))]),
  stringsAsFactors = FALSE)

for (trait_name in colnames(module_trait_cor)) {
  summary_df[[paste0("r_", trait_name)]] <- round(module_trait_cor[, trait_name], 4)
  summary_df[[paste0("p_", trait_name)]] <- signif(moduleTraitPvalue[, trait_name], 3)
  summary_df[[paste0("n_", trait_name)]] <- n_per_trait[trait_name]
  summary_df[[paste0("sig_", trait_name)]] <- ifelse(
    abs(module_trait_cor[, trait_name]) > r_threshold &
      moduleTraitPvalue[, trait_name] < p_threshold, "YES", "")
}

addWorksheet(wb, "Summary")
writeData(wb, "Summary", summary_df)

# One sheet per module (sorted by size, excluding grey)
modules <- unique(gene_info$module)
modules <- modules[modules != "grey"]
modules <- modules[order(sapply(modules, function(m) sum(gene_info$module == m)), decreasing = TRUE)]

for (mod in modules) {
  mod_genes <- gene_info[gene_info$module == mod, ]
  mod_genes <- mod_genes[order(-mod_genes$hub_score), ]
  me_col <- paste0("ME", mod)
  if (me_col %in% rownames(module_trait_cor)) {
    for (trait_name in colnames(module_trait_cor)) {
      mod_genes[[paste0("r_", trait_name)]] <- round(module_trait_cor[me_col, trait_name], 4)
      mod_genes[[paste0("p_", trait_name)]] <- signif(moduleTraitPvalue[me_col, trait_name], 3)
    }
  }
  sheet_name <- substr(mod, 1, 31)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, mod_genes)
}

# # GO enrichment sheet (exclude grey terms)
# if (length(all_enrichment) > 0) {
#   enrichment_clean <- enrichment_combined[enrichment_combined$module != "grey", ]
#   addWorksheet(wb, "GO_Enrichment")
#   writeData(wb, "GO_Enrichment", enrichment_clean)
# }
# 
# # Core genes sheet
# if (length(all_core_genes) > 0) {
#   addWorksheet(wb, "Core_Genes")
#   writeData(wb, "Core_Genes", core_all)
# }

# Save Excel (write to /workspace/ first, then cp to /mnt/results/)
tmp_xlsx <- "~/OSD152/OSD152_WGCNA/wgcna_gene_modules_by_module.xlsx"
saveWorkbook(wb, tmp_xlsx, overwrite = TRUE)
# system2("cp", args = c(tmp_xlsx, file.path(out_dir, "wgcna_gene_modules_by_module.xlsx")))

logit(sprintf("  Saved: wgcna_gene_modules_by_module.xlsx (Summary + %d module sheets)",
              length(modules)))

# --- Save RDS objects ---
logit("")
logit("  --- Save RDS objects ---")

trait_results <- list(
  cor = module_trait_cor,
  pvalue = moduleTraitPvalue,
  n_per_trait = n_per_trait,
  MEs = MEs)
saveRDS(trait_results, file.path(out_dir, "wgcna_trait_results.rds"))

full_results <- list(
  net = net,
  module_colors = module_colors,
  datExpr = datExpr,
  traits = traits,
  trait_results = trait_results,
  gene_info = gene_info,
  datKME = datKME,
  # enrichment = if (length(all_enrichment) > 0) enrichment_combined else NULL,
  # enrichment_excl_grey = if (length(all_enrichment) > 0) enrichment_clean else NULL,
  # core_genes = if (length(all_core_genes) > 0) core_all else NULL,
  parameters = list(
    power = soft_power,
    networkType = network_type,
    corType = cor_type,
    kme_threshold = kme_threshold,
    r_threshold = r_threshold,
    p_threshold = p_threshold,
    normalization = "LIMMA_Norm_TRUE"
  ))
saveRDS(full_results, file.path(out_dir, "wgcna_full_results.rds"))
logit("  Saved: wgcna_trait_results.rds, wgcna_full_results.rds")

###############################################################################
# Save consolidated log
###############################################################################
writeLines(log_lines, file.path(out_dir, "wgcna_phase2_log.txt"))
logit("")
logit("=== WGCNA Phase 2 COMPLETE ===")
