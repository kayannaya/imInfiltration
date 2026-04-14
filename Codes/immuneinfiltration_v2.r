# Immune Infiltration Analysis
# Modified from Google Colab (immuneinfiltration_v2)

# This uses data from Xenabrowser

# path configuration

# root folder
LOCAL_WD <- "/Users/kayannaya/Documents/Work/TEEP/imInfiltration2/BRCA"
DATASET_PATH <- "/Users/kayannaya/Documents/Work/TEEP/imInfiltration2/dataset"

# output csv directory
OUTPUT_DIR <- file.path(LOCAL_WD, "/outputs")

# path to cnv file
CNV_FILE <- file.path(DATASET_PATH, "TCGA.BRCA.sampleMap_Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz")

# path to survival file
SURVIVAL_FILE <- file.path(DATASET_PATH, "/01_cox/BRCA.csv")

# path to GSVA file
GENESETS_FILE <- file.path(DATASET_PATH, "geneSets.csv")

# path to BRCA exp file (for the deconvolution)
BRCA_FILE <- file.path(DATASET_PATH, "/brca_tcga_gdc/TCGA-BRCA.star_tpm.tsv")

# path to BRCA counts exp file (for DESeq2)
BRCA_COUNTS_FILE <- file.path(DATASET_PATH, "/brca_tcga_gdc/TCGA-BRCA.star_counts.tsv")

# path to BRCA RDS exp file
TCGA_FILE <- file.path(LOCAL_WD, "/TCGA_BRCA_se.rds")

# path to plot folder
PLOT_DIR <- file.path(OUTPUT_DIR, "plots")

# local directories

dir.create(LOCAL_WD,    showWarnings = FALSE, recursive = TRUE)
dir.create(OUTPUT_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(LOCAL_WD, "GDCdata"), showWarnings = FALSE, recursive = TRUE)
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

setwd(LOCAL_WD)
cat("Working directory:", getwd(), "\n")

# save to local output folder (replaces save_and_upload)

save_local <- function(df, filename) {
  out_path <- file.path(OUTPUT_DIR, filename)
  write.csv(df, out_path)
  cat("Saved:", out_path, "\n")
}

# install the required packages

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

bioc_pkgs <- c("TCGAbiolinks", "SummarizedExperiment", "SingleCellExperiment",
               "DESeq2", "limma", "GSVA", "clusterProfiler", "org.Hs.eg.db",
               "survminer")

cran_pkgs <- c("tidyverse", "survival", "dplyr", "purrr", "readxl", "sva", "remotes")

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, ask = FALSE)
}
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

if (!requireNamespace("immunedeconv", quietly = TRUE)) {
  remotes::install_github("omnideconv/immunedeconv", upgrade = "never")
}

# load the library

library(dplyr)
library(tidyverse)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(DESeq2)
library(limma)
library(GSVA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(survival)
library(survminer)
library(purrr)
library(readxl)
library(immunedeconv)


# CHANGE
# i'm changing the data, so i don't need to load it from GDC
# 1. download the TCGA data using GDC API

# delete the incomplete tar.gz chunks

# check if the TCGA_FILE is present or not
# if (!file.exists(TCGA_FILE)) {
#   # file not found → download from GDC
#   file.remove(list.files(pattern = "*.tar.gz", full.names = TRUE))
#   
#   options(TCGA_downloaded_data = file.path(LOCAL_WD, "GDCdata"))
#   
#   query_expr <- GDCquery(project       = "TCGA-BRCA",
#                          data.category = "Transcriptome Profiling",
#                          data.type     = "Gene Expression Quantification",
#                          workflow.type = "STAR - Counts")
#   
#   GDCdownload(query_expr, method = "api", files.per.chunk = 10)
#   
#   se <- GDCprepare(query_expr)
#   
#   saveRDS(se, file = TCGA_FILE)
#   cat("SE object saved locally.\n")
#   
# } else {
#   # immediately access the downloaded data from the previous session
#   se <- readRDS(TCGA_FILE)
#   cat("SE object loaded from:", TCGA_FILE, "\n")
# }


# 1. load the exp data
raw <- read.table(BRCA_FILE, header = TRUE, sep = "\t",
                  row.names = 1, check.names = FALSE)

expr_tpm_mrna <- as.matrix(raw)

# strip version numbers and deduplicate at the source
rownames(expr_tpm_mrna) <- sub("\\..*$", "", rownames(expr_tpm_mrna))
expr_tpm_mrna <- expr_tpm_mrna[order(-rowMeans(expr_tpm_mrna)), ]
expr_tpm_mrna <- expr_tpm_mrna[!duplicated(rownames(expr_tpm_mrna)), ]

symbol_mrna      <- rownames(expr_tpm_mrna)

cat("Loaded:", nrow(expr_tpm_mrna), "genes x", ncol(expr_tpm_mrna), "samples\n")

# 2. data pre-processing 

expr_tpm_mrna_symbol <- cbind(data.frame(symbol_mrna), as.data.frame(expr_tpm_mrna))

BRCA_tpm <- expr_tpm_mrna_symbol %>%
  as_tibble() %>%
  mutate(meanrow = rowMeans(.[, -1]), .before = 2) %>%
  filter(meanrow >= 1) %>%
  arrange(desc(meanrow)) %>%
  dplyr::distinct(symbol_mrna, .keep_all = TRUE) %>%
  dplyr::select(-meanrow) %>%
  column_to_rownames(var = "symbol_mrna") %>%
  as.data.frame()

colnames(BRCA_tpm) <- substr(colnames(BRCA_tpm), 1, 15)
BRCA_tpm           <- BRCA_tpm[, !duplicated(colnames(BRCA_tpm))]

selected_cols  <- colnames(BRCA_tpm)[substr(colnames(BRCA_tpm), 14, 15) == "01"]
clean_BRCA_tpm <- BRCA_tpm[, selected_cols]

cat("clean_BRCA_tpm dimensions:", dim(clean_BRCA_tpm), "\n")
cat("Sample type codes:\n")
print(table(substr(colnames(clean_BRCA_tpm), 14, 15)))

# 2.1 convert clean_BRCA_tpm rownames from Ensembl to symbols
id_map <- AnnotationDbi::select(org.Hs.eg.db,
                                keys    = rownames(clean_BRCA_tpm),
                                columns = "SYMBOL",
                                keytype = "ENSEMBL")

# drop NAs and duplicated Ensembl IDs
id_map <- id_map[!is.na(id_map$SYMBOL), ]
id_map <- id_map[!duplicated(id_map$ENSEMBL), ]

clean_BRCA_tpm           <- as.data.frame(clean_BRCA_tpm)
clean_BRCA_tpm$symbol    <- id_map$SYMBOL[match(rownames(clean_BRCA_tpm), id_map$ENSEMBL)]
clean_BRCA_tpm$meanexpr  <- rowMeans(clean_BRCA_tpm[, -ncol(clean_BRCA_tpm)], na.rm = TRUE)

clean_BRCA_tpm <- clean_BRCA_tpm[!is.na(clean_BRCA_tpm$symbol), ]
clean_BRCA_tpm <- clean_BRCA_tpm[order(-clean_BRCA_tpm$meanexpr), ]
clean_BRCA_tpm <- clean_BRCA_tpm[!duplicated(clean_BRCA_tpm$symbol), ]

rownames(clean_BRCA_tpm) <- clean_BRCA_tpm$symbol
clean_BRCA_tpm$symbol   <- NULL
clean_BRCA_tpm$meanexpr <- NULL

cat("clean_BRCA_tpm after symbol conversion:", dim(clean_BRCA_tpm), "\n")
cat("First 5 rownames:", head(rownames(clean_BRCA_tpm), 5), "\n")

# 3. batch effect removal
# skipping - barcodes are truncated to 16 chars (sample level only),
# no plate/batch info available. STAR FPKM-UQ is uniformly processed so this is fine.

cat("Skipping batch correction - single unified pipeline detected.\n")

dat        <- expr_tpm_mrna
dat_log2   <- log2(dat + 1)
expr_limma <- as.matrix(2^dat_log2 - 1)   # passthrough, no correction applied

rownames(expr_limma) <- symbol_mrna

expr_limma_df       <- as.data.frame(expr_limma)
expr_limma_df$genes <- rownames(expr_limma_df)

expr_limma_sym <- expr_limma_df %>%
  group_by(genes) %>%
  summarise_all(mean) %>%
  as.data.frame()

rownames(expr_limma_sym) <- expr_limma_sym$genes
expr_limma_sym           <- expr_limma_sym[, colnames(expr_limma_sym) != "genes"]
expr_limma_sym           <- na.omit(expr_limma_sym)

matrix_in_cols <- colnames(expr_limma_sym)[substr(colnames(expr_limma_sym), 14, 15) == "01"]
matrix_in      <- expr_limma_sym[, matrix_in_cols]

# subset matrix to mappable genes
matrix_in_sym <- matrix_in[rownames(matrix_in) %in% id_map$ENSEMBL, ]

# map symbols via temporary column (avoids duplicate rowname error)
matrix_in_sym           <- as.data.frame(matrix_in_sym)
matrix_in_sym$symbol    <- id_map$SYMBOL[match(rownames(matrix_in_sym), id_map$ENSEMBL)]
matrix_in_sym$meanexpr  <- rowMeans(matrix_in_sym[, -ncol(matrix_in_sym)])

# deduplicate symbols — keep highest mean expression
matrix_in_sym <- matrix_in_sym[order(-matrix_in_sym$meanexpr), ]
matrix_in_sym <- matrix_in_sym[!duplicated(matrix_in_sym$symbol), ]

# set symbols as rownames and clean up helper columns
rownames(matrix_in_sym) <- matrix_in_sym$symbol
matrix_in_sym$symbol    <- NULL
matrix_in_sym$meanexpr  <- NULL
matrix_in_sym           <- as.matrix(matrix_in_sym)

cat("matrix_in_sym dimensions:", dim(matrix_in_sym), "\n")
cat("First 5 gene names:", head(rownames(matrix_in_sym), 5), "\n")

# 4. TME deconvolution

TME <- function(expr, method) {
  imm           <- immunedeconv::deconvolute(gene_expression = expr, method)
  imm           <- t(imm)
  colnames(imm) <- imm[1, ]
  imm           <- imm[-1, ]
  imm           <- as.data.frame(imm)
  samples       <- rownames(imm)
  imm           <- sapply(imm, as.numeric)
  rownames(imm) <- samples
  return(imm)
}

frac_norm <- function(x) {
  x / rowSums(x, na.rm = TRUE)
}

cat("Running deconvolution... this will take a while.\n")

deconvoluted <- list()
for (method in c("quantiseq", "mcp_counter", "abis", "epic", "estimate")) {
  cat("  Running:", method, "\n")
  tryCatch({
    deconvoluted[[method]] <- TME(matrix_in_sym, method)
    cat("  Done:", method, "\n")
  }, error = function(e) cat("  Error in", method, ":", e$message, "\n"))
}

deconvoluted[["quantiseq_norm"]] <- frac_norm(deconvoluted[["quantiseq"]])
deconvoluted[["epic_norm"]]      <- frac_norm(deconvoluted[["epic"]])
deconvoluted[["quantiseq"]]      <- NULL
deconvoluted[["epic"]]           <- NULL

cat("\nABIS cell type columns:\n")
print(colnames(deconvoluted[["abis"]]))

for (method_name in names(deconvoluted)) {
  save_local(deconvoluted[[method_name]],
             paste0("BRCA_TIL_", method_name, ".csv"))
}

# reload fresh from the deconvoluted object
head(rownames(deconvoluted[["abis"]]), 5)
nchar(rownames(deconvoluted[["abis"]])[1])

# 5. load CNV data from local directory

# NOTE: Download this file from UCSC Xena and place it at CNV_FILE (set above)
# https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz

BRCA_cnv              <- read.delim(CNV_FILE, header = TRUE, check.names = FALSE)
rownames(BRCA_cnv)    <- BRCA_cnv$`Gene Symbol`
BRCA_cnv$`Gene Symbol` <- NULL
BRCA_cnv              <- BRCA_cnv[, !duplicated(colnames(BRCA_cnv))]

selected_cols  <- colnames(BRCA_cnv)[substr(colnames(BRCA_cnv), 14, 15) == "01"]
clean_BRCA_cnv <- BRCA_cnv[, selected_cols]

cat("clean_BRCA_cnv dimensions:", dim(clean_BRCA_cnv), "\n")
cat("TPM barcode format:", colnames(clean_BRCA_tpm)[1], "\n")
cat("CNV barcode format:", colnames(clean_BRCA_cnv)[1], "\n")
cat("Common samples:", length(intersect(colnames(clean_BRCA_tpm), colnames(clean_BRCA_cnv))), "\n")

# 6. cnv + expression grouping

common_cols <- intersect(colnames(clean_BRCA_tpm), colnames(clean_BRCA_cnv))
common_rows <- intersect(rownames(clean_BRCA_tpm), rownames(clean_BRCA_cnv))

clean_BRCA_tpm_common <- clean_BRCA_tpm[common_rows, common_cols]
clean_BRCA_cnv_common <- clean_BRCA_cnv[common_rows, common_cols]

cnv_data    <- clean_BRCA_cnv_common
expr_data   <- clean_BRCA_tpm_common
cnv_results <- list()

cat("Grouping", nrow(cnv_data), "genes by CNV and expression...\n")

for (gene in rownames(cnv_data)) {
  samples_neg_cnv <- which(cnv_data[gene, ] %in% c(-2, -1, 0))
  samples_pos_cnv <- which(cnv_data[gene, ] %in% c(1, 2))
  
  expr_values    <- expr_data[gene, ]
  ranked_samples <- rank(expr_values, na.last = "keep", ties.method = "average")
  
  n_samples           <- length(ranked_samples)
  bottom_30_threshold <- 0.3 * n_samples
  top_30_threshold    <- 0.7 * n_samples
  
  samples_neg_and_bottom_30 <- names(ranked_samples[samples_neg_cnv][ranked_samples[samples_neg_cnv] <= bottom_30_threshold])
  samples_pos_and_top_30    <- names(ranked_samples[samples_pos_cnv][ranked_samples[samples_pos_cnv] > top_30_threshold])
  
  cnv_results[[gene]] <- list(neg_cnv_bottom_30 = samples_neg_and_bottom_30,
                              pos_cnv_top_30    = samples_pos_and_top_30)
}

filtered_cnv_results   <- keep(cnv_results, ~ length(.$neg_cnv_bottom_30) >= 30 & length(.$pos_cnv_top_30) >= 30)
final_BRCA_cnv_results <- discard(filtered_cnv_results, ~ any(is.na(unlist(.x))))

cat("Genes passing CNV filter:", length(final_BRCA_cnv_results), "\n")

# 7. correlation with TIL abundance (ABIS / CD8+)

BRCA_TIL_abis           <- deconvoluted[["abis"]]
rownames(BRCA_TIL_abis) <- substr(rownames(BRCA_TIL_abis), 1, 15)

cd8_col_names <- grep("CD8", colnames(BRCA_TIL_abis), ignore.case = TRUE, value = TRUE)
cat("CD8-related ABIS columns:\n")
print(cd8_col_names)

# merge TIL and tpm data
clean_BRCA_tpm_filtered   <- subset(clean_BRCA_tpm,
                                    rownames(clean_BRCA_tpm) %in% names(final_BRCA_cnv_results))
clean_BRCA_tpm_filtered_t <- as.data.frame(t(clean_BRCA_tpm_filtered))

rownames(clean_BRCA_tpm_filtered_t) <- trimws(rownames(clean_BRCA_tpm_filtered_t))
rownames(BRCA_TIL_abis)             <- trimws(rownames(BRCA_TIL_abis))

# bug handling 
# note this part

# confirm TIL columns are numeric
cat("TIL col class:", class(BRCA_TIL_abis[, "T cell CD8+ memory"]), "\n")
cat("TIL col NAs:", sum(is.na(BRCA_TIL_abis[, "T cell CD8+ memory"])), "\n")
cat("BRCA_tpm_TIL dims:", dim(BRCA_TIL_abis), "\n")                                   

# overlapping barcodes check
sum(rownames(clean_BRCA_tpm_filtered_t) %in% rownames(BRCA_TIL_abis))

# check for hidden characters
nchar(rownames(clean_BRCA_tpm_filtered_t)[1])
nchar(rownames(BRCA_TIL_abis)[1])

# convert to ASCII format
rownames(clean_BRCA_tpm_filtered_t) <- iconv(trimws(rownames(clean_BRCA_tpm_filtered_t)),
                                             to = "ASCII", sub = "-")
rownames(BRCA_TIL_abis)             <- iconv(trimws(rownames(BRCA_TIL_abis)),
                                             to = "ASCII", sub = "-")

chartr("", "", rownames(clean_BRCA_tpm_filtered_t)[1])  # TPM
chartr("", "", rownames(BRCA_TIL_abis)[1])              # ABIS

# or more directly
utf8ToInt(rownames(clean_BRCA_tpm_filtered_t)[1])
utf8ToInt(rownames(BRCA_TIL_abis)[1])

# try a tiny manual merge with just 3 rows each
test_tpm  <- clean_BRCA_tpm_filtered_t[1:3, 1:3]
test_abis <- BRCA_TIL_abis[1:3, 1:3]

rownames(test_tpm)
rownames(test_abis)

merge(test_tpm, test_abis, by = "row.names", all = FALSE)

# check for duplicate rownames in either table
sum(duplicated(rownames(clean_BRCA_tpm_filtered_t)))
sum(duplicated(rownames(BRCA_TIL_abis)))

# THE KEY is to deduplicate the BRCA_TIL_abis

# deduplicate the 16 samples
BRCA_TIL_abis <- BRCA_TIL_abis[!duplicated(rownames(BRCA_TIL_abis)), ]
cat("BRCA_TIL_abis rows after dedup:", nrow(BRCA_TIL_abis), "\n")

# merge
BRCA_tpm_TIL <- merge(clean_BRCA_tpm_filtered_t, BRCA_TIL_abis,
                      by = "row.names", all = FALSE)

rownames(BRCA_tpm_TIL) <- BRCA_tpm_TIL$Row.names
BRCA_tpm_TIL$Row.names <- NULL

gene_col_idx <- 1:ncol(clean_BRCA_tpm_filtered_t)

# check dimensions before and after merge
cat("clean_BRCA_tpm_filtered_t rows:", nrow(clean_BRCA_tpm_filtered_t), "\n")
cat("BRCA_TIL_abis rows:", nrow(BRCA_TIL_abis), "\n")
cat("BRCA_tpm_TIL rows:", nrow(BRCA_tpm_TIL), "\n")

# save to local
save_local(BRCA_tpm_TIL, "BRCA_tpm_TIL.csv")

# if still error
# inspect barcode formats
cat("TPM barcodes (first 3):\n")
print(head(rownames(clean_BRCA_tpm_filtered_t), 3))
cat("ABIS barcodes (first 3):\n")
print(head(rownames(BRCA_TIL_abis), 3))

cor_results_all <- lapply(cd8_col_names, function(til_col) {
  til_vec <- as.numeric(BRCA_tpm_TIL[, til_col])  # force numeric
  res <- sapply(BRCA_tpm_TIL[, gene_col_idx], function(gene_vec) {
    gene_vec <- as.numeric(gene_vec)
    valid <- is.finite(gene_vec) & is.finite(til_vec)
    if (sum(valid) < 3) return(c(correlation = NA, p.value = NA))
    test <- cor.test(gene_vec[valid], til_vec[valid], method = "spearman")
    c(correlation = test$estimate, p.value = test$p.value)
  })
  df            <- as.data.frame(t(res))
  df$adjusted.p <- p.adjust(df$p.value, method = "BH")
  df$til_column <- til_col
  df
})
names(cor_results_all) <- cd8_col_names

passing_genes_TIL <- unique(unlist(lapply(cor_results_all, function(df) {
  rownames(df)[df$correlation.rho < -0.2 & df$adjusted.p < 0.01]
})))

til_vec <- as.numeric(BRCA_tpm_TIL[, "T cell CD8+ memory"])

cat("Genes passing TIL correlation filter:", length(passing_genes_TIL), "\n")

primary_cd8_col <- cd8_col_names[grep("activ", cd8_col_names, ignore.case = TRUE)]
if (length(primary_cd8_col) == 0) primary_cd8_col <- cd8_col_names[1]

BRCA_expr_T_activated          <- cor_results_all[[primary_cd8_col]]
BRCA_expr_T_activated_filtered <- BRCA_expr_T_activated[
  rownames(BRCA_expr_T_activated) %in% passing_genes_TIL, ]

save_local(BRCA_expr_T_activated,          "BRCA_expr_T_activated.csv")
save_local(BRCA_expr_T_activated_filtered, "BRCA_expr_T_activated_filtered.csv")

# 8. correlation with immune markers (GSVA)

# place geneSets.csv at GENESETS_FILE (set above in the configurations)
# genesets can be found in the original paper
gene_list_gs <- read.csv(GENESETS_FILE, header = TRUE, stringsAsFactors = FALSE)
geneSets     <- list(gene_list_gs$genes)

expr_matrix <- as.matrix(clean_BRCA_tpm)
gsvaPar     <- gsvaParam(expr_matrix, geneSets, kcdf = "Gaussian")
gsva.es     <- gsva(gsvaPar, verbose = FALSE)
rownames(gsva.es) <- "gsva.es"

clean_BRCA_tpm_filtered2 <- subset(clean_BRCA_tpm,
                                   rownames(clean_BRCA_tpm) %in% passing_genes_TIL)
BRCA_tpm_ssGSEA   <- rbind(gsva.es, clean_BRCA_tpm_filtered2)
BRCA_tpm_ssGSEA_t <- t(BRCA_tpm_ssGSEA)

results_list_gsva <- list()
for (i in 2:ncol(BRCA_tpm_ssGSEA_t)) {
  test <- cor.test(BRCA_tpm_ssGSEA_t[, i], BRCA_tpm_ssGSEA_t[, 1], method = "spearman")
  results_list_gsva[[colnames(BRCA_tpm_ssGSEA_t)[i]]] <-
    c(cor_coefficient = test$estimate, p_value = test$p.value)
}

BRCA_expr_immune_markers            <- as.data.frame(do.call(rbind, results_list_gsva))
BRCA_expr_immune_markers$adjusted.p <- p.adjust(BRCA_expr_immune_markers$p_value, method = "BH")
BRCA_expr_immune_markers_filtered   <- BRCA_expr_immune_markers[
  BRCA_expr_immune_markers$cor_coefficient.rho < -0.20 &
    BRCA_expr_immune_markers$adjusted.p < 0.01, ]

cat("Genes passing immune marker filter:", nrow(BRCA_expr_immune_markers_filtered), "\n")

save_local(BRCA_expr_immune_markers,          "BRCA_expr_immune_markers.csv")
save_local(BRCA_expr_immune_markers_filtered, "BRCA_expr_immune_markers_filtered.csv")

# 9. multivariate survival analysis

BRCA_expr_survival_coef <- read_csv(SURVIVAL_FILE, show_col_types = FALSE)
BRCA_expr_survival_coef <- as.data.frame(BRCA_expr_survival_coef)
rownames(BRCA_expr_survival_coef) <- BRCA_expr_survival_coef$gene
BRCA_expr_survival_coef <- dplyr::select(BRCA_expr_survival_coef, coef, p)

BRCA_expr_survival_coef <- subset(BRCA_expr_survival_coef,
                                  rownames(BRCA_expr_survival_coef) %in% rownames(BRCA_expr_immune_markers_filtered))
BRCA_expr_survival_coef$adjusted.p <- p.adjust(BRCA_expr_survival_coef$p, method = "BH")

BRCA_expr_survival_multivariate_filtered <- subset(BRCA_expr_survival_coef, coef > 0.15)

cat("Genes passing multivariate survival filter:", nrow(BRCA_expr_survival_multivariate_filtered), "\n")

save_local(BRCA_expr_survival_coef,                   "BRCA_expr_survival_coef.csv")
save_local(BRCA_expr_survival_multivariate_filtered,  "BRCA_expr_survival_multivariate_filtered.csv")

# 10. kaplan-meier survival analysis

BRCA_clinical <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")
survival_data <- BRCA_clinical[, c("bcr_patient_barcode", "days_to_last_follow_up", "vital_status")]
colnames(survival_data) <- c("patient_id", "OS.time", "OS.status")
survival_data$OS.status <- ifelse(survival_data$OS.status == "Dead", 1, 0)

clean_BRCA_tpm_temp <- subset(clean_BRCA_tpm,
                              rownames(clean_BRCA_tpm) %in% rownames(BRCA_expr_survival_multivariate_filtered))
colnames(clean_BRCA_tpm_temp) <- substr(colnames(clean_BRCA_tpm_temp), 1, 12)

BRCA_surv_results <- data.frame(Gene = character(), P_value = numeric(), stringsAsFactors = FALSE)

for (gene in rownames(clean_BRCA_tpm_temp)) {
  tryCatch({
    gene.df            <- as.data.frame(t(subset(clean_BRCA_tpm_temp,
                                                 rownames(clean_BRCA_tpm_temp) == gene)))
    gene.df$patient_id <- rownames(gene.df)
    merged.data        <- na.omit(merge(survival_data, gene.df, by = "patient_id"))
    cutpoint           <- surv_cutpoint(merged.data, time = "OS.time",
                                        event = "OS.status", variables = gene)
    if (is.na(cutpoint$cutpoint[[1]])) { next }
    merged.data$group  <- ifelse(merged.data[, gene] > cutpoint$cutpoint[[1]], "High", "Low")
    logrank_test       <- survdiff(Surv(OS.time, OS.status) ~ group, data = merged.data)
    BRCA_surv_results  <- rbind(BRCA_surv_results,
                                data.frame(Gene = gene, P_value = logrank_test$pval))
  }, error = function(e) cat("Error for", gene, ":", e$message, "\n"))
}

BRCA_surv_results$adjusted_P_value <- p.adjust(BRCA_surv_results$P_value, method = "BH")
BRCA_expr_survival_KM              <- BRCA_surv_results
rownames(BRCA_expr_survival_KM)    <- BRCA_expr_survival_KM[, 1]
BRCA_expr_survival_KM$Gene        <- NULL
BRCA_expr_survival_KM_filtered     <- BRCA_expr_survival_KM[
  BRCA_expr_survival_KM$adjusted_P_value < 0.05, ]

cat("Genes passing KM filter:", nrow(BRCA_expr_survival_KM_filtered), "\n")

save_local(BRCA_expr_survival_KM,          "BRCA_expr_survival_KM.csv")
save_local(BRCA_expr_survival_KM_filtered, "BRCA_expr_survival_KM_filtered.csv")

# prep clinical + expression data (same as step 10 in main pipeline)
BRCA_clinical <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")
survival_data <- BRCA_clinical[, c("bcr_patient_barcode", "days_to_last_follow_up", "vital_status")]
colnames(survival_data) <- c("patient_id", "OS.time", "OS.status")
survival_data$OS.status <- ifelse(survival_data$OS.status == "Dead", 1, 0)

# 11. merge all filters

BRCA_expr_T_activated_filtered <- BRCA_expr_T_activated_filtered %>%
  rename_with(~ paste0("expr_T_activated_", .), everything())
BRCA_expr_immune_markers_filtered <- BRCA_expr_immune_markers_filtered %>%
  rename_with(~ paste0("expr_immune_markers_", .), everything())
BRCA_expr_survival_multivariate_filtered <- BRCA_expr_survival_multivariate_filtered %>%
  rename_with(~ paste0("expr_survival_multivariate_", .), everything())
BRCA_expr_survival_KM_filtered <- BRCA_expr_survival_KM_filtered %>%
  rename_with(~ paste0("expr_survival_KM_", .), everything())

surviving_genes <- rownames(BRCA_expr_survival_KM_filtered)
for (df_name in c("BRCA_expr_T_activated_filtered",
                  "BRCA_expr_immune_markers_filtered",
                  "BRCA_expr_survival_multivariate_filtered")) {
  df <- get(df_name)
  assign(df_name, subset(df, rownames(df) %in% surviving_genes))
}

BRCA_gene_list <- Reduce(function(x, y) {
  m           <- merge(x, y, by = "row.names", all = FALSE)
  rownames(m) <- m$Row.names
  m$Row.names <- NULL
  m
}, list(BRCA_expr_T_activated_filtered,
        BRCA_expr_immune_markers_filtered,
        BRCA_expr_survival_multivariate_filtered,
        BRCA_expr_survival_KM_filtered))

cat("Final gene list size:", nrow(BRCA_gene_list), "\n")
save_local(BRCA_gene_list, "BRCA_gene_list.csv")

# 12. DESeq2 + GSEA

# load raw counts (separate file)
raw_counts <- read.table(BRCA_COUNTS_FILE, header = TRUE, sep = "\t",
                         row.names = 1, check.names = FALSE)

# strip version numbers and deduplicate
# the before code just removed the version numbers and make the corrected rows to have the same name
# example: ENSG00000002586.1 and ENSG00000002586.2 both become ENSG00000002586
# the new code keeps the version with the highest expressed version
raw_counts$gene_id <- sub("\\..*$", "", rownames(raw_counts))
raw_counts$meanrow <- rowMeans(raw_counts[, !colnames(raw_counts) %in% c("gene_id", "meanrow")])
raw_counts         <- raw_counts[order(-raw_counts$meanrow), ]
raw_counts         <- raw_counts[!duplicated(raw_counts$gene_id), ]
rownames(raw_counts) <- raw_counts$gene_id
raw_counts$gene_id <- NULL
raw_counts$meanrow <- NULL

# convert Ensembl to symbols using id_map from section 3
raw_counts_sym          <- as.data.frame(raw_counts)
raw_counts_sym$symbol   <- id_map$SYMBOL[match(rownames(raw_counts_sym), id_map$ENSEMBL)]
raw_counts_sym$meanexpr <- rowMeans(raw_counts_sym[, !colnames(raw_counts_sym) %in% "symbol"],
                                    na.rm = TRUE)
raw_counts_sym <- raw_counts_sym[!is.na(raw_counts_sym$symbol), ]
raw_counts_sym <- raw_counts_sym[order(-raw_counts_sym$meanexpr), ]
raw_counts_sym <- raw_counts_sym[!duplicated(raw_counts_sym$symbol), ]
rownames(raw_counts_sym) <- raw_counts_sym$symbol
raw_counts_sym$symbol    <- NULL
raw_counts_sym$meanexpr  <- NULL

# filter, subset to tumour samples, fix barcodes
BRCA_counts <- raw_counts_sym[rowMeans(raw_counts_sym) >= 10, ]
colnames(BRCA_counts) <- substr(colnames(BRCA_counts), 1, 15)
BRCA_counts           <- BRCA_counts[, !duplicated(colnames(BRCA_counts))]
selected_cols         <- colnames(BRCA_counts)[substr(colnames(BRCA_counts), 14, 15) == "01"]
clean_BRCA_counts     <- BRCA_counts[, selected_cols]

# convert to integer for DESeq2
clean_BRCA_counts <- round(clean_BRCA_counts)
rn                <- rownames(clean_BRCA_counts)
clean_BRCA_counts <- as.data.frame(lapply(clean_BRCA_counts, as.integer))
rownames(clean_BRCA_counts) <- rn

# standardize separators to match final_BRCA_cnv_results (hyphens)
colnames(clean_BRCA_counts) <- gsub("\\.", "-", colnames(clean_BRCA_counts))

cat("clean_BRCA_counts dimensions:", dim(clean_BRCA_counts), "\n")
# expect thousands of genes x ~1000 samples

BRCA_screened_list  <- rownames(BRCA_gene_list)
immune_pathway_rows <- list()

for (gene in BRCA_screened_list) {
  cat("Processing:", gene, "\n")
  tryCatch({
    group1_samples <- final_BRCA_cnv_results[[gene]][[1]]
    group2_samples <- final_BRCA_cnv_results[[gene]][[2]]
    
    subset_counts <- clean_BRCA_counts[,
                                       colnames(clean_BRCA_counts) %in% c(group1_samples, group2_samples)]
    
    sample_info       <- data.frame(
      sample_name = colnames(subset_counts),
      group       = ifelse(colnames(subset_counts) %in% group1_samples, "Group1", "Group2")
    )
    sample_info$group <- relevel(as.factor(sample_info$group), ref = "Group2")
    
    dds <- DESeqDataSetFromMatrix(countData = subset_counts,
                                  colData   = sample_info,
                                  design    = ~ group)
    dds <- DESeq(dds)
    res <- results(dds)
    
    gene_log2FC        <- data.frame(gene           = rownames(res),
                                     log2FoldChange = res$log2FoldChange)
    gene_log2FC        <- gene_log2FC[order(-gene_log2FC$log2FoldChange), ]
    gene_log2FC        <- gene_log2FC[!is.na(gene_log2FC$log2FoldChange), ]
    vec                <- setNames(gene_log2FC$log2FoldChange, gene_log2FC$gene)
    
    gsea_result <- gseGO(geneList      = vec,
                         OrgDb         = org.Hs.eg.db,
                         ont           = "ALL",
                         keyType       = "SYMBOL",
                         minGSSize     = 10,
                         maxGSSize     = 500,
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         verbose       = FALSE)
    
    immune_term_count      <- 0
    antigen_MHC_term_count <- 0
    
    for (ont in c("BP", "CC", "MF")) {
      filtered_result <- gsea_result[gsea_result$ONTOLOGY == ont &
                                       gsea_result$pvalue < 0.05 &
                                       gsea_result$p.adjust < 0.25, ]
      immune_hits     <- grep("cytokine|interferon|interleukin|chemokine",
                              filtered_result$Description, ignore.case = TRUE)
      immune_hits     <- immune_hits[!grepl("cytokinesis|cytokinetic",
                                            filtered_result$Description[immune_hits],
                                            ignore.case = TRUE)]
      immune_term_count      <- immune_term_count + length(immune_hits)
      antigen_MHC_term_count <- antigen_MHC_term_count +
        length(grep("antigen|MHC", filtered_result$Description, ignore.case = TRUE))
    }
    
    immune_pathway_rows[[gene]] <- data.frame(
      gene                   = gene,
      immune_term_count      = immune_term_count,
      antigen_MHC_term_count = antigen_MHC_term_count
    )
    
  }, error = function(e) cat("  Error for", gene, ":", e$message, "\n"))
}

immune_pathway_df           <- do.call(rbind, immune_pathway_rows)
rownames(immune_pathway_df) <- immune_pathway_df$gene
immune_pathway_df           <- immune_pathway_df[, -1]
immune_pathway_df           <- immune_pathway_df %>%
  rename_with(~ paste0("immune_pathway_", .), everything())

save_local(immune_pathway_df, "BRCA_immune_pathway_df.csv")

BRCA_gene_list_final <- merge(BRCA_gene_list, immune_pathway_df, by = "row.names")
rownames(BRCA_gene_list_final) <- BRCA_gene_list_final$Row.names  # fix rownames
BRCA_gene_list_final$Row.names <- NULL                            # clean up column
save_local(BRCA_gene_list_final, "BRCA_gene_list_final.csv")

# the KM plot should be here as the gene list are defined latter in the code
# previous running works bc i added the km plot maker AFTER i ran the gene list part as well
# my mistake was adding it into the survival analysis step (still in the filter part) instead of adding it to the result step

# subset tpm to final candidate genes only
clean_BRCA_tpm_plot <- subset(clean_BRCA_tpm,
                              rownames(clean_BRCA_tpm) %in% rownames(BRCA_gene_list_final))
colnames(clean_BRCA_tpm_plot) <- substr(colnames(clean_BRCA_tpm_plot), 1, 12)

# plot one KM curve per gene
km_plots <- list()

for (gene in rownames(clean_BRCA_tpm_plot)) {
  tryCatch({
    # prepare gene expression + survival data
    gene.df            <- as.data.frame(t(clean_BRCA_tpm_plot[gene, , drop = FALSE]))
    colnames(gene.df)  <- gene
    gene.df$patient_id <- rownames(gene.df)
    merged.data        <- na.omit(merge(survival_data, gene.df, by = "patient_id"))
    
    # find optimal cutpoint
    cutpoint          <- surv_cutpoint(merged.data, time = "OS.time",
                                       event = "OS.status", variables = gene)
    if (is.na(cutpoint$cutpoint[[1]])) next
    merged.data$group <- ifelse(merged.data[, gene] > cutpoint$cutpoint[[1]], "High", "Low")
    
    # log-rank p-value
    logrank_test <- survdiff(Surv(OS.time, OS.status) ~ group, data = merged.data)
    p_val        <- signif(logrank_test$pval, 3)
    
    # KM fit
    fit <- survfit(Surv(OS.time, OS.status) ~ group, data = merged.data)
    
    # plot
    km_plots[[gene]] <- ggsurvplot(
      fit,
      data          = merged.data,
      title         = gene,
      pval          = TRUE,
      pval.method   = TRUE,
      conf.int      = FALSE,
      risk.table    = TRUE,
      risk.table.height = 0.25,
      palette       = c("#E7524A", "#4A90D9"),
      legend.labs   = c("High", "Low"),
      legend.title  = "Expression",
      xlab          = "Days",
      ylab          = "Overall survival probability",
      ggtheme       = theme_classic(base_size = 12)
    )
    
    cat("KM plot done:", gene, "\n")
    
  }, error = function(e) cat("Error for", gene, ":", e$message, "\n"))
}

# save individual KM plots
for (gene in names(km_plots)) {
  out_path <- file.path(PLOT_DIR, paste0("KM_", gene, ".pdf"))
  pdf(out_path, width = 7, height = 7)
  print(km_plots[[gene]])
  dev.off()
  cat("Saved:", out_path, "\n")
}

# save all KM plots in one combined PDF
combined_path <- file.path(PLOT_DIR, "KM_all_genes.pdf")
pdf(combined_path, width = 7, height = 7)
for (gene in names(km_plots)) print(km_plots[[gene]])
dev.off()
cat("Saved combined KM PDF:", combined_path, "\n")

# 13. gene ranking
# scoring logic:
#   candidate pool = union of all genes that passed any filter
#   each filter contributes a normalized 0-1 score based on effect strength
#   missing = gene didn't pass that filter, treated as 0
#   final score = mean of all normalized scores (higher = better candidate)

# helper: normalize a vector to 0-1 range
norm01 <- function(x) {
  if (all(is.na(x)) || diff(range(x, na.rm = TRUE)) == 0) return(rep(0, length(x)))
  (x - min(x, na.rm = TRUE)) / diff(range(x, na.rm = TRUE))
}

# load filter result CSVs
T_activated  <- read.csv(file.path(OUTPUT_DIR, "BRCA_expr_T_activated.csv"),         row.names = 1)
immune_mark  <- read.csv(file.path(OUTPUT_DIR, "BRCA_expr_immune_markers.csv"),       row.names = 1)
surv_coef    <- read.csv(file.path(OUTPUT_DIR, "BRCA_expr_survival_coef.csv"),        row.names = 1)
surv_KM      <- read.csv(file.path(OUTPUT_DIR, "BRCA_expr_survival_KM.csv"),          row.names = 1)
immune_path  <- read.csv(file.path(OUTPUT_DIR, "BRCA_immune_pathway_df.csv"),         row.names = 1)

# candidate pool = union of all genes across all filters
candidate_genes <- unique(c(
  rownames(T_activated),
  rownames(immune_mark),
  rownames(surv_coef),
  rownames(surv_KM),
  rownames(immune_path)
))
cat("Total candidate genes for ranking:", length(candidate_genes), "\n")

# build scoring table, one row per candidate gene
ranking           <- data.frame(gene = candidate_genes, stringsAsFactors = FALSE)
rownames(ranking) <- candidate_genes

# score 1: TIL correlation strength (more negative rho = stronger immune exclusion)
ranking$score_TIL <- ifelse(candidate_genes %in% rownames(T_activated),
                            -T_activated[candidate_genes, "correlation.rho"],
                            NA)

# score 2: GSVA immune marker correlation strength (more negative = better)
ranking$score_GSVA <- ifelse(candidate_genes %in% rownames(immune_mark),
                             -immune_mark[candidate_genes, "cor_coefficient.rho"],
                             NA)

# score 3: multivariate survival coefficient (higher coef = worse prognosis = more relevant)
ranking$score_surv_coef <- ifelse(candidate_genes %in% rownames(surv_coef),
                                  surv_coef[candidate_genes, "coef"],
                                  NA)

# score 4: KM significance (-log10 adjusted p, higher = more significant)
ranking$score_KM <- ifelse(candidate_genes %in% rownames(surv_KM),
                           -log10(surv_KM[candidate_genes, "adjusted_P_value"] + 1e-10),
                           NA)

# score 5: immune pathway hits from GSEA (more hits = more immune-relevant)
ranking$score_immune_pathway <- ifelse(candidate_genes %in% rownames(immune_path),
                                       immune_path[candidate_genes, "immune_pathway_immune_term_count"] +
                                         immune_path[candidate_genes, "immune_pathway_antigen_MHC_term_count"],
                                       NA)

# normalize each score to 0-1 across all candidates (NA stays NA, becomes 0 after)
score_cols <- c("score_TIL", "score_GSVA", "score_surv_coef", "score_KM", "score_immune_pathway")
for (col in score_cols) {
  norm_col                <- norm01(ranking[[col]])
  norm_col[is.na(ranking[[col]])] <- 0  # genes that didn't reach this filter score 0
  ranking[[paste0(col, "_norm")]] <- norm_col
}

# final score = mean of normalized scores
norm_cols           <- paste0(score_cols, "_norm")
ranking$final_score <- rowMeans(ranking[, norm_cols], na.rm = TRUE)

# track how many filters each gene passed
ranking$filters_passed <- rowSums(!is.na(ranking[, score_cols]))

# sort by final score descending
ranking <- ranking[order(-ranking$final_score), ]
ranking$rank <- seq_len(nrow(ranking))

cat("\nTop 10 ranked genes:\n")
print(ranking[1:10, c("rank", "gene", "final_score", "filters_passed", norm_cols)])

# save ranking table
save_local(ranking, "BRCA_gene_ranking.csv")

# ranking visualization

# reshape to long format for ggplot
score_labels <- c(
  score_TIL_norm            = "TIL correlation",
  score_GSVA_norm           = "GSVA immune score",
  score_surv_coef_norm      = "Survival coefficient",
  score_KM_norm             = "KM significance",
  score_immune_pathway_norm = "Immune pathway hits"
)

ranking_long <- ranking %>%
  dplyr::select(gene, rank, all_of(norm_cols)) %>%
  pivot_longer(cols = all_of(norm_cols), names_to = "filter", values_to = "score") %>%
  mutate(
    filter   = recode(filter, !!!score_labels),
    gene     = factor(gene, levels = rev(ranking$gene))  # top ranked at top
  )

p_ranking <- ggplot(ranking_long, aes(x = filter, y = gene, size = score, color = score)) +
  geom_point() +
  scale_size_continuous(range = c(1, 8), name = "Normalized\nscore") +
  scale_color_gradient(low = "#D0E8F5", high = "#1A5276", name = "Normalized\nscore") +
  labs(
    title    = "Gene ranking by filter scores",
    subtitle = "Each dot = normalized score for that filter (larger/darker = stronger)",
    x        = NULL,
    y        = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x  = element_text(angle = 35, hjust = 1),
    axis.text.y  = element_text(size = 9),
    plot.title   = element_text(face = "bold"),
    legend.position = "right"
  )

# save ranking plot
ranking_plot_path <- file.path(PLOT_DIR, "gene_ranking_dotplot.pdf")
ggsave(ranking_plot_path, p_ranking,
       width  = min(20, max(6, nrow(ranking) * 0.35)),
       height = 5)
cat("Saved:", ranking_plot_path, "\n")

# also save a simple bar chart of final scores
p_bar <- ggplot(ranking, aes(x = reorder(gene, final_score), y = final_score)) +
  geom_col(fill = "#1A5276", alpha = 0.85) +
  coord_flip() +
  labs(
    title = "Overall gene ranking",
    x     = NULL,
    y     = "Final score (mean of normalized filter scores)"
  ) +
  theme_classic(base_size = 11)
theme(plot.title = element_text(face = "bold"))

bar_plot_path <- file.path(PLOT_DIR, "gene_ranking_barplot.pdf")
ggsave(bar_plot_path, p_bar,
       width  = min(20, max(6, nrow(ranking) * 0.35)),
       height = 5)
cat("Saved:", bar_plot_path, "\n")

cat("\nAll plots saved to:", PLOT_DIR, "\n")

cat("\nDone! All results saved to:", OUTPUT_DIR, "\n")

# checking where did LINC00707 got erased
target_ensembl <- "ENSG00000238266"
target_symbol  <- "LINC00707"

cat("1. In raw expr_tpm_mrna:", target_ensembl %in% rownames(expr_tpm_mrna), "\n")
cat("2. In BRCA_tpm (after meanrow>=1):", target_ensembl %in% rownames(BRCA_tpm), "\n")
cat("3. In clean_BRCA_tpm (after symbol map):", target_symbol %in% rownames(clean_BRCA_tpm), "\n")
cat("4. In final_BRCA_cnv_results (CNV filter):", target_symbol %in% names(final_BRCA_cnv_results), "\n")
cat("5. In passing_genes_TIL (TIL filter):", target_symbol %in% passing_genes_TIL, "\n")
cat("6. In immune_markers_filtered (GSVA filter):", target_symbol %in% rownames(BRCA_expr_immune_markers_filtered), "\n")
cat("7. In survival_multivariate_filtered:", target_symbol %in% rownames(BRCA_expr_survival_multivariate_filtered), "\n")
cat("8. In KM_filtered:", target_symbol %in% rownames(BRCA_expr_survival_KM_filtered), "\n")
cat("9. In final gene list:", target_symbol %in% rownames(BRCA_gene_list_final), "\n")