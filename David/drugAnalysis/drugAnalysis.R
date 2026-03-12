# gene Expression vs drug response box plot (TCGA)

# don't forget to install these packages first
library(tidyverse)
library(dplyr)
library(ggplot2)
library(org.Hs.eg.db)

select    <- dplyr::select
filter    <- dplyr::filter
mutate    <- dplyr::mutate
summarise <- dplyr::summarise

# configurations

# you need to download clinical treatment and gene expression data. i got both from cBioPortal -> TCGA
TREATMENT_FILE        <- "/Users/kayannaya/Documents/Work/TEEP/David/data/dataColon/coad_tcga_gdc/data_timeline_treatment.txt"
EXPRESSION_FILE       <- "/Users/kayannaya/Documents/Work/TEEP/David/data/dataColon/coad_tcga_gdc/data_mrna_seq_read_counts.txt"
GENES_OF_INTEREST     <- c("TEX36-AS1")
MIN_SAMPLES_PER_GROUP <- 3
OUTPUT_DIR            <- "/Users/kayannaya/Documents/Work/TEEP/David/drugAnalysis/output"

dir.create(OUTPUT_DIR, showWarnings = FALSE)

# load treatment file

treatment_raw <- read_tsv(TREATMENT_FILE, col_types = cols(.default = "c"))

# pre-processing treatment data

treatment <- treatment_raw %>%
  # drop other columns except these columns
  select(PATIENT_ID, THERAPEUTIC_AGENT, TREATMENT_OUTCOME) %>%
  
  # removes any accidental leading or trailing whitespace from every column 
  # — so for example if a cell in your data contains " complete response" (with a space at the start) or "Fluorouracil " (with a space at the end), 
  # it gets cleaned to "complete response" and "Fluorouracil".
  mutate(across(everything(), str_trim)) %>%
  
  # filter the rows that has both columns' value present
  filter(!is.na(THERAPEUTIC_AGENT), !is.na(TREATMENT_OUTCOME), TREATMENT_OUTCOME != "") %>%
  
  mutate(
    
    # shorten the PATIENT_ID barcode into 12 characters only
    barcode = str_sub(PATIENT_ID, 1, 12),
    
    # groups the TREATMENT_OUTCOME into two groups
    RESPONSE_GROUP = case_when(
      str_to_lower(TREATMENT_OUTCOME) %in% c("complete response", "partial response")                        ~ "Response",
      str_to_lower(TREATMENT_OUTCOME) %in% c("stable disease", "progressive disease",
                                             "clinical progressive disease",
                                             "radiographic progressive disease")                            ~ "Non-Response",
      TRUE ~ NA_character_
    )
  ) %>%
  
  # drops the row that doesn't match either group requirement
  filter(!is.na(RESPONSE_GROUP)) %>%
  
  # removes fully duplicated rows
  distinct()

cat("Treatment rows after filtering:", nrow(treatment), "\n")
cat("Response group breakdown:\n")
print(table(treatment$RESPONSE_GROUP))

# save the file 
write_tsv(treatment,
          file.path(OUTPUT_DIR, "treatment_response.txt"))

treatment_raw %>%
  select(THERAPEUTIC_AGENT, TREATMENT_OUTCOME) %>%
  distinct(TREATMENT_OUTCOME) %>%
  print(n = 50)

# load expression matrix

expr_raw <- read_tsv(EXPRESSION_FILE,
                     col_types = cols(Entrez_Gene_Id = col_character(),
                                      .default       = col_double()))

# map Entrez ID to Hugo gene symbol
entrez_to_symbol <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = na.omit(expr_raw$Entrez_Gene_Id),
  columns = "SYMBOL",
  keytype = "ENTREZID"
)
colnames(entrez_to_symbol) <- c("Entrez_Gene_Id", "Hugo_Symbol")
entrez_to_symbol <- entrez_to_symbol[!is.na(entrez_to_symbol$Hugo_Symbol), ]

# get the Entrez ID for your gene of interest
target_entrez <- entrez_to_symbol$Entrez_Gene_Id[entrez_to_symbol$Hugo_Symbol %in% GENES_OF_INTEREST]
cat("Entrez IDs for", GENES_OF_INTEREST, ":", target_entrez, "\n")

# subset expression matrix using Entrez ID directly — no Hugo_Symbol column needed
expr_long <- expr_raw[expr_raw$Entrez_Gene_Id %in% target_entrez, ] %>%
  pivot_longer(
    cols      = -Entrez_Gene_Id,
    names_to  = "barcode_full",
    values_to = "expression"
  ) %>%
  mutate(
    gene   = GENES_OF_INTEREST[match(Entrez_Gene_Id, target_entrez)],
    barcode = str_sub(barcode_full, 1, 12)
  )

cat("Expression rows for genes of interest:", nrow(expr_long), "\n")

write_tsv(expr_long,
          file.path(OUTPUT_DIR, "expression_response.txt"))

# merge both data

merged <- expr_long %>%
  inner_join(treatment, by = "barcode", relationship = "many-to-many") %>%
  filter(!is.na(expression))

cat("Merged rows:", nrow(merged), "\n")
cat("Unique patients:", n_distinct(merged$barcode), "\n")
cat("Response group breakdown in merged data:\n")
print(table(merged$RESPONSE_GROUP))

if (nrow(merged) == 0) {
  stop("Join returned 0 rows -- check that barcodes match between files.")
}

# deduplication and filter
# one expression value per patient per gene per response group

merged_filtered <- merged %>%
  distinct(gene, barcode, expression, RESPONSE_GROUP) %>%
  group_by(gene, RESPONSE_GROUP) %>%
  filter(n() >= MIN_SAMPLES_PER_GROUP) %>%
  ungroup()

cat("Rows after deduplication and filtering:", nrow(merged_filtered), "\n")

# save the merged data as txt and csv files

write_tsv(merged_filtered, file.path(OUTPUT_DIR, "merged_expression_response.txt"))
write_csv(merged_filtered, file.path(OUTPUT_DIR, "merged_expression_response.csv"))
cat("Data saved.\n")

# create the box plot

plot_gene <- function(gene_name, data) {
  
  df <- data %>%
    filter(gene == gene_name) %>%
    mutate(RESPONSE_GROUP = factor(RESPONSE_GROUP, levels = c("Response", "Non-Response")))
  
  n_labels <- df %>%
    group_by(RESPONSE_GROUP) %>%
    summarise(n = n(), .groups = "drop")
  
  p <- ggplot(df, aes(x = RESPONSE_GROUP, y = expression, fill = RESPONSE_GROUP)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.75, width = 0.5) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.5, color = "gray30") +
    geom_text(
      data = n_labels,
      aes(x = RESPONSE_GROUP, y = -Inf, label = paste0("n=", n)),
      inherit.aes = FALSE,
      vjust = -0.4, size = 3.5, color = "gray40"
    ) +
    scale_fill_manual(
      values = c("Response" = "#4E79A7", "Non-Response" = "#E15759"),
      guide  = "none"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.12, 0.05))) +
    labs(
      title    = paste(gene_name, "Expression by Drug Response"),
      subtitle = "Response = complete/partial  |  Non-Response = stable/progressive disease",
      x        = "Response Group",
      y        = paste(gene_name, "Expression (z-score)")
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title         = element_text(face = "bold"),
      panel.grid.major.x = element_blank()
    )
  
  out_path <- file.path(OUTPUT_DIR, paste0(gene_name, "_response_boxplot.png"))
  ggsave(out_path, p, width = 6, height = 6, dpi = 150)
  cat("Plot saved:", out_path, "\n")
}

walk(GENES_OF_INTEREST, plot_gene, data = merged_filtered)