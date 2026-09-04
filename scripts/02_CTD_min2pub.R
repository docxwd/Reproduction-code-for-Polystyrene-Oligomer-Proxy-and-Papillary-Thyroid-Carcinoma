# CTD chemical-gene interaction screening
# Rule:
# 1) same chemical search object: Chemical ID == D011137 (Polystyrenes)
# 2) same gene: merged by Gene Symbol and Gene ID
# 3) retain genes supported by at least 2 distinct PMID values in References
# 4) mRNA/protein/other interaction records are all acceptable evidence, but
#    repeated rows supported by the same PMID are counted as one publication.

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

cli_args <- commandArgs(trailingOnly = TRUE)
work_dir <- if (length(cli_args) >= 1L) cli_args[[1L]] else getwd()
if (!dir.exists(work_dir)) stop("Analysis directory not found: ", work_dir)
setwd(work_dir)
input_file <- if (length(cli_args) >= 2L) cli_args[[2L]] else "CTD_D011137_ixns_export.xlsx"
chemical_id_keep <- "D011137"
min_publications <- 2

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file,
       "\nPlease run this script in the CTD folder.")
}

raw_ctd <- readxl::read_excel(input_file)

required_cols <- c(
  "Chemical Name", "Chemical ID", "Gene Symbol", "Gene ID",
  "Interaction", "Interaction Actions", "References", "Organisms"
)
missing_cols <- setdiff(required_cols, names(raw_ctd))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

ctd <- raw_ctd %>%
  mutate(across(everything(), as.character)) %>%
  filter(`Chemical ID` == chemical_id_keep)

expanded_pmids <- ctd %>%
  mutate(References = if_else(is.na(References), "", References)) %>%
  separate_rows(References, sep = "\\|") %>%
  mutate(PMID = str_squish(References)) %>%
  filter(PMID != "", str_detect(PMID, "^[0-9]+$"))

gene_summary <- expanded_pmids %>%
  group_by(`Chemical Name`, `Chemical ID`, `Gene Symbol`, `Gene ID`) %>%
  summarise(
    PMID_count = n_distinct(PMID),
    PMIDs = paste(sort(unique(PMID)), collapse = "|"),
    CTD_record_count = n_distinct(Interaction),
    Interaction_actions = paste(sort(unique(`Interaction Actions`)), collapse = "; "),
    Organisms = paste(sort(unique(unlist(str_split(Organisms, "\\|")))), collapse = "|"),
    .groups = "drop"
  ) %>%
  arrange(desc(PMID_count), `Gene Symbol`)

retained_genes <- gene_summary %>%
  filter(PMID_count >= min_publications)

excluded_genes <- gene_summary %>%
  filter(PMID_count < min_publications)

filtered_records <- ctd %>%
  semi_join(
    retained_genes,
    by = c("Chemical ID", "Gene Symbol", "Gene ID")
  ) %>%
  arrange(`Gene Symbol`)

retained_gene_list <- retained_genes %>%
  arrange(`Gene Symbol`) %>%
  pull(`Gene Symbol`)

write.csv(gene_summary, "CTD_gene_summary_by_distinct_PMID.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(retained_genes, "CTD_genes_retained_min2PMID_summary.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(excluded_genes, "CTD_genes_excluded_lt2PMID_summary.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(filtered_records, "CTD_records_retained_min2PMID.csv", row.names = FALSE, fileEncoding = "UTF-8")
writeLines(retained_gene_list, "CTD_genes_retained_min2PMID.txt", useBytes = TRUE)

if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(
    list(
      gene_summary = gene_summary,
      retained_genes = retained_genes,
      excluded_genes = excluded_genes,
      retained_records = filtered_records
    ),
    "CTD_filtered_min2PMID.xlsx"
  )
} else {
  message("Package 'writexl' is not installed; xlsx output was skipped. CSV/TXT outputs were created.")
}

message("CTD screening completed.")
message("Input CTD records: ", nrow(ctd))
message("Unique genes before filtering: ", n_distinct(ctd$`Gene Symbol`))
message("Retained genes with >= ", min_publications, " distinct PMID(s): ", nrow(retained_genes))
message("Retained CTD records: ", nrow(filtered_records))
