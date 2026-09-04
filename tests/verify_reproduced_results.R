# Run with R 4.3.3; paths are relative to this package.
if (.Platform$OS.type == "windows") {
  invisible(suppressWarnings(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8")))
}
if (as.character(getRversion()) != "4.3.3") stop("Use R 4.3.3")
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
source_files <- vapply(sys.frames(), function(x) if (is.null(x$ofile)) "" else x$ofile, character(1))
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else tail(source_files[nzchar(source_files)], 1)
if (!length(this_file)) stop("Run this file with Rscript or source().")
root <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts/lib/workflow_helpers.R"))
args <- commandArgs(trailingOnly = TRUE)

# Check numerical outputs against the publication reference tables.
require_packages("jsonlite")
if (!length(args)) stop("Supply the run output directory.")
out <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
skip_s13 <- "--skip-s13" %in% args
ref <- file.path(root, "derived_outputs")
checked <- character()
canonical <- function(x) {
  x <- as.character(x); x[is.na(x)] <- ""
  bool <- tolower(x) %in% c("true", "false")
  x[bool] <- toupper(x[bool])
  x
}
compare <- function(actual, expected, keys, omit = character()) {
  stopifnot(nrow(actual) == nrow(expected), setequal(names(actual), names(expected)))
  a <- actual[do.call(order, c(actual[keys], list(method = "radix"))), , drop = FALSE]
  b <- expected[do.call(order, c(expected[keys], list(method = "radix"))), , drop = FALSE]
  for (col in setdiff(names(b), omit)) {
    if (is.numeric(b[[col]])) {
      if (!identical(is.na(a[[col]]), is.na(b[[col]]))) stop(col, ": NA pattern differs")
      ok <- !is.na(b[[col]])
      if (any(abs(a[[col]][ok]-b[[col]][ok]) > 1e-12 + 1e-10*abs(b[[col]][ok]))) stop(col, ": values differ")
    } else if (!identical(canonical(a[[col]]), canonical(b[[col]]))) stop(col, ": labels differ")
  }
  invisible(TRUE)
}
check_table <- function(folder, name, keys, omit = character()) {
  compare(read_data(file.path(out, folder, name)), read_data(file.path(ref, name)), keys, omit)
  checked <<- c(checked, name)
}
# R target/sensitivity calculations: exact gene sets and integer counts.
shared <- read_data(file.path(out, "targets/shared_targets_75_reproduced.csv"))
stopifnot(setequal(shared$Gene, read_data(file.path(root, "data/processed/shared_targets_75.csv"))$Gene))
threshold <- read_data(file.path(out, "targets/genecards_threshold_sensitivity_reproduced.csv"))
stopifnot(identical(threshold$Shared_gene_count, c(39L,53L,68L,75L,84L)))
checked <- c(checked, "shared_targets_75_reproduced.csv", "genecards_threshold_sensitivity_reproduced.csv")
if (dir.exists(file.path(out, "deg"))) {
  a <- read_data(file.path(out, "deg/GEO2R_DEG_gene_summary_adjP0.05_logFC1_official.csv"))
  compare(a, read_data(file.path(ref, "GSE33630_DEG_gene_summary.csv")), "Gene.symbol")
  probes <- read_data(file.path(out, "deg/GEO2R_DEG_probes_adjP0.05_logFC1_official.csv"))
  stopifnot(nrow(probes)==1901, sum(probes$DEG_status=="Up")==1060, sum(probes$DEG_status=="Down")==841)
  checked <- c(checked, "GSE33630_DEG_gene_summary.csv")
}
check_table("candidates", "seven_gene_candidates.csv", "Gene.symbol")
if (!skip_s13) {
  check_table("nested_ml", "ml_nested_resampling_split_results.csv", "Split")
  check_table("nested_ml", "ml_nested_resampling_selection_stability.csv", "Gene")
  check_table("nested_ml", "ml_nested_resampling_performance.csv", "Method", "Evaluation_note")
}
check_table("immune_groups", "immune_group_comparisons_with_BH_FDR.csv", "Immune_cell")
name <- "gene_immune_correlations_full_and_tumor_only_with_BH_FDR.csv"
actual <- read_data(file.path(out, "gene_immune", name))
expected <- read_data(file.path(ref, name))
expected$PTC_only_significant_FDR_lt_0.05[is.na(expected$PTC_only_P)] <- "Not_estimable"
compare(actual, expected, c("Gene", "Immune_cell"))
counts <- read_data(file.path(out, "gene_immune/gene_immune_test_counts.csv"))
stopifnot(sum(is.na(actual$PTC_only_P))==3L, identical(counts$BH_denominator, c(66L,63L)))
checked <- c(checked, name)
check_table("ihc_statistics", "ihc_paired_statistics.csv", "Gene")
# CI values fluctuate because the 1,000 bootstrap resamples are unseeded.
roc_reference <- read_data(file.path(ref, "single_gene_ROC_results.csv"))
roc_delta <- list()
for (dataset in c("GSE33630", "GSE60542")) {
  folder <- if (dataset=="GSE33630") "roc_discovery" else "roc_external"
  a <- read_data(file.path(out, folder, "gene_analysis_results.csv"))
  b <- roc_reference[roc_reference$Dataset==dataset, names(a), drop=FALSE]
  compare(a,b,"Gene",c("AUC_Lower_CI","AUC_Upper_CI"))
  b <- b[match(a$Gene,b$Gene),,drop=FALSE]
  stopifnot(all(is.finite(a$AUC_Lower_CI)),all(is.finite(a$AUC_Upper_CI)),
            all(a$AUC_Lower_CI>=0 & a$AUC_Upper_CI<=1 &
                a$AUC_Lower_CI<=a$AUC & a$AUC<=a$AUC_Upper_CI))
  delta <- pmax(abs(a$AUC_Lower_CI-b$AUC_Lower_CI),abs(a$AUC_Upper_CI-b$AUC_Upper_CI))
  if (any(delta>0.02)) warning("Inspect bootstrap CI endpoint differences exceeding 0.02: ",dataset)
  roc_delta[[dataset]] <- data.frame(Dataset=dataset,Gene=a$Gene,max_CI_endpoint_difference=delta)
}
read_genes <- function(path) gsub('"','',trimws(readLines(path,warn=FALSE)))
stopifnot(setequal(read_genes(file.path(out,"lasso/LASSO_genes_modified.txt")),c("BAX","BCL2","COL1A1","FN1")),
          setequal(read_genes(file.path(out,"svm/SVM-RFE.gene.txt")),c("BAX","FN1","BCL2","CXCL8")))
panels <- list.files(file.path(out,"immune_figures"),pattern="[.]png$",full.names=TRUE)
stopifnot(length(panels)==5L,all(file.info(panels)$size>10000),
          length(list.files(file.path(out,"ppi"),pattern="[.]png$"))>0)
if (file.exists(file.path(out,"go/GO_results.txt"))) {
  a <- read.delim(file.path(out,"go/GO_results.txt"),check.names=FALSE)
  b <- read.delim(file.path(ref,"submitted_75_GO_results.txt"),check.names=FALSE)
  stopifnot(nrow(a)==1899L,setequal(a$ID,b$ID))
  cols <- c("ID","pvalue","p.adjust","qvalue","Count")
  compare(a[cols],b[cols],"ID"); checked <- c(checked,"GO_results.txt")
}
report <- list(validated_tables=checked, S13_checked=!skip_s13, roc_AUC_means_P_match=TRUE,
               roc_unseeded_CI_differences=do.call(rbind,roc_delta),
               PTC_estimable_tests=63L,PTC_undefined_tests=3L)
jsonlite::write_json(report,file.path(out,"verification_summary.json"),auto_unbox=TRUE,pretty=TRUE,digits=NA)
cat("Recomputed statistics agree with reference values. ROC CI differences recorded.\n")

