# Run with R 4.3.3; paths are relative to this package.
if (.Platform$OS.type == "windows") {
  invisible(suppressWarnings(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8")))
}
if (as.character(getRversion()) != "4.3.3") stop("Use R 4.3.3")
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
source_files <- vapply(sys.frames(), function(x) if (is.null(x$ofile)) "" else x$ofile, character(1))
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else tail(source_files[nzchar(source_files)], 1)
if (!length(this_file)) stop("Run this file with Rscript or source().")
root <- normalizePath(dirname(this_file), winslash = "/", mustWork = TRUE)
source(file.path(root, "scripts/lib/workflow_helpers.R"))
args <- commandArgs(trailingOnly = TRUE)

# Run the compact-data analyses. Use --output to preserve separate runs.
require_packages("jsonlite")
output <- file.path(root, "reproduced_outputs/complete_run")
include_go <- FALSE
geo2r_table <- NULL
skip_s13 <- FALSE
i <- 1L
while (i <= length(args)) {
  flag <- args[i]
  if (flag %in% c("--include-go", "--skip-s13")) {
    if (flag == "--include-go") include_go <- TRUE else skip_s13 <- TRUE
    i <- i + 1L
  } else if (flag %in% c("--output", "--geo2r-table")) {
    if (i == length(args)) stop("Missing value for ", flag)
    value <- args[i + 1L]
    if (flag == "--output") output <- value
    if (flag == "--geo2r-table") geo2r_table <- normalizePath(value, winslash = "/", mustWork = TRUE)
    i <- i + 2L
  } else stop("Unknown option: ", flag)
}
out <- safe_output(output, root, fresh = TRUE)
logs <- file.path(out, "logs"); dir.create(logs)
rscript <- file.path(R.home("bin"), "Rscript")
jobs <- list()
run_job <- function(name, executable, arguments) {
  old_locale <- Sys.getenv(c("LC_ALL", "LANG"), unset = NA_character_)
  on.exit({
    for (key in names(old_locale)) {
      if (is.na(old_locale[[key]])) Sys.unsetenv(key) else do.call(Sys.setenv, setNames(list(old_locale[[key]]), key))
    }
  }, add = TRUE)
  if (.Platform$OS.type == "windows") Sys.setenv(LC_ALL = "Chinese (Simplified)_China.utf8", LANG = "Chinese (Simplified)_China.utf8")
  cat("Running ", name, "...\n", sep = ""); flush.console()
  start <- Sys.time()
  status <- system2(executable, args = vapply(arguments, shQuote, character(1)),
                    stdout = file.path(logs, paste0(name, ".log")), stderr = file.path(logs, paste0(name, ".log")))
  jobs[[length(jobs) + 1L]] <<- data.frame(Analysis = name, Exit_status = status,
                                         Seconds = as.numeric(difftime(Sys.time(), start, units = "secs")))
  write_data(do.call(rbind, jobs), file.path(out, "run_summary.csv"))
  if (status != 0L) stop(name, " failed; inspect ", name, ".log in the chosen output/logs directory.")
}
run_r <- function(name, script, arguments = character()) run_job(name, rscript, c(file.path(root, script), arguments))
run_r("targets", "scripts/01_target_screening_and_sensitivity.R", file.path(out, "targets"))
deg_summary <- file.path(root, "derived_outputs/GSE33630_DEG_gene_summary.csv")
if (!is.null(geo2r_table)) {
  dir.create(file.path(out, "deg"))
  run_r("deg", "scripts/06_GEO2R_DEG.R", c(file.path(out, "deg"), geo2r_table))
  deg_summary <- file.path(out, "deg/GEO2R_DEG_gene_summary_adjP0.05_logFC1_official.csv")
}
run_r("candidates", "scripts/06b_candidate_intersection.R", c(deg_summary, file.path(out, "candidates")))
run_r("ppi", "scripts/03b_STRING_PPI_concentric_figure.R", c(file.path(root, "derived_outputs"), file.path(out, "ppi")))
inputs <- list(
  lasso = c("07_LASSO_original.R", "discovery_7gene_expression.csv"),
  svm = c("08_SVM_RFE_original.R", "discovery_7gene_expression.csv"),
  roc_discovery = c("09_ROC_discovery.R", "discovery_7gene_expression.csv"),
  roc_external = c("10_ROC_external.R", "external_core_expression.csv")
)
for (name in names(inputs)) {
  work <- file.path(out, name); dir.create(work)
  file.copy(file.path(root, "data/processed", inputs[[name]][2]), file.path(work, "geneexp.csv"))
  if (startsWith(name, "roc")) file.copy(file.path(root, "data/processed/core_genes_no_header.csv"),
                                       file.path(work, "IntersectionGenes.csv"))
  run_r(name, file.path("scripts", inputs[[name]][1]), work)
}
if (!skip_s13) {
  run_r("nested_ml", "scripts/08b_nested_resampling.R", file.path(out, "nested_ml"))
}
run_r("immune_groups", "scripts/12_immune_group_comparison.R",
      c(file.path(out, "immune_groups"), file.path(out, "immune_figures")))
run_r("gene_immune", "scripts/13_gene_immune_correlations.R",
      c(file.path(out, "gene_immune"), file.path(out, "immune_figures")))
run_r("ihc_statistics", "scripts/17b_IHC_paired_statistics.R", file.path(out, "ihc_statistics"))
if (include_go) run_r("go", "scripts/04_GO_enrichment.R",
                     c(file.path(out, "go"), file.path(root, "data/processed/shared_targets_75.csv")))
run_r("verification", "tests/verify_reproduced_results.R", c(out, if (skip_s13) "--skip-s13"))
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
cat("Completed. All generated files are in the chosen output directory.\n")
