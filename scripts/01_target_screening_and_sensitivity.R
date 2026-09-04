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

# 1. Load source membership and disease-gene lists.
out <- safe_output(if (length(args)) args[1] else file.path(root, "reproduced_outputs/targets"), root)
processed <- file.path(root, "data/processed")
source_table <- read_data(file.path(processed, "target_source_membership_144.csv"))
disease <- read_data(file.path(processed, "thyroid_cancer_genes_2139.csv"))
genecards <- read_data(file.path(processed, "genecards_top3000.csv"))
shared_reference <- read_data(file.path(processed, "shared_targets_75.csv"))
clean_symbols <- function(x) sort(unique(trimws(as.character(x[!is.na(x) & nzchar(trimws(x))]))), method = "radix")
yes <- function(x) tolower(trimws(as.character(x))) %in% c("yes", "true", "1")
candidates <- c("BAX", "BCL2", "CDKN1A", "COL1A1", "CXCL8", "FN1", "TGFB1")
core <- c("BAX", "BCL2", "FN1")
exposure <- clean_symbols(source_table$Gene)
disease_genes <- clean_symbols(disease$Gene)
shared <- sort(intersect(exposure, disease_genes), method = "radix")
stopifnot(length(exposure) == 144L, length(disease_genes) == 2139L,
          length(shared) == 75L, setequal(shared, shared_reference$Gene),
          all(candidates %in% shared), all(core %in% shared))
ctd <- clean_symbols(source_table$Gene[yes(source_table$CTD_literature_curated)])
ctd_shared <- intersect(ctd, disease_genes)

# 2. Repeat intersections over the specified GeneCards thresholds.
omim <- clean_symbols(disease$Gene[yes(disease$In_OMIM)])
ranks <- suppressWarnings(as.numeric(genecards$Rank))
thresholds <- do.call(rbind, lapply(c(500L, 1000L, 1500L, 2000L, 3000L), function(cutoff) {
  disease_set <- union(clean_symbols(genecards$Gene[!is.na(ranks) & ranks <= cutoff]), omim)
  shared_set <- intersect(exposure, disease_set)
  kept <- sort(intersect(shared_set, candidates), method = "radix")
  data.frame(GeneCards_cutoff = cutoff, Disease_gene_count = length(disease_set),
             Shared_gene_count = length(shared_set), PTC_DEG_candidate_count = length(kept),
             PTC_DEG_candidates = paste(kept, collapse = "; "),
             Core_genes_retained = all(core %in% shared_set))
}))
write_data(data.frame(Gene = shared), file.path(out, "shared_targets_75_reproduced.csv"))
write_data(thresholds, file.path(out, "genecards_threshold_sensitivity_reproduced.csv"))
require_packages("jsonlite")
summary <- list(exposure_target_union = length(exposure), thyroid_cancer_gene_union = length(disease_genes),
                submitted_shared_targets = length(shared), ptc_deg_candidates = candidates, final_core_genes = core,
                ctd_curated_targets = length(ctd), ctd_only_shared_targets = length(ctd_shared),
                ctd_only_candidates = sort(intersect(ctd_shared, candidates), method = "radix"),
                genecards_thresholds = thresholds)
jsonlite::write_json(summary, file.path(out, "target_screening_summary.json"),
                     auto_unbox = TRUE, pretty = TRUE, digits = NA)
cat("Target and threshold intersections completed.\n")

