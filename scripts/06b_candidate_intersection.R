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

# Intersect the shared targets with the thresholded gene-level DEG summary.
deg_path <- if (length(args)) args[1] else file.path(root, "derived_outputs/GSE33630_DEG_gene_summary.csv")
out <- safe_output(if (length(args) >= 2L) args[2] else file.path(root, "reproduced_outputs/candidates"), root)
deg <- read_data(deg_path)
shared <- read_data(file.path(root, "data/processed/shared_targets_75.csv"))
stopifnot(all(c("Gene.symbol", "logFC", "adj.P.Val") %in% names(deg)),
          !anyDuplicated(deg$Gene.symbol), all(deg$adj.P.Val < 0.05 & abs(deg$logFC) > 1))
selected <- deg[deg$Gene.symbol %in% shared$Gene, , drop = FALSE]
selected <- selected[order(selected$Gene.symbol, method = "radix"), , drop = FALSE]
selected$In_GEO2R_DEGs <- TRUE
selected <- selected[c(names(deg)[1], "In_GEO2R_DEGs", names(deg)[-1])]
write_data(selected, file.path(out, "seven_gene_candidates.csv"))
writeLines(selected$Gene.symbol, file.path(out, "candidate_7genes.txt"))
cat(nrow(selected), "candidate genes derived from the shared-target/DEG intersection.\n")

