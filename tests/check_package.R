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

# Parse R sources, check reference data, and scan the distributable tree.
require_packages("digest")
files <- list.files(root, recursive=TRUE, full.names=TRUE, all.files=TRUE, no..=TRUE)
files <- files[!file.info(files)$isdir]
relative <- substring(gsub("\\\\","/",files),nchar(root)+2L)
keep <- !grepl("^(reproduced_outputs|results|logs|external_data|[.]git)/|(^|/)node_modules/",relative)
files <- files[keep]; relative <- relative[keep]
r_files <- files[grepl("[.]R$",files,ignore.case=TRUE)]
for (path in r_files) parse(file=path,keep.source=FALSE)
stopifnot(!any(grepl("[.]py$|[.](png|jpe?g|tiff?|pdf|zip|rds|gz|h5)$",relative,ignore.case=TRUE)))
patterns <- c(
  windows_path="(?<![A-Za-z])[A-Za-z]:[/\\\\]",
  home_path="(?:^|[/\\\\])(?:[U]sers|[h]ome)[/\\\\][^/\\\\[:space:]]+[/\\\\]",
  identifier="\\b(?:patient[_]name|medical[_]record[_](?:id|number))\\b",
  credential="gh[pousr]_[A-Za-z0-9]{30,}|sk-[A-Za-z0-9]{32,}"
)
text_files <- files[grepl("[.](R|md|txt|csv|tsv|json|yml|yaml|cff)$",files,ignore.case=TRUE)]
findings <- character()
for (path in text_files) {
  text <- paste(readLines(path,warn=FALSE,encoding="UTF-8"),collapse="\n")
  for (p in names(patterns)) {
    if (grepl(patterns[[p]],text,perl=TRUE,ignore.case=TRUE)) findings <- c(findings,paste(basename(path),p))
  }
}
if (length(findings)) stop("Review sensitive text: ",paste(findings,collapse="; "))
p <- function(name) read_data(file.path(root,"data/processed",name))
r <- function(name) read_data(file.path(root,"derived_outputs",name))
stopifnot(nrow(p("target_source_membership_144.csv"))==144L,
          nrow(p("thyroid_cancer_genes_2139.csv"))==2139L,nrow(p("shared_targets_75.csv"))==75L,
          nrow(r("submitted_75_PPI_nodes.csv"))==75L,nrow(r("submitted_75_PPI_edges.csv"))==508L)
stopifnot(setequal(intersect(p("target_source_membership_144.csv")$Gene,p("thyroid_cancer_genes_2139.csv")$Gene),
                    p("shared_targets_75.csv")$Gene))
go <- read.delim(file.path(root,"derived_outputs/submitted_75_GO_results.txt"))
kegg <- read.delim(file.path(root,"derived_outputs/submitted_75_KEGG_results.txt"))
stopifnot(nrow(go)==1899L,nrow(kegg)==169L,sum(kegg$p.adjust<0.05)==166L)
immune <- r("immune_group_comparisons_with_BH_FDR.csv")
correlation <- r("gene_immune_correlations_full_and_tumor_only_with_BH_FDR.csv")
stopifnot(sum(immune$Significant_BH_FDR_lt_0.05=="Yes")==1L,
          sum(correlation$All_samples_significant_FDR_lt_0.05=="Yes")==7L,
          sum(correlation$PTC_only_significant_FDR_lt_0.05=="Yes")==0L,
          sum(is.na(correlation$PTC_only_P))==3L)
clinical <- p("ihc_patient_level_clinical_deidentified.csv")
hscore <- p("ihc_hscore_deidentified.csv")
hpa <- p("hpa_panel_provenance.csv")
stopifnot(nrow(clinical)==15L,length(unique(clinical$Pair_ID))==15L,nrow(hscore)==45L,
          setequal(clinical$Pair_ID,hscore$Pair_ID),nrow(hpa)==3L,
          all(hpa$Normal_antibody_ID==hpa$PTC_antibody_ID))
performance <- r("ml_nested_resampling_performance.csv")
stopifnot(all(performance$Outer_split_count==50L),
          abs(performance$Mean_outer_AUC[1]-0.9699506172839507)<1e-12)
# The optional manifest identifies exact source/reference files for a release.
checksum <- file.path(root,"docs/checksums_sha256.txt")
hash_paths <- files[files != checksum]
hash_names <- relative[files != checksum]
ord <- order(hash_names,method="radix"); hash_paths <- hash_paths[ord]; hash_names <- hash_names[ord]
hashes <- vapply(hash_paths,function(p) digest::digest(file=p,algo="sha256",serialize=FALSE),character(1))
lines <- paste(hashes,hash_names,sep="  ")
if ("--write-checksums" %in% args) writeLines(lines,checksum,useBytes=TRUE)
if (file.exists(checksum)) {
  saved <- readLines(checksum,warn=FALSE,encoding="UTF-8")
  if (!identical(saved,lines)) stop("Checksum manifest differs; inspect edits before updating it.")
}
cat(length(r_files),"R files parsed;",length(text_files),"text files scanned; package checks passed.\n")

