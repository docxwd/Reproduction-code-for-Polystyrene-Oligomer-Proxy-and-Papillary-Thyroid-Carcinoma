# FN1 entry point for the shared H-score analysis workflow.

if (.Platform$OS.type == "windows") {
  invisible(
    suppressWarnings(try(Sys.setlocale("LC_CTYPE", ".UTF-8"), silent = TRUE))
  )
}

get_entry_script_dir <- function() {
  command_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", command_args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(
      sub("^--file=", "", file_arg[[1L]]),
      winslash = "/",
      mustWork = TRUE
    )))
  }
  sourced_file <- tryCatch(sys.frame(1L)$ofile, error = function(e) NULL)
  if (!is.null(sourced_file)) {
    return(dirname(normalizePath(
      sourced_file,
      winslash = "/",
      mustWork = TRUE
    )))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

options(
  hscore.marker_name = "FN1",
  hscore.script_index = "03",
  hscore.primary_test = "paired_t",
  hscore.input_filename = file.path("data", "processed", "FN1_Hscore_deidentified.csv"),
  hscore.annotation_mode = "both",
  hscore.export_bar_only = TRUE
)
source(
  file.path(get_entry_script_dir(), "17_IHC_Hscore_core.R"),
  encoding = "UTF-8",
  chdir = FALSE
)
