#!/usr/bin/env Rscript

if (.Platform$OS.type == "windows") {
  invisible(suppressWarnings(try(Sys.setlocale("LC_CTYPE", ".UTF-8"), silent = TRUE)))
}
trailing <- commandArgs(trailingOnly = TRUE)
if (length(trailing) >= 1L) {
  root <- normalizePath(trailing[[1L]], winslash = "/", mustWork = TRUE)
} else {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/", mustWork = TRUE)
  root <- dirname(dirname(script_path))
}
files <- list.files(file.path(root, "scripts"), pattern = "[.]R$", full.names = TRUE)
for (path in files) {
  tryCatch(
    parse(file = path, keep.source = FALSE),
    error = function(e) stop("Parse failure in ", basename(path), ": ", conditionMessage(e))
  )
}
cat("R PARSE TEST: PASS (", length(files), " scripts)\n", sep = "")
