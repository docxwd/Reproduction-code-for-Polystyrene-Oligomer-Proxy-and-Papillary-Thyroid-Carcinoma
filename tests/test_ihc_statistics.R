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

expected <- data.frame(
  marker = c("BAX", "BCL2", "FN1"),
  mean_difference = c(20.351333333333333, -23.43, 41.964),
  t_p = c(0.001220, 0.001297, 0.000434),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(expected))) {
  marker <- expected$marker[[i]]
  path <- file.path(root, "data", "processed", paste0(marker, "_Hscore_deidentified.csv"))
  dat <- read.csv(path, check.names = FALSE)
  stopifnot(nrow(dat) == 15L)
  difference <- dat$tumor_hscore - dat$adjacent_normal_hscore
  stopifnot(max(abs(difference - dat$tumor_minus_adjacent_difference)) <= 0.011)
  stopifnot(abs(mean(difference) - expected$mean_difference[[i]]) < 1e-12)
  result <- t.test(dat$tumor_hscore, dat$adjacent_normal_hscore, paired = TRUE)
  stopifnot(abs(result$p.value - expected$t_p[[i]]) < 1e-6)
  wilcox <- suppressWarnings(wilcox.test(
    dat$tumor_hscore, dat$adjacent_normal_hscore,
    paired = TRUE, exact = TRUE, conf.int = TRUE
  ))
  stopifnot(wilcox$p.value < 0.001)
}

cat("IHC BASE-R STATISTICS TEST: PASS\n")
