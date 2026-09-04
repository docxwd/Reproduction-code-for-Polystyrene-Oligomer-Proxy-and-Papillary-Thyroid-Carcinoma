# Shared input, output and validation utilities. No analysis is run on source().
read_data <- function(path) read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                                    fileEncoding = "UTF-8-BOM", na.strings = c("", "NA"))
write_data <- function(x, path) write.csv(x, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))
}
safe_output <- function(path, root, fresh = FALSE) {
  # Create the parent first, then resolve it before writing any analysis files.
  path <- path.expand(path)
  if (!dir.exists(dirname(path))) dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  path <- file.path(normalizePath(dirname(path), winslash = "/", mustWork = TRUE), basename(path))
  resolved <- if (dir.exists(path)) normalizePath(path, winslash = "/", mustWork = TRUE) else path
  key <- tolower(resolved)
  protected <- tolower(c(root, file.path(root, c("data", "derived_outputs", "scripts", "tests",
                                                "figures", "docs", "config", "environment"))))
  if (key == protected[1] || any(vapply(protected[-1], function(p)
      key == p || startsWith(key, paste0(p, "/")), logical(1)))) stop("Output overlaps package inputs")
  if (fresh && file.exists(resolved)) stop("Output already exists; choose a new directory.")
  dir.create(resolved, recursive = TRUE, showWarnings = FALSE)
  normalizePath(resolved, winslash = "/", mustWork = TRUE)
}
bh_valid <- function(p) {
  result <- rep(NA_real_, length(p))
  valid <- is.finite(p)
  result[valid] <- p.adjust(p[valid], method = "BH")
  result
}
significance <- function(q) ifelse(is.na(q), "Not_estimable", ifelse(q < 0.05, "Yes", "No"))
load_fractions <- function(path) {
  dat <- read_data(path)
  stopifnot(all(c("id", "Group") %in% names(dat)))
  cells <- setdiff(names(dat), c("id", "Group"))
  stopifnot(length(cells) == 22L, nrow(dat) == 42L, !anyDuplicated(dat$id),
            sum(dat$Group == "Control") == 15L, sum(dat$Group == "Tumor") == 27L,
            all(is.finite(as.matrix(dat[cells]))), all(as.matrix(dat[cells]) >= 0),
            max(abs(rowSums(dat[cells]) - 1)) < 1e-12)
  dat
}

