# =============================================================================
# IHC H-score: shared paired statistics and publication-style figure workflow
#
# Purpose:
#   1) Validate paired H-score data without altering the source CSV.
#   2) Run a pre-specified primary paired test plus a sensitivity test.
#   3) Export a paired-data figure, a reference-style paired bar figure,
#      a QQ QC figure, a machine-readable statistics table, and session information.
#
# Required package: ggplot2
# Tested with R 4.3.3.
# =============================================================================

# ---- 0. User configuration ---------------------------------------------------
marker_name <- getOption("hscore.marker_name", "BAX")
script_index <- getOption("hscore.script_index", "00")
input_filename <- getOption("hscore.input_filename", NULL)

# Choose exactly one primary test through the marker entry script. Both paired
# tests are saved; annotation_mode controls whether one or both P values appear.
primary_test <- getOption("hscore.primary_test", "paired_t")
annotation_mode <- getOption("hscore.annotation_mode", "primary")
export_bar_only <- isTRUE(getOption("hscore.export_bar_only", FALSE))
alternative <- "two.sided"
confidence_level <- 0.95

# NPG-inspired colors commonly used in biomedical figures.
npg_colors <- c("#4DBBD5", "#E64B35")

# X-axis labels for the SCI figure.
plot_group_labels <- c(adjacent_normal = "Con", tumor = "Tumor")

# Publication export settings.
figure_width_mm <- 100
figure_height_mm <- 92
png_dpi <- 300
tiff_dpi <- 600
analysis_date <- format(Sys.Date(), "%Y-%m-%d")

# H-score is conventionally bounded between 0 and 300.
hscore_lower_bound <- 0
hscore_upper_bound <- 300
difference_check_tolerance <- 0.011  # allows a two-decimal supplied difference

# Optional explicit module directory. Keep NULL for automatic resolution.
analysis_dir_override <- NULL

# ---- 1. Locale, package, and path setup -------------------------------------
# Some Windows sessions inherit the Linux locale name C.UTF-8, which prevents
# R for Windows from enumerating Chinese paths. This switches to Windows UTF-8.
if (.Platform$OS.type == "windows") {
  invisible(
    suppressWarnings(try(Sys.setlocale("LC_CTYPE", ".UTF-8"), silent = TRUE))
  )
}
options(stringsAsFactors = FALSE)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop(
    "Package 'ggplot2' is required; see environment/requirements-r.txt.",
    call. = FALSE
  )
}

get_script_path <- function() {
  command_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", command_args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(normalizePath(
      sub("^--file=", "", file_arg[[1L]]),
      winslash = "/",
      mustWork = TRUE
    ))
  }

  sourced_file <- tryCatch(sys.frame(1L)$ofile, error = function(e) NULL)
  if (!is.null(sourced_file)) {
    return(normalizePath(sourced_file, winslash = "/", mustWork = TRUE))
  }
  NA_character_
}

script_path <- get_script_path()

if (!is.null(analysis_dir_override)) {
  analysis_dir <- normalizePath(
    analysis_dir_override,
    winslash = "/",
    mustWork = TRUE
  )
} else if (!is.na(script_path)) {
  analysis_dir <- dirname(dirname(script_path))
} else {
  # Fallback for interactive execution: run from either the module root or
  # its scripts directory.
  working_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  analysis_dir <- if (basename(working_dir) == "scripts") {
    dirname(working_dir)
  } else {
    working_dir
  }
}

if (!is.null(input_filename)) {
  input_candidate <- file.path(analysis_dir, input_filename)
  if (!file.exists(input_candidate)) {
    stop(
      "Configured source CSV does not exist: ", input_candidate,
      call. = FALSE
    )
  }
  input_file <- normalizePath(
    input_candidate,
    winslash = "/",
    mustWork = TRUE
  )
} else {
  source_pattern <- paste0("^", marker_name, "_Hscore_.*[.]csv$")
  input_candidates <- list.files(
    analysis_dir,
    pattern = source_pattern,
    full.names = TRUE,
    recursive = FALSE
  )

  if (length(input_candidates) != 1L) {
    stop(
      "Expected exactly one source CSV in: ", analysis_dir,
      "\nPattern: ", source_pattern,
      "\nFound: ", length(input_candidates),
      "\nSet option 'hscore.input_filename' to select the intended file explicitly.",
      call. = FALSE
    )
  }

  input_file <- normalizePath(
    input_candidates[[1L]],
    winslash = "/",
    mustWork = TRUE
  )
}

output_root <- getOption("hscore.output_dir", file.path(analysis_dir, "reproduced_outputs", "ihc_figures"))
results_dir <- file.path(output_root, "results")
figures_dir <- file.path(output_root, "figures")
qc_dir <- file.path(figures_dir, "qc")
logs_dir <- file.path(output_root, "logs")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 2. Import as text first, then validate and convert ----------------------
# Reading every source field as character preserves identifiers exactly as they
# appear in the CSV (including any future leading zeroes).
raw_data <- utils::read.csv(
  input_file,
  fileEncoding = "UTF-8-BOM",
  check.names = FALSE,
  colClasses = "character",
  na.strings = c("", "NA", "N/A")
)

# Spreadsheet exports can contain completely empty trailing columns. They are
# ignored explicitly and counted in the QC output without altering the source.
raw_column_count <- ncol(raw_data)
all_empty_columns <- vapply(
  raw_data,
  function(column) {
    normalized <- ifelse(is.na(column), "", trimws(column))
    all(!nzchar(normalized))
  },
  logical(1)
)
ignored_empty_column_count <- sum(all_empty_columns)
if (ignored_empty_column_count > 0L) {
  raw_data <- raw_data[, !all_empty_columns, drop = FALSE]
}

expected_columns <- c(
  "pair_id",
  "adjacent_normal_hscore",
  "tumor_hscore"
)
missing_columns <- setdiff(expected_columns, names(raw_data))
if (length(missing_columns) > 0L) {
  stop(
    "Missing required column(s): ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

id_column <- expected_columns[[1L]]
adjacent_column <- expected_columns[[2L]]
tumor_column <- expected_columns[[3L]]
difference_column <- "tumor_minus_adjacent_difference"

to_numeric_strict <- function(x, column_name) {
  source_text <- trimws(x)
  converted <- suppressWarnings(as.numeric(source_text))
  bad <- is.na(converted) & !is.na(source_text) & nzchar(source_text)
  if (any(bad)) {
    stop(
      "Non-numeric value(s) in column '", column_name,
      "' at row(s): ", paste(which(bad) + 1L, collapse = ", "),
      call. = FALSE
    )
  }
  converted
}

sample_id <- trimws(raw_data[[id_column]])
adjacent <- to_numeric_strict(raw_data[[adjacent_column]], adjacent_column)
tumor <- to_numeric_strict(raw_data[[tumor_column]], tumor_column)

if (nrow(raw_data) < 3L) {
  stop("At least three complete biological pairs are required.", call. = FALSE)
}
if (anyNA(sample_id) || any(!nzchar(sample_id))) {
  stop("Sample ID contains missing or blank values.", call. = FALSE)
}
if (anyDuplicated(sample_id)) {
  stop(
    "Duplicated sample ID(s): ",
    paste(unique(sample_id[duplicated(sample_id)]), collapse = ", "),
    call. = FALSE
  )
}
if (anyNA(adjacent) || anyNA(tumor)) {
  stop("H-score columns contain missing values.", call. = FALSE)
}
if (any(c(adjacent, tumor) < hscore_lower_bound) ||
    any(c(adjacent, tumor) > hscore_upper_bound)) {
  stop(
    "H-score outside the expected range ", hscore_lower_bound,
    " to ", hscore_upper_bound, ".",
    call. = FALSE
  )
}

paired_difference <- round(tumor - adjacent, digits = 10L)
n_pairs <- length(paired_difference)

# Validate the supplied difference column when present; downstream statistics
# always use the recomputed tumor-minus-adjacent difference.
difference_check_available <- difference_column %in% names(raw_data)
difference_check_max_error <- NA_real_
difference_check_reverse_max_error <- NA_real_
difference_check_pass <- NA
difference_check_reversed <- NA
if (difference_check_available) {
  supplied_difference <- to_numeric_strict(
    raw_data[[difference_column]],
    difference_column
  )
  if (anyNA(supplied_difference)) {
    stop("The supplied difference column contains missing values.", call. = FALSE)
  }
  difference_check_max_error <- max(
    abs(supplied_difference - paired_difference)
  )
  difference_check_reverse_max_error <- max(
    abs(supplied_difference + paired_difference)
  )
  difference_check_pass <-
    difference_check_max_error <= difference_check_tolerance
  difference_check_reversed <-
    difference_check_reverse_max_error <= difference_check_tolerance
  if (!difference_check_pass) {
    if (difference_check_reversed) {
      warning(
        "The supplied difference column matches adjacent-minus-tumor, which ",
        "is opposite to its tumor-minus-adjacent header. Recomputed ",
        "tumor-minus-adjacent values will be used."
      )
    } else {
      warning(
        "The supplied difference column does not agree with the recomputed ",
        "tumor-minus-adjacent difference. Recomputed values will be used."
      )
    }
  }
}

data_wide <- data.frame(
  sample_id = sample_id,
  adjacent = adjacent,
  tumor = tumor,
  difference = paired_difference,
  stringsAsFactors = FALSE
)

# Retain source identifiers exactly. This note is logged instead of silently
# padding numeric-looking identifiers with leading zeroes.
numeric_like_ids <- grepl("^[0-9]+$", data_wide$sample_id)
possible_leading_zero_note <- if (
  all(numeric_like_ids) && any(nchar(data_wide$sample_id) < 8L)
) {
  paste(
    "Source IDs were retained verbatim; some contain fewer than eight digits.",
    "Confirm against the specimen registry before publication if leading",
    "zeroes are meaningful."
  )
} else {
  "Source IDs were retained verbatim."
}

# ---- 3. Paired statistical analysis -----------------------------------------
if (n_pairs >= 3L && n_pairs <= 5000L) {
  shapiro_result <- stats::shapiro.test(data_wide$difference)
} else {
  shapiro_result <- NULL
}

paired_t_result <- stats::t.test(
  data_wide$tumor,
  data_wide$adjacent,
  paired = TRUE,
  alternative = alternative,
  conf.level = confidence_level
)

# Base R can calculate an exact Wilcoxon signed-rank P value when there are no
# zero differences or tied absolute differences. Make that decision explicit.
nonzero_difference <- data_wide$difference[data_wide$difference != 0]
can_use_exact_wilcoxon <-
  n_pairs < 50L &&
  length(nonzero_difference) == n_pairs &&
  !anyDuplicated(abs(nonzero_difference))

wilcoxon_warnings <- character()
wilcoxon_result <- withCallingHandlers(
  tryCatch(
    stats::wilcox.test(
      data_wide$tumor,
      data_wide$adjacent,
      paired = TRUE,
      alternative = alternative,
      exact = can_use_exact_wilcoxon,
      correct = !can_use_exact_wilcoxon,
      conf.int = TRUE,
      conf.level = confidence_level
    ),
    error = function(e) {
      stats::wilcox.test(
        data_wide$tumor,
        data_wide$adjacent,
        paired = TRUE,
        alternative = alternative,
        exact = can_use_exact_wilcoxon,
        correct = !can_use_exact_wilcoxon,
        conf.int = FALSE
      )
    }
  ),
  warning = function(w) {
    wilcoxon_warnings <<- c(wilcoxon_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

absolute_ranks <- rank(abs(nonzero_difference), ties.method = "average")
positive_rank_sum <- sum(absolute_ranks[nonzero_difference > 0])
negative_rank_sum <- sum(absolute_ranks[nonzero_difference < 0])
rank_biserial <-
  (positive_rank_sum - negative_rank_sum) /
  (positive_rank_sum + negative_rank_sum)

cohen_dz <- mean(data_wide$difference) / stats::sd(data_wide$difference)

if (!primary_test %in% c("paired_t", "wilcox")) {
  stop("primary_test must be either 'paired_t' or 'wilcox'.", call. = FALSE)
}
if (!annotation_mode %in% c("primary", "both")) {
  stop("annotation_mode must be either 'primary' or 'both'.", call. = FALSE)
}

primary_result <- if (primary_test == "paired_t") {
  paired_t_result
} else {
  wilcoxon_result
}
primary_p_value <- unname(primary_result$p.value)
primary_method_label <- if (primary_test == "paired_t") {
  "Paired t-test"
} else {
  "Wilcoxon signed-rank test"
}

format_p_value <- function(p_value) {
  if (is.na(p_value)) return("P = NA")
  if (p_value < 0.0001) return("P < 0.0001")
  if (p_value < 0.01) return(sprintf("P = %.4f", p_value))
  sprintf("P = %.3f", p_value)
}

wilcoxon_figure_label <- if (can_use_exact_wilcoxon) {
  "Exact Wilcoxon signed-rank test"
} else {
  "Wilcoxon signed-rank test"
}
primary_annotation <- if (annotation_mode == "both") {
  paste0(
    "Paired t-test: ", format_p_value(paired_t_result$p.value),
    "\n", wilcoxon_figure_label, ": ",
    format_p_value(wilcoxon_result$p.value)
  )
} else {
  paste(
    primary_method_label,
    format_p_value(primary_p_value),
    sep = "\n"
  )
}
annotation_text_size <- if (annotation_mode == "both") 2.55 else 3.0

# ---- 4. Build a machine-readable statistics table ---------------------------
result_rows <- list()
add_result <- function(
  section,
  group = "",
  metric,
  estimate = NA_real_,
  conf_low = NA_real_,
  conf_high = NA_real_,
  statistic = NA_real_,
  df = NA_real_,
  p_value = NA_real_,
  method = "",
  note = ""
) {
  result_rows[[length(result_rows) + 1L]] <<- data.frame(
    section = section,
    group = group,
    metric = metric,
    estimate = estimate,
    conf_low = conf_low,
    conf_high = conf_high,
    statistic = statistic,
    df = df,
    p_value = p_value,
    method = method,
    note = note,
    stringsAsFactors = FALSE
  )
}

add_result("data_qc", metric = "n_pairs", estimate = n_pairs)
add_result(
  "data_qc",
  metric = "source_column_count",
  estimate = raw_column_count
)
add_result(
  "data_qc",
  metric = "ignored_completely_empty_columns",
  estimate = ignored_empty_column_count
)
add_result(
  "data_qc",
  metric = "difference_check_max_absolute_error",
  estimate = difference_check_max_error,
  note = if (isTRUE(difference_check_pass)) "pass" else if (
    identical(difference_check_pass, FALSE)
  ) "fail; recomputed difference used" else "not available"
)
add_result(
  "data_qc",
  metric = "difference_check_reverse_max_absolute_error",
  estimate = difference_check_reverse_max_error,
  note = if (isTRUE(difference_check_reversed)) {
    "supplied difference matches adjacent-minus-tumor; recomputed tumor-minus-adjacent used"
  } else if (identical(difference_check_reversed, FALSE)) {
    "not a complete sign reversal"
  } else {
    "not available"
  }
)
add_result(
  "data_qc",
  metric = "pairs_with_tumor_higher",
  estimate = sum(data_wide$difference > 0)
)
add_result(
  "data_qc",
  metric = "pairs_with_tumor_lower",
  estimate = sum(data_wide$difference < 0)
)

group_values <- list(
  adjacent_normal = data_wide$adjacent,
  tumor = data_wide$tumor,
  paired_difference_tumor_minus_adjacent = data_wide$difference
)

for (group_name in names(group_values)) {
  values <- group_values[[group_name]]
  add_result("descriptive", group_name, "n", length(values))
  add_result("descriptive", group_name, "mean", mean(values))
  add_result("descriptive", group_name, "sd", stats::sd(values))
  add_result("descriptive", group_name, "median", stats::median(values))
  add_result(
    "descriptive", group_name, "q1",
    unname(stats::quantile(values, 0.25, names = FALSE))
  )
  add_result(
    "descriptive", group_name, "q3",
    unname(stats::quantile(values, 0.75, names = FALSE))
  )
  add_result("descriptive", group_name, "minimum", min(values))
  add_result("descriptive", group_name, "maximum", max(values))
}

if (!is.null(shapiro_result)) {
  add_result(
    "assumption_check",
    "paired_difference_tumor_minus_adjacent",
    "Shapiro-Wilk",
    statistic = unname(shapiro_result$statistic),
    p_value = shapiro_result$p.value,
    method = shapiro_result$method,
    note = "Normality check applies to paired differences, not raw groups."
  )
}

add_result(
  "paired_test",
  "tumor_vs_adjacent_normal",
  "mean_difference",
  estimate = unname(paired_t_result$estimate),
  conf_low = unname(paired_t_result$conf.int[[1L]]),
  conf_high = unname(paired_t_result$conf.int[[2L]]),
  statistic = unname(paired_t_result$statistic),
  df = unname(paired_t_result$parameter),
  p_value = paired_t_result$p.value,
  method = paired_t_result$method,
  note = paste0("Cohen dz = ", formatC(cohen_dz, digits = 6L, format = "fg"))
)

wilcoxon_estimate <- if (!is.null(wilcoxon_result$estimate)) {
  unname(wilcoxon_result$estimate)
} else {
  NA_real_
}
wilcoxon_conf_low <- if (!is.null(wilcoxon_result$conf.int)) {
  unname(wilcoxon_result$conf.int[[1L]])
} else {
  NA_real_
}
wilcoxon_conf_high <- if (!is.null(wilcoxon_result$conf.int)) {
  unname(wilcoxon_result$conf.int[[2L]])
} else {
  NA_real_
}

add_result(
  "paired_test",
  "tumor_vs_adjacent_normal",
  "Hodges-Lehmann_pseudomedian_difference",
  estimate = wilcoxon_estimate,
  conf_low = wilcoxon_conf_low,
  conf_high = wilcoxon_conf_high,
  statistic = unname(wilcoxon_result$statistic),
  p_value = wilcoxon_result$p.value,
  method = wilcoxon_result$method,
  note = paste0(
    if (can_use_exact_wilcoxon) "Exact P value; " else "Asymptotic P value; ",
    "matched rank-biserial = ",
    trimws(formatC(rank_biserial, digits = 6L, format = "fg")),
    if (length(wilcoxon_warnings) > 0L) {
      paste0("; warning: ", paste(unique(wilcoxon_warnings), collapse = " | "))
    } else {
      ""
    }
  )
)

add_result(
  "primary_analysis",
  "tumor_vs_adjacent_normal",
  primary_test,
  p_value = primary_p_value,
  method = primary_method_label,
  note = if (annotation_mode == "both") {
    paste(
      "Primary test is identified here;",
      "figures display paired t-test and Wilcoxon signed-rank P values",
      "at the author's request."
    )
  } else {
    "Only this pre-specified P value is displayed on the main figure."
  }
)

statistics_table <- do.call(rbind, result_rows)
statistics_file <- file.path(
  results_dir,
  paste0(marker_name, "_Hscore_paired_statistics_", analysis_date, ".csv")
)

write_csv_utf8_bom <- function(x, output_file) {
  temp_file <- tempfile(fileext = ".csv")
  on.exit(unlink(temp_file), add = TRUE)
  utils::write.csv(
    x,
    temp_file,
    row.names = FALSE,
    na = "",
    fileEncoding = "UTF-8"
  )
  csv_size <- file.info(temp_file)$size
  csv_bytes <- readBin(temp_file, what = "raw", n = csv_size)
  output_connection <- file(output_file, open = "wb")
  on.exit(close(output_connection), add = TRUE)
  writeBin(as.raw(c(0xEF, 0xBB, 0xBF)), output_connection)
  writeBin(csv_bytes, output_connection)
}
write_csv_utf8_bom(statistics_table, statistics_file)

# ---- 5. Prepare plotting data and style -------------------------------------
adjacent_label <- unname(plot_group_labels[["adjacent_normal"]])
tumor_label <- unname(plot_group_labels[["tumor"]])
if (any(!nzchar(c(adjacent_label, tumor_label))) || adjacent_label == tumor_label) {
  stop("plot_group_labels must contain two distinct, non-empty labels.")
}
group_levels <- c(adjacent_label, tumor_label)
group_colors <- stats::setNames(npg_colors, group_levels)

plot_data <- rbind(
  data.frame(
    sample_id = data_wide$sample_id,
    group = adjacent_label,
    hscore = data_wide$adjacent,
    stringsAsFactors = FALSE
  ),
  data.frame(
    sample_id = data_wide$sample_id,
    group = tumor_label,
    hscore = data_wide$tumor,
    stringsAsFactors = FALSE
  )
)
plot_data$group <- factor(plot_data$group, levels = group_levels)

if (.Platform$OS.type == "windows") {
  grDevices::windowsFonts(
    `Microsoft YaHei` = grDevices::windowsFont("Microsoft YaHei")
  )
  figure_font <- "Microsoft YaHei"
} else {
  figure_font <- "sans"
}

theme_publication <- function() {
  ggplot2::theme_classic(base_size = 9, base_family = figure_font) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(color = "#202020", linewidth = 0.45),
      axis.ticks = ggplot2::element_line(color = "#202020", linewidth = 0.40),
      axis.ticks.length = grid::unit(1.8, "mm"),
      axis.text = ggplot2::element_text(color = "#202020", size = 8.5),
      axis.text.x = ggplot2::element_text(size = 9),
      axis.title.y = ggplot2::element_text(size = 9.5, margin = ggplot2::margin(r = 6)),
      plot.title = ggplot2::element_text(
        face = "bold", size = 11, hjust = 0.5, margin = ggplot2::margin(b = 2)
      ),
      plot.subtitle = ggplot2::element_text(
        size = 8.5, hjust = 0.5, color = "#404040", margin = ggplot2::margin(b = 5)
      ),
      plot.caption = ggplot2::element_text(
        size = 7.2, hjust = 0, color = "#555555", margin = ggplot2::margin(t = 4)
      ),
      legend.position = "none",
      plot.margin = ggplot2::margin(5, 8, 5, 6)
    )
}

data_maximum <- max(plot_data$hscore)
vertical_scale <- max(data_maximum, 20)
bracket_y <- data_maximum + 0.075 * vertical_scale
bracket_tick <- 0.025 * vertical_scale
annotation_y <- bracket_y + 0.060 * vertical_scale
y_axis_upper <- annotation_y + 0.065 * vertical_scale
y_axis_breaks <- pretty(c(0, y_axis_upper), n = 6)
y_axis_breaks <- y_axis_breaks[y_axis_breaks >= 0 & y_axis_breaks <= y_axis_upper]

add_significance_bracket <- function(plot_object) {
  plot_object +
    ggplot2::annotate(
      "segment", x = 1, xend = 2, y = bracket_y, yend = bracket_y,
      linewidth = 0.5, color = "#202020"
    ) +
    ggplot2::annotate(
      "segment", x = 1, xend = 1, y = bracket_y,
      yend = bracket_y - bracket_tick, linewidth = 0.5, color = "#202020"
    ) +
    ggplot2::annotate(
      "segment", x = 2, xend = 2, y = bracket_y,
      yend = bracket_y - bracket_tick, linewidth = 0.5, color = "#202020"
    ) +
    ggplot2::annotate(
      "text",
      x = 1.5,
      y = annotation_y,
      label = primary_annotation,
      family = figure_font,
      size = annotation_text_size,
      lineheight = 0.95,
      color = "#202020"
    )
}

mean_sd <- function(x) {
  center <- mean(x)
  spread <- stats::sd(x)
  data.frame(y = center, ymin = max(0, center - spread), ymax = center + spread)
}

median_iqr <- function(x) {
  data.frame(
    y = stats::median(x),
    ymin = unname(stats::quantile(x, 0.25, names = FALSE)),
    ymax = unname(stats::quantile(x, 0.75, names = FALSE))
  )
}

summary_function <- if (primary_test == "paired_t") mean_sd else median_iqr
summary_center <- if (primary_test == "paired_t") mean else stats::median
base_layers <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(x = group, y = hscore)
) +
  ggplot2::geom_line(
    ggplot2::aes(group = sample_id),
    color = "#9AA0A6",
    linewidth = 0.45,
    alpha = 0.72
  ) +
  ggplot2::geom_point(
    ggplot2::aes(fill = group),
    shape = 21,
    color = "#202020",
    size = 2.6,
    stroke = 0.55,
    alpha = 0.96
  ) +
  ggplot2::stat_summary(
    fun.data = summary_function,
    geom = "errorbar",
    width = 0.13,
    linewidth = 0.65,
    color = "#202020"
  ) +
  ggplot2::stat_summary(
    fun = summary_center,
    geom = "point",
    shape = 95,
    size = 7,
    color = "#202020"
  ) +
  ggplot2::scale_fill_manual(values = group_colors, guide = "none") +
  ggplot2::scale_y_continuous(
    limits = c(0, y_axis_upper),
    breaks = y_axis_breaks,
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::labs(
    title = marker_name,
    subtitle = paste0("n = ", n_pairs, " biological pairs"),
    x = NULL,
    y = "H-score"
  ) +
  theme_publication()

paired_plot <- add_significance_bracket(base_layers)

# Reference-style variant retaining the user's preferred mean +/- SD bars.
bar_summary <- data.frame(
  group = factor(group_levels, levels = group_levels),
  mean = c(mean(data_wide$adjacent), mean(data_wide$tumor)),
  sd = c(stats::sd(data_wide$adjacent), stats::sd(data_wide$tumor))
)
bar_summary$lower <- pmax(0, bar_summary$mean - bar_summary$sd)
bar_summary$upper <- bar_summary$mean + bar_summary$sd

bar_plot <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(x = group, y = hscore)
) +
  ggplot2::geom_col(
    data = bar_summary,
    ggplot2::aes(x = group, y = mean, fill = group),
    inherit.aes = FALSE,
    width = 0.58,
    color = "#202020",
    linewidth = 0.5,
    alpha = 0.66
  ) +
  ggplot2::geom_line(
    ggplot2::aes(group = sample_id),
    color = "#9AA0A6",
    linewidth = 0.45,
    alpha = 0.72
  ) +
  ggplot2::geom_errorbar(
    data = bar_summary,
    ggplot2::aes(x = group, ymin = lower, ymax = upper),
    inherit.aes = FALSE,
    width = 0.16,
    linewidth = 0.68,
    color = "#202020"
  ) +
  ggplot2::geom_point(
    shape = 21,
    fill = "white",
    color = "#202020",
    size = 2.6,
    stroke = 0.55,
    alpha = 0.96
  ) +
  ggplot2::scale_fill_manual(values = group_colors, guide = "none") +
  ggplot2::scale_y_continuous(
    limits = c(0, y_axis_upper),
    breaks = y_axis_breaks,
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::labs(
    title = marker_name,
    subtitle = paste0("n = ", n_pairs, " biological pairs"),
    x = NULL,
    y = "H-score"
  ) +
  theme_publication()
bar_plot <- add_significance_bracket(bar_plot)

# Optional pure-bar comparison requested by the author. It deliberately omits
# all individual points, pairing lines, error bars, and the significance bracket;
# bar height is the mean. Its vertical scale is based on the displayed means so
# that the two-line test annotation sits in balanced whitespace above the bars.
bar_only_mean_maximum <- max(bar_summary$mean)
bar_only_vertical_scale <- max(bar_only_mean_maximum, 20)
bar_only_annotation_y <- bar_only_mean_maximum + 0.30 * bar_only_vertical_scale
bar_only_y_axis_upper <- bar_only_mean_maximum + 0.48 * bar_only_vertical_scale
bar_only_y_axis_breaks <- pretty(c(0, bar_only_y_axis_upper), n = 6)
bar_only_y_axis_breaks <- bar_only_y_axis_breaks[
  bar_only_y_axis_breaks >= 0 & bar_only_y_axis_breaks <= bar_only_y_axis_upper
]

add_bar_only_annotation <- function(plot_object) {
  plot_object +
    ggplot2::annotate(
      "text",
      x = 1.5,
      y = bar_only_annotation_y,
      label = primary_annotation,
      family = figure_font,
      size = annotation_text_size,
      lineheight = 1.05,
      color = "#202020"
    )
}

bar_only_plot <- ggplot2::ggplot(
  bar_summary,
  ggplot2::aes(x = group, y = mean, fill = group)
) +
  ggplot2::geom_col(
    width = 0.58,
    color = "#202020",
    linewidth = 0.5,
    alpha = 0.72
  ) +
  ggplot2::scale_fill_manual(values = group_colors, guide = "none") +
  ggplot2::scale_y_continuous(
    limits = c(0, bar_only_y_axis_upper),
    breaks = bar_only_y_axis_breaks,
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::labs(
    title = marker_name,
    subtitle = paste0("n = ", n_pairs, " biological pairs"),
    x = NULL,
    y = "H-score"
  ) +
  theme_publication()
bar_only_plot <- add_bar_only_annotation(bar_only_plot)

# ---- 6. Export publication and QC figures ----------------------------------
save_figure_set <- function(plot_object, stem) {
  png_file <- file.path(figures_dir, paste0(stem, ".png"))
  tiff_file <- file.path(figures_dir, paste0(stem, ".tiff"))
  pdf_file <- file.path(figures_dir, paste0(stem, ".pdf"))

  png_args <- list(
    filename = png_file,
    plot = plot_object,
    device = "png",
    width = figure_width_mm,
    height = figure_height_mm,
    units = "mm",
    dpi = png_dpi,
    bg = "white"
  )
  if (isTRUE(capabilities("cairo"))) png_args$type <- "cairo-png"
  do.call(ggplot2::ggsave, png_args)

  tiff_args <- list(
    filename = tiff_file,
    plot = plot_object,
    device = "tiff",
    width = figure_width_mm,
    height = figure_height_mm,
    units = "mm",
    dpi = tiff_dpi,
    compression = "lzw",
    bg = "white"
  )
  if (isTRUE(capabilities("cairo"))) tiff_args$type <- "cairo"
  do.call(ggplot2::ggsave, tiff_args)

  ggplot2::ggsave(
    filename = pdf_file,
    plot = plot_object,
    device = grDevices::cairo_pdf,
    width = figure_width_mm,
    height = figure_height_mm,
    units = "mm",
    bg = "white"
  )

  invisible(c(png = png_file, tiff = tiff_file, pdf = pdf_file))
}

paired_stem <- paste0(
  marker_name,
  "_Hscore_paired_plot_publication_",
  analysis_date
)
bar_stem <- paste0(
  marker_name,
  "_Hscore_paired_bar_reference_",
  analysis_date
)
paired_files <- save_figure_set(paired_plot, paired_stem)
bar_files <- save_figure_set(bar_plot, bar_stem)
bar_only_files <- NULL
if (export_bar_only) {
  annotation_suffix <- if (annotation_mode == "both") {
    "both_tests"
  } else {
    "primary_test"
  }
  bar_only_stem <- paste0(
    marker_name,
    "_Hscore_paired_bar_only_",
    annotation_suffix,
    "_",
    analysis_date
  )
  bar_only_files <- save_figure_set(bar_only_plot, bar_only_stem)
}

# QQ plot of paired differences: assumption/QC output, not a manuscript panel.
qq_data <- data.frame(difference = data_wide$difference)
qq_label <- if (is.null(shapiro_result)) {
  "Shapiro-Wilk not run"
} else {
  paste0(
    "Shapiro-Wilk W = ",
    formatC(unname(shapiro_result$statistic), digits = 4, format = "f"),
    ", ", format_p_value(shapiro_result$p.value)
  )
}
qq_plot <- ggplot2::ggplot(
  qq_data,
  ggplot2::aes(sample = difference)
) +
  ggplot2::stat_qq(
    shape = 21,
    size = 2.4,
    stroke = 0.5,
    fill = npg_colors[[2L]],
    color = "#202020"
  ) +
  ggplot2::stat_qq_line(color = "#4B5563", linewidth = 0.55) +
  ggplot2::labs(
    title = paste0(marker_name, " paired-difference QQ plot"),
    subtitle = qq_label,
    x = "Theoretical quantiles",
    y = "Sample quantiles"
  ) +
  theme_publication() +
  ggplot2::theme(legend.position = "none")

qq_file <- file.path(
  qc_dir,
  paste0(marker_name, "_Hscore_paired_difference_QQ_QC_", analysis_date, ".png")
)
ggplot2::ggsave(
  filename = qq_file,
  plot = qq_plot,
  device = "png",
  type = if (isTRUE(capabilities("cairo"))) "cairo-png" else NULL,
  width = 90,
  height = 78,
  units = "mm",
  dpi = png_dpi,
  bg = "white"
)

# ---- 7. Reproducibility log and console summary -----------------------------
session_file <- file.path(
  logs_dir,
  paste0(script_index, "_", marker_name, "_Hscore_sessionInfo_", analysis_date, ".txt")
)

input_md5 <- unname(tools::md5sum(input_file))
log_lines <- c(
  paste0("Run timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
  paste0("Script: ", if (is.na(script_path)) "interactive" else basename(script_path)),
  paste0("Input: ", basename(input_file)),
  paste0("Input MD5: ", input_md5),
  paste0("Pairs: ", n_pairs),
  paste0("Primary test: ", primary_method_label),
  paste0("Primary P value: ", format(primary_p_value, digits = 12)),
  paste0("Figure annotation mode: ", annotation_mode),
  paste0("Pure bar-only figure exported: ", export_bar_only),
  paste0("Exact Wilcoxon used: ", can_use_exact_wilcoxon),
  paste0("Difference-column check pass: ", difference_check_pass),
  paste0("Difference-column maximum absolute error: ", difference_check_max_error),
  paste0("Difference-column reversed direction: ", difference_check_reversed),
  paste0("Difference-column reverse maximum absolute error: ", difference_check_reverse_max_error),
  paste0("Source column count: ", raw_column_count),
  paste0("Ignored completely empty columns: ", ignored_empty_column_count),
  paste0("Identifier note: ", possible_leading_zero_note),
  "",
  "Session information:",
  capture.output(utils::sessionInfo())
)
writeLines(log_lines, con = session_file, useBytes = TRUE)

message("Analysis completed successfully.")
message("Input: ", input_file)
message("Pairs: ", n_pairs)
message(
  "Primary result: ", primary_method_label, ", ",
  format_p_value(primary_p_value)
)
message("Statistics: ", statistics_file)
message("Main paired figure: ", paired_files[["png"]])
message("Reference bar figure: ", bar_files[["png"]])
if (!is.null(bar_only_files)) {
  message("Pure bar-only figure: ", bar_only_files[["png"]])
}
message("QQ QC figure: ", qq_file)
message("Session log: ", session_file)
