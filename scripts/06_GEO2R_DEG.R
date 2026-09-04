# Official extraction and visualization of GEO2R differential expression results
#
# Input:
#   GSE33630_GEO2R_top_table.tsv
#
# DEG criteria:
#   adj.P.Val < 0.05 and |logFC| > 1
#
# Notes:
#   This script does not rerun differential expression analysis. It uses the
#   GEO2R-exported table as the only input and performs downstream filtering,
#   gene-symbol parsing, summary statistics, and figure generation.

required_packages <- c("ggplot2", "dplyr", "ggrepel", "writexl")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Missing R package(s): ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before execution as documented in environment/requirements-r.txt."
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
})

cli_args <- commandArgs(trailingOnly = TRUE)
work_dir <- if (length(cli_args) >= 1L) cli_args[[1L]] else getwd()
if (!dir.exists(work_dir)) stop("Analysis directory not found: ", work_dir)
setwd(work_dir)
input_file <- if (length(cli_args) >= 2L) cli_args[[2L]] else "GSE33630_GEO2R_top_table.tsv"
adj_p_cutoff <- 0.05
logfc_cutoff <- 1
top_label_n <- 8

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

geo2r <- read.delim(
  input_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_cols <- c("ID", "adj.P.Val", "P.Value", "t", "B", "logFC", "Gene.symbol", "Gene.title")
missing_cols <- setdiff(required_cols, colnames(geo2r))
if (length(missing_cols) > 0) {
  stop("Missing required column(s): ", paste(missing_cols, collapse = ", "))
}

geo2r <- geo2r %>%
  mutate(
    adj.P.Val = as.numeric(adj.P.Val),
    P.Value = as.numeric(P.Value),
    logFC = as.numeric(logFC),
    neg_log10_adjP = -log10(pmax(adj.P.Val, .Machine$double.xmin)),
    Gene.symbol = trimws(as.character(Gene.symbol)),
    Gene.title = trimws(as.character(Gene.title)),
    DEG_status = case_when(
      adj.P.Val < adj_p_cutoff & logFC > logfc_cutoff ~ "Up",
      adj.P.Val < adj_p_cutoff & logFC < -logfc_cutoff ~ "Down",
      TRUE ~ "Not significant"
    )
  )

deg_probes <- geo2r %>%
  filter(DEG_status %in% c("Up", "Down")) %>%
  arrange(adj.P.Val, desc(abs(logFC)))

raw_gene_symbols <- deg_probes$Gene.symbol
raw_gene_symbols <- raw_gene_symbols[!is.na(raw_gene_symbols) & raw_gene_symbols != ""]
raw_gene_symbols <- sort(unique(raw_gene_symbols))

split_gene_symbols <- unlist(strsplit(raw_gene_symbols, "///|//|;|,|\\|", perl = TRUE))
split_gene_symbols <- trimws(split_gene_symbols)
split_gene_symbols <- split_gene_symbols[!is.na(split_gene_symbols) & split_gene_symbols != ""]
split_gene_symbols <- sort(unique(split_gene_symbols))

split_probe_rows <- lapply(seq_len(nrow(deg_probes)), function(i) {
  symbols <- unlist(strsplit(deg_probes$Gene.symbol[[i]], "///|//|;|,|\\|", perl = TRUE))
  symbols <- trimws(symbols)
  symbols <- symbols[!is.na(symbols) & symbols != ""]
  if (length(symbols) == 0) {
    return(NULL)
  }
  data.frame(
    Gene.symbol = symbols,
    Probe.ID = deg_probes$ID[[i]],
    logFC = deg_probes$logFC[[i]],
    adj.P.Val = deg_probes$adj.P.Val[[i]],
    P.Value = deg_probes$P.Value[[i]],
    DEG_status = deg_probes$DEG_status[[i]],
    Gene.title = deg_probes$Gene.title[[i]],
    stringsAsFactors = FALSE
  )
})

deg_gene_long <- bind_rows(split_probe_rows)

deg_gene_summary <- deg_gene_long %>%
  group_by(Gene.symbol) %>%
  arrange(adj.P.Val, desc(abs(logFC)), .by_group = TRUE) %>%
  summarise(
    Probe_count = n_distinct(Probe.ID),
    Representative_probe = Probe.ID[1],
    logFC = logFC[1],
    adj.P.Val = adj.P.Val[1],
    P.Value = P.Value[1],
    DEG_status = DEG_status[1],
    Direction_consistency = if (dplyr::n_distinct(DEG_status) == 1) DEG_status[1] else "Mixed",
    All_probe_IDs = paste(unique(Probe.ID), collapse = "|"),
    .groups = "drop"
  ) %>%
  arrange(adj.P.Val, desc(abs(logFC)))

deg_counts <- data.frame(
  Category = c("Up", "Down", "Significant probes", "Raw unique Gene.symbol", "Split unique genes"),
  Count = c(
    sum(deg_probes$DEG_status == "Up"),
    sum(deg_probes$DEG_status == "Down"),
    nrow(deg_probes),
    length(raw_gene_symbols),
    length(split_gene_symbols)
  )
)

write.csv(geo2r, "GEO2R_all_results_with_DEG_status.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(deg_probes, "GEO2R_DEG_probes_adjP0.05_logFC1_official.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(deg_gene_long, "GEO2R_DEG_gene_probe_long_adjP0.05_logFC1_official.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(deg_gene_summary, "GEO2R_DEG_gene_summary_adjP0.05_logFC1_official.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(deg_counts, "GEO2R_DEG_counts_adjP0.05_logFC1_official.csv", row.names = FALSE, fileEncoding = "UTF-8")
writeLines(raw_gene_symbols, "GEO2R_DEG_gene_symbols_raw_unique_adjP0.05_logFC1_official.txt", useBytes = TRUE)
writeLines(split_gene_symbols, "GEO2R_DEG_gene_symbols_split_unique_adjP0.05_logFC1_official.txt", useBytes = TRUE)

writexl::write_xlsx(
  list(
    DEG_counts = deg_counts,
    DEG_probes = deg_probes,
    DEG_gene_summary = deg_gene_summary,
    DEG_gene_probe_long = deg_gene_long
  ),
  "GEO2R_DEG_extract_summary_official.xlsx"
)

figure_dir <- "GEO2R_official_figures"
if (!dir.exists(figure_dir)) {
  dir.create(figure_dir, recursive = TRUE)
}

nature_blue <- "#4E79A7"
nature_red <- "#C44E52"
neutral_gray <- "#B8BFC7"
text_black <- "#222222"

top_labels <- bind_rows(
  geo2r %>%
    filter(DEG_status == "Up", !is.na(Gene.symbol), Gene.symbol != "") %>%
    mutate(Label = sub("///.*$", "", Gene.symbol)) %>%
    arrange(adj.P.Val, desc(abs(logFC))) %>%
    distinct(Label, .keep_all = TRUE) %>%
    arrange(adj.P.Val, desc(abs(logFC))) %>%
    slice_head(n = top_label_n),
  geo2r %>%
    filter(DEG_status == "Down", !is.na(Gene.symbol), Gene.symbol != "") %>%
    mutate(Label = sub("///.*$", "", Gene.symbol)) %>%
    arrange(adj.P.Val, desc(abs(logFC))) %>%
    distinct(Label, .keep_all = TRUE) %>%
    arrange(adj.P.Val, desc(abs(logFC))) %>%
    slice_head(n = top_label_n)
)

volcano_plot <- ggplot(
  geo2r,
  aes(x = logFC, y = neg_log10_adjP, color = DEG_status)
) +
  geom_point(size = 0.75, alpha = 0.70, stroke = 0) +
  geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = "dashed", linewidth = 0.35, color = "#666666") +
  geom_hline(yintercept = -log10(adj_p_cutoff), linetype = "dashed", linewidth = 0.35, color = "#666666") +
  ggrepel::geom_text_repel(
    data = top_labels,
    aes(label = Label),
    color = text_black,
    size = 2.15,
    family = "sans",
    box.padding = 0.25,
    point.padding = 0.12,
    min.segment.length = 0,
    segment.size = 0.20,
    segment.color = "#777777",
    max.overlaps = Inf,
    seed = 20260524,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c("Down" = nature_blue, "Not significant" = neutral_gray, "Up" = nature_red),
    breaks = c("Up", "Down", "Not significant")
  ) +
  labs(
    x = "log2 fold change",
    y = expression(-log[10]("adjusted P value")),
    color = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_family = "sans", base_size = 8.5) +
  theme(
    axis.line = element_line(linewidth = 0.35, color = text_black),
    axis.ticks = element_line(linewidth = 0.30, color = text_black),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.text = element_text(size = 7.2, color = text_black),
    legend.key.width = unit(4.5, "mm"),
    legend.key.height = unit(3.0, "mm"),
    legend.spacing.x = unit(1.5, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(5.5, 12, 5.5, 5.5),
    plot.background = element_rect(fill = "white", color = NA)
  )

count_plot_data <- data.frame(
  DEG_status = factor(c("Down", "Up"), levels = c("Down", "Up")),
  Count = c(sum(deg_probes$DEG_status == "Down"), sum(deg_probes$DEG_status == "Up"))
)

count_plot <- ggplot(count_plot_data, aes(x = DEG_status, y = Count, fill = DEG_status)) +
  geom_col(width = 0.62, color = text_black, linewidth = 0.25) +
  geom_text(aes(label = Count), vjust = -0.45, family = "sans", size = 3.0, color = text_black) +
  scale_fill_manual(values = c("Down" = nature_blue, "Up" = nature_red)) +
  labs(x = NULL, y = "Number of significant probes") +
  theme_classic(base_family = "sans", base_size = 8.5) +
  theme(
    legend.position = "none",
    axis.line = element_line(linewidth = 0.35, color = text_black),
    axis.ticks.x = element_blank(),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  expand_limits(y = max(count_plot_data$Count) * 1.12)

save_plot <- function(plot, filename, width_mm, height_mm) {
  ggplot2::ggsave(
    file.path(figure_dir, paste0(filename, ".jpeg")),
    plot = plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = 600,
    device = "jpeg",
    quality = 95,
    bg = "white"
  )
  ggplot2::ggsave(
    file.path(figure_dir, paste0(filename, ".tiff")),
    plot = plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = 600,
    compression = "lzw"
  )
  ggplot2::ggsave(
    file.path(figure_dir, paste0(filename, ".png")),
    plot = plot,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = 600
  )
}

save_plot(volcano_plot, "GEO2R_volcano_nature_official", 95, 78)
save_plot(count_plot, "GEO2R_DEG_count_bar_nature_official", 62, 58)

message("Official GEO2R DEG extraction and visualization completed.")
message("Input rows: ", nrow(geo2r))
message("DEG probe rows: ", nrow(deg_probes))
message("Up-regulated probe rows: ", sum(deg_probes$DEG_status == "Up"))
message("Down-regulated probe rows: ", sum(deg_probes$DEG_status == "Down"))
message("Raw unique Gene.symbol: ", length(raw_gene_symbols))
message("Split unique genes: ", length(split_gene_symbols))
