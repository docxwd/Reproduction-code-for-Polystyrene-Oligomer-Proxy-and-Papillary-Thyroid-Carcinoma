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
# Immune-cell fractions, group comparisons (S16) and Figure 8A-C.
options(stringsAsFactors = FALSE)
require_packages(c("reshape2", "ggplot2", "ggpubr", "dplyr", "RColorBrewer", "linkET", "pheatmap"))
suppressPackageStartupMessages({
  library(reshape2); library(ggplot2); library(ggpubr); library(dplyr)
  library(RColorBrewer); library(linkET); library(pheatmap)
})
source(file.path(root, "scripts/lib/plot_helpers.R"))

out <- safe_output(if (length(args)) args[1] else file.path(root, "reproduced_outputs/immune_groups"), root)
output_dir <- safe_output(if (length(args) >= 2L) args[2] else file.path(out, "figures"), root)
fraction_path <- if (length(args) >= 3L) args[3] else file.path(root, "derived_outputs/cibersort_retained_samples.csv")
retained <- load_fractions(fraction_path)
dat <- retained
cell_names <- setdiff(names(retained), c("id", "Group"))
cells <- cell_names
fraction_matrix <- as.matrix(retained[cell_names])
storage.mode(fraction_matrix) <- "double"
rownames(fraction_matrix) <- retained$id
group <- do.call(rbind, lapply(sort(cells), function(cell) {
  x <- dat[dat$Group=="Control", cell]; y <- dat[dat$Group=="Tumor", cell]
  data.frame(Immune_cell=cell, Control_n=length(x), PTC_n=length(y),
    Control_mean=mean(x), PTC_mean=mean(y), Control_median=median(x), PTC_median=median(y),
    Control_SD=sd(x), PTC_SD=sd(y), Raw_P=suppressWarnings(wilcox.test(x,y)$p.value))
}))
group$BH_FDR_22_cells <- bh_valid(group$Raw_P)
group$Direction_in_PTC <- ifelse(group$PTC_mean > group$Control_mean, "Increased", "Decreased")
group$Significant_raw_P_lt_0.05 <- significance(group$Raw_P)
group$Significant_BH_FDR_lt_0.05 <- significance(group$BH_FDR_22_cells)

write_data(group, file.path(out, "immune_group_comparisons_with_BH_FDR.csv"))
group_stats <- group[match(cell_names, group$Immune_cell), , drop = FALSE]
# -----------------------------------------------------------------------------
# Figure 8A: estimated immune-cell fractions
# -----------------------------------------------------------------------------
plot_data <- retained
plot_data$Sample <- plot_data$id
plot_data$Group <- factor(plot_data$Group, levels = c("Control", "Tumor"))
control_samples <- plot_data$id[plot_data$Group == "Control"]
tumor_samples <- plot_data$id[plot_data$Group == "Tumor"]
all_samples_ordered <- c(control_samples, "gap", tumor_samples)

data_long <- reshape2::melt(
  plot_data[, c("Sample", "Group", cell_names)],
  id.vars = c("Sample", "Group"),
  variable.name = "Immune",
  value.name = "Fraction"
)
data_long$Sample <- factor(data_long$Sample, levels = all_samples_ordered)
data_long$Immune <- factor(data_long$Immune, levels = cell_names)

my_colors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(length(cell_names))
names(my_colors) <- cell_names

figure8a <- ggplot(data_long, aes(x = Sample, y = Fraction, fill = Immune)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_colors) +
  scale_x_discrete(drop = FALSE) +
  theme_minimal(base_size = 18) +
  theme(
    text = element_text(face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    x = NULL,
    y = "Relative Percent",
    fill = "Immune\nCell Type",
    title = "Immune Cell Distribution",
    subtitle = "Control vs. Tumor"
  ) +
  coord_cartesian(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.05))) +
  annotate(
    "segment", x = 0.5, xend = length(control_samples) + 0.5,
    y = -0.04, yend = -0.04, color = "#D65DB1", linewidth = 5
  ) +
  annotate(
    "text", x = length(control_samples) / 2 + 0.5, y = -0.08,
    label = "Control", color = "#D65DB1", size = 7, fontface = "bold"
  ) +
  annotate(
    "segment", x = length(control_samples) + 1.5,
    xend = length(control_samples) + length(tumor_samples) + 1.5,
    y = -0.04, yend = -0.04, color = "#0089BA", linewidth = 5
  ) +
  annotate(
    "text",
    x = length(control_samples) + length(tumor_samples) / 2 + 1.5,
    y = -0.08, label = "Tumor", color = "#0089BA",
    size = 7, fontface = "bold"
  )

figure8a_files <- save_gg_panel(
  figure8a,
  "Figure8A_immune_distribution",
  width = 12,
  height = 7.5
)

# -----------------------------------------------------------------------------
# Figure 8B: group comparisons with BH-adjusted significance
# -----------------------------------------------------------------------------
label_y <- max(data_long$Fraction, na.rm = TRUE) * 1.08
figure8b_labels <- data.frame(
  Immune = factor(group_stats$Immune_cell, levels = cell_names),
  q_label = significance_symbol(group_stats$BH_FDR_22_cells),
  y.position = label_y
)

my_cute_colors <- c("Control" = "#FFC0CB", "Tumor" = "#87CEFA")
figure8b <- ggpubr::ggboxplot(
  data_long,
  x = "Immune",
  y = "Fraction",
  fill = "Group",
  palette = my_cute_colors,
  xlab = "",
  ylab = "Fraction",
  legend.title = "Group",
  notch = FALSE,
  width = 0.8,
  outlier.shape = NA
) +
  geom_point(
    aes(group = Group),
    color = "#222222",
    alpha = 0.55,
    size = 1.1,
    position = position_jitterdodge(
      jitter.width = 0.12,
      jitter.height = 0,
      dodge.width = 0.8,
      seed = 20260127
    ),
    show.legend = FALSE
  ) +
  geom_text(
    data = figure8b_labels,
    aes(x = Immune, y = y.position, label = q_label),
    inherit.aes = FALSE,
    size = 4.0
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1),
    plot.margin = margin(8, 18, 18, 18)
  ) +
  labs(
    title = "Immune Cell Comparison",
    subtitle = paste0(
      "(Control n=", length(control_samples),
      ", Tumor n=", length(tumor_samples), ")"
    )
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.16))) +
  coord_cartesian(clip = "off")

figure8b_files <- save_gg_panel(
  figure8b,
  "Figure8B_group_comparison_BH_q",
  width = 12,
  height = 7.5
)

# -----------------------------------------------------------------------------
# Figure 8C: descriptive Pearson correlation matrix without significance tests
# -----------------------------------------------------------------------------
cell_cor_matrix <- cor(fraction_matrix, method = "pearson")
figure8c_labels <- matrix(
  sprintf("%.2f", cell_cor_matrix),
  nrow = nrow(cell_cor_matrix),
  ncol = ncol(cell_cor_matrix),
  dimnames = dimnames(cell_cor_matrix)
)
figure8c_colors <- colorRampPalette(c("skyblue", "white", "#ba6262"))(100)

draw_figure8c <- function() {
  pheatmap::pheatmap(
    cell_cor_matrix,
    color = figure8c_colors,
    display_numbers = figure8c_labels,
    number_color = "black",
    fontsize_number = 8,
    main = "Immune Cells Correlation Heatmap",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    fontsize = 12,
    labels_row = rownames(cell_cor_matrix),
    labels_col = colnames(cell_cor_matrix),
    border_color = NA
  )
}

figure8c_files <- save_grid_panel(
  draw_figure8c,
  "Figure8C_immune_cell_correlation_descriptive",
  width = 12,
  height = 12
)


capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
cat("Immune group statistics and three Figure 8 panels completed.\n")

