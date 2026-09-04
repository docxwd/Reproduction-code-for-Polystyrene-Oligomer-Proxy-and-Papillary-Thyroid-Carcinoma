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
# Pooled/PTC-only gene-fraction correlations (S17) and Figure 9A-B.
options(stringsAsFactors = FALSE)
require_packages(c("reshape2", "ggplot2", "ggpubr", "dplyr", "RColorBrewer", "linkET", "pheatmap"))
suppressPackageStartupMessages({
  library(reshape2); library(ggplot2); library(ggpubr); library(dplyr)
  library(RColorBrewer); library(linkET); library(pheatmap)
})
source(file.path(root, "scripts/lib/plot_helpers.R"))

out <- safe_output(if (length(args)) args[1] else file.path(root, "reproduced_outputs/gene_immune"), root)
output_dir <- safe_output(if (length(args) >= 2L) args[2] else file.path(out, "figures"), root)
fraction_path <- if (length(args) >= 3L) args[3] else file.path(root, "derived_outputs/cibersort_retained_samples.csv")
retained <- load_fractions(fraction_path)
dat <- retained
cell_names <- setdiff(names(retained), c("id", "Group"))
cells <- cell_names
fraction_matrix <- as.matrix(retained[cell_names])
storage.mode(fraction_matrix) <- "double"
rownames(fraction_matrix) <- retained$id

expression_path <- if (length(args) >= 4L) args[4] else file.path(root, "data/processed/immune_core_expression.csv")
expr_data <- read_data(expression_path)
stopifnot(!anyDuplicated(expr_data[[1]]))
expr <- as.matrix(expr_data[-1])
rownames(expr) <- expr_data[[1]]
storage.mode(expr) <- "double"
genes <- c("BCL2", "BAX", "FN1")
gene_names <- genes
stopifnot(all(genes %in% rownames(expr)), all(dat$id %in% colnames(expr)),
          all(is.finite(expr)))
spearman <- function(x,y) {
  if (length(unique(x))<2L || length(unique(y))<2L) return(c(rho=NA_real_, p=NA_real_))
  test <- cor.test(x,y,method="spearman",exact=FALSE)
  c(rho=unname(test$estimate), p=test$p.value)
}
tumor <- dat$Group=="Tumor"
rows <- list()
for (gene in genes) for (cell in cells) {
  x <- as.numeric(expr[gene, dat$id]); y <- dat[[cell]]
  pooled <- spearman(x,y); ptc <- spearman(x[tumor],y[tumor])
  rows[[length(rows)+1L]] <- data.frame(Gene=gene, Immune_cell=cell,
    All_samples_n=length(x), All_samples_rho=pooled[1], All_samples_P=pooled[2],
    PTC_only_n=sum(tumor), PTC_only_rho=ptc[1], PTC_only_P=ptc[2])
}
corr <- do.call(rbind, rows); rownames(corr) <- NULL
corr$All_samples_BH_FDR_66_tests <- bh_valid(corr$All_samples_P)
# Compatibility column name from the reference table; actual finite family is 63.
# Explicit test counts and estimability are exported below. No missing P is filled.
corr$PTC_only_BH_FDR_66_tests <- bh_valid(corr$PTC_only_P)
corr$All_samples_significant_FDR_lt_0.05 <- significance(corr$All_samples_BH_FDR_66_tests)
corr$PTC_only_significant_FDR_lt_0.05 <- significance(corr$PTC_only_BH_FDR_66_tests)
corr <- corr[c("Gene","Immune_cell","All_samples_n","All_samples_rho","All_samples_P",
 "All_samples_BH_FDR_66_tests","All_samples_significant_FDR_lt_0.05","PTC_only_n",
 "PTC_only_rho","PTC_only_P","PTC_only_BH_FDR_66_tests","PTC_only_significant_FDR_lt_0.05")]
metadata <- data.frame(Analysis=c("Pooled", "PTC_only"), Planned_pairs=66L,
 Estimable_pairs=c(sum(is.finite(corr$All_samples_P)),sum(is.finite(corr$PTC_only_P))),
 BH_denominator=c(sum(is.finite(corr$All_samples_P)),sum(is.finite(corr$PTC_only_P))))
status <- data.frame(Gene=corr$Gene, Immune_cell=corr$Immune_cell,
 Pooled_status=ifelse(is.na(corr$All_samples_P),"Not_estimable_constant_input","Estimable"),
 PTC_only_status=ifelse(is.na(corr$PTC_only_P),"Not_estimable_constant_input","Estimable"))

write_data(corr, file.path(out, "gene_immune_correlations_full_and_tumor_only_with_BH_FDR.csv"))
write_data(metadata, file.path(out, "gene_immune_test_counts.csv"))
write_data(status, file.path(out, "gene_immune_estimability.csv"))
corr_stats <- corr
# -----------------------------------------------------------------------------
# Figure 9A: cell-cell matrix and gene-fraction links
# -----------------------------------------------------------------------------
figure9a_links <- data.frame(
  spec = corr_stats$Gene,
  env = corr_stats$Immune_cell,
  signed_r = corr_stats$All_samples_rho,
  q = corr_stats$All_samples_BH_FDR_66_tests,
  stringsAsFactors = FALSE
)
figure9a_links$pd <- ifelse(
  !is.na(figure9a_links$q) & figure9a_links$q < 0.05,
  ifelse(
    figure9a_links$signed_r > 0,
    "positive correlation",
    "negative correlation"
  ),
  "not significant"
)
figure9a_links$r <- abs(figure9a_links$signed_r)
figure9a_links$rd <- cut(
  figure9a_links$r,
  breaks = c(-Inf, 0.2, 0.4, 0.6, Inf),
  labels = c("< 0.2", "0.2 - 0.4", "0.4 - 0.6", ">= 0.6")
)
figure9a_links$pd <- factor(
  figure9a_links$pd,
  levels = c("negative correlation", "not significant", "positive correlation")
)

immune_matrix_for_linket <- as.data.frame(fraction_matrix, check.names = FALSE)
figure9a <- linkET::qcorrplot(
  linkET::correlate(immune_matrix_for_linket, method = "spearman"),
  type = "lower",
  diag = FALSE
) +
  linkET::geom_square(size = 1.0, colour = "#222222") +
  linkET::geom_couple(
    aes(colour = pd, size = rd),
    data = figure9a_links,
    curvature = linkET::nice_curvature()
  ) +
  scale_fill_gradientn(
    colours = rev(RColorBrewer::brewer.pal(11, "Spectral"))
  ) +
  scale_size_manual(
    values = c("< 0.2" = 0.5, "0.2 - 0.4" = 1, "0.4 - 0.6" = 2, ">= 0.6" = 3),
    drop = FALSE
  ) +
  scale_colour_manual(
    values = c(
      "positive correlation" = "#FF1493",
      "negative correlation" = "#00CED1",
      "not significant" = "#999999"
    ),
    drop = FALSE
  ) +
  guides(
    size = guide_legend(
      title = "abs(Cor)",
      override.aes = list(colour = "grey35"),
      order = 2
    ),
    colour = guide_legend(
      title = "BH q-value",
      override.aes = list(size = 3),
      order = 1
    ),
    fill = guide_colorbar(title = "Cell-cell cor", order = 3)
  ) +
  labs(
    x = "immune infiltrating cells",
    y = "immune infiltrating cells",
    title = "Gene-immune infiltrating cells Correlation"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.x = element_text(size = 14, face = "bold", color = "black"),
    axis.title.y = element_text(size = 14, face = "bold", color = "black"),
    axis.text.x = element_text(
      size = 12, face = "bold", color = "black", angle = 45, hjust = 1
    ),
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )

figure9a_files <- save_gg_panel(
  figure9a,
  "Figure9A_gene_immune_linkET_BH_q",
  width = 12,
  height = 7
)

# -----------------------------------------------------------------------------
# Figure 9B: pooled Spearman coefficients with BH-adjusted significance
# -----------------------------------------------------------------------------
gene_cor_matrix <- matrix(
  NA_real_,
  nrow = length(gene_names),
  ncol = length(cell_names),
  dimnames = list(gene_names, cell_names)
)
gene_q_matrix <- gene_cor_matrix

for (gene in gene_names) {
  for (cell in cell_names) {
    one_row <- corr_stats[
      corr_stats$Gene == gene & corr_stats$Immune_cell == cell,
      , drop = FALSE
    ]
    if (nrow(one_row) != 1L) stop("Missing or duplicated pair: ", gene, " / ", cell)
    gene_cor_matrix[gene, cell] <- one_row$All_samples_rho
    gene_q_matrix[gene, cell] <- one_row$All_samples_BH_FDR_66_tests
  }
}

star_only <- function(q_value) {
  ifelse(
    is.na(q_value), "",
    ifelse(q_value < 0.001, "***",
      ifelse(q_value < 0.01, "**", ifelse(q_value < 0.05, "*", ""))
    )
  )
}

figure9b_labels <- matrix(
  paste0(
    sprintf("%.2f", as.vector(gene_cor_matrix)),
    "\n",
    star_only(as.vector(gene_q_matrix))
  ),
  nrow = nrow(gene_cor_matrix),
  ncol = ncol(gene_cor_matrix),
  dimnames = dimnames(gene_cor_matrix)
)
figure9b_colors <- colorRampPalette(c("#6b8e23", "white", "#ba6262"))(100)

draw_figure9b <- function() {
  pheatmap::pheatmap(
    gene_cor_matrix,
    color = figure9b_colors,
    display_numbers = figure9b_labels,
    number_color = "black",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    fontsize_number = 7,
    fontsize_row = 11,
    fontsize_col = 10,
    border_color = NA,
    main = "Gene-immune Correlation Heatmap"
  )
}

figure9b_files <- save_grid_panel(
  draw_figure9b,
  "Figure9B_gene_immune_heatmap_BH_q",
  width = 15,
  height = 4
)


capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
print(metadata)
cat("Gene-fraction statistics and two Figure 9 panels completed.\n")

