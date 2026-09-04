## Compare core-gene cell-type localization between GSE191288 and GSE184362
## Run with: Rscript scripts/16_scRNA_cross_dataset_comparison.R [PROJECT_DIR]

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(openxlsx)
})

cli_args <- commandArgs(trailingOnly = TRUE)
project_dir <- normalizePath(
  if (length(cli_args) >= 1L) cli_args[[1L]] else getwd(),
  winslash = "/",
  mustWork = TRUE
)
gse191_dir <- file.path(project_dir, "28.单细胞分析_GSE191288_2026-06-12")
gse184_dir <- file.path(project_dir, "29.单细胞分析_GSE184362_2026-06-12")
result_dir <- file.path(gse184_dir, "results")
figure_dir <- file.path(gse184_dir, "figures")
log_dir <- file.path(gse184_dir, "logs")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

sink(file.path(log_dir, "02_compare_GSE191288_GSE184362_core_gene_patterns.log"), split = TRUE)
cat("Compare GSE191288 and GSE184362 core-gene cell-type patterns\n")
cat("R version:", R.version.string, "\n")
cat("Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")

core_genes <- c("BAX", "BCL2", "FN1")

gse191 <- fread(file.path(gse191_dir, "results", "13_core_gene_expression_by_celltype.csv"))
gse191 <- gse191[condition == "PTC"]
gse191[, dataset := "GSE191288"]
gse191[, tissue_context := "PTC"]
setnames(gse191, "cluster_annotation", "cell_type")
gse191 <- gse191[, .(dataset, tissue_context, gene, cell_type, n_cells, pct_expressing, avg_log_expression, median_log_expression)]

gse184 <- fread(file.path(gse184_dir, "results", "12_core_gene_expression_primary_vs_paratumor_by_celltype.csv"))
gse184 <- gse184[tissue_type == "Primary tumor"]
gse184[, dataset := "GSE184362"]
gse184[, tissue_context := "Primary tumor"]
setnames(gse184, "marker_based_cell_type", "cell_type")
gse184 <- gse184[, .(dataset, tissue_context, gene, cell_type, n_cells, pct_expressing, avg_log_expression, median_log_expression)]

combined <- rbindlist(list(gse191, gse184), fill = TRUE)
combined <- combined[gene %in% core_genes]
combined[, rank_avg_expr := frank(-avg_log_expression, ties.method = "min"), by = .(dataset, gene)]
combined[, rank_pct_expr := frank(-pct_expressing, ties.method = "min"), by = .(dataset, gene)]
combined[, rank_integrated := (rank_avg_expr + rank_pct_expr) / 2]
setorder(combined, dataset, gene, rank_integrated, -avg_log_expression)
fwrite(combined, file.path(result_dir, "15_GSE191288_GSE184362_core_gene_celltype_comparison_all.csv"))

top5 <- combined[rank_integrated <= 5]
fwrite(top5, file.path(result_dir, "16_GSE191288_GSE184362_core_gene_top5_celltypes.csv"))

overlap_rows <- rbindlist(lapply(core_genes, function(g) {
  a <- top5[dataset == "GSE191288" & gene == g, cell_type]
  b <- top5[dataset == "GSE184362" & gene == g, cell_type]
  data.table(
    gene = g,
    GSE191288_top5 = paste(a, collapse = "; "),
    GSE184362_top5 = paste(b, collapse = "; "),
    shared_top_celltypes = paste(intersect(a, b), collapse = "; "),
    n_shared_top_celltypes = length(intersect(a, b))
  )
}))
fwrite(overlap_rows, file.path(result_dir, "17_GSE191288_GSE184362_core_gene_top_celltype_overlap.csv"))
print(overlap_rows)

wb <- createWorkbook()
addWorksheet(wb, "all_ranked_celltypes")
writeData(wb, "all_ranked_celltypes", combined)
addWorksheet(wb, "top5")
writeData(wb, "top5", top5)
addWorksheet(wb, "overlap")
writeData(wb, "overlap", overlap_rows)
saveWorkbook(wb, file.path(result_dir, "GSE191288_GSE184362_core_gene_pattern_comparison.xlsx"), overwrite = TRUE)

plot_dt <- combined[cell_type %in% unique(top5$cell_type)]
plot_dt[, gene := factor(gene, levels = rev(core_genes))]
plot_dt[, dataset := factor(dataset, levels = c("GSE191288", "GSE184362"))]

p <- ggplot(plot_dt, aes(x = cell_type, y = gene)) +
  geom_point(aes(size = pct_expressing, color = avg_log_expression)) +
  facet_grid(dataset ~ ., scales = "free_x", space = "free_x") +
  scale_color_gradient(low = "#D8DEE9", high = "#BF3F3F") +
  scale_size(range = c(0.5, 7)) +
  labs(
    x = NULL,
    y = NULL,
    title = "Core-gene cell-type localization consistency across GSE191288 and GSE184362",
    size = "% expressing",
    color = "Avg log expr"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "grey92")
  )

ggsave(file.path(figure_dir, "07_GSE191288_GSE184362_core_gene_celltype_consistency_dotplot.png"), p, width = 12, height = 7.5, dpi = 320)
ggsave(file.path(figure_dir, "07_GSE191288_GSE184362_core_gene_celltype_consistency_dotplot.pdf"), p, width = 12, height = 7.5)

cat("\nSession info:\n")
print(sessionInfo())
sink()
