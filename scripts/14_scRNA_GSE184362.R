## GSE184362 scRNA-seq analysis for PS-MPs and thyroid cancer project
## Main goal: validate whether BAX/BCL2/FN1 cell-type localization is consistent
## with the smaller GSE191288 single-cell analysis.
## Run with: Rscript scripts/14_scRNA_GSE184362.R [PROJECT_DIR]

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 12 * 1024^3)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(openxlsx)
})

set.seed(20260612)

cli_args <- commandArgs(trailingOnly = TRUE)
project_dir <- normalizePath(
  if (length(cli_args) >= 1L) cli_args[[1L]] else getwd(),
  winslash = "/",
  mustWork = TRUE
)
analysis_dir <- file.path(project_dir, "29.单细胞分析_GSE184362_2026-06-12")
mtx_dir <- file.path(analysis_dir, "raw", "mtx_files")
input_dir <- file.path(analysis_dir, "inputs")
result_dir <- file.path(analysis_dir, "results")
figure_dir <- file.path(analysis_dir, "figures")
log_dir <- file.path(analysis_dir, "logs")

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

sink(file.path(log_dir, "01_seurat_GSE184362_core_gene_analysis.log"), split = TRUE)
cat("GSE184362 Seurat core gene analysis\n")
cat("Dataset: GSE184362\n")
cat("Analysis time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("R version:", R.version.string, "\n")
cat("Seurat version:", as.character(packageVersion("Seurat")), "\n\n")

read_gene_list <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("txt", "list")) {
    genes <- readLines(path, warn = FALSE)
  } else {
    x <- fread(path, header = TRUE)
    if ("Gene" %in% names(x)) {
      genes <- x$Gene
    } else {
      genes <- unlist(x[[1]], use.names = FALSE)
    }
  }
  genes <- trimws(genes)
  genes <- genes[genes != "" & toupper(genes) != "GENE"]
  sort(unique(toupper(genes)))
}

parse_sample_name <- function(prefix) {
  sample_part <- sub("^GSM[0-9]+_", "", prefix)
  gsm <- sub("_.*$", "", prefix)
  patient <- sub("^(PTC[0-9]+).*", "\\1", sample_part)
  tissue_code <- sub("^PTC[0-9]+_", "", sample_part)
  tissue_type <- ifelse(
    tissue_code == "T", "Primary tumor",
    ifelse(tissue_code == "P", "Paratumor",
           ifelse(grepl("LN$", tissue_code), "Lymph node metastasis",
                  ifelse(tissue_code == "SC", "Subcutaneous metastasis", tissue_code)))
  )
  tissue_label <- ifelse(
    tissue_code == "T", "Tumor",
    ifelse(tissue_code == "P", "Paratumor",
           ifelse(tissue_code == "LeftLN", "Left LN met.",
                  ifelse(tissue_code == "RightLN", "Right LN met.",
                         ifelse(tissue_code == "SC", "SC met.", tissue_code))))
  )
  analysis_scope <- ifelse(tissue_type %in% c("Primary tumor", "Paratumor"),
                           "Primary_vs_paratumor", "Metastatic_lesion")
  data.table(
    prefix = prefix,
    gsm = gsm,
    sample_id = sample_part,
    sample_label = paste(patient, tissue_label),
    patient = patient,
    tissue_code = tissue_code,
    tissue_type = tissue_type,
    analysis_scope = analysis_scope
  )
}

make_unique_features <- function(feature_names) {
  feature_names <- trimws(feature_names)
  feature_names[feature_names == ""] <- paste0("unnamed_", which(feature_names == ""))
  make.unique(feature_names)
}

read_10x_prefixed <- function(prefix) {
  matrix_path <- file.path(mtx_dir, paste0(prefix, "_matrix.mtx.gz"))
  features_path <- file.path(mtx_dir, paste0(prefix, "_features.tsv.gz"))
  barcodes_path <- file.path(mtx_dir, paste0(prefix, "_barcodes.tsv.gz"))
  if (!file.exists(matrix_path) || !file.exists(features_path) || !file.exists(barcodes_path)) {
    stop("Missing 10X files for prefix: ", prefix)
  }
  counts <- readMM(gzfile(matrix_path))
  features <- fread(features_path, header = FALSE)
  barcodes <- fread(barcodes_path, header = FALSE)
  gene_col <- if (ncol(features) >= 2) 2 else 1
  rownames(counts) <- make_unique_features(features[[gene_col]])
  colnames(counts) <- paste(prefix, barcodes[[1]], sep = "_")
  counts
}

candidate_genes <- read_gene_list(file.path(input_dir, "candidate_7genes.txt"))
core_genes <- read_gene_list(file.path(input_dir, "final_core_genes_BCL2_BAX_FN1.csv"))
genes_of_interest <- unique(c(core_genes, candidate_genes))
cat("Candidate genes:", paste(candidate_genes, collapse = ", "), "\n")
cat("Final core genes:", paste(core_genes, collapse = ", "), "\n\n")

matrix_files <- list.files(mtx_dir, pattern = "_matrix\\.mtx\\.gz$", full.names = TRUE)
prefixes <- sub("_matrix\\.mtx\\.gz$", "", basename(matrix_files))
if (length(prefixes) == 0L) stop("No matrix files found in: ", mtx_dir)
sample_meta <- rbindlist(lapply(prefixes, parse_sample_name), fill = TRUE)
setorder(sample_meta, patient, tissue_type, sample_id)
fwrite(sample_meta, file.path(result_dir, "01_sample_metadata.csv"))
print(sample_meta)

make_seurat_object <- function(meta_row) {
  cat("Reading", meta_row$prefix, "...\n")
  counts <- read_10x_prefixed(meta_row$prefix)
  obj <- CreateSeuratObject(
    counts = counts,
    project = meta_row$sample_id,
    min.cells = 3,
    min.features = 200
  )
  obj$gsm <- meta_row$gsm
  obj$sample_id <- meta_row$sample_id
  obj$patient <- meta_row$patient
  obj$tissue_code <- meta_row$tissue_code
  obj$tissue_type <- meta_row$tissue_type
  obj$analysis_scope <- meta_row$analysis_scope
  obj
}

objs <- lapply(seq_len(nrow(sample_meta)), function(i) make_seurat_object(sample_meta[i]))
names(objs) <- sample_meta$sample_id

cat("\nMerging samples...\n")
seu <- merge(objs[[1]], y = objs[-1], add.cell.ids = names(objs), project = "GSE184362_PTC")
seu <- JoinLayers(seu)
rm(objs)
gc()

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
qc_before <- as.data.table(seu@meta.data, keep.rownames = "cell")
qc_summary_before <- qc_before[, .(
  n_cells = .N,
  median_nCount_RNA = median(nCount_RNA),
  median_nFeature_RNA = median(nFeature_RNA),
  median_percent_mt = median(percent.mt)
), by = .(sample_id, tissue_type, analysis_scope)]
setorder(qc_summary_before, analysis_scope, tissue_type, sample_id)
fwrite(qc_summary_before, file.path(result_dir, "02_QC_summary_before_filtering.csv"))
print(qc_summary_before)

cat("\nFiltering cells: nFeature_RNA >= 200, nCount_RNA >= 500, percent.mt <= 25\n")
seu <- subset(seu, subset = nFeature_RNA >= 200 & nCount_RNA >= 500 & percent.mt <= 25)

qc_after <- as.data.table(seu@meta.data, keep.rownames = "cell")
qc_summary_after <- qc_after[, .(
  n_cells = .N,
  median_nCount_RNA = median(nCount_RNA),
  median_nFeature_RNA = median(nFeature_RNA),
  median_percent_mt = median(percent.mt)
), by = .(sample_id, tissue_type, analysis_scope)]
setorder(qc_summary_after, analysis_scope, tissue_type, sample_id)
fwrite(qc_summary_after, file.path(result_dir, "03_QC_summary_after_filtering.csv"))
print(qc_summary_after)
cat("\nCells retained after QC:", ncol(seu), "\n")

cat("\nNormalizing full object for marker-score annotation and gene-expression summaries...\n")
seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)

marker_sets <- list(
  Thyroid_epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "TG", "TPO", "SLC5A5", "PAX8", "FOXE1", "CLDN4"),
  T_NK_cells = c("PTPRC", "CD3D", "CD3E", "TRAC", "CD2", "CD7", "NKG7", "GNLY", "KLRD1", "GZMB"),
  B_cells = c("PTPRC", "MS4A1", "CD79A", "CD79B", "CD74", "BANK1", "CD19"),
  Plasma_cells = c("MZB1", "JCHAIN", "SDC1", "XBP1", "IGHG1", "IGKC"),
  Myeloid_cells = c("PTPRC", "LYZ", "LST1", "TYROBP", "CD68", "C1QA", "C1QB", "FCGR3A", "CSF1R"),
  Dendritic_cells = c("PTPRC", "FCER1A", "CLEC10A", "ITGAX", "LILRA4", "IRF7"),
  Mast_cells = c("TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2"),
  Endothelial_cells = c("PECAM1", "VWF", "KDR", "CDH5", "EMCN", "CLDN5", "RAMP2"),
  Fibroblasts = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FAP", "PDGFRA", "THY1"),
  Pericyte_smooth_muscle = c("ACTA2", "MYH11", "RGS5", "MCAM", "CSPG4", "PDGFRB", "TAGLN"),
  Proliferating_cells = c("MKI67", "TOP2A", "UBE2C", "STMN1", "HMGB2", "TYMS")
)
marker_sets <- lapply(marker_sets, function(x) intersect(x, rownames(seu)))
marker_sets <- marker_sets[lengths(marker_sets) > 0]
cat("\nMarker sets used for annotation:\n")
print(lapply(marker_sets, length))

seu <- AddModuleScore(seu, features = marker_sets, name = "marker_score_", seed = 20260612, search = FALSE)
score_cols <- paste0("marker_score_", seq_along(marker_sets))
score_name_map <- data.table(score_col = score_cols, marker_set = names(marker_sets))
fwrite(score_name_map, file.path(result_dir, "04_marker_score_column_map.csv"))

score_dt <- as.data.table(seu@meta.data[, c(score_cols), drop = FALSE], keep.rownames = "cell")
score_mat <- as.matrix(score_dt[, ..score_cols])
max_idx <- max.col(score_mat, ties.method = "first")
second_score <- apply(score_mat, 1, function(x) sort(x, decreasing = TRUE)[min(2, length(x))])
score_dt[, marker_based_cell_type := names(marker_sets)[max_idx]]
score_dt[, marker_score_max := score_mat[cbind(seq_len(nrow(score_mat)), max_idx)]]
score_dt[, marker_score_second := second_score]
score_dt[, marker_score_margin := marker_score_max - marker_score_second]
score_dt[, marker_annotation_confidence := ifelse(marker_score_max < 0.05 | marker_score_margin < 0.02, "low_confidence", "confident")]

seu$marker_based_cell_type <- score_dt$marker_based_cell_type[match(colnames(seu), score_dt$cell)]
seu$marker_annotation_confidence <- score_dt$marker_annotation_confidence[match(colnames(seu), score_dt$cell)]
seu$marker_based_cell_type_display <- ifelse(
  seu$marker_annotation_confidence == "low_confidence",
  paste0(seu$marker_based_cell_type, "_low_confidence"),
  seu$marker_based_cell_type
)

celltype_counts <- as.data.table(seu@meta.data)[, .N, by = .(sample_id, tissue_type, analysis_scope, marker_based_cell_type, marker_annotation_confidence)]
setorder(celltype_counts, analysis_scope, tissue_type, sample_id, -N)
fwrite(celltype_counts, file.path(result_dir, "05_celltype_counts_by_sample.csv"))

genes_present <- intersect(genes_of_interest, rownames(seu))
genes_missing <- setdiff(genes_of_interest, rownames(seu))
fwrite(data.table(gene = genes_present), file.path(result_dir, "06_genes_of_interest_present.csv"))
fwrite(data.table(gene = genes_missing), file.path(result_dir, "07_genes_of_interest_missing.csv"))
cat("\nGenes of interest present:", paste(genes_present, collapse = ", "), "\n")
if (length(genes_missing) > 0) cat("Genes missing:", paste(genes_missing, collapse = ", "), "\n")

fetch_vars <- unique(c("sample_id", "patient", "tissue_type", "analysis_scope",
                       "marker_based_cell_type", "marker_annotation_confidence",
                       "marker_based_cell_type_display", genes_present))
expr_df <- FetchData(seu, vars = fetch_vars)
expr_dt <- as.data.table(expr_df, keep.rownames = "cell")
expr_long <- melt(
  expr_dt,
  id.vars = c("cell", "sample_id", "patient", "tissue_type", "analysis_scope",
              "marker_based_cell_type", "marker_annotation_confidence", "marker_based_cell_type_display"),
  measure.vars = genes_present,
  variable.name = "gene",
  value.name = "log_normalized_expression"
)
expr_long[, expressed := log_normalized_expression > 0]

summarize_expr <- function(dt, by_cols) {
  dt[, .(
    n_cells = .N,
    pct_expressing = round(mean(expressed) * 100, 3),
    avg_log_expression = round(mean(log_normalized_expression), 5),
    median_log_expression = round(median(log_normalized_expression), 5)
  ), by = by_cols]
}

summary_by_tissue_celltype <- summarize_expr(
  expr_long,
  c("gene", "tissue_type", "analysis_scope", "marker_based_cell_type")
)
setorder(summary_by_tissue_celltype, gene, analysis_scope, tissue_type, marker_based_cell_type)
fwrite(summary_by_tissue_celltype, file.path(result_dir, "08_gene_expression_by_tissue_celltype_all_cells.csv"))

summary_confident <- summarize_expr(
  expr_long[marker_annotation_confidence == "confident"],
  c("gene", "tissue_type", "analysis_scope", "marker_based_cell_type")
)
setorder(summary_confident, gene, analysis_scope, tissue_type, marker_based_cell_type)
fwrite(summary_confident, file.path(result_dir, "09_gene_expression_by_tissue_celltype_confident_cells.csv"))

summary_by_sample_celltype <- summarize_expr(
  expr_long,
  c("gene", "sample_id", "patient", "tissue_type", "analysis_scope", "marker_based_cell_type")
)
setorder(summary_by_sample_celltype, gene, analysis_scope, tissue_type, sample_id, marker_based_cell_type)
fwrite(summary_by_sample_celltype, file.path(result_dir, "10_gene_expression_by_sample_celltype.csv"))

overall_by_tissue <- summarize_expr(expr_long, c("gene", "tissue_type", "analysis_scope"))
setorder(overall_by_tissue, gene, analysis_scope, tissue_type)
fwrite(overall_by_tissue, file.path(result_dir, "11_gene_expression_overall_by_tissue.csv"))

primary_paratumor_summary <- summary_by_tissue_celltype[analysis_scope == "Primary_vs_paratumor"]
core_primary_paratumor <- primary_paratumor_summary[gene %in% core_genes]
candidate_primary_paratumor <- primary_paratumor_summary[gene %in% candidate_genes]
fwrite(core_primary_paratumor, file.path(result_dir, "12_core_gene_expression_primary_vs_paratumor_by_celltype.csv"))
fwrite(candidate_primary_paratumor, file.path(result_dir, "13_candidate_7gene_expression_primary_vs_paratumor_by_celltype.csv"))

all_lesion_core <- summary_by_tissue_celltype[gene %in% core_genes]
fwrite(all_lesion_core, file.path(result_dir, "14_core_gene_expression_all_tissues_by_celltype.csv"))

cat("\nCreating balanced downsample object for UMAP visualization...\n")
meta_dt <- as.data.table(seu@meta.data, keep.rownames = "cell")
set.seed(20260612)
max_cells_per_sample_for_umap <- 4000L
umap_cells <- meta_dt[, {
  n_take <- min(.N, max_cells_per_sample_for_umap)
  .(cell = sample(cell, n_take))
}, by = sample_id]$cell
cat("UMAP downsample cells:", length(umap_cells), "\n")
seu_umap <- subset(seu, cells = umap_cells)
seu_umap <- FindVariableFeatures(seu_umap, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
seu_umap <- ScaleData(seu_umap, features = VariableFeatures(seu_umap), verbose = FALSE)
seu_umap <- RunPCA(seu_umap, features = VariableFeatures(seu_umap), npcs = 30, verbose = FALSE)
seu_umap <- FindNeighbors(seu_umap, dims = 1:20, verbose = FALSE)
seu_umap <- FindClusters(seu_umap, resolution = 0.5, verbose = FALSE)
seu_umap <- RunUMAP(seu_umap, dims = 1:20, seed.use = 20260612, verbose = FALSE)

cat("\nSaving figures...\n")
png(file.path(figure_dir, "01_UMAP_downsample_by_tissue_type.png"), width = 2600, height = 1900, res = 300)
print(DimPlot(seu_umap, reduction = "umap", group.by = "tissue_type", pt.size = 0.12) + ggtitle("GSE184362 downsampled UMAP by tissue type"))
dev.off()

png(file.path(figure_dir, "02_UMAP_downsample_by_marker_cell_type.png"), width = 2800, height = 2100, res = 300)
print(DimPlot(seu_umap, reduction = "umap", group.by = "marker_based_cell_type", label = TRUE, repel = TRUE, pt.size = 0.12) + ggtitle("GSE184362 marker-based cell type annotation"))
dev.off()

present_core <- intersect(core_genes, rownames(seu_umap))
if (length(present_core) > 0) {
  png(file.path(figure_dir, "03_UMAP_featureplots_core_genes_downsample.png"), width = 3200, height = 2100, res = 300)
  print(FeaturePlot(seu_umap, features = present_core, reduction = "umap", ncol = length(present_core), pt.size = 0.08, order = TRUE))
  dev.off()
}

plot_dot <- function(dat, genes, file_stub, title_text, tissue_filter = NULL) {
  dat <- copy(dat[gene %in% genes])
  if (!is.null(tissue_filter)) dat <- dat[tissue_type %in% tissue_filter]
  if (nrow(dat) == 0) return(invisible(NULL))
  dat[, gene := factor(gene, levels = rev(genes))]
  tissue_levels <- c("Paratumor", "Primary tumor", "Lymph node metastasis", "Subcutaneous metastasis")
  dat[, tissue_type := factor(tissue_type, levels = intersect(tissue_levels, unique(dat$tissue_type)))]
  dat[, marker_based_cell_type := factor(marker_based_cell_type, levels = unique(marker_based_cell_type))]
  p <- ggplot(dat, aes(x = marker_based_cell_type, y = gene)) +
    geom_point(aes(size = pct_expressing, color = avg_log_expression)) +
    facet_wrap(~ tissue_type, ncol = 1) +
    scale_color_gradient(low = "#D8DEE9", high = "#BF3F3F") +
    scale_size(range = c(0.5, 7)) +
    labs(x = NULL, y = NULL, title = title_text, size = "% expressing", color = "Avg log expr") +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_line(color = "grey92")
    )
  ggsave(file.path(figure_dir, paste0(file_stub, ".png")), p, width = 11.5, height = 8.5, dpi = 320)
  ggsave(file.path(figure_dir, paste0(file_stub, ".pdf")), p, width = 11.5, height = 8.5)
}

plot_dot(
  summary_by_tissue_celltype,
  genes_present,
  "04_dotplot_candidate_and_core_genes_by_celltype_tissue_all_cells",
  "GSE184362 candidate/core gene expression across cell compartments"
)
plot_dot(
  summary_by_tissue_celltype,
  intersect(core_genes, genes_present),
  "05_dotplot_core_genes_primary_vs_paratumor_by_celltype",
  "GSE184362 BAX/BCL2/FN1 expression in primary tumor and paratumor",
  tissue_filter = c("Primary tumor", "Paratumor")
)

composition <- as.data.table(seu@meta.data)[, .N, by = .(sample_id, tissue_type, marker_based_cell_type)]
composition[, fraction := N / sum(N), by = sample_id]
composition <- merge(composition, sample_meta[, .(sample_id, sample_label)], by = "sample_id", all.x = TRUE)
composition[, sample_label := factor(sample_label, levels = sample_meta$sample_label)]
p_comp <- ggplot(composition, aes(x = sample_label, y = fraction, fill = marker_based_cell_type)) +
  geom_col(width = 0.82) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = NULL, y = "Cell fraction", fill = "Cell type", title = "GSE184362 marker-based cell type composition by sample") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.text.x = element_text(angle = 55, hjust = 1, size = 8))
ggsave(file.path(figure_dir, "06_celltype_composition_by_sample.png"), p_comp, width = 13, height = 6.8, dpi = 320)
ggsave(file.path(figure_dir, "06_celltype_composition_by_sample.pdf"), p_comp, width = 13, height = 6.8)

wb <- createWorkbook()
addWorksheet(wb, "sample_metadata")
writeData(wb, "sample_metadata", sample_meta)
addWorksheet(wb, "QC_after")
writeData(wb, "QC_after", qc_summary_after)
addWorksheet(wb, "celltype_counts")
writeData(wb, "celltype_counts", celltype_counts)
addWorksheet(wb, "core_primary_paratumor")
writeData(wb, "core_primary_paratumor", core_primary_paratumor)
addWorksheet(wb, "candidate_primary_paratumor")
writeData(wb, "candidate_primary_paratumor", candidate_primary_paratumor)
addWorksheet(wb, "core_all_tissues")
writeData(wb, "core_all_tissues", all_lesion_core)
addWorksheet(wb, "overall_by_tissue")
writeData(wb, "overall_by_tissue", overall_by_tissue)
saveWorkbook(wb, file.path(result_dir, "GSE184362_core_gene_single_cell_summary.xlsx"), overwrite = TRUE)

cat("\nSaving RDS objects...\n")
saveRDS(seu_umap, file.path(result_dir, "GSE184362_downsampled_umap_marker_annotated.rds"))
saveRDS(seu@meta.data, file.path(result_dir, "GSE184362_full_metadata_marker_annotated.rds"))

cat("\nAnalysis complete.\n")
cat("Full cells retained:", ncol(seu), "\n")
cat("UMAP downsample cells:", ncol(seu_umap), "\n")
cat("Marker-based cell compartments:", paste(sort(unique(seu$marker_based_cell_type)), collapse = ", "), "\n")

cat("\nSession info:\n")
print(sessionInfo())
sink()
