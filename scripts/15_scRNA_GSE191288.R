## GSE191288 scRNA-seq analysis for PS-MPs and thyroid cancer project
## Main goal: locate candidate/core genes across PTC single-cell compartments
## Data: GSE191288 processed 10X H5 files from GEO
## Run with: Rscript scripts/15_scRNA_GSE191288.R [PROJECT_DIR]

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 8 * 1024^3)

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
analysis_dir <- file.path(project_dir, "28.单细胞分析_GSE191288_2026-06-12")
h5_dir <- file.path(analysis_dir, "raw", "h5_files")
input_dir <- file.path(analysis_dir, "inputs")
result_dir <- file.path(analysis_dir, "results")
figure_dir <- file.path(analysis_dir, "figures")
log_dir <- file.path(analysis_dir, "logs")

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

sink(file.path(log_dir, "01_seurat_GSE191288_core_gene_analysis.log"), split = TRUE)
cat("GSE191288 Seurat core gene analysis\n")
cat("Project directory:", project_dir, "\n")
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

parse_sample <- function(file) {
  base <- tools::file_path_sans_ext(basename(file))
  gsm <- sub("_.*$", "", base)
  sample_code <- sub("^GSM[0-9]+_", "", base)
  condition <- ifelse(sample_code == "NT", "Non-tumor thyroid", "PTC")
  patient <- ifelse(sample_code == "NT", "P4", paste0("P", sub("^T([0-9]+).*", "\\1", sample_code)))
  lobe <- ifelse(
    grepl("L$", sample_code), "Left",
    ifelse(grepl("R$", sample_code), "Right", "Non-tumor")
  )
  sample_label <- ifelse(
    sample_code == "NT",
    "Non-tumor\nthyroid",
    paste0("PTC-", sub("^P", "", patient), "\n", lobe, " tumor")
  )
  data.table(
    file = file,
    gsm = gsm,
    sample_id = sample_code,
    sample_label = sample_label,
    condition = condition,
    patient = patient,
    lobe = lobe
  )
}

candidate_genes <- read_gene_list(file.path(input_dir, "candidate_7genes.txt"))
core_genes <- read_gene_list(file.path(input_dir, "final_core_genes_BCL2_BAX_FN1.csv"))
genes_of_interest <- unique(c(core_genes, candidate_genes))

cat("Candidate genes:", paste(candidate_genes, collapse = ", "), "\n")
cat("Final core genes:", paste(core_genes, collapse = ", "), "\n\n")

h5_files <- list.files(h5_dir, pattern = "\\.h5$", full.names = TRUE)
if (length(h5_files) == 0) stop("No H5 files found in: ", h5_dir)
sample_meta <- rbindlist(lapply(h5_files, parse_sample))
setorder(sample_meta, sample_id)
fwrite(sample_meta[, !"file"], file.path(result_dir, "01_sample_metadata.csv"))
print(sample_meta)

make_seurat_object <- function(meta_row) {
  cat("Reading", basename(meta_row$file), "...\n")
  counts <- Read10X_h5(meta_row$file, use.names = TRUE, unique.features = TRUE)
  if (is.list(counts)) {
    if ("Gene Expression" %in% names(counts)) {
      counts <- counts[["Gene Expression"]]
    } else {
      counts <- counts[[1]]
    }
  }
  obj <- CreateSeuratObject(
    counts = counts,
    project = meta_row$sample_id,
    min.cells = 3,
    min.features = 200
  )
  obj$gsm <- meta_row$gsm
  obj$sample_id <- meta_row$sample_id
  obj$condition <- meta_row$condition
  obj$patient <- meta_row$patient
  obj$lobe <- meta_row$lobe
  obj
}

objs <- lapply(seq_len(nrow(sample_meta)), function(i) make_seurat_object(sample_meta[i]))
names(objs) <- sample_meta$sample_id

cat("\nMerging samples...\n")
seu <- merge(objs[[1]], y = objs[-1], add.cell.ids = names(objs), project = "GSE191288_PTC")
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
), by = .(sample_id, condition)]
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
), by = .(sample_id, condition)]
fwrite(qc_summary_after, file.path(result_dir, "03_QC_summary_after_filtering.csv"))
print(qc_summary_after)

cat("\nRunning standard Seurat workflow...\n")
seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
seu <- ScaleData(seu, features = VariableFeatures(seu), verbose = FALSE)
seu <- RunPCA(seu, features = VariableFeatures(seu), npcs = 30, verbose = FALSE)
seu <- FindNeighbors(seu, dims = 1:20, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.5, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:20, seed.use = 20260612, verbose = FALSE)

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

score_dt <- as.data.table(seu@meta.data[, c("seurat_clusters", score_cols), drop = FALSE], keep.rownames = "cell")
score_mat <- as.matrix(score_dt[, ..score_cols])
max_idx <- max.col(score_mat, ties.method = "first")
second_score <- apply(score_mat, 1, function(x) sort(x, decreasing = TRUE)[min(2, length(x))])
score_dt[, marker_based_cell_type := names(marker_sets)[max_idx]]
score_dt[, marker_score_max := score_mat[cbind(seq_len(nrow(score_mat)), max_idx)]]
score_dt[, marker_score_second := second_score]
score_dt[, marker_score_margin := marker_score_max - marker_score_second]
score_dt[marker_score_max < 0.05 | marker_score_margin < 0.02, marker_based_cell_type := paste0(marker_based_cell_type, "_low_confidence")]

seu$marker_based_cell_type <- score_dt$marker_based_cell_type[match(colnames(seu), score_dt$cell)]

celltype_counts <- as.data.table(seu@meta.data)[, .N, by = .(sample_id, condition, marker_based_cell_type)]
setorder(celltype_counts, sample_id, -N)
fwrite(celltype_counts, file.path(result_dir, "05_celltype_counts_by_sample.csv"))

cluster_celltype <- as.data.table(seu@meta.data)[, .N, by = .(seurat_clusters, marker_based_cell_type)]
cluster_celltype[, cluster_total := sum(N), by = seurat_clusters]
cluster_celltype[, fraction := N / cluster_total]
setorder(cluster_celltype, seurat_clusters, -N)
cluster_annotation <- cluster_celltype[, .SD[1], by = seurat_clusters]
setnames(cluster_annotation, "marker_based_cell_type", "cluster_annotation")
fwrite(cluster_celltype, file.path(result_dir, "06_cluster_celltype_composition.csv"))
fwrite(cluster_annotation, file.path(result_dir, "07_cluster_annotation_majority_marker_score.csv"))

seu$cluster_annotation <- cluster_annotation$cluster_annotation[match(as.character(seu$seurat_clusters), as.character(cluster_annotation$seurat_clusters))]
Idents(seu) <- "cluster_annotation"

genes_present <- intersect(genes_of_interest, rownames(seu))
genes_missing <- setdiff(genes_of_interest, rownames(seu))
fwrite(data.table(gene = genes_present), file.path(result_dir, "08_genes_of_interest_present.csv"))
fwrite(data.table(gene = genes_missing), file.path(result_dir, "09_genes_of_interest_missing.csv"))
cat("\nGenes of interest present:", paste(genes_present, collapse = ", "), "\n")
if (length(genes_missing) > 0) cat("Genes missing:", paste(genes_missing, collapse = ", "), "\n")

fetch_vars <- unique(c("sample_id", "condition", "patient", "lobe", "seurat_clusters", "marker_based_cell_type", "cluster_annotation", genes_present))
expr_df <- FetchData(seu, vars = fetch_vars)
expr_dt <- as.data.table(expr_df, keep.rownames = "cell")
expr_long <- melt(
  expr_dt,
  id.vars = c("cell", "sample_id", "condition", "patient", "lobe", "seurat_clusters", "marker_based_cell_type", "cluster_annotation"),
  measure.vars = genes_present,
  variable.name = "gene",
  value.name = "log_normalized_expression"
)
expr_long[, expressed := log_normalized_expression > 0]

summary_by_celltype <- expr_long[, .(
  n_cells = .N,
  pct_expressing = round(mean(expressed) * 100, 3),
  avg_log_expression = round(mean(log_normalized_expression), 5),
  median_log_expression = round(median(log_normalized_expression), 5)
), by = .(gene, condition, cluster_annotation)]
setorder(summary_by_celltype, gene, condition, cluster_annotation)
fwrite(summary_by_celltype, file.path(result_dir, "10_gene_expression_by_condition_celltype.csv"))

summary_by_sample_celltype <- expr_long[, .(
  n_cells = .N,
  pct_expressing = round(mean(expressed) * 100, 3),
  avg_log_expression = round(mean(log_normalized_expression), 5),
  median_log_expression = round(median(log_normalized_expression), 5)
), by = .(gene, sample_id, condition, patient, lobe, cluster_annotation)]
setorder(summary_by_sample_celltype, gene, sample_id, cluster_annotation)
fwrite(summary_by_sample_celltype, file.path(result_dir, "11_gene_expression_by_sample_celltype.csv"))

overall_gene_summary <- expr_long[, .(
  n_cells = .N,
  pct_expressing = round(mean(expressed) * 100, 3),
  avg_log_expression = round(mean(log_normalized_expression), 5),
  median_log_expression = round(median(log_normalized_expression), 5)
), by = .(gene, condition)]
fwrite(overall_gene_summary, file.path(result_dir, "12_gene_expression_overall_by_condition.csv"))

core_summary <- summary_by_celltype[gene %in% core_genes]
candidate_summary <- summary_by_celltype[gene %in% candidate_genes]
fwrite(core_summary, file.path(result_dir, "13_core_gene_expression_by_celltype.csv"))
fwrite(candidate_summary, file.path(result_dir, "14_candidate_7gene_expression_by_celltype.csv"))

wb <- createWorkbook()
addWorksheet(wb, "sample_metadata")
writeData(wb, "sample_metadata", sample_meta[, !"file"])
addWorksheet(wb, "QC_after")
writeData(wb, "QC_after", qc_summary_after)
addWorksheet(wb, "celltype_counts")
writeData(wb, "celltype_counts", celltype_counts)
addWorksheet(wb, "cluster_annotation")
writeData(wb, "cluster_annotation", cluster_annotation)
addWorksheet(wb, "core_gene_celltype")
writeData(wb, "core_gene_celltype", core_summary)
addWorksheet(wb, "candidate_gene_celltype")
writeData(wb, "candidate_gene_celltype", candidate_summary)
saveWorkbook(wb, file.path(result_dir, "GSE191288_core_gene_single_cell_summary.xlsx"), overwrite = TRUE)

cat("\nSaving figures...\n")
png(file.path(figure_dir, "01_UMAP_by_condition.png"), width = 2400, height = 1800, res = 300)
print(DimPlot(seu, reduction = "umap", group.by = "condition", pt.size = 0.2) + ggtitle("GSE191288 UMAP by condition"))
dev.off()

png(file.path(figure_dir, "02_UMAP_by_sample.png"), width = 2600, height = 1900, res = 300)
print(DimPlot(seu, reduction = "umap", group.by = "sample_id", pt.size = 0.2) + ggtitle("GSE191288 UMAP by sample"))
dev.off()

png(file.path(figure_dir, "03_UMAP_by_cluster_annotation.png"), width = 2800, height = 2100, res = 300)
print(DimPlot(seu, reduction = "umap", group.by = "cluster_annotation", label = TRUE, repel = TRUE, pt.size = 0.18) + ggtitle("Marker-based cell type annotation"))
dev.off()

if (length(core_genes) > 0) {
  present_core <- intersect(core_genes, rownames(seu))
  if (length(present_core) > 0) {
    png(file.path(figure_dir, "04_UMAP_featureplots_core_genes.png"), width = 3200, height = 2100, res = 300)
    print(FeaturePlot(seu, features = present_core, reduction = "umap", ncol = length(present_core), pt.size = 0.15, order = TRUE))
    dev.off()
  }
}

plot_dot <- function(dat, genes, file_stub, title_text) {
  dat <- copy(dat[gene %in% genes])
  if (nrow(dat) == 0) return(invisible(NULL))
  dat[, gene := factor(gene, levels = rev(genes))]
  dat[, cluster_annotation := factor(cluster_annotation, levels = unique(cluster_annotation))]
  p <- ggplot(dat, aes(x = cluster_annotation, y = gene)) +
    geom_point(aes(size = pct_expressing, color = avg_log_expression)) +
    facet_wrap(~ condition, ncol = 1) +
    scale_color_gradient(low = "#D8DEE9", high = "#BF3F3F") +
    scale_size(range = c(0.5, 7)) +
    labs(x = NULL, y = NULL, title = title_text, size = "% expressing", color = "Avg log expr") +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_line(color = "grey92")
    )
  ggsave(file.path(figure_dir, paste0(file_stub, ".png")), p, width = 11, height = 7.5, dpi = 320)
  ggsave(file.path(figure_dir, paste0(file_stub, ".pdf")), p, width = 11, height = 7.5)
}

plot_dot(summary_by_celltype, genes_present, "05_dotplot_candidate_and_core_genes_by_celltype_condition", "Candidate/core gene expression across single-cell compartments")
plot_dot(summary_by_celltype, intersect(core_genes, genes_present), "06_dotplot_core_genes_by_celltype_condition", "BAX/BCL2/FN1 expression across single-cell compartments")

composition <- as.data.table(seu@meta.data)[, .N, by = .(sample_id, condition, cluster_annotation)]
composition[, fraction := N / sum(N), by = sample_id]
composition <- merge(composition, sample_meta[, .(sample_id, sample_label)], by = "sample_id", all.x = TRUE)
composition[, sample_label := factor(sample_label, levels = sample_meta$sample_label)]
p_comp <- ggplot(composition, aes(x = sample_label, y = fraction, fill = cluster_annotation)) +
  geom_col(width = 0.82) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = NULL, y = "Cell fraction", fill = "Cell type", title = "Cell type composition by sample") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.text.x = element_text(angle = 0, hjust = 0.5, lineheight = 0.95))
ggsave(file.path(figure_dir, "07_celltype_composition_by_sample.png"), p_comp, width = 10, height = 6, dpi = 320)
ggsave(file.path(figure_dir, "07_celltype_composition_by_sample.pdf"), p_comp, width = 10, height = 6)

present_core <- intersect(core_genes, rownames(seu))
if (length(present_core) > 0) {
  p_vln <- VlnPlot(seu, features = present_core, group.by = "cluster_annotation", split.by = "condition", pt.size = 0, ncol = 1) &
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(figure_dir, "08_violin_core_genes_by_celltype_condition.png"), p_vln, width = 12, height = 10, dpi = 320)
  ggsave(file.path(figure_dir, "08_violin_core_genes_by_celltype_condition.pdf"), p_vln, width = 12, height = 10)
}

cat("\nFinding cluster markers on downsampled cells for annotation audit...\n")
Idents(seu) <- "cluster_annotation"
markers <- tryCatch({
  FindAllMarkers(
    seu,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    max.cells.per.ident = 500,
    random.seed = 20260612,
    verbose = FALSE
  )
}, error = function(e) {
  warning("FindAllMarkers failed: ", conditionMessage(e))
  data.frame()
})
if (nrow(markers) > 0) {
  fwrite(as.data.table(markers), file.path(result_dir, "15_cluster_markers_downsampled_for_annotation_audit.csv"))
}

saveRDS(seu, file.path(result_dir, "GSE191288_seurat_processed_marker_annotated.rds"))

cat("\nAnalysis complete.\n")
cat("Cells retained:", ncol(seu), "\n")
cat("Clusters:", length(unique(seu$seurat_clusters)), "\n")
cat("Annotated cell compartments:", paste(sort(unique(seu$cluster_annotation)), collapse = ", "), "\n")

cat("\nSession info:\n")
print(sessionInfo())
sink()
